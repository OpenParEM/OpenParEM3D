////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//    OpenParEM3D - A fullwave 3D electromagnetic simulator.                  //
//    Copyright (C) 2022 Brian Young                                          //
//                                                                            //
//    This program is free software: you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation, either version 3 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <complex>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <chrono>
#include <unistd.h>
#include <lapacke.h>
#include <HYPRE.h>
#include "port.hpp"
#include "fem3D.hpp"
#include "results.hpp"
#include "project.h"
#include "jobrelated.hpp"
#include "OpenParEMmaterials.hpp"
#include "frequencyPlan.hpp"
#include "license.hpp"
#include "mesh.hpp"
#include "petscErrorHandler.hpp"

using namespace std;
using namespace mfem;

extern "C" void init_project (struct projectData *);
extern "C" int load_project_file(const char*, projectData*, const char*);
extern "C" void print_project (struct projectData *, struct projectData *, const char *);
extern "C" void free_project(projectData*);
extern "C" char* get_project_name (const char *filename);

double edgeLength (double *a, double *b)
{
   return sqrt((a[0]-b[0])*(a[0]-b[0])+
               (a[1]-b[1])*(a[1]-b[1])+
               (a[2]-b[2])*(a[2]-b[2]));
}

void get_edgeLengthsPerWavelength (ParMesh *pmesh, double frequency, Vector er, Vector mur, double *shortest, double *longest)
{
   double eps0=8.8541878176e-12;

   *shortest=DBL_MAX;
   *longest=0;

   int i=0;
   while (i < pmesh->GetNE()) {

      int attribute=pmesh->GetAttribute(i);
      double velocity=1/sqrt(eps0*er[attribute-1]*4e-7*M_PI*mur[attribute-1]);
      double wavelength=velocity/frequency;

      Array<int> vertices;
      pmesh->GetElementVertices(i,vertices);

      int j=0;
      while (j < vertices.Size()) {

         double *coordsj=pmesh->GetVertex(vertices[j]);

         int k=0;
         while (k < vertices.Size()) {
            if (k != j) {
               double *coordsk=pmesh->GetVertex(vertices[k]);
               double length=edgeLength(coordsj,coordsk)/wavelength;
               if (length > *longest) *longest=length;
               if (length < *shortest) *shortest=length;
            }
            k++;
         }
         j++;
      }
      i++;
   }
}

void help () {
   PetscPrintf(PETSC_COMM_WORLD,"usage: OpenParEM3D [-h] filename\n");
   PetscPrintf(PETSC_COMM_WORLD,"       -h          : Print this help text\n");
   PetscPrintf(PETSC_COMM_WORLD,"       filename    : Filename of an OpenParEM setup file.\n");
   PetscPrintf(PETSC_COMM_WORLD,"\nOpenParEM3D is a full-wave 3D electromagnetic solver.\n");
}

bool print_mesh_quality_message (ParMesh *pmesh, struct projectData *projData)
{
   double h_min,h_max,kappa_min,kappa_max;
   pmesh->GetCharacteristics(h_min, h_max, kappa_min, kappa_max);
   PetscPrintf(PETSC_COMM_WORLD,"|      mesh worst element aspect ratio: %g",kappa_max);
   if (kappa_max < 5) PetscPrintf(PETSC_COMM_WORLD," < target: 5\n");
   else PetscPrintf(PETSC_COMM_WORLD," > target: 5\n");

   if (kappa_max > projData->mesh_quality_limit) {
      PetscPrintf(PETSC_COMM_WORLD,"|         ERROR3005: aspect ratio > limit: %g\n",projData->mesh_quality_limit);
      return true;
   }

   return false;
}


void isOpenParEM2Dreachable ()
{
   PetscMPIInt rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   if (rank == 0) {
      if (system("which OpenParEM2D > /dev/null")) {
         PetscPrintf(PETSC_COMM_WORLD,"ERROR3008: OpenParEM2D not found.\n");
         PetscFinalize();
         exit(1);
      }
   }
   MPI_Barrier(PETSC_COMM_WORLD);
}

void load_project_file (const char *projFile, struct projectData *defaultData, struct projectData *projData, char *lockfile, chrono::system_clock::time_point job_start_time)
{
   PetscPrintf(PETSC_COMM_WORLD,"   loading project file \"%s\"\n",projFile);

   init_project (defaultData);
   init_project (projData);

   if (load_project_file (projFile,projData,"   ")) {
      if (projData->debug_show_project) {print_project (projData,defaultData,"      ");}
      exit_job_on_error (job_start_time,lockfile,true);
   }
   if (projData->debug_show_project) {print_project (projData,defaultData,"      ");}
}

void delete_stale_files (const char *baseName)
{
   PetscMPIInt rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   if (rank == 0) {
      delete_file(baseName,"","_results.csv");
      delete_file(baseName,"","_results.txt");
      delete_file(baseName,"","_iterations.txt");
      delete_file(baseName,"","_protoype_test_cases.csv");
      delete_file(baseName,"","_fields.csv");
      delete_file(baseName,"","_attributes.csv");
      delete_file(baseName,"temp_","");
      delete_file(baseName,"ParaView_","");
      delete_file(baseName,"ParaView_2D_port_","");
      delete_file(baseName,"ParaView_3D_port_","");
      delete_file(baseName,"ParaView_modal_2D_","");
      delete_file(baseName,"ParaView_solution_2D_","");
   }
   MPI_Barrier(PETSC_COMM_WORLD);
}

bool saveSerialMesh (struct projectData *projData, meshMaterialList *meshMaterials, ParMesh *pmesh)
{
   stringstream ssMesh;
   ssMesh << "refined_" << projData->mesh_file;
   ofstream outMesh;
   outMesh.open(ssMesh.str().c_str(),ofstream::out);
   if (outMesh.is_open()) {
      pmesh->PrintAsSerial(outMesh);
      outMesh.close();
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3177: Failed to open file \"%s\" for writing.\n",ssMesh.str().c_str());
      return true;
   }
   
   stringstream ssRegions;
   ssRegions << "materials_for_refined_" << projData->mesh_file;
   bool retval=meshMaterials->saveRegionsFile (ssRegions.str().c_str());
   if (retval) return true;

   return false;
}

int main(int argc, char *argv[])
{
   double eps0=8.8541878176e-12;
   const char *projFile;
   struct projectData projData,defaultData;
   BoundaryDatabase boundaryDatabase;
   FrequencyPlan frequencyPlan;
   meshMaterialList meshMaterials;
   MaterialDatabase materialDatabase;
   vector<ParSubMesh> parSubMeshes;
   ResultDatabase resultDatabase;
   GammaDatabase gammaDatabase;

   // Initialize Petsc and MPI
   PetscInitializeNoArguments();
   PetscMPIInt size,rank;
   MPI_Comm_size(PETSC_COMM_WORLD, &size);
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   MPI_Barrier(PETSC_COMM_WORLD);

   // trap PETSc errors to enable graceful exit, primarily for out-of-memory errors
   struct applicationContext appCtx;
   PetscPushErrorHandler(errorHandler,(struct applicationContext *) &appCtx);

   chrono::system_clock::time_point job_start_time=chrono::system_clock::now();

   // parse inputs
   int printHelp=0;
   if (argc <= 1) printHelp=1;
   else {
      if (strcmp(argv[1],"-h") == 0) printHelp=1;
      else projFile=argv[1];
   }
   if (printHelp) {help(); PetscFinalize(); exit(1);}
   
   print_copyright_notice ("OpenParEM3D");
   char *baseName=get_project_name(projFile);
   char *lockfile=create_lock_file(baseName);
   isOpenParEM2Dreachable();

   appCtx.job_start_time=job_start_time;
   appCtx.lockfile=lockfile;

   PetscPrintf(PETSC_COMM_WORLD,"Setting up ...\n");

   // project
   load_project_file(projFile,&defaultData,&projData,lockfile,job_start_time);
   delete_stale_files(baseName);
   if (projData.output_show_license) {print_license(); exit_job_on_error (job_start_time,lockfile,true);}
   show_memory (projData.debug_show_memory, "   ");

   // materials
   if (materialDatabase.load_materials(projData.materials_global_path,projData.materials_global_name,
                                       projData.materials_local_path,projData.materials_local_name,
                                       projData.materials_check_limits)) exit_job_on_error (job_start_time,lockfile,true);
   if (projData.debug_show_materials) {materialDatabase.print("   ");}

   // boundaries
   if (boundaryDatabase.load(projData.port_definition_file,projData.solution_check_closed_loop)) exit_job_on_error (job_start_time,lockfile,true);
   boundaryDatabase.markPortModes();
   if (boundaryDatabase.checkSportNumbering()) exit_job_on_error (job_start_time,lockfile,true);
   boundaryDatabase.assignAttributes();
   if (projData.debug_show_port_definitions) boundaryDatabase.print();
   if (boundaryDatabase.check_overlaps()) exit_job_on_error (job_start_time,lockfile,true);
   boundaryDatabase.set2DModeNumbers();
   boundaryDatabase.fillUnused2DModeNumbers();
   boundaryDatabase.fillIntegrationPoints();

   // mesh
   PetscPrintf(PETSC_COMM_WORLD,"Loading mesh and assigning materials ...\n");
   if (!projData.materials_check_limits) PetscPrintf(PETSC_COMM_WORLD,"   Skipping limit checks on material values\n");
   if (meshMaterials.load(projData.mesh_file,3)) exit_job_on_error (job_start_time,lockfile,true);
   //meshMaterials.print();
   Mesh mesh(projData.mesh_file, 1, 1);
   //   mesh.ScaleElements(0.001);    // hard coded for now to convert from mm to m.  - ToDo - Generalize, but ScaleElements may be broken.
   int dim=mesh.Dimension();
   if (! (dim == 3)) {PetscPrintf(PETSC_COMM_WORLD,"ERROR3011: Mesh must be 3-dimensional.\n"); exit_job_on_error (job_start_time,lockfile,true);}

   // scale checks for basic error detection
   ParMesh *pmesh=nullptr;
   if (check_field_points (projFile,&mesh,pmesh,projData.mesh_order,3,
                           projData.field_points_count,projData.field_points_x,projData.field_points_y,projData.field_points_z)) {
      exit_job_on_error (job_start_time,lockfile,true);
   }
   if (boundaryDatabase.check_scale(&mesh,projData.mesh_order)) exit_job_on_error (job_start_time,lockfile,true);

   // boundary and mesh linkage
   boundaryDatabase.markMeshBoundaries (&mesh);
   if (boundaryDatabase.setDefaultMeshBoundary(&projData,&mesh,&materialDatabase,&boundaryDatabase)) exit_job_on_error (job_start_time,lockfile,true);
   mesh.SetAttributes(); // recalculates the support data structures

   // mesh and mesh materials consistency check
   if (mesh.attributes.Max() != meshMaterials.size()) {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3012: Mesh file \"%s\" does not include the correct number of regions for material definitions.\n",projData.mesh_file);
      PetscPrintf(PETSC_COMM_WORLD,"       The $PhysicalNames block should have %d entries, but only %d were found\n.",
                                           mesh.attributes.Max()+1,meshMaterials.size()+1);
      exit_job_on_error (job_start_time,lockfile,true);
   }

   // vectors to hold material properties per distinct material region
   Vector Negko2ReRelativePermittivity(mesh.attributes.Max());
   Vector Negko2ImRelativePermittivity(mesh.attributes.Max());
   Vector InvMur(mesh.attributes.Max());
   Vector InvOmegaMu(mesh.attributes.Max());
   Vector mur(mesh.attributes.Max());
   Vector er(mesh.attributes.Max());

   show_memory (projData.debug_show_memory,"");

   // set up the frequency plan
   if (frequencyPlan.assemble(projData.refinement_frequency,projData.inputFrequencyPlansCount,projData.inputFrequencyPlans))
       exit_job_on_error (job_start_time,lockfile,true);
   if (projData.debug_show_frequency_plan) frequencyPlan.print();
   double lastFrequency=0;

   // set up pmesh when not using adaptive mesh refinement
   if (! frequencyPlan.is_refining()) pmesh=new ParMesh(PETSC_COMM_WORLD, mesh);

   // loop over frequencies

   double frequency=-1;
   bool refineMesh=false;
   bool restartMesh=false;
   bool saveFieldsHeader=true;
   while (frequencyPlan.get_frequency(projData.refinement_frequency,&frequency,&refineMesh,&restartMesh)) {

      PetscPrintf(PETSC_COMM_WORLD,"|*******************************************************************************************************************************************************\n");
      PetscPrintf(PETSC_COMM_WORLD,"| Frequency: %g\n",frequency);
      PetscPrintf(PETSC_COMM_WORLD,"|*******************************************************************************************************************************************************\n");

      double ko=2*M_PI*frequency*sqrt(4.0e-7*M_PI*eps0);

      if (lastFrequency == 0) lastFrequency=frequency;

      // set up pmesh when using adaptive mesh refinement
      if (refineMesh && restartMesh) {
         if (pmesh != NULL) delete pmesh;
         pmesh=new ParMesh(PETSC_COMM_WORLD,mesh);
         boundaryDatabase.resetElementNumbers();
      }

      if (refineMesh) {
         if (restartMesh) {PetscPrintf(PETSC_COMM_WORLD,"|   Refining mesh starting with the initial mesh.\n");}
         else {PetscPrintf(PETSC_COMM_WORLD,"|   Refining mesh starting with the last mesh.\n");}
      }

      if (! pmesh) {
         PetscPrintf(PETSC_COMM_WORLD,"ASSERT: pmesh is not defined.\n");
         exit_job_on_error (job_start_time,lockfile,true);
      }

      // create frequency-dependent material constants for integration across the finite-element space

      int j=0;
      while (j < pmesh->attributes.Max()) {
         Material *useMaterial=materialDatabase.get(meshMaterials.get_name(meshMaterials.get_index(j)));

         if (useMaterial != NULL) {

            // permittivity
            complex<double> e=useMaterial->get_eps(projData.solution_temperature,frequency,materialDatabase.get_tol(),materialDatabase.get_indent());
            if (e == complex<double>(-DBL_MAX,0)) exit_job_on_error (job_start_time,lockfile,true);

            Negko2ReRelativePermittivity[j]=-ko*ko*real(e)/eps0;
            Negko2ImRelativePermittivity[j]=-ko*ko*imag(e)/eps0;

            er[j]=real(e)/eps0;

            // permeability
            double mu=useMaterial->get_mu(projData.solution_temperature,frequency,materialDatabase.get_tol(),materialDatabase.get_indent());
            if (mu == -DBL_MAX) exit_job_on_error (job_start_time,lockfile,true);

            InvMur[j]=1.0/(mu/(4e-7*M_PI));
            InvOmegaMu[j]=1.0/(2*M_PI*frequency*mu);    // 1/(w*permeability)

            mur[j]=mu/(4e-7*M_PI);

         } else {
            PetscPrintf(PETSC_COMM_WORLD,"ERROR3013: Material \"%s\" for region %d is not present in the material database.\n",
                                         meshMaterials.get_name(j).c_str(),j+1);
            exit_job_on_error (job_start_time,lockfile,true);
         }

         j++;
      }

      PWConstCoefficient neg_ko2_Re_er(Negko2ReRelativePermittivity);
      PWConstCoefficient neg_ko2_Im_er(Negko2ImRelativePermittivity);
      PWConstCoefficient Inv_mur(InvMur);
      PWConstCoefficient Inv_w_mu(InvOmegaMu);

      // solve 

      int iteration=resultDatabase.get_lastIteration(frequency);
      bool iterate=true;
      while (iterate) {

         iterate=refineMesh;
 
         if (iterate) {PetscPrintf(PETSC_COMM_WORLD,"|   Iteration %d ...\n",iteration+1);}
         else          PetscPrintf(PETSC_COMM_WORLD,"|   Using existing mesh ...\n");

         chrono::system_clock::time_point fem_setup_start_time=chrono::system_clock::now();

         if (print_mesh_quality_message (pmesh,&projData)) exit_job_on_error (job_start_time,lockfile,true);
         PetscPrintf(PETSC_COMM_WORLD,"|      mesh size: %d\n",getGlobalNE(pmesh));

         double shortestPerWavelength=-1;
         double longestPerWavelength=-1;
         get_edgeLengthsPerWavelength(pmesh,frequency,er,mur,&shortestPerWavelength,&longestPerWavelength);
         PetscPrintf(PETSC_COMM_WORLD,"|      mesh edge length/wavelength in element: shortest=%g longest=%g\n",shortestPerWavelength,longestPerWavelength);

         // solve the 2D ports
         PetscPrintf(PETSC_COMM_WORLD,"|      solving 2D ports ...\n");
         if (boundaryDatabase.solve2Dports(pmesh,&parSubMeshes,&projData,frequency,&meshMaterials,&gammaDatabase)) {
            exit_job_on_error (job_start_time,lockfile,true);
         }

         // save the propgation constant to use as an initial guess in OpenParEM2D
         if (boundaryDatabase.getGamma(&gammaDatabase,frequency)) exit_job_on_error (job_start_time,lockfile,true);

         PetscPrintf(PETSC_COMM_WORLD,"|      building finite element spaces ...\n");
         fem3D *fem=new fem3D();
         fem->set_data(&pmesh,&projData,frequency);
         fem->build_fe_spaces();
         fem->build_portEssTdofLists(&boundaryDatabase);
         boundaryDatabase.buildGrids(fem);
         boundaryDatabase.save2DModalParaView(&parSubMeshes,&projData,frequency,false);
         boundaryDatabase.save3DParaView(pmesh,&projData,frequency,false);

         chrono::system_clock::time_point fem_setup_end_time=chrono::system_clock::now();

         // loop through the ports to drive each one in sequence

         int drivingSport=1;
         int SportCount=boundaryDatabase.get_SportCount();
         resultDatabase.set_SportCount(SportCount);
         while (drivingSport <= SportCount) {

            chrono::system_clock::time_point solve_start_time=chrono::system_clock::now();

            PetscPrintf(PETSC_COMM_WORLD,"|      Driving S-port %d ...\n",drivingSport);

            PetscPrintf(PETSC_COMM_WORLD,"|         solving E field in 3D volume ...\n");
            if (fem->solve(&boundaryDatabase,&materialDatabase,projData.solution_temperature,&neg_ko2_Re_er,&neg_ko2_Im_er,&Inv_mur,&Inv_w_mu,drivingSport)) {
               exit_job_on_error (job_start_time,lockfile,true);
            }

            Result *result=new Result();
            result->set(drivingSport,frequency,shortestPerWavelength,longestPerWavelength,iteration,projData.reference_impedance,fem,pmesh);

            boundaryDatabase.transfer_3Dsolution_3Dgrids_to_2Dgrids(fem);
            boundaryDatabase.save2DSolutionParaView(&parSubMeshes,&projData,frequency,drivingSport,false);

            PetscPrintf(PETSC_COMM_WORLD,"|         calculating S-parameters ...\n");

            Mode *drivingMode=boundaryDatabase.getDrivingMode(drivingSport);
            boundaryDatabase.set_V3Dsize ();
            boundaryDatabase.calculateVoltages(pmesh,fem->get_gridReE(),fem->get_gridImE(),drivingMode);
            boundaryDatabase.calculateWeights(drivingMode);
            boundaryDatabase.calculateSparameters(result,drivingMode);

            chrono::system_clock::time_point solve_end_time=chrono::system_clock::now();
            result->set_solve_time(elapsed_time(solve_start_time,solve_end_time));
            result->set_fem_setup_time(elapsed_time(fem_setup_start_time,fem_setup_end_time));

            if (iterate) {
               chrono::system_clock::time_point mesh_error_start_time=chrono::system_clock::now();
               PetscPrintf(PETSC_COMM_WORLD,"|         calculating mesh errors ...\n");
               fem->calculateMeshErrors(&projData,&boundaryDatabase,&Inv_mur);
               chrono::system_clock::time_point mesh_error_end_time=chrono::system_clock::now();
               result->set_mesh_error_time(elapsed_time(mesh_error_start_time,mesh_error_end_time));
            }

            resultDatabase.push(result);

            if (saveFieldsHeader) {fem->saveFieldValuesHeader(&projData); saveFieldsHeader=false;}
            fem->saveFieldValues(&projData,pmesh,iteration,drivingSport);

            drivingSport++;
         }

         // reflections as a quality measure
         drivingSport=1;
         while (drivingSport <= SportCount) {
            Result *result=resultDatabase.get_Result(frequency,drivingSport,iteration+1);
            double maxReflection=boundaryDatabase.get_maxReflection(drivingSport);
            result->set_maxReflection(maxReflection);
            drivingSport++;
         }
         boundaryDatabase.printPortReflections();

         if (iterate) {
            chrono::system_clock::time_point refine_start_time=chrono::system_clock::now();

            // mark as having iterations
            resultDatabase.set_hasIterations(true,frequency,iteration+1);

            // save the relative error
            double maxRelativeError=resultDatabase.calculate_maxRelativeError(frequency,iteration+1);
            resultDatabase.set_maxRelativeError(maxRelativeError,frequency,iteration+1);

            // check if converged - preliminary
            bool isConverged=false;
            if (!(maxRelativeError < 0)  && maxRelativeError <= projData.refinement_tolerance) {isConverged=true;}

            // check for minimum iterations
            if (iteration < projData.refinement_iteration_min-1) {isConverged=false;}

            // save the convergence status
            resultDatabase.set_isConverged(isConverged,frequency,iteration+1);

            // check for convergence
            if (resultDatabase.isSequentialConverged(&projData,frequency)) iterate=false;

            // check for maximum iterations
            if (iteration == projData.refinement_iteration_max-1) iterate=false;

            // refine
            if (iterate) {
               PetscPrintf(PETSC_COMM_WORLD,"|         refining mesh ...\n");
               fem->refineMesh();

               if (projData.mesh_save_refined) {
                  PetscPrintf(PETSC_COMM_WORLD,"|         saving refined mesh ...\n");
                  if (saveSerialMesh(&projData,&meshMaterials,pmesh)) exit_job_on_error (job_start_time,lockfile,true);
               }
            }

            chrono::system_clock::time_point refine_end_time=chrono::system_clock::now();
            resultDatabase.set_refine_time(elapsed_time(refine_start_time,refine_end_time),frequency,iteration+1);
         }

         if (projData.reference_impedance > 0) {
            if (resultDatabase.renormalize(iteration+1,frequency)) {
               cout << "ASSERT: renomalization failed." << endl;
            }
         }
         resultDatabase.save(&projData);
         resultDatabase.saveFormatted(&projData);
         resultDatabase.saveCSV (&projData,boundaryDatabase.get_SportCount(),true);
         resultDatabase.saveTouchstone (&projData,boundaryDatabase.get_SportCount());

         delete fem; fem=nullptr;
 
         boundaryDatabase.reset();

         ++iteration;
      }

      lastFrequency=frequency;
   }

   // save the results as test cases
   if (projData.test_create_cases && rank == 0) {
      resultDatabase.save_as_test(&projData);
   }

   // print convergence status
   if (resultDatabase.hasIterations()) {
      if (resultDatabase.isAllConverged()) PetscPrintf(PETSC_COMM_WORLD,"|Converged\n");
      else PetscPrintf(PETSC_COMM_WORLD,"|NOT CONVERGED\n");
   }

   PetscPrintf(PETSC_COMM_WORLD,"|Job Complete\n");

   if (rank == 0) {
      stringstream ss;
      ss << "temp_" << projData.project_name;
      if (! projData.debug_tempfiles_keep && std::filesystem::exists(ss.str().c_str())) {
        std::filesystem::remove_all(ss.str().c_str());
      }
   }

   write_attributes (baseName,pmesh);
   delete pmesh;

   resultDatabase.save(&projData);
   resultDatabase.saveCSV (&projData,boundaryDatabase.get_SportCount(),true);

   free_project (&projData);

   show_memory (projData.debug_show_memory, "");

   remove_lock_file(lockfile);

   chrono::system_clock::time_point job_end_time=chrono::system_clock::now();
   resultDatabase.set_solve_time(elapsed_time(job_start_time,job_end_time));
   PetscPrintf(PETSC_COMM_WORLD,"|Elapsed time: %g s\n",resultDatabase.get_solve_time());

   PetscFinalize();

   exit(0);
}

