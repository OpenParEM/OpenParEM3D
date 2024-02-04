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

#ifndef PORT_H
#define PORT_H

#include "mfem.hpp"
#include "petscsys.h"
#include <petsc.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <limits>
#include <filesystem>
#include <unistd.h>
#include "project.h"
#include "path.hpp"
#include "OpenParEMmaterials.hpp"
#include "mesh.hpp"
#include "sourcefile.hpp"

using namespace std;
using namespace mfem;

class RotatedMesh : public Mesh
{
   public:
      ~RotatedMesh();
      void set_spaceDim (int);
      bool rotate (Path *, bool);
};

struct mpi_complex_int {
   double real;
   double imag;
   int location;
};

class BoundaryDatabase;
class Result;
class fem3D;

double elapsed_time (chrono::system_clock::time_point, chrono::system_clock::time_point);
bool isClose (double, double);

class Gamma
{
   private:
      double alpha;
      double beta;
      double frequency;
   public:
      void set (double alpha_, double beta_, double frequency_) {alpha=alpha_; beta=beta_; frequency=frequency_;}
      double get_alpha () {return alpha;}
      double get_beta () {return beta;}
      double get_frequency () {return frequency;}
};

class GammaDatabase
{
   private:
      vector<Gamma *> gammaList;
   public:
      ~GammaDatabase();
      void push (Gamma *gamma) {gammaList.push_back(gamma);}
      Gamma* getGamma (long unsigned int);
      void reset ();
};

// Boundary and Mode share a structure in OpenParEM2D.
// Here they are split to accommodate ports more easily.

class Boundary
{
   private:
      int startLine;
      int endLine;
      keywordPair name;
      keywordPair type;                    // surface_impedance | perfect_electric_conductor (PEC) | perfect_magnetic_conductor (PMC) | radiation
      keywordPair material;                // surface impedance boundary only
      keywordPair wave_impedance;          // radiation only
      vector<keywordPair *> pathNameList;
      vector<long unsigned int> pathIndexList;
      vector<bool> reverseList;
      int attribute=-1;                    // attribute assigned to the mesh indicating this boundary
      Path *rotated=nullptr;               // rotated path boundary for ports and some geometrical operations
   public:
      Boundary (int,int);
      ~Boundary();
      bool load(string *, inputFile *);
      bool inBlock (int);
      bool check(string *, vector<Path *>);
      bool assignPathIndices(vector<Path *> *);
      bool checkBoundingBox(Vector *, Vector *, string *, double, vector<Path *> *);
      int get_startLine() {return startLine;}
      int get_endLine() {return endLine;}
      string get_name() {return name.get_value();}
      int get_attribute() {return attribute;}
      int get_pathIndex(int i) {return pathIndexList[i];}
      bool name_is_loaded() {return name.is_loaded();}
      int get_name_lineNumber() {return name.get_lineNumber();}
      void set_name (string name_) {name.set_value(name_);}
      void set_type (string type_) {type.set_value(type_);}
      void set_material (string material_) {material.set_value(material_);}
      void set_attribute(int attribute_) {attribute=attribute_;};
      string get_type() {return type.get_value();}
      string get_material() {return material.get_value();}
      double get_wave_impedance() {return wave_impedance.get_dbl_value();}
      Path* get_rotated() {return rotated;}
      string get_pathName(long unsigned int i) {return pathNameList[i]->get_value();}
      int get_pathName_lineNumber(long unsigned int i) {return pathNameList[i]->get_lineNumber();}
      bool get_reverse(long unsigned int i) {return reverseList[i];}
      long unsigned int get_path_size() {return pathIndexList.size();}
      long unsigned int get_path(long unsigned int i) {return pathIndexList[i];}
      void push(long unsigned int a) {pathIndexList.push_back(a);}
      bool is_surface_impedance();
      bool is_perfect_electric_conductor();
      bool is_perfect_magnetic_conductor();
      bool is_radiation();
      bool is_modal();
      bool is_line();
      bool has_attribute(int attribute_) {if (attribute == attribute_) return true; return false;}
      bool merge(vector<Path *> *);
      bool createRotated(vector<Path *> *, string);
      bool is_point_inside (double, double, double);
      bool is_triangleInside (DenseMatrix *);
      bool is_overlapPath (vector<Path *> *, Path *);
      Boundary* get_matchBoundary (double, double, double, double, double, double);
      void addMassImpedanceIntegrator (double, double, ParMesh *, ParMixedBilinearForm *,
                                           MaterialDatabase *, vector<Array<int> *> &,
                                           vector<ConstantCoefficient *> &, bool);
      void print();
};

class OPEMIntegrationPoint
{
   private:
      int pointNumber;
      int rank;
      int initialized;
      DenseMatrix point;
      int elementNumber;
      IntegrationPoint integrationPoint;
      complex<double> fieldX,fieldY,fieldZ;
      Vector pt;  // working space
   public:
      OPEMIntegrationPoint(int, double, double, double);
      void update (ParMesh *);
      void get_location (double *x, double *y, double *z);
      void set (double, double, double, double, double, double);
      void get_fields (complex<double> *, complex<double> *, complex<double> *);
      void get_fieldValue (ParGridFunction *, ParGridFunction *);  // from grids to local value
      void get_field (complex<double> *fieldX_, complex<double> *fieldY_, complex<double> *fieldZ_) {
         *fieldX_=fieldX; *fieldY_=fieldY, *fieldZ_=fieldZ;
      }
      void resetElementNumber ();
      void send (int);
      void print ();
};

class OPEMIntegrationPointList
{
   private:
      vector<OPEMIntegrationPoint *> points;
      complex<double> integratedValue;
   public:
      ~OPEMIntegrationPointList ();
      void push (OPEMIntegrationPoint *a) {points.push_back(a);}
      long unsigned int get_size() {return points.size();}
      OPEMIntegrationPoint* get_point (long unsigned int i) {return points[i];}
      void update (ParMesh *);
      void get_fieldValues (ParGridFunction *, ParGridFunction *);
      void resetElementNumbers ();
      void assemble ();
      void integrate ();
      complex<double> get_integratedValue () {return integratedValue;}
      void print ();
};

class Weights
{
   private:
      int drivenSport;
      complex<double> e0,e1,e2,e3;
      complex<double> h0,h1,h2,h3;
      complex<double> voltage2D;
      complex<double> Zo2D;
      complex<double> voltage3D;
   public:
      Weights (int, complex<double>, complex<double>, complex<double>, complex<double>,
                    complex<double>, complex<double>, complex<double>, complex<double>,
                    complex<double>, complex<double>, complex<double>);
      int get_drivenSport () {return drivenSport;}
      bool isPortMatch (int);
      complex<double> calculateSii ();
      complex<double> calculateSij (Weights *);
      void print (string);
};

// Mode or Line - equivalent from the spec
class Mode
{
   private:
      int startLine;
      int endLine;
      keywordPair Sport;                 // integer value for the S-parameter port number
      keywordPair type;                  // voltage | current
      vector<keywordPair *> pathNameList;
      vector<long unsigned int> pathIndexList;
      vector<bool> reverseList;
      string calculation;                // modal | line - for output formatting
      int modeNumber2D;                  // mode number used for the 2D solution

      vector<OPEMIntegrationPointList *> pointsList;  // integration points for the voltage line

      bool isUsed;
      bool isSolutionLoaded;
      double *eVecReE=nullptr;
      double *eVecImE=nullptr;
      double *eVecReH=nullptr;
      double *eVecImH=nullptr;
      Vector *eigenVecReEt=nullptr;
      Vector *eigenVecImEt=nullptr;
      Vector *eigenVecReEz=nullptr;
      Vector *eigenVecImEz=nullptr;
      Vector *eigenVecReHt=nullptr;
      Vector *eigenVecImHt=nullptr;
      Vector *eigenVecReHz=nullptr;
      Vector *eigenVecImHz=nullptr;
      double alpha,beta;

      // grid functions to hold the 2D modal fields from the 2D solution
      ParGridFunction *grid2DReEt=nullptr;
      ParGridFunction *grid2DImEt=nullptr;
      ParGridFunction *grid2DReEz=nullptr;
      ParGridFunction *grid2DImEz=nullptr;
      ParGridFunction *grid2DReHt=nullptr;
      ParGridFunction *grid2DImHt=nullptr;
      ParGridFunction *grid2DReHz=nullptr;
      ParGridFunction *grid2DImHz=nullptr;

      // grid functions to hold the 2D modal fields projected onto the 3D space
      ParGridFunction *grid3DReEt=nullptr;
      ParGridFunction *grid3DImEt=nullptr;
      ParGridFunction *grid3DReEz=nullptr;
      ParGridFunction *grid3DImEz=nullptr;
      ParGridFunction *grid3DReHt=nullptr;
      ParGridFunction *grid3DImHt=nullptr;
      ParGridFunction *grid3DReHz=nullptr;
      ParGridFunction *grid3DImHz=nullptr;

      // grid functions to hold the 2D modal solutions projected back onto 2D spaces
      // (to align with the grid functions grid2Dsolution* on the ports)
      ParGridFunction *grid2DmodalReEt=nullptr;
      ParGridFunction *grid2DmodalImEt=nullptr;
      ParGridFunction *grid2DmodalReEz=nullptr;
      ParGridFunction *grid2DmodalImEz=nullptr;
      ParGridFunction *grid2DmodalReHt=nullptr;
      ParGridFunction *grid2DmodalImHt=nullptr;
      ParGridFunction *grid2DmodalReHz=nullptr;
      ParGridFunction *grid2DmodalImHz=nullptr;

      // for S-parameter calculation
      vector<Weights *> weightsList;       // one set of weights per driven port (Sport)

   public:
      Mode(int,int,string);
      ~Mode();
      int get_startLine() {return startLine;}
      int get_endLine() {return endLine;}
      int get_Sport() {return Sport.get_int_value();}
      int get_Sport_lineNumber() {return Sport.get_lineNumber();}
      bool get_isUsed() {return isUsed;}
      double get_alpha() {return alpha;}
      double get_beta() {return beta;}
      int get_modeNumber2D() {return modeNumber2D;}
      void set_isUsed() {isUsed=true;}
      void set_modeNumber2D (int modeNumber2D_) {modeNumber2D=modeNumber2D_;}
      void set_alpha(double alpha_) {alpha=alpha_;}
      void set_beta(double beta_) {beta=beta_;}
      bool is_modal() {if (calculation.compare("modal") == 0) return true; return false;}
      bool is_line() {if (calculation.compare("line") == 0) return true; return false;}
      string get_type() {return type.get_value();}
      bool load(string *, inputFile *);
      bool inModeBlock(int);
      bool is_voltage() {if (type.get_value().compare("voltage") == 0) return true; return false;}
      bool is_current() {if (type.get_value().compare("current") == 0) return true; return false;}
      bool check(string *, vector<Path *>);
      bool check_current_paths (string *, vector<Path *> *, bool);
      bool checkBoundingBox(Vector *, Vector *, string *, double, vector<Path *> *);
      bool assignPathIndices(vector<Path *> *);
      bool is_enclosedByPath (vector<Path *> *, Path *);
      void print(string);
      void output (ofstream *, vector<Path *> *, Path *, bool);
      bool loadSolution(string *, string, size_t, size_t);
      bool scaleSolution ();
      void printSolution ();
      void setSign (double *, double *);
      void build2Dgrids (ParFiniteElementSpace *, ParFiniteElementSpace *);
      void build3Dgrids (ParFiniteElementSpace *, ParFiniteElementSpace *);
      void build2DModalGrids (ParFiniteElementSpace *, ParFiniteElementSpace *);
      void flipModalHsign (bool);
      void fillX (Vec *, Vec *, Array<int> *, HYPRE_BigInt *);
      void fillIntegrationPoints (vector<Path *> *);
      complex<double> calculateLineIntegral (ParMesh *, vector<Path *> *, ParGridFunction *, ParGridFunction *);
      void calculateWeights (ParFiniteElementSpace *, ParGridFunction *, ParGridFunction *, ParGridFunction *,
                             ParGridFunction *, ParGridFunction *, ParGridFunction *, ParGridFunction *, ParGridFunction *, 
                             Vector, int, complex<double>, complex<double>, complex<double>);
      void clearWeights ();
      double get_maxReflection ();
      double get_maxReflection (int);
      void printPortReflections ();
      complex<double> calculateSparameter (Mode *);
      void transfer_2Dsolution_2Dgrids_to_3Dgrids ();
      void transfer_2Dsolution_3Dgrids_to_2Dgrids ();
      void save2DParaView (ParSubMesh *, struct projectData *, double, bool);
      void save3DParaView (ParMesh *, struct projectData *, double, bool);
      void save2DModalParaView (ParSubMesh *, struct projectData *, double, bool);
      void resetElementNumbers ();
      void reset ();
};

class Port
{
   private:
      int startLine;
      int endLine;
      keywordPair name;                            // alphanumeric name 
      vector<keywordPair *> pathNameList;
      vector<long unsigned int> pathIndexList;
      vector<bool> reverseList;
      keywordPair impedance_definition;            // VI, PV, or PI
      keywordPair impedance_calculation;           // modal | line
      vector<Mode *> modeList;
      int attribute=-1;

      ND_FECollection *fec2D_ND=nullptr;           // Et
      ParFiniteElementSpace *fes2D_ND=nullptr;     // on 2D ParMesh

      H1_FECollection *fec2D_H1=nullptr;           // Ez
      ParFiniteElementSpace *fes2D_H1=nullptr;     // on 2D ParMesh

      L2_FECollection *fec2D_L2=nullptr;           // for S-parameter calculations using modal projections
      ParFiniteElementSpace *fes2D_L2=nullptr;     // on 2D ParMesh

      bool spin180degrees;
      size_t t_size, z_size;
      string meshFilename="";
      string modesFilename="";
      Path *rotated=nullptr;                       // rotated path boundary for the port
      Vector normal;                               // normal facing outward from the 3D volume
      Vector rotated_normal;                       // outward facing normal after rotation for the rotated 2D mesh
      DenseMatrix ReZ2D,ImZ2D;                     // from 2D modal simulations (OpenParEM2D)
      DenseMatrix ReV2D,ImV2D;                     // from 2D modal simulations (OpenParEM2D)
      DenseMatrix ReV3D,ImV3D;                     // from 3D simulation (OpenParEM3D)

      Array<int> *ess_tdof_list=nullptr;           // on 3D mesh and space
      HYPRE_BigInt *offset;

      // grid functions to hold the 3D solutions on the 2D ports
      ParGridFunction *grid2DsolutionReEt=nullptr;
      ParGridFunction *grid2DsolutionImEt=nullptr;
      ParGridFunction *grid2DsolutionReEz=nullptr;
      ParGridFunction *grid2DsolutionImEz=nullptr;
      ParGridFunction *grid2DsolutionReHt=nullptr;
      ParGridFunction *grid2DsolutionImHt=nullptr;
      ParGridFunction *grid2DsolutionReHz=nullptr;
      ParGridFunction *grid2DsolutionImHz=nullptr;

   public:
      Port(int,int);
      ~Port();
      int get_startLine() {return startLine;}
      int get_endLine() {return endLine;}
      string get_name() {return name.get_value();}
      int get_name_lineNumber() {return name.get_lineNumber();}
      int get_modeCount();
      int get_SportCount();
      int get_attribute() {return attribute;}
      Path* get_rotated() {return rotated;}
      int get_pathIndex(int i) {return pathIndexList[i];}
      void set_attribute(int attribute_) {attribute=attribute_;};
      bool is_modal() {if (impedance_calculation.get_value().compare("modal") == 0) return true; return false;}
      bool has_attribute(int attribute_) {if (attribute == attribute_) return true; return false;}
      bool load(string *, inputFile *);
      bool inBlock(int);
      bool inModeBlocks(int);
      bool findModeBlocks(inputFile *);
      bool findLineBlocks(inputFile *);
      bool check(string *, vector<Path *>, bool);
      bool assignPathIndices (vector<Path *> *);
      bool checkBoundingBox(Vector *, Vector *, string *, double, vector<Path *> *);
      vector<int> get_SportList ();
      void push(long unsigned int a) {pathIndexList.push_back(a);}
      bool merge(vector<Path *> *);
      bool is_point_inside (double, double, double);
      bool is_triangleInside (DenseMatrix *);
      bool is_modePathInside (string *, vector<Path *> *);
      bool createRotated(vector<Path *> *, string);
      void create2Dmesh (int, ParMesh *, vector<ParSubMesh> *, long unsigned int, double);
      void saveMesh (meshMaterialList *, string *, ParSubMesh *);
      bool postProcessMesh (string, string);
      void save2Dsetup(struct projectData *, string *, double, Gamma *);
      void saveModeFile (struct projectData *, vector<Path *> *, BoundaryDatabase *);
      void set_filenames();
      bool is_overlapPath (Path *);
      Mode* get_mode (long unsigned int);
      Vector& get_normal() {return normal;}
      Vector get_rotated_normal() {return rotated_normal;}
      bool createDirectory(string *);
      void set2DModeNumbers();
      void fillUnused2DModeNumbers();
      void print();
      void printSolution (string);
      void printPaths(vector<Path *> *);
      bool solve(string *, MPI_Comm *);
      bool loadSolution(string *,double);
      bool getGamma (Gamma *, double);
      bool loadSizes_tz (string *);
      bool uses_current();
      bool uses_voltage();
      void markUsedModes ();
      bool set_alphaBeta (double, double, int);
      bool load_gammaZ (string *, double);
      void build2Dgrids ();
      void build3Dgrids (ParFiniteElementSpace *, ParFiniteElementSpace *);
      void build2DSolutionGrids ();
      void build2DModalGrids ();
      void build_essTdofList (ParFiniteElementSpace *, ParMesh *);
      bool addMassPortIntegrators (ParMesh *, ParMixedBilinearForm *, PWConstCoefficient *, vector<Array<int> *> &, vector<ConstantCoefficient *> &, vector<ConstantCoefficient *> &, bool);
      void fillX (Vec *, Vec *, int);
      void flipModalHsign(long unsigned int);
      void extract2Dmesh(ParMesh *, vector<ParSubMesh> *);
      void set_V3Dsize (int);
      void calculateWeights(Mode *);
      double get_maxReflection();
      double get_maxReflection(int);
      void printPortReflections();
      Mode* getDrivingMode (int);
      void fillIntegrationPoints (vector<Path *> *);
      void calculateVoltages (ParMesh *, vector<Path *> *, ParGridFunction *, ParGridFunction *, Mode *);
      void calculateSparameter (Result *, Mode *);
      void transfer_2Dsolution_2Dgrids_to_3Dgrids();
      void transfer_2Dsolution_3Dgrids_to_2Dgrids();
      void transfer_3Dsolution_3Dgrids_to_2Dgrids(fem3D *);
      void save2DParaView(ParSubMesh *, struct projectData *, double, bool);
      void save3DParaView(ParMesh *, struct projectData *, double, bool);
      void save2DSolutionParaView(ParSubMesh *, struct projectData *, double, int, bool);
      void save2DModalParaView(ParSubMesh *, struct projectData *, double, bool);
      void resetElementNumbers ();
      void reset();
};

// boundary conditions and ports
class BoundaryDatabase
{
   private:
      inputFile inputs;
      vector<SourceFile *> sourceFileList;
      vector<Path *> pathList;
      vector<Boundary *> boundaryList;
      vector<Port *> portList;
      double tol=1e-14;
      string tempDirectory="";
      string indent="   ";
      string version_name="#OpenParEMports";
      string version_value="1.0";
   public:
      ~BoundaryDatabase();
      void set_tempDirectory(string tempDirectory_) {tempDirectory=tempDirectory_;}
      string get_tempDirectory() {return tempDirectory;}
      void assignAttributes ();
      void markMeshBoundaries (Mesh *mesh);
      bool setDefaultMeshBoundary (struct projectData *, Mesh *, MaterialDatabase *, BoundaryDatabase *);
      bool inBlocks(int);
      bool findSourceFileBlocks();
      bool findPathBlocks();
      bool findBoundaryBlocks();
      bool findPortBlocks();
      bool load(const char *, bool);
      bool check(bool);
      bool checkSportNumbering();
      bool check_scale (Mesh *, int);
      bool check_overlaps ();
      void subdivide_paths ();
      void print();
      void create2Dmeshes(int, ParMesh *, vector<ParSubMesh> *);
      vector<Port *> get_portList() {return portList;}
      vector<Path *> get_pathList() {return pathList;}
      Path* get_path (long unsigned int i) {return pathList[i];}
      Port* get_port (long unsigned int i) {return portList[i];}
      Port* get_port (int);
      int get_portAttribute (long unsigned int i) {return portList[i]->get_attribute();}
      int getLastAttribute ();
      void savePortMeshes(meshMaterialList *, vector<ParSubMesh> *);
      void save2Dsetups (struct projectData *, double, GammaDatabase *);
      void saveModeFiles (struct projectData *);
      Boundary* get_matchBoundary (double, double, double, double, double, double);
      bool createPortDirectories ();
      bool solvePorts(int, ParMesh *, vector<ParSubMesh> *, double, meshMaterialList *, struct projectData *, GammaDatabase *);
      bool loadPortSolutions (double);
      bool getGamma (GammaDatabase *, double);
      void markPortModes ();
      void printPortSolutions();
      int get_totalModeCount();
      int get_SportCount();
      int get_maxDofCount ();
      void set2DModeNumbers();
      void fillUnused2DModeNumbers();
      void build_portEssTdofLists (ParFiniteElementSpace *, ParMesh *);
      bool addMassPortIntegrators (ParMesh *, ParMixedBilinearForm *, PWConstCoefficient *,
                                   vector<Array<int> *> &, vector<ConstantCoefficient *> &,
                                   vector<ConstantCoefficient *> &, bool); 
      void addMassImpedanceIntegrators (double, double, ParMesh *, ParMixedBilinearForm *, MaterialDatabase *,
                                        vector<Array<int> *> &, vector<ConstantCoefficient *> &, bool);
      void fillX (Vec *, Vec *, int);
      void showPortDofCounts ();
      void extract2Dmesh(ParMesh *, vector<ParSubMesh> *);
      void set_V3Dsize ();
      void calculateWeights (Mode *);
      Mode* getDrivingMode (int);
      void fillIntegrationPoints ();
      void calculateVoltages (ParMesh *, ParGridFunction *, ParGridFunction *, Mode *);
      void calculateSparameters (Result *, Mode *);
      double get_maxReflection();
      double get_maxReflection(int);
      void printPortReflections();
      void build2Dgrids ();
      void build2DModalGrids ();
      void build3Dgrids (ParFiniteElementSpace *, ParFiniteElementSpace *);
      void build2DSolutionGrids ();
      void buildGrids (fem3D *);
      void transfer_2Dsolution_2Dgrids_to_3Dgrids ();
      void transfer_2Dsolution_3Dgrids_to_2Dgrids ();
      void transfer_3Dsolution_3Dgrids_to_2Dgrids (fem3D *);
      void save2DParaView (vector<ParSubMesh> *, struct projectData *, double, bool);
      void save3DParaView (ParMesh *, struct projectData *, double, bool);
      void save2DSolutionParaView (vector<ParSubMesh> *, struct projectData *, double, int, bool);
      void save2DModalParaView (vector<ParSubMesh> *, struct projectData *, double, bool);
      bool solve2Dports (ParMesh *, vector<ParSubMesh> *, struct projectData *, double, meshMaterialList *, GammaDatabase *);
      void resetElementNumbers ();
      void reset();
};

#endif
