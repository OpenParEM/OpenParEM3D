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

#include "port.hpp"
#include "fem3D.hpp"
#include "results.hpp"

string convertLogic (int a)
{
   string retval;
   if (a) retval="true";
   else retval="false";
   return retval;
}

string convertLogic (bool a)
{
   string retval;
   if (a) retval="true";
   else retval="false";
   return retval;
}

///////////////////////////////////////////////////////////////////////////////////////////
// RotatedMesh
///////////////////////////////////////////////////////////////////////////////////////////

void RotatedMesh::set_spaceDim (int spaceDim_)
{
   spaceDim=spaceDim_;
}

bool RotatedMesh::rotate (Path *rotated, bool spin180degrees)
{
   double coord[3];

   if (rotated == nullptr) {
      cout << "ASSERT: RotatedMesh::rotate passed invalid Path." << endl;
      return true;
   }

   int i=0;
   while (i < vertices.Size()) {
      coord[0]=vertices[i](0);
      coord[1]=vertices[i](1);
      coord[2]=vertices[i](2);

      rotated->rotatePoint(&coord[0],&coord[1],&coord[2],spin180degrees);

      // set z to 0
      coord[2]=0;

      vertices[i].SetCoords(3,coord);

      i++;
   }

   return false;
}

RotatedMesh::~RotatedMesh ()
{
   //cout << "~RotatedMesh: " << this << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Gamma
///////////////////////////////////////////////////////////////////////////////////////////

void GammaDatabase::reset ()
{
   long unsigned int i=0;
   while (i < gammaList.size()) {
      delete gammaList[i];
      i++;
   }

   gammaList.clear();
}

Gamma* GammaDatabase::getGamma (long unsigned int i)
{
   if (i < gammaList.size()) return gammaList[i];
   return nullptr;
}

GammaDatabase::~GammaDatabase ()
{
   long unsigned int i=0;
   while (i < gammaList.size()) {
      delete gammaList[i];
      i++;
   }
}

///////////////////////////////////////////////////////////////////////////////////////////
// OPEMIntegrationPoint
///////////////////////////////////////////////////////////////////////////////////////////

OPEMIntegrationPoint::OPEMIntegrationPoint (int pointNumber_, double x, double y, double z)
{
   pointNumber=pointNumber_;
   rank=-1;
   initialized=0;
   point.SetSize(3,1);
   point(0,0)=x;
   point(1,0)=y;
   point(2,0)=z;
   elementNumber=-1;
   fieldX=complex<double>(DBL_MAX,DBL_MAX);
   fieldY=complex<double>(DBL_MAX,DBL_MAX);
   fieldZ=complex<double>(DBL_MAX,DBL_MAX);
   pt.SetSize(3);
}

void OPEMIntegrationPoint::get_location (double *x, double *y, double *z)
{
   *x=point(0,0);
   *y=point(1,0);
   *z=point(2,0);
}

void OPEMIntegrationPoint::update (ParMesh *pmesh)
{
   PetscMPIInt rank_;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);

   bool search=false;
   ElementTransformation *eltransf;
   InverseElementTransformation *inv_tr=new InverseElementTransformation;

   inv_tr->SetPrintLevel(-1);         // -1 - never print
   inv_tr->SetInitialGuessType(InverseElementTransformation::Center);
   inv_tr->SetSolverType(InverseElementTransformation::Newton);

   pt(0)=point(0,0);
   pt(1)=point(1,0);
   pt(2)=point(2,0);

   if (elementNumber < 0) {
      if (!initialized) search=true;
   } else {
      if (rank == rank_) {
         eltransf=pmesh->GetElementTransformation(elementNumber);
         inv_tr->SetTransformation(*eltransf);
         int res=inv_tr->Transform(pt,integrationPoint);
         if (res != InverseElementTransformation::Inside) search=true;
      }
   }

   if (search) {
      int i=0;
      while (i < pmesh->GetNE()) {
         eltransf=pmesh->GetElementTransformation(i);
         inv_tr->SetTransformation(*eltransf);
         int res=inv_tr->Transform(pt,integrationPoint);
         if (res == InverseElementTransformation::Inside) {
            elementNumber=i;
            rank=rank_;
            break;
         }
         i++;
      }
   }

   initialized=1;

   delete inv_tr;
}

void OPEMIntegrationPoint::set (double ReFieldX, double ImFieldX, double ReFieldY, double ImFieldY, double ReFieldZ, double ImFieldZ)
{
   fieldX=complex<double>(ReFieldX,ImFieldX);
   fieldY=complex<double>(ReFieldY,ImFieldY);
   fieldZ=complex<double>(ReFieldZ,ImFieldZ);
}

void OPEMIntegrationPoint::get_fields (complex<double> *fieldX_, complex<double> *fieldY_, complex<double> *fieldZ_)
{
   *fieldX_=fieldX;
   *fieldY_=fieldY;
   *fieldZ_=fieldZ;
}

void OPEMIntegrationPoint::get_fieldValue (ParGridFunction *grid_re, ParGridFunction *grid_im)
{
   if (elementNumber >= 0) {
      Vector ReValue,ImValue;
      grid_re->GetVectorValue(elementNumber,integrationPoint,ReValue);
      grid_im->GetVectorValue(elementNumber,integrationPoint,ImValue);
      fieldX=complex<double>(ReValue.Elem(0),ImValue.Elem(0));
      fieldY=complex<double>(ReValue.Elem(1),ImValue.Elem(1));
      fieldZ=complex<double>(ReValue.Elem(2),ImValue.Elem(2));
   }
}

void OPEMIntegrationPoint::resetElementNumber ()
{
   elementNumber=-1;
}

void OPEMIntegrationPoint::send(int destination)
{
   MPI_Send(&pointNumber,1,MPI_INT,destination,100,PETSC_COMM_WORLD);
   MPI_Send(&rank,1,MPI_INT,destination,101,PETSC_COMM_WORLD);
   MPI_Send(&elementNumber,1,MPI_INT,destination,103,PETSC_COMM_WORLD);

   double valueReX=real(fieldX);
   MPI_Send(&valueReX,1,MPI_DOUBLE,destination,104,PETSC_COMM_WORLD);

   double valueImX=imag(fieldX);
   MPI_Send(&valueImX,1,MPI_DOUBLE,destination,105,PETSC_COMM_WORLD);

   double valueReY=real(fieldY);
   MPI_Send(&valueReY,1,MPI_DOUBLE,destination,106,PETSC_COMM_WORLD);

   double valueImY=imag(fieldY);
   MPI_Send(&valueImY,1,MPI_DOUBLE,destination,107,PETSC_COMM_WORLD);

   double valueReZ=real(fieldZ);
   MPI_Send(&valueReZ,1,MPI_DOUBLE,destination,108,PETSC_COMM_WORLD);

   double valueImZ=imag(fieldZ);
   MPI_Send(&valueImZ,1,MPI_DOUBLE,destination,109,PETSC_COMM_WORLD);
}

void OPEMIntegrationPoint::print()
{
   cout << "OPEMIntegrationPoint: this=" << this << endl;
   cout << "   pointNumber=" << pointNumber << endl;
   cout << "   rank=" << rank << endl;
   cout << "   initialized=" << initialized << endl;
   cout << "   point=(" << point(0,0) << "," << point(1,0) << "," << point(2,0) << ")" << endl;
   cout << "   elementNumber=" << elementNumber << endl;
   cout << "   integrationPoint=(" << integrationPoint.x << "," << integrationPoint.y << "," << integrationPoint.z << ")" << endl;
   cout << "   fieldX=" << fieldX << endl;
   cout << "   fieldY=" << fieldY << endl;
   cout << "   fieldZ=" << fieldZ << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// OPEMIntegrationPointList
///////////////////////////////////////////////////////////////////////////////////////////

void OPEMIntegrationPointList::update (ParMesh *pmesh)
{
   long unsigned int i=0;
   while (i < points.size()) {
      points[i]->update(pmesh);
      i++;
   }
}

void OPEMIntegrationPointList::get_fieldValues (ParGridFunction *grid_re, ParGridFunction *grid_im)
{
   long unsigned int i=0;
   while (i < points.size()) {
      points[i]->get_fieldValue(grid_re,grid_im);
      i++;
   }
}

void OPEMIntegrationPointList::assemble()
{
   PetscMPIInt rank,size;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   MPI_Comm_size(PETSC_COMM_WORLD, &size);

   // collect to 0
   if (rank == 0) {
     int i=1;
     while (i < size) {
        long unsigned int j=0;
        while (j < points.size()) {
           int pointNumber_;
           int rank_;
           int elementNumber_;
           double ReFieldX,ImFieldX,ReFieldY,ImFieldY,ReFieldZ,ImFieldZ;

           MPI_Recv(&pointNumber_,1,MPI_INT,i,100,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
           MPI_Recv(&rank_,1,MPI_INT,i,101,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
           MPI_Recv(&elementNumber_,1,MPI_INT,i,103,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
           MPI_Recv(&ReFieldX,1,MPI_DOUBLE,i,104,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
           MPI_Recv(&ImFieldX,1,MPI_DOUBLE,i,105,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
           MPI_Recv(&ReFieldY,1,MPI_DOUBLE,i,106,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
           MPI_Recv(&ImFieldY,1,MPI_DOUBLE,i,107,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
           MPI_Recv(&ReFieldZ,1,MPI_DOUBLE,i,108,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
           MPI_Recv(&ImFieldZ,1,MPI_DOUBLE,i,109,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

           if (elementNumber_ >= 0) {
              points[pointNumber_]->set(ReFieldX,ImFieldX,ReFieldY,ImFieldY,ReFieldZ,ImFieldZ);
           }

           j++;
        }
        i++;
     }
   } else {
      long unsigned int i=0;
      while (i < points.size()) {
         points[i]->send(0);
         i++;
      }
   }

   // send the complete set back out
   if (rank == 0) {
      int i=1;
      while (i < size) {
         long unsigned int j=0;
         while (j < points.size()) {
           points[j]->send(i);
           j++;
         }
         i++;
      }
   } else {
      long unsigned int j=0;
      while (j < points.size()) {
         int pointNumber_;
         int rank_;
         int elementNumber_;
         double ReFieldX,ImFieldX,ReFieldY,ImFieldY,ReFieldZ,ImFieldZ;

         MPI_Recv(&pointNumber_,1,MPI_INT,0,100,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         MPI_Recv(&rank_,1,MPI_INT,0,101,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         MPI_Recv(&elementNumber_,1,MPI_INT,0,103,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         MPI_Recv(&ReFieldX,1,MPI_DOUBLE,0,104,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         MPI_Recv(&ImFieldX,1,MPI_DOUBLE,0,105,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         MPI_Recv(&ReFieldY,1,MPI_DOUBLE,0,106,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         MPI_Recv(&ImFieldY,1,MPI_DOUBLE,0,107,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         MPI_Recv(&ReFieldZ,1,MPI_DOUBLE,0,108,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         MPI_Recv(&ImFieldZ,1,MPI_DOUBLE,0,109,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

         points[pointNumber_]->set(ReFieldX,ImFieldX,ReFieldY,ImFieldY,ReFieldZ,ImFieldZ);

         j++;
      }
   }
}

void OPEMIntegrationPointList::integrate ()
{
   double x1,y1,z1;
   double x2,y2,z2;
   double length,nhatx,nhaty,nhatz;
   double segmentLength;
   complex<double> fieldX1,fieldY1,fieldZ1;
   complex<double> fieldX2,fieldY2,fieldZ2; 
   complex<double> valueX,valueY,valueZ;

   integratedValue=complex<double>(0,0);

   nhatx=1; nhaty=1; nhatz=1;
   if (points.size() > 2) {
      points[0]->get_location(&x1,&y1,&z1);
      points[1]->get_location(&x2,&y2,&z2);

      length=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
      nhatx=(x2-x1)/length;
      nhaty=(y2-y1)/length;
      nhatz=(z2-z1)/length;
   }

   long unsigned int i=0;
   while (i < points.size()-1) {
      points[i]->get_location(&x1,&y1,&z1);
      points[i+1]->get_location(&x2,&y2,&z2);

      segmentLength=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));

      points[i]->get_field(&fieldX1,&fieldY1,&fieldZ1);
      points[i+1]->get_field(&fieldX2,&fieldY2,&fieldZ2);

      // enable the calculation to continue in case a point lies just outside the mesh
      if (fieldX1 == complex<double>(DBL_MAX,DBL_MAX) || fieldX2 == complex<double>(DBL_MAX,DBL_MAX) ||
          fieldY1 == complex<double>(DBL_MAX,DBL_MAX) || fieldY2 == complex<double>(DBL_MAX,DBL_MAX) ||
          fieldZ1 == complex<double>(DBL_MAX,DBL_MAX) || fieldZ2 == complex<double>(DBL_MAX,DBL_MAX)) {i++; continue;}

      valueX=(fieldX1+fieldX2)/2;
      valueY=(fieldY1+fieldY2)/2;
      valueZ=(fieldZ1+fieldZ2)/2;

      integratedValue+=(valueX*nhatx+valueY*nhaty+valueZ*nhatz)*segmentLength;

      i++;
   }
}

void OPEMIntegrationPointList::resetElementNumbers ()
{
   long unsigned int i=0;
   while (i < points.size()) {
      points[i]->resetElementNumber();
      i++;
   }
}

void OPEMIntegrationPointList::print()
{
   long unsigned int i=0;
   while (i < points.size()) {
      points[i]->print();
      i++;
   }
}

OPEMIntegrationPointList::~OPEMIntegrationPointList ()
{
   long unsigned int i=0;
   while (i < points.size()) {
      if (points[i]) delete points[i];
      i++;
   }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Boundary
///////////////////////////////////////////////////////////////////////////////////////////

Boundary::Boundary(int startLine_, int endLine_)
{
   startLine=startLine_;
   endLine=endLine_;

   // name
   name.push_alias("name");
   name.set_loaded(false);
   name.set_positive_required(false);
   name.set_non_negative_required(false);
   name.set_lowerLimit(0);
   name.set_upperLimit(0);
   name.set_checkLimits(false);

   // type
   type.push_alias("type");
   type.set_loaded(false);
   type.set_positive_required(false);
   type.set_non_negative_required(false);
   type.set_lowerLimit(0);
   type.set_upperLimit(0);
   type.set_checkLimits(false);

   // material
   material.push_alias("material");
   material.set_loaded(false);
   material.set_positive_required(false);
   material.set_non_negative_required(false);
   material.set_lowerLimit(0);
   material.set_upperLimit(0);
   material.set_checkLimits(false);

   // wave_impedance
   wave_impedance.push_alias("wave_impedance");
   wave_impedance.set_loaded(false);
   wave_impedance.set_positive_required(true);
   wave_impedance.set_non_negative_required(false);
   wave_impedance.set_lowerLimit(1e-6);
   wave_impedance.set_upperLimit(1e9);
   wave_impedance.set_checkLimits(true);
}

bool Boundary::load(string *indent, inputFile *inputs)
{
   bool fail=false;
   bool found_first_path=false;

   int lineNumber=inputs->get_next_lineNumber(startLine);
   int stopLineNumber=inputs->get_previous_lineNumber(endLine);
   while (lineNumber <= stopLineNumber) {

      string token,value,line;
      line=inputs->get_line(lineNumber);
      get_token_pair(&line,&token,&value,&lineNumber,*indent);

      int recognized=0;

      if (name.match_alias(&token)) {
         if (name.is_loaded()) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3014: Duplicate entry at line %d for previous entry at line %d.\n",
                                         indent->c_str(),indent->c_str(),lineNumber,name.get_lineNumber());
            fail=true;
         } else {
            name.set_keyword(token);
            name.set_value(value);
            name.set_lineNumber(lineNumber);
            name.set_loaded(true);
         }
         recognized++;
      }

      if (type.match_alias(&token)) {
         if (type.is_loaded()) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3015: Duplicate entry at line %d for previous entry at line %d.\n",
                                         indent->c_str(),indent->c_str(),lineNumber,type.get_lineNumber());
            fail=true;
         } else {
            type.set_keyword(token);
            type.set_value(value);
            type.set_lineNumber(lineNumber);
            type.set_loaded(true);
         }
         recognized++;
      }

      if (material.match_alias(&token)) {
         if (material.is_loaded()) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3016: Duplicate entry at line %d for previous entry at line %d.\n",
                                         indent->c_str(),indent->c_str(),lineNumber,material.get_lineNumber());
            fail=true;
         } else {
            material.set_keyword(token);
            material.set_value(value);
            material.set_lineNumber(lineNumber);
            material.set_loaded(true);
         }
         recognized++;
      }

      if (wave_impedance.match_alias(&token)) {
         recognized++;
         if (wave_impedance.loadDouble(&token, &value, lineNumber)) fail=true;
      }

      if (token.compare("path") == 0) {
         if (found_first_path) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3017: Extraneous path= statement at line %d.\n",
                                          indent->c_str(),indent->c_str(),lineNumber);
            fail=true;
         } else {
            bool reverse=false;
            if (value.substr(0,1).compare("+") == 0) value=value.substr(1);
            else if (value.substr(0,1).compare("-") == 0) {value=value.substr(1); reverse=true;}

            keywordPair *path=new keywordPair();
            path->push_alias("path");
            path->set_keyword(token);
            path->set_value(value);
            path->set_lineNumber(lineNumber);
            path->set_positive_required(false);
            path->set_non_negative_required(false);
            path->set_lowerLimit(0);
            path->set_upperLimit(0);
            path->set_checkLimits(false);
            path->set_loaded(true);

            pathNameList.push_back(path);
            reverseList.push_back(reverse);
         }

         found_first_path=true;
         recognized++;
      }

      if (token.compare("path+") == 0) {
         if (found_first_path) {
            if (value.substr(0,1).compare("+")  == 0 || value.substr(0,1).compare("-") == 0) {
               PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3018: Misformatted path at line %d.\n",
                                            indent->c_str(),indent->c_str(),lineNumber);
               fail=true;
            } else {

               keywordPair *path=new keywordPair();
               path->push_alias("path");
               path->set_keyword(token);
               path->set_value(value);
               path->set_lineNumber(lineNumber);
               path->set_positive_required(false);
               path->set_non_negative_required(false);
               path->set_lowerLimit(0);
               path->set_upperLimit(0);
               path->set_checkLimits(false);
               path->set_loaded(true);

               pathNameList.push_back(path);
               reverseList.push_back(false);
            }
         } else {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3019: Missing path= statement before line %d.\n",
                                         indent->c_str(),indent->c_str(),lineNumber);
            fail=true;
         }
         recognized++;
      }

      if (token.compare("path-") == 0) {
         if (found_first_path) {
            if (value.substr(0,1).compare("+") == 0 || value.substr(0,1).compare("-") == 0) {
               PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3020: Misformatted path at line %d.\n",
                                            indent->c_str(),indent->c_str(),lineNumber);
               fail=true;
            } else {

               keywordPair *path=new keywordPair();
               path->push_alias("path");
               path->set_keyword(token);
               path->set_value(value);
               path->set_lineNumber(lineNumber);
               path->set_positive_required(false);
               path->set_non_negative_required(false);
               path->set_lowerLimit(0);
               path->set_upperLimit(0);
               path->set_checkLimits(false);
               path->set_loaded(true);

               pathNameList.push_back(path);
               reverseList.push_back(true);
            }
         } else {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3021: Missing path= statement before line %d.\n",
                                          indent->c_str(),indent->c_str(),lineNumber);
            fail=true;
         }
         recognized++;
      }

      // should recognize one keyword
      if (recognized != 1) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3022: Unrecognized keyword at line %d.\n",
                                       indent->c_str(),indent->c_str(),lineNumber);
         fail=true;
      }
      lineNumber=inputs->get_next_lineNumber(lineNumber);
   }

   // apply defaults

   if (is_radiation() && !wave_impedance.is_loaded()) {
      wave_impedance.set_dbl_value(sqrt(M_PI*4e-7/8.8541878176e-12));
      wave_impedance.set_loaded(true);
   }

   return fail;
}

bool Boundary::is_surface_impedance()
{
   if (type.get_value().compare("surface_impedance") == 0) return true;
   return false;
}

bool Boundary::is_perfect_electric_conductor()
{
   if (type.get_value().compare("perfect_electric_conductor") == 0) return true;
   return false;
}

bool Boundary::is_perfect_magnetic_conductor()
{
   if (type.get_value().compare("perfect_magnetic_conductor") == 0) return true;
   return false;
}

bool Boundary::is_radiation()
{
   if (type.get_value().compare("radiation") == 0) return true;
   return false;
}

void Boundary::print()
{
   PetscPrintf(PETSC_COMM_WORLD,"Boundary\n");
   PetscPrintf(PETSC_COMM_WORLD,"   name=%s\n",get_name().c_str());
   PetscPrintf(PETSC_COMM_WORLD,"   type=%s\n",get_type().c_str());
   if (is_surface_impedance() && type.is_loaded()) PetscPrintf(PETSC_COMM_WORLD,"   material=%s\n",get_material().c_str());
   if (is_radiation()) PetscPrintf(PETSC_COMM_WORLD,"   wave_impedance=%g\n",get_wave_impedance());
   long unsigned int i=0;
   while (i < pathNameList.size()) {
      if (i == 0) {
         if (reverseList[i]) PetscPrintf(PETSC_COMM_WORLD,"   path=-%s\n",pathNameList[i]->get_value().c_str());
         else PetscPrintf(PETSC_COMM_WORLD,"   path=%s\n",pathNameList[i]->get_value().c_str());
      } else {
         if (reverseList[i]) PetscPrintf(PETSC_COMM_WORLD,"   path-=%s\n",pathNameList[i]->get_value().c_str());
         else PetscPrintf(PETSC_COMM_WORLD,"   path+=%s\n",pathNameList[i]->get_value().c_str());
      }
      i++;
   }
   PetscPrintf(PETSC_COMM_WORLD,"   attribute=%d\n",attribute);
   PetscPrintf(PETSC_COMM_WORLD,"   rotated=%p\n",rotated);
   PetscPrintf(PETSC_COMM_WORLD,"EndBoundary\n");

   return;
}

bool Boundary::inBlock (int lineNumber)
{
   if (lineNumber >= startLine && lineNumber <= endLine) return true;
   return false;
}

bool Boundary::check(string *indent, vector<Path *> pathList)
{
   bool fail=false;

   // name
   if (!name.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3023: Boundary block at line %d must specify a name.\n",
                                   indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   // type
   if (! type.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3024: Block at line %d must specify a type.\n",
                                   indent->c_str(),indent->c_str(),startLine);
      return true;
   }

   // must have a path
   if (pathNameList.size() == 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3025: Boundary block at line %d must specify a path.\n",
                                   indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   // material
   if (is_surface_impedance()) {
      if (!material.is_loaded()) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3026: Boundary block at line %d must specify a material.\n",
                                      indent->c_str(),indent->c_str(),startLine);
         fail=true;
      }
   } else if (is_perfect_electric_conductor() || is_perfect_magnetic_conductor() || is_radiation()) {
      if (material.is_loaded()) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3027: Boundary block at line %d must not specify a material.\n",
                                      indent->c_str(),indent->c_str(),startLine);
         fail=true;
      }
   }

   // wave_impedance
   if (is_surface_impedance() || is_perfect_electric_conductor() || is_perfect_magnetic_conductor()) {
      if (wave_impedance.is_loaded()) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3010: Boundary block at line %d must not specify a wave impedance.\n",
                                      indent->c_str(),indent->c_str(),startLine);
         fail=true;
      }
   }

   // paths exist
   long unsigned int i=0;
   while (i < pathNameList.size()) {
      bool found=false;
      long unsigned int j=0;
      while (j < pathList.size()) {
         if (pathNameList[i]->get_value().compare(pathList[j]->get_name()) == 0) {found=true; break;}
         j++;
      }
      if (! found) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3028: Boundary block at line %d specifies a non-existent path.\n",
                                      indent->c_str(),indent->c_str(),startLine);
         fail=true;
      }
      i++;
   }

   // paths are not duplicated
   i=0;
   while (i < pathNameList.size()-1) {
      long unsigned int j=i+1;
      while (j < pathNameList.size()) {
         if (pathNameList[i]->get_value().compare(pathNameList[j]->get_value()) == 0) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3029: Boundary block at line %d duplicates path \"%s\".\n",
                                         indent->c_str(),indent->c_str(),startLine,pathNameList[j]->get_value().c_str());
            fail=true;
         }
         j++;
      }
      i++;
   }

   return fail;
}

bool Boundary::assignPathIndices (vector<Path *> *pathList)
{
   bool fail=false;

   long unsigned int i=0;
   while (i < pathNameList.size()) {
      bool found=false;
      long unsigned int j=0;
      while (j < (*pathList).size()) {
         if ((*pathList)[j]->get_name().compare(pathNameList[i]->get_value()) == 0) {
            pathIndexList.push_back(j);
            found=true;
            break;
         }
         j++;
      }
      if (! found) {
         // errors previously reported
         fail=true;
      }
      i++;
   }

   return fail;
}

bool Boundary::checkBoundingBox (Vector *lowerLeft, Vector *upperRight, string *indent, double tol, vector<Path *> *pathList)
{
   bool fail=false;

   long unsigned int i=0;
   while (i < pathIndexList.size()) {
      if ((*pathList)[pathIndexList[i]]->checkBoundingBox(lowerLeft, upperRight, indent, tol)) fail=true;
      i++;
   }
   return fail;
}

bool Boundary::merge(vector<Path *> *pathList)
{
   Path *mergedPath=nullptr;

   // merge

   bool fail=mergePaths(pathList,&pathIndexList,&reverseList,"Boundary",get_name(),&mergedPath);
   if (fail) return fail;

   if (! mergedPath) return fail;

   stringstream pathName;
   pathName << "B" << get_name() << "_OpenParEM3D_generated";
   mergedPath->set_name(pathName.str());

   pathList->push_back(mergedPath);

   // update the tracking information

   pathNameList.clear();
   keywordPair *newName=new keywordPair();
   newName->set_value(pathName.str());
   pathNameList.push_back(newName);

   pathIndexList.clear();
   pathIndexList.push_back(pathList->size()-1);

   reverseList.clear();
   reverseList.push_back(false);

   return fail;
}

bool Boundary::createRotated (vector<Path *> *pathList, string indent)
{
   if (pathIndexList.size() != 1) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3178: Boundary at line %d is incorrectly formatted.\n",
                                    indent.c_str(),indent.c_str(),startLine);
      return true;
   }

   if (rotated != nullptr) delete rotated; 
   rotated=(*pathList)[pathIndexList[0]]->rotateToXYplane();

   if (! rotated) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3179: Boundary at line %d does not form a closed polygon with nonzero area.\n",
                                    indent.c_str(),indent.c_str(),startLine);
      return true;
   }
   return false;
}

bool Boundary::is_point_inside (double x, double y, double z)
{
   if (! rotated) return false;  // not a closed boundary because it failed the rotation operation, so no point can be inside
   if (rotated->is_point_inside (x,y,z)) return true;
   return false;
}

bool Boundary::is_triangleInside (DenseMatrix *pointMat)
{
   if (is_point_inside (pointMat->Elem(0,0),pointMat->Elem(1,0),pointMat->Elem(2,0)) &&
       is_point_inside (pointMat->Elem(0,1),pointMat->Elem(1,1),pointMat->Elem(2,1)) &&
       is_point_inside (pointMat->Elem(0,2),pointMat->Elem(1,2),pointMat->Elem(2,2))) {
      return true;
   }
   return false;
}

bool Boundary::is_overlapPath (vector<Path *> *pathList, Path *testPath)
{
   if (pathIndexList.size() != 1) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Boundary::is_overlapPath operation on a Boundary with an invalid path definition.\n");
      return false;
   }

   // any point interior constitutes an overlap
   long unsigned int i=0;
   while (i < testPath->get_points_size()) {
      if (rotated->is_point_interior(testPath->get_point_x(i),testPath->get_point_y(i),testPath->get_point_z(i))) return true;
      i++;
   }

   // all points inside constitutes an overlap
   bool outside=false;
   i=0;
   while (i < testPath->get_points_size()) {
      if (! rotated->is_point_inside(testPath->get_point_x(i),testPath->get_point_y(i),testPath->get_point_z(i))) {
         outside=true;
         break;
      }
      i++;
   }
   if (! outside) return true;

   // crossing lines constitutes an overlap
   if (rotated->is_path_overlap(testPath)) return true;

   return false;
}

void Boundary::addMassImpedanceIntegrator (double frequency, double temperature, ParMesh *pmesh, ParMixedBilinearForm *pmblf, 
                                           MaterialDatabase *materialDatabase, vector<Array<int> *> &borderAttributesList,
                                           vector<ConstantCoefficient *> &ZconstList, bool isReal)
{
   if (is_perfect_electric_conductor() || is_perfect_magnetic_conductor()) return;
   if (isReal) return;

   Array<int> *border_attributes=new Array<int>;
   border_attributes->SetSize(pmesh->bdr_attributes.Max());
   (*border_attributes)=0;                // default to nothing
   (*border_attributes)[attribute-1]=1;   // enable this boundary

   double Rs=DBL_MAX;

   if (is_surface_impedance()) {
      Material *material=materialDatabase->get(get_material());
      string indent="    ";
      Rs=material->get_Rs(temperature,frequency,1e-12,indent);
   }

   if (is_radiation()) {
      Rs=get_wave_impedance();
   }

//xxx
//PetscMPIInt rank;
//MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//cout << rank << ": name=" << get_name() << "  type=" << get_type() << "  attribute=" << get_attribute() << "  Rs=" << Rs << endl;

   double coef=2*M_PI*frequency*4e-7*M_PI/Rs;
   ConstantCoefficient *Zconst=new ConstantCoefficient(coef);

   pmblf->AddBoundaryIntegrator(new VectorFEMassIntegrator(*Zconst),*border_attributes);

   // save for later deleting
   borderAttributesList.push_back(border_attributes);
   ZconstList.push_back(Zconst);
}

Boundary::~Boundary()
{
   if (rotated != nullptr) {delete rotated; rotated=nullptr;}
}

///////////////////////////////////////////////////////////////////////////////////////////
// Weights
///////////////////////////////////////////////////////////////////////////////////////////

Weights::Weights (int Sport_, complex<double> e0_, complex<double> e1_, complex<double> e2_, complex<double> e3_,
                              complex<double> h0_, complex<double> h1_, complex<double> h2_, complex<double> h3_,
                              complex<double> voltage2D_, complex<double> Zo2D_, complex<double> voltage3D_)
{
   drivenSport=Sport_;
   e0=e0_; e1=e1_; e2=e2_; e3=e3_;
   h0=h0_; h1=h1_; h2=h2_; h3=h3_;
   voltage2D=voltage2D_;
   Zo2D=Zo2D_;
   voltage3D=voltage3D_;
}

bool Weights::isPortMatch (int Sport)
{
   if (drivenSport == Sport) return true;
   return false;
}

complex<double> Weights::calculateSii ()
{
   return (e0/e2+h0/h2)/(e0/e2-h0/h2);
}

complex<double> Weights::calculateSij (Weights *driving) 
{
   complex<double> driving_c2=0.5*(driving->e0/driving->e2-driving->h0/driving->h2);

   // method using voltage integrated from the 3D fields
   return voltage3D/(driving_c2*driving->voltage3D)*sqrt(driving->Zo2D/Zo2D);

   // methods using 2D modal voltage - induces energy conservation errors
   //S=(e0+e1)/(e2+e3)*voltage2D/(driving_c1*driving->voltage2D)*sqrt(driving->Zo2D/Zo2D);
   //S=e0/e2*voltage2D/(driving_c1*driving->voltage2D)*sqrt(driving->Zo2D/Zo2D);
}

void Weights::print (string indent)
{
   PetscPrintf(PETSC_COMM_WORLD,"%s   Weights:\n",indent.c_str());
   PetscPrintf(PETSC_COMM_WORLD,"%s      drivenSport=%d\n",indent.c_str(),drivenSport);
   PetscPrintf(PETSC_COMM_WORLD,"%s      e0=(%g,%g)\n",indent.c_str(),real(e0),imag(e0));
   PetscPrintf(PETSC_COMM_WORLD,"%s      e1=(%g,%g)\n",indent.c_str(),real(e1),imag(e1));
   PetscPrintf(PETSC_COMM_WORLD,"%s      e2=(%g,%g)\n",indent.c_str(),real(e2),imag(e2));
   PetscPrintf(PETSC_COMM_WORLD,"%s      e3=(%g,%g)\n",indent.c_str(),real(e3),imag(e3));
   PetscPrintf(PETSC_COMM_WORLD,"%s      h0=(%g,%g)\n",indent.c_str(),real(h0),imag(h0));
   PetscPrintf(PETSC_COMM_WORLD,"%s      h1=(%g,%g)\n",indent.c_str(),real(h1),imag(h1));
   PetscPrintf(PETSC_COMM_WORLD,"%s      h2=(%g,%g)\n",indent.c_str(),real(h2),imag(h2));
   PetscPrintf(PETSC_COMM_WORLD,"%s      h3=(%g,%g)\n",indent.c_str(),real(h3),imag(h3));
   PetscPrintf(PETSC_COMM_WORLD,"%s      voltage2D=(%g,%g)\n",indent.c_str(),real(voltage2D),imag(voltage2D));
   PetscPrintf(PETSC_COMM_WORLD,"%s      Zo2D=(%g,%g)\n",indent.c_str(),real(Zo2D),imag(Zo2D));
   PetscPrintf(PETSC_COMM_WORLD,"%s      voltage3D=(%g,%g)\n",indent.c_str(),real(voltage3D),imag(voltage3D));
}

///////////////////////////////////////////////////////////////////////////////////////////
// Mode
///////////////////////////////////////////////////////////////////////////////////////////

Mode::Mode(int startLine_, int endLine_, string calculation_)
{
   startLine=startLine_;
   endLine=endLine_;

   // mode
   Sport.push_alias("Sport");
   Sport.set_loaded(false);
   Sport.set_positive_required(true);
   Sport.set_non_negative_required(false);
   Sport.set_lowerLimit(1);
   Sport.set_upperLimit(100);
   Sport.set_checkLimits(true);

   // type
   type.push_alias("type");
   type.set_loaded(false);
   type.set_positive_required(false);
   type.set_non_negative_required(false);
   type.set_lowerLimit(0);
   type.set_upperLimit(0);
   type.set_checkLimits(false);

   calculation=calculation_;
   isUsed=false;
   isSolutionLoaded=false;
}

bool Mode::load(string *indent, inputFile *inputs)
{
   bool fail=false;
   bool found_first_path=false;

   int lineNumber=inputs->get_next_lineNumber(startLine);
   int stopLineNumber=inputs->get_previous_lineNumber(endLine);
   while (lineNumber <= stopLineNumber) {

      string token,value,line;
      line=inputs->get_line(lineNumber);
      get_token_pair(&line,&token,&value,&lineNumber,*indent);

      int recognized=0;

      if (Sport.match_alias(&token)) {
         recognized++;
         if (Sport.loadInt(&token, &value, lineNumber)) fail=true;
      }

      if (type.match_alias(&token)) {
         if (type.is_loaded()) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3030: Duplicate entry at line %d for previous entry at line %d.\n",
                                         indent->c_str(),indent->c_str(),lineNumber,type.get_lineNumber());
            fail=true;
         } else {
            type.set_keyword(token);
            type.set_value(value);
            type.set_lineNumber(lineNumber);
            type.set_loaded(true);
         }
         recognized++;
      }

      if (token.compare("path") == 0) {
         if (found_first_path) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3031: Extraneous path= statement at line %d.\n",
                                         indent->c_str(),indent->c_str(),lineNumber);
            fail=true;
         } else {
            bool reverse=false;
            if (value.substr(0,1).compare("+") == 0) value=value.substr(1);
            else if (value.substr(0,1).compare("-") == 0) {value=value.substr(1); reverse=true;}

            keywordPair *path=new keywordPair();
            path->push_alias("path");
            path->set_keyword(token);
            path->set_value(value);
            path->set_lineNumber(lineNumber);
            path->set_positive_required(false);
            path->set_non_negative_required(false);
            path->set_lowerLimit(0);
            path->set_upperLimit(0);
            path->set_checkLimits(false);
            path->set_loaded(true);

            pathNameList.push_back(path);
            reverseList.push_back(reverse);
         }

         found_first_path=true;
         recognized++;
      }

      if (token.compare("path+") == 0) {
         if (found_first_path) {
            if (value.substr(0,1).compare("+")  == 0 || value.substr(0,1).compare("-") == 0) {
               PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3032: Misformatted path at line %d.\n",
                                            indent->c_str(),indent->c_str(),lineNumber);
               fail=true;
            } else {

               keywordPair *path=new keywordPair();
               path->push_alias("path");
               path->set_keyword(token);
               path->set_value(value);
               path->set_lineNumber(lineNumber);
               path->set_positive_required(false);
               path->set_non_negative_required(false);
               path->set_lowerLimit(0);
               path->set_upperLimit(0);
               path->set_checkLimits(false);
               path->set_loaded(true);

               pathNameList.push_back(path);
               reverseList.push_back(false);
            }
         } else {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3033: Missing path= statement before line %d.\n",
                                         indent->c_str(),indent->c_str(),lineNumber);
            fail=true;
         }
         recognized++;
      }

      if (token.compare("path-") == 0) {
         if (found_first_path) {
            if (value.substr(0,1).compare("+") == 0 || value.substr(0,1).compare("-") == 0) {
               PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3034: Misformatted path at line %d.\n",
                                            indent->c_str(),indent->c_str(),lineNumber);
               fail=true;
            } else {

               keywordPair *path=new keywordPair();
               path->push_alias("path");
               path->set_keyword(token);
               path->set_value(value);
               path->set_lineNumber(lineNumber);
               path->set_positive_required(false);
               path->set_non_negative_required(false);
               path->set_lowerLimit(0);
               path->set_upperLimit(0);
               path->set_checkLimits(false);
               path->set_loaded(true);

               pathNameList.push_back(path);
               reverseList.push_back(true);
            }
         } else {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3035: Missing path= statement before line %d.\n",
                                         indent->c_str(),indent->c_str(),lineNumber);
            fail=true;
         }
         recognized++;
      }

      // should recognize one keyword
      if (recognized != 1) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3036: Unrecognized keyword at line %d.\n",
                                      indent->c_str(),indent->c_str(),lineNumber);
         fail=true;
      }
      lineNumber=inputs->get_next_lineNumber(lineNumber);
   }

   return fail;
}

bool Mode::inModeBlock (int lineNumber)
{
   if (lineNumber >= startLine && lineNumber <= endLine) return true;
   return false;
}

bool Mode::check(string *indent, vector<Path *> pathList)
{
   bool fail=false;

   // Sport
   if (!Sport.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3037: Mode block at line %d must specify an Sport number.\n",
                                   indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   // type
   if (type.is_loaded()) {
      if (!is_voltage() && !is_current()) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3038: Input at line %d must be \"voltage\" or \"current\".\n",
                                      indent->c_str(),indent->c_str(),type.get_lineNumber());
         fail=true;
      }
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3039: Mode block at line %d must specify a type.\n",
                                   indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   // at least one path is specified
   if (pathNameList.size() == 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3040: Mode block at line %d must specify a path.\n",
                                   indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   // paths exist
   long unsigned int i=0;
   while (i < pathNameList.size()) {
      bool found=false;
      long unsigned int j=0;
      while (j < pathList.size()) {
         if (pathNameList[i]->get_value().compare(pathList[j]->get_name()) == 0) {found=true; break;}
         j++;
      }
      if (! found) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3041: Mode block at line %d specifies a non-existent path.\n",
                                      indent->c_str(),indent->c_str(),startLine);
         fail=true;
      }
      i++;
   }

   // paths are not duplicated
   i=0;
   while (pathNameList.size() > 0 && i < pathNameList.size()-1) {
      long unsigned int j=i+1;
      while (j < pathNameList.size()) {
         if (pathNameList[i]->get_value().compare(pathNameList[j]->get_value()) == 0) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3042: Mode block at line %d duplicates path \"%s\".\n",
                                         indent->c_str(),indent->c_str(),startLine,pathNameList[j]->get_value().c_str());
            fail=true;
         }
         j++;
      }
      i++;
   }

   return fail;
}

bool Mode::checkBoundingBox (Vector *lowerLeft, Vector *upperRight, string *indent, double tol, vector<Path *> *pathList)
{
   bool fail=false;

   long unsigned int i=0;
   while (i < pathIndexList.size()) {
      if ((*pathList)[pathIndexList[i]]->checkBoundingBox(lowerLeft, upperRight, indent, tol)) fail=true;
      i++;
   }
   return fail;
}

// for currents, paths must form 1 or more closed loops
bool Mode::check_current_paths (string *indent, vector<Path *> *pathList, bool check_closed_loop)
{
   bool fail=false;

   if (!is_current()) return fail;

   vector<bool> closed;
   vector<bool> connectedStart;
   vector<bool> connectedEnd;

   // to keep track of what has been looked at
   long unsigned int i=0;
   while (i < pathIndexList.size()) {
      if ((*pathList)[pathIndexList[i]]->is_closed()) {
         closed.push_back(true);
         connectedStart.push_back(true);
         connectedEnd.push_back(true);
      } else {
         closed.push_back(false);
         connectedStart.push_back(false);
         connectedEnd.push_back(false);
      }
      i++;
   }

   // line up the ends of the open sections
   i=0;
   while (pathIndexList.size() > 0 && i < pathIndexList.size()-1) {
      if (! closed[i]) {

         long unsigned int j=i+1;
         while (j < pathIndexList.size()) {
            if (! closed[j]) {

               // start to start
               if ((*pathList)[pathIndexList[i]]->get_startPoint()->point_compare((*pathList)[pathIndexList[j]]->get_startPoint())) {
                  if (! connectedStart[j]) {
                     connectedStart[i]=true;
                     connectedStart[j]=true;
                  } else {
                     PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3043: Mode block at line %d topology error at (%g,%g,%g).\n",
                                                  indent->c_str(),indent->c_str(),startLine,
                                                  (*pathList)[pathIndexList[j]]->get_startPoint()->get_point_value_x(),
                                                  (*pathList)[pathIndexList[j]]->get_startPoint()->get_point_value_y(),
                                                  (*pathList)[pathIndexList[j]]->get_startPoint()->get_point_value_z());
                     fail=true;
                  }
               }

               // start to end
               if ((*pathList)[pathIndexList[i]]->get_startPoint()->point_compare((*pathList)[pathIndexList[j]]->get_endPoint())) {
                  if (! connectedEnd[j]) {
                     connectedStart[i]=true;
                     connectedEnd[j]=true;
                  } else {
                     PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3044: Mode block at line %d topology error at (%g,%g,%g).\n",
                                                  indent->c_str(),indent->c_str(),startLine,
                                                  (*pathList)[pathIndexList[j]]->get_endPoint()->get_point_value_x(),
                                                  (*pathList)[pathIndexList[j]]->get_endPoint()->get_point_value_y(),
                                                  (*pathList)[pathIndexList[j]]->get_endPoint()->get_point_value_z());
                     fail=true;
                  }
               }

               // end to start
               if ((*pathList)[pathIndexList[i]]->get_endPoint()->point_compare((*pathList)[pathIndexList[j]]->get_startPoint())) {
                  if (! connectedStart[j]) {
                     connectedEnd[i]=true;
                     connectedStart[j]=true;
                  } else {
                     PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3045: Mode block at line %d topology error at (%g,%g,%g).\n",
                                                  indent->c_str(),indent->c_str(),startLine,
                                                  (*pathList)[pathIndexList[j]]->get_startPoint()->get_point_value_x(),
                                                  (*pathList)[pathIndexList[j]]->get_startPoint()->get_point_value_y(),
                                                  (*pathList)[pathIndexList[j]]->get_startPoint()->get_point_value_z());
                     fail=true;
                  }
               }

               // end to end
               if ((*pathList)[pathIndexList[i]]->get_endPoint()->point_compare((*pathList)[pathIndexList[j]]->get_endPoint())) {
                  if (! connectedEnd[j]) {
                     connectedEnd[i]=true;
                     connectedEnd[j]=true;
                  } else {
                     PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3046: Mode block at line %d topology error at (%g,%g,%g).\n",
                                                  indent->c_str(),indent->c_str(),startLine,
                                                  (*pathList)[pathIndexList[j]]->get_endPoint()->get_point_value_x(),
                                                  (*pathList)[pathIndexList[j]]->get_endPoint()->get_point_value_y(),
                                                  (*pathList)[pathIndexList[j]]->get_endPoint()->get_point_value_z());
                     fail=true;
                  }
               }

            }
            j++;
         }
      }
      i++;
   }

   // check for dangling ends
   if (check_closed_loop) {
      i=0;
      while (i < pathIndexList.size()) {
         if (! closed[i]) {
            if (! connectedStart[i]) {
               PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3047: Mode block at line %d topology error with dangling point at (%g,%g,%g).\n",
                                            indent->c_str(),indent->c_str(),startLine,
                                            (*pathList)[pathIndexList[i]]->get_startPoint()->get_point_value_x(),
                                            (*pathList)[pathIndexList[i]]->get_startPoint()->get_point_value_y(),
                                            (*pathList)[pathIndexList[i]]->get_startPoint()->get_point_value_z());
               fail=true;
            }
            if (! connectedEnd[i]) {
               PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3048: Mode block at line %d topology error with dangling point at (%g,%g,%g).\n",
                                            indent->c_str(),indent->c_str(),startLine,
                                            (*pathList)[pathIndexList[i]]->get_endPoint()->get_point_value_x(),
                                            (*pathList)[pathIndexList[i]]->get_endPoint()->get_point_value_y(),
                                            (*pathList)[pathIndexList[i]]->get_endPoint()->get_point_value_z());
               fail=true;
            }
         }
         i++;
      }
   }

   return fail;
}

void Mode::print(string indent)
{
   if (is_modal()) PetscPrintf(PETSC_COMM_WORLD,"%sMode\n",indent.c_str());
   if (is_line()) PetscPrintf(PETSC_COMM_WORLD,"%sLine\n",indent.c_str());

   PetscPrintf(PETSC_COMM_WORLD,"%s   type=%s\n",indent.c_str(),type.get_value().c_str());
   long unsigned int i=0;
   while (i < pathNameList.size()) {
      if (i == 0) {
         if (reverseList[i]) PetscPrintf(PETSC_COMM_WORLD,"%s   path=-%s\n",indent.c_str(),pathNameList[i]->get_value().c_str());
         else PetscPrintf(PETSC_COMM_WORLD,"%s   path=%s\n",indent.c_str(),pathNameList[i]->get_value().c_str());
      } else {
         if (reverseList[i]) PetscPrintf(PETSC_COMM_WORLD,"%s   path-=%s\n",indent.c_str(),pathNameList[i]->get_value().c_str());
         else PetscPrintf(PETSC_COMM_WORLD,"%s   path+=%s\n",indent.c_str(),pathNameList[i]->get_value().c_str());
      }
      i++;
   }

   PetscPrintf(PETSC_COMM_WORLD,"%s   isUsed=%s\n",indent.c_str(),convertLogic(isUsed).c_str());
   PetscPrintf(PETSC_COMM_WORLD,"%s   isSolutionLoaded=%s\n",indent.c_str(),convertLogic(isSolutionLoaded).c_str());
   PetscPrintf(PETSC_COMM_WORLD,"%s   modeNumber2D=%d\n",indent.c_str(),modeNumber2D);
   i=0;
   while (i < weightsList.size()) {
      weightsList[i]->print(indent);
      i++;
   }
   if (is_modal()) PetscPrintf(PETSC_COMM_WORLD,"%sEndMode\n",indent.c_str());
   if (is_line()) PetscPrintf(PETSC_COMM_WORLD,"%sEndLine\n",indent.c_str());

   return;
}

void printMatlabComplex(double re, double im)
{
   cout << setprecision(16) << scientific;
   cout << re;
   if (im >= 0) cout << " + ";
   else cout << " - ";
   cout << abs(im);
   cout << "i";
   cout << endl;
}

void Mode::printSolution()
{
   PetscMPIInt rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   if (rank == 0) {
      cout << "Mode S-port=" << get_Sport() << endl;
      cout << "   alpha=" << alpha << endl;
      cout << "   beta=" << beta << endl;

      cout << "   E:" << endl;
      int i=0;
      while (i < eigenVecReEt->Size()) {
         printMatlabComplex((*eigenVecReEt)(i),(*eigenVecImEt)(i));
         i++;
      }

      i=0;
      while (i < eigenVecReEz->Size()) {
         printMatlabComplex((*eigenVecReEz)(i),(*eigenVecImEz)(i));
         i++;
      }

      cout << "   H:" << endl;
      i=0;
      while (i < eigenVecReHt->Size()) {
         printMatlabComplex((*eigenVecReHt)(i),(*eigenVecImHt)(i));
         i++;
      }

      i=0;
      while (i < eigenVecReHz->Size()) {
         printMatlabComplex((*eigenVecReHz)(i),(*eigenVecImHz)(i));
         i++;
      }

   }
}

bool Mode::assignPathIndices(vector<Path *> *pathList)
{
   bool fail=false;

   long unsigned int i=0;
   while (i < pathNameList.size()) {
      bool found=false;
      long unsigned int j=0;
      while (j < (*pathList).size()) {
         if ((*pathList)[j]->get_name().compare(pathNameList[i]->get_value()) == 0) {
            pathIndexList.push_back(j);
            found=true;
            break;
         }
         j++;
      }
      if (! found) {
         // errors previously reported
         fail=true;
      }
      i++;
   }

   return fail;
}

bool Mode::is_enclosedByPath (vector<Path *> *pathList, Path *testPath)
{
   long unsigned int i=0;
   while (i < pathIndexList.size()) {
      if (! testPath->is_path_inside((*pathList)[pathIndexList[i]])) return false;
      i++;
   }
   return true;
}

void Mode::output (ofstream *out, vector<Path *> *pathList, Path *rotatedPath, bool spin180degrees)
{
   Path *path;

   // paths
   long unsigned int i=0;
   while (i < pathIndexList.size()) {
      path=(*pathList)[pathIndexList[i]]->clone();
      path->rotateToPath(rotatedPath,spin180degrees);
      path->output(out,2);  // drop the z component since the path is rotated
      delete path;
      *out << endl;
      i++;
   }

   // mode
   if (is_modal()) *out << "Mode" << endl;
   if (is_line()) *out << "Line" << endl;
   if (is_modal()) *out << "   mode=" << modeNumber2D << endl;
   if (is_line()) *out << "   line=" << modeNumber2D << endl;
   *out << "   type=" << get_type() << endl;
   i=0;
   while (i < pathNameList.size()) {
      if (i == 0) {
         if (reverseList[0]) *out << "   path=-" << pathNameList[0]->get_value() << endl;
         else *out << "   path=" << pathNameList[0]->get_value() << endl;
      } else {
         if (reverseList[i]) *out << "   path-=" << pathNameList[i]->get_value() << endl;
         else *out << "   path+=" << pathNameList[i]->get_value() << endl;
      }
      i++;
   }
   if (is_modal()) *out << "EndMode" << endl;
   if (is_line()) *out << "EndLine" << endl << endl;
}

bool Mode::loadSolution(string *directory, string portName, size_t t_size, size_t z_size)
{
   char filename[128];
   size_t vecSize;

   if (!isUsed) return false;

   // cd to the project directory

   stringstream projDirectory;
   projDirectory << *directory << "/S" << portName;

   try {
      std::filesystem::current_path(projDirectory.str().c_str());
   } catch (std::filesystem::filesystem_error const& ex) {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3049: Missing project directory for the 2D solution of Port %s.\n",portName.c_str());
      return true;
   }

   // cd to the temp directory

   stringstream tempDirectory;
   tempDirectory << "temp_S" << portName;

   try {
      std::filesystem::current_path(tempDirectory.str().c_str());
   } catch (std::filesystem::filesystem_error const& ex) {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3050: Missing temporary directory for the 2D solution of Port %s.\n",portName.c_str());
      std::filesystem::current_path("../../");
      return true;
   }

   // Efield

   sprintf(filename,"Efield_mode_%d.dat",modeNumber2D);
   ifstream ssEigenVecE;
   ssEigenVecE.open(filename,ifstream::in);
   if (ssEigenVecE.is_open()) {
      if (eVecReE) free(eVecReE);
      if (eVecImE) free(eVecImE);
      bool fail=loadData (&ssEigenVecE,&eVecReE,&eVecImE,&vecSize,filename);
      if (! fail) {

         if (vecSize == t_size+z_size) {

            // for building 2D grids

            // Et

            if (eigenVecReEt) delete eigenVecReEt;
            eigenVecReEt=new Vector(eVecReE,t_size);

            if (eigenVecImEt) delete eigenVecImEt;
            eigenVecImEt=new Vector (eVecImE,t_size);
            // Ez

            if (eigenVecReEz) delete eigenVecReEz;
            eigenVecReEz=new Vector;
            eigenVecReEz->SetDataAndSize(eVecReE+t_size,z_size);

            if (eigenVecImEz) delete eigenVecImEz;
            eigenVecImEz=new Vector;
            eigenVecImEz->SetDataAndSize(eVecImE+t_size,z_size);

         } else {
            PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Mode::loadSolution: Mismatched data sizes.\n");
            std::filesystem::current_path("../../../");
            fail=true;
         }
      }
      ssEigenVecE.close();
      if (fail) return true;
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3051: Unable to open file \"%s\" for reading.\n",filename);
      std::filesystem::current_path("../../../");
      return true;
   }

   // Hfield

   sprintf(filename,"Hfield_mode_%d.dat",modeNumber2D);
   ifstream ssEigenVecH;
   ssEigenVecH.open(filename,ifstream::in);
   if (ssEigenVecH.is_open()) {
      if (eVecReH) free(eVecReH);
      if (eVecImH) free(eVecImH);
      bool fail=loadData (&ssEigenVecH,&eVecReH,&eVecImH,&vecSize,filename);
      if (! fail) {

         if (vecSize == t_size+z_size) {

            // Ht

            if (eigenVecReHt) delete eigenVecReHt;
            eigenVecReHt=new Vector (eVecReH,t_size);

            if (eigenVecImHt) delete eigenVecImHt;
            eigenVecImHt=new Vector (eVecImH,t_size);

            // Hz

            if (eigenVecReHz) delete eigenVecReHz;
            eigenVecReHz=new Vector;
            eigenVecReHz->SetDataAndSize(eVecReH+t_size,z_size);

            if (eigenVecImHz) delete eigenVecImHz;
            eigenVecImHz=new Vector;
            eigenVecImHz->SetDataAndSize(eVecImH+t_size,z_size);

         } else {
            PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Mode::loadSolution: Mismatched data sizes.\n");
            std::filesystem::current_path("../../../");
            fail=true;
         }
      }
      ssEigenVecE.close();
      if (fail) return true;
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3052: Unable to open file \"%s\" for reading.\n",filename);
      std::filesystem::current_path("../../../");
      return true;
   }

   // cd back to the 3D project directory
   std::filesystem::current_path("../../../");

   isSolutionLoaded=true;

   return false;
}

// scale by the largest Et component
bool Mode::scaleSolution ()
{
   double mag2,magMax=0;
   double indexMax=0;

   if (!isUsed) return false;

   // find the max Et value
   int i=0;
   while (i < eigenVecReEt->Size()) {
      mag2=(*eigenVecReEt)[i]*(*eigenVecReEt)[i]+(*eigenVecImEt)[i]*(*eigenVecImEt)[i];
      if (mag2 > magMax) {indexMax=i; magMax=mag2;}
      i++;
   }
   complex<double> maxEt=complex<double>((*eigenVecReEt)(indexMax),(*eigenVecImEt)(indexMax));

   // scale Et, Ht

   i=0;
   while (i < eigenVecReEt->Size()) {
      complex<double> E=complex<double>((*eigenVecReEt)(i),(*eigenVecImEt)(i))/maxEt;
      (*eigenVecReEt)(i)=real(E);
      (*eigenVecImEt)(i)=imag(E);

      complex<double> H=complex<double>((*eigenVecReHt)(i),(*eigenVecImHt)(i))/maxEt;
      (*eigenVecReHt)(i)=real(H);
      (*eigenVecImHt)(i)=imag(H);

      i++;
   }

   // scale Ez, Hz
   i=0;
   while (i < eigenVecReEz->Size()) {
      complex<double> E=complex<double>((*eigenVecReEz)(i),(*eigenVecImEz)(i))/maxEt;
      (*eigenVecReEz)(i)=real(E);
      (*eigenVecImEz)(i)=imag(E);

      complex<double> H=complex<double>((*eigenVecReHz)(i),(*eigenVecImHz)(i))/maxEt;
      (*eigenVecReHz)(i)=real(H);
      (*eigenVecImHz)(i)=imag(H);

      i++;
   }

   return false;
}

// Set the voltage sign so that it aligns with the user's definition.
// This drives the requirement that all ports must have a voltage line defined
// even if using the PI definition for impedance.
void Mode::setSign (double *realV, double *imagV)
{
   PetscMPIInt rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   if (!isUsed) return;

   complex<double> scale=complex<double>(1,0);
   if (*realV < 0) {
      scale=-scale;
      *realV=-*realV;
      *imagV=-*imagV;
   }

   // scale Et, Ht

   int i=0;
   while (i < eigenVecReEt->Size()) {
      complex<double> E=complex<double>((*eigenVecReEt)(i),(*eigenVecImEt)(i))/scale;
      (*eigenVecReEt)(i)=real(E);
      (*eigenVecImEt)(i)=imag(E);

      complex<double> H=complex<double>((*eigenVecReHt)(i),(*eigenVecImHt)(i))/scale;
      (*eigenVecReHt)(i)=real(H);
      (*eigenVecImHt)(i)=imag(H);

      i++;
   }

   // scale Ez, Hz
   i=0;
   while (i < eigenVecReEz->Size()) {
      complex<double> E=complex<double>((*eigenVecReEz)(i),(*eigenVecImEz)(i))/scale;
      (*eigenVecReEz)(i)=real(E);
      (*eigenVecImEz)(i)=imag(E);

      complex<double> H=complex<double>((*eigenVecReHz)(i),(*eigenVecImHz)(i))/scale;
      (*eigenVecReHz)(i)=real(H);
      (*eigenVecImHz)(i)=imag(H);

      i++;
   }
}

void Mode::fillX (Vec *X, Vec *Xdofs, Array<int> *ess_tdof_port_list, HYPRE_BigInt *offset_port)
{
   PetscScalar value;
   HypreParVector *hypreRe=grid3DReEt->GetTrueDofs();
   HypreParVector *hypreIm=grid3DImEt->GetTrueDofs();

   HYPRE_BigInt *partitioningRe=hypreRe->Partitioning();

   // Assumes partioning is the same for grid3DReEt and grid3DImEt

   int i=0;
   while (i < ess_tdof_port_list->Size()) {
      int k=offset_port[0]+(*ess_tdof_port_list)[i];

      if (k >= partitioningRe[0] && k < partitioningRe[1]) {
         value=hypreRe->Elem(k-partitioningRe[0])+PETSC_i*hypreIm->Elem(k-partitioningRe[0]);
         VecSetValue(*X,k,value,INSERT_VALUES);
         VecSetValue(*Xdofs,k,1,INSERT_VALUES);
      }

      i++;
   }

   VecAssemblyBegin(*X);
   VecAssemblyBegin(*Xdofs);
   VecAssemblyEnd(*X);
   VecAssemblyEnd(*Xdofs);

   delete hypreRe;
   delete hypreIm;
}

void Mode::build2Dgrids (ParFiniteElementSpace *fes_ND, ParFiniteElementSpace *fes_H1)
{
   PetscMPIInt rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   HYPRE_BigInt *offset_ND=fes_ND->GetTrueDofOffsets();
   HYPRE_BigInt *offset_H1=fes_H1->GetTrueDofOffsets();

   // E

   grid2DReEt=new ParGridFunction(fes_ND);
   *grid2DReEt=0.0;
   if (offset_ND[0] < offset_ND[1]) {
      Vector local(offset_ND[1]-offset_ND[0]);
      int i=0;
      while (i < offset_ND[1]-offset_ND[0]) {
         local.Elem(i)=eigenVecReEt->Elem(offset_ND[0]+i);
         i++;
      }
      grid2DReEt->Distribute(&local);
   }

   grid2DImEt=new ParGridFunction(fes_ND);
   *grid2DImEt=0.0;
   if (offset_ND[0] < offset_ND[1]) {
      Vector local(offset_ND[1]-offset_ND[0]);
      int i=0;
      while (i < offset_ND[1]-offset_ND[0]) {
         local.Elem(i)=eigenVecImEt->Elem(offset_ND[0]+i);
         i++;
      }
      grid2DImEt->Distribute(&local);
   }

   grid2DReEz=new ParGridFunction(fes_H1);
   *grid2DReEz=0.0;
   if (offset_H1[0] < offset_H1[1]) {
      Vector local(offset_H1[1]-offset_H1[0]);
      int i=0;
      while (i < offset_H1[1]-offset_H1[0]) {
         local.Elem(i)=eigenVecReEz->Elem(offset_H1[0]+i);
         i++;
      }
      grid2DReEz->Distribute(&local);
   }

   grid2DImEz=new ParGridFunction(fes_H1);
   *grid2DImEz=0.0;
   if (offset_H1[0] < offset_H1[1]) {
      Vector local(offset_H1[1]-offset_H1[0]);
      int i=0;
      while (i < offset_H1[1]-offset_H1[0]) {
         local.Elem(i)=eigenVecImEz->Elem(offset_H1[0]+i);
         i++;
      }
      grid2DImEz->Distribute(&local);
   }

   // H

   grid2DReHt=new ParGridFunction(fes_ND);
   *grid2DReHt=0.0;
   if (offset_ND[0] < offset_ND[1]) {
      Vector local(offset_ND[1]-offset_ND[0]);
      int i=0;
      while (i < offset_ND[1]-offset_ND[0]) {
         local.Elem(i)=eigenVecReHt->Elem(offset_ND[0]+i);
         i++;
      }
      grid2DReHt->Distribute(&local);
   }

   grid2DImHt=new ParGridFunction(fes_ND);
   *grid2DImHt=0.0;
   if (offset_ND[0] < offset_ND[1]) {
      Vector local(offset_ND[1]-offset_ND[0]);
      int i=0;
      while (i < offset_ND[1]-offset_ND[0]) {
         local.Elem(i)=eigenVecImHt->Elem(offset_ND[0]+i);
         i++;
      }
      grid2DImHt->Distribute(&local);
   }

   grid2DReHz=new ParGridFunction(fes_H1);
   *grid2DReHz=0.0;
   if (offset_H1[0] < offset_H1[1]) {
      Vector local(offset_H1[1]-offset_H1[0]);
      int i=0;
      while (i < offset_H1[1]-offset_H1[0]) {
         local.Elem(i)=eigenVecReHz->Elem(offset_H1[0]+i);
         i++;
      }
      grid2DReHz->Distribute(&local);
   }

   grid2DImHz=new ParGridFunction(fes_H1);
   *grid2DImHz=0.0;
   if (offset_H1[0] < offset_H1[1]) {
      Vector local(offset_H1[1]-offset_H1[0]);
      int i=0;
      while (i < offset_H1[1]-offset_H1[0]) {
         local.Elem(i)=eigenVecImHz->Elem(offset_H1[0]+i);
         i++;
      }
      grid2DImHz->Distribute(&local);
   }
}

void Mode::build3Dgrids (ParFiniteElementSpace *fes_ND, ParFiniteElementSpace *fes_H1)
{
   // E

   grid3DReEt=new ParGridFunction(fes_ND);
   *grid3DReEt=0.0;

   grid3DImEt=new ParGridFunction(fes_ND);
   *grid3DImEt=0.0;

   grid3DReEz=new ParGridFunction(fes_H1);
   *grid3DReEz=0.0;

   grid3DImEz=new ParGridFunction(fes_H1);
   *grid3DImEz=0.0;

   // H

   grid3DReHt=new ParGridFunction(fes_ND);
   *grid3DReHt=0.0;

   grid3DImHt=new ParGridFunction(fes_ND);
   *grid3DImHt=0.0;

   grid3DReHz=new ParGridFunction(fes_H1);
   *grid3DReHz=0.0;

   grid3DImHz=new ParGridFunction(fes_H1);
   *grid3DImHz=0.0;
}

void Mode::build2DModalGrids (ParFiniteElementSpace *fes_ND, ParFiniteElementSpace *fes_H1)
{
   // E

   grid2DmodalReEt=new ParGridFunction(fes_ND);
   *grid2DmodalReEt=0.0;

   grid2DmodalImEt=new ParGridFunction(fes_ND);
   *grid2DmodalImEt=0.0;

   grid2DmodalReEz=new ParGridFunction(fes_H1);
   *grid2DmodalReEz=0.0;

   grid2DmodalImEz=new ParGridFunction(fes_H1);
   *grid2DmodalImEz=0.0;

   // H

   grid2DmodalReHt=new ParGridFunction(fes_ND);
   *grid2DmodalReHt=0.0;

   grid2DmodalImHt=new ParGridFunction(fes_ND);
   *grid2DmodalImHt=0.0;

   grid2DmodalReHz=new ParGridFunction(fes_H1);
   *grid2DmodalReHz=0.0;

   grid2DmodalImHz=new ParGridFunction(fes_H1);
   *grid2DmodalImHz=0.0;
}

void Mode::flipModalHsign(bool flip)
{
   if (flip) {
      *grid2DReHt*=-1;
      *grid2DImHt*=-1;
      *grid2DReHz*=-1;
      *grid2DImHz*=-1;
   }
}

void Mode::fillIntegrationPoints (vector<Path *> *pathList)
{
   PetscMPIInt rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   double x1,y1,z1,x2,y2,z2;

   int pointsCount=100;     // 100, number of points along each line segment for integration

   // loop through the paths
   long unsigned int iPath=0;
   while (iPath < pathIndexList.size()) {

      Path *path=(*pathList)[pathIndexList[iPath]];

      // loop through the path segments

      long unsigned int limit=path->get_points_size();
      if (path->is_closed()) limit++;

      x2=path->get_point_x(0);
      y2=path->get_point_y(0);
      z2=path->get_point_z(0);

      long unsigned int j=1;
      while (j < limit) {

         x1=x2;
         y1=y2;
         z1=z2;

         if (path->is_closed()) {
            if (j < limit-1) {
               x2=path->get_point_x(j);
               y2=path->get_point_y(j);
               z2=path->get_point_z(j);
            } else {
               x2=path->get_point_x(0);
               y2=path->get_point_y(0);
               z2=path->get_point_z(0);
            }
         } else {
            x2=path->get_point_x(j);
            y2=path->get_point_y(j);
            z2=path->get_point_z(j);
         }

         // get points along the integration line
         OPEMIntegrationPointList *points=new OPEMIntegrationPointList();
         int i=0;
         while (i < pointsCount) {
            OPEMIntegrationPoint *point=new OPEMIntegrationPoint(i,x1+i*(x2-x1)/(pointsCount-1),
                                                                   y1+i*(y2-y1)/(pointsCount-1),
                                                                   z1+i*(z2-z1)/(pointsCount-1));
            points->push(point);
            i++;
         }

         pointsList.push_back(points);

         j++;
      }

      iPath++;
   }
}

// numerically integrate
complex<double> Mode::calculateLineIntegral (ParMesh *pmesh, vector<Path *> *pathList, ParGridFunction *grid_re, ParGridFunction *grid_im)
{
   complex<double> integral=complex<double>(0,0);

   long unsigned int iPath=0;
   while (iPath < pathIndexList.size()) {
      OPEMIntegrationPointList *points=pointsList[iPath];
      points->update(pmesh);
      points->get_fieldValues(grid_re,grid_im);
      points->assemble();
      points->integrate();

      if (reverseList[iPath]) integral-=points->get_integratedValue();
      else                    integral+=points->get_integratedValue();

      iPath++;
   }

   return integral;
}

// use n x tangential field
complex<double> tangentialInnerProduct(VectorCoefficient *ReA, VectorCoefficient *ImA,
                                       VectorCoefficient *ReB, VectorCoefficient *ImB,
                                       ParGridFunction grid_fes2D_L2,
                                       ParLinearForm *lf, VectorConstantCoefficient *normal)
{
   VectorCrossProductCoefficient ReAt(*normal,*ReA);
   VectorCrossProductCoefficient ImAt(*normal,*ImA);

   VectorCrossProductCoefficient ReBt(*normal,*ReB);
   VectorCrossProductCoefficient ImBt(*normal,*ImB);

   InnerProductCoefficient ReRe(ReAt,ReBt);
   InnerProductCoefficient ReIm(ReAt,ImBt);
   InnerProductCoefficient ImRe(ImAt,ReBt);
   InnerProductCoefficient ImIm(ImAt,ImBt);

   SumCoefficient Re(ReRe,ImIm,1,+1);
   SumCoefficient Im(ReIm,ImRe,1,-1);

   grid_fes2D_L2.ProjectCoefficient(Re);
   double r=(*lf)(grid_fes2D_L2);

   grid_fes2D_L2.ProjectCoefficient(Im);
   double i=(*lf)(grid_fes2D_L2);

   return complex<double>(r,i);
}

complex<double> normalProduct(Coefficient *ReA, Coefficient *ImA,
                              Coefficient *ReB, Coefficient *ImB,
                              ParGridFunction grid_fes2D_L2,
                              ParLinearForm *lf)
{
   ProductCoefficient ReRe(*ReA,*ReB);
   ProductCoefficient ReIm(*ReA,*ImB);
   ProductCoefficient ImRe(*ImA,*ReB);
   ProductCoefficient ImIm(*ImA,*ImB);

   SumCoefficient Re(ReRe,ImIm,1,+1);
   SumCoefficient Im(ReIm,ImRe,1,-1);

   grid_fes2D_L2.ProjectCoefficient(Re);
   double r=(*lf)(grid_fes2D_L2);

   grid_fes2D_L2.ProjectCoefficient(Im);
   double i=(*lf)(grid_fes2D_L2);

   return complex<double>(r,i);
}

complex<double> normalProduct(Coefficient *ReA, Coefficient *ImA,
                              VectorCoefficient *ReB, VectorCoefficient *ImB,
                              ParGridFunction grid_fes2D_L2,
                              ParLinearForm *lf, VectorConstantCoefficient *normal)
{
   InnerProductCoefficient ReBn(*normal,*ReB);
   InnerProductCoefficient ImBn(*normal,*ImB);

   return normalProduct (ReA,ImA,&ReBn,&ImBn,grid_fes2D_L2,lf);
}

void Mode::calculateWeights (ParFiniteElementSpace *fes2D_L2,
                             ParGridFunction *grid2DsolutionReEt, ParGridFunction *grid2DsolutionImEt,
                             ParGridFunction *grid2DsolutionReEz, ParGridFunction *grid2DsolutionImEz,
                             ParGridFunction *grid2DsolutionReHt, ParGridFunction *grid2DsolutionImHt,
                             ParGridFunction *grid2DsolutionReHz, ParGridFunction *grid2DsolutionImHz,
                             Vector normal, int Sport, complex<double> voltage2D, complex<double> Zo2D, complex<double> voltage3D)
{
   VectorConstantCoefficient vccNormal(normal);

   ParGridFunction fes2D_L2_grid=ParGridFunction(fes2D_L2);

   ConstantCoefficient one(1.0);
   ParLinearForm lf(fes2D_L2);
   lf.AddDomainIntegrator(new DomainLFIntegrator(one));
   lf.Assemble();

   // mode

   VectorGridFunctionCoefficient ReEmt(grid2DmodalReEt);
   VectorGridFunctionCoefficient ImEmt(grid2DmodalImEt);
   GridFunctionCoefficient ReEmz(grid2DmodalReEz);
   GridFunctionCoefficient ImEmz(grid2DmodalImEz);

   VectorGridFunctionCoefficient ReHmt(grid2DmodalReHt);
   VectorGridFunctionCoefficient ImHmt(grid2DmodalImHt);
   GridFunctionCoefficient ReHmz(grid2DmodalReHz);
   GridFunctionCoefficient ImHmz(grid2DmodalImHz);

   // solution 

   VectorGridFunctionCoefficient ReEst(grid2DsolutionReEt);
   VectorGridFunctionCoefficient ImEst(grid2DsolutionImEt);

   GridFunctionCoefficient ReEsz(grid2DsolutionReEz);
   GridFunctionCoefficient ImEsz(grid2DsolutionImEz);

   VectorGridFunctionCoefficient ReHst(grid2DsolutionReHt);
   VectorGridFunctionCoefficient ImHst(grid2DsolutionImHt);

   GridFunctionCoefficient ReHsz(grid2DsolutionReHz);
   GridFunctionCoefficient ImHsz(grid2DsolutionImHz);

   // weight calculations

   complex<double> e0=tangentialInnerProduct(&ReEmt,&ImEmt,&ReEst,&ImEst,fes2D_L2_grid,&lf,&vccNormal);
   complex<double> e1=normalProduct(&ReEmz,&ImEmz,&ReEsz,&ImEsz,fes2D_L2_grid,&lf);
   complex<double> e2=tangentialInnerProduct(&ReEmt,&ImEmt,&ReEmt,&ImEmt,fes2D_L2_grid,&lf,&vccNormal);
   complex<double> e3=normalProduct(&ReEmz,&ImEmz,&ReEmz,&ImEmz,fes2D_L2_grid,&lf);

   complex<double> h0=tangentialInnerProduct(&ReHmt,&ImHmt,&ReHst,&ImHst,fes2D_L2_grid,&lf,&vccNormal);
   complex<double> h1=normalProduct(&ReHmz,&ImHmz,&ReHsz,&ImHsz,fes2D_L2_grid,&lf);
   complex<double> h2=tangentialInnerProduct(&ReHmt,&ImHmt,&ReHmt,&ImHmt,fes2D_L2_grid,&lf,&vccNormal);
   complex<double> h3=normalProduct(&ReHmz,&ImHmz,&ReHmz,&ImHmz,fes2D_L2_grid,&lf);

   Weights *weights=new Weights (Sport,e0,e1,e2,e3,h0,h1,h2,h3,voltage2D,Zo2D,voltage3D);
   weightsList.push_back(weights);
}

void Mode::clearWeights ()
{
   long unsigned int i=0;
   while (i < weightsList.size()) {
      if (weightsList[i]) delete weightsList[i];
      i++;
   }
   weightsList.clear();
}

complex<double> Mode::calculateSparameter (Mode *drivingMode)
{
   if (this == drivingMode) {
      long unsigned int i=0;
      while (i < weightsList.size()) {
         if (weightsList[i]->isPortMatch(get_Sport())) return weightsList[i]->calculateSii();
         i++;
      }
   } else {

      Weights *from=nullptr;
      Weights *to=nullptr;

      long unsigned int i=0;
      while (i < weightsList.size()) {
         if (weightsList[i]->isPortMatch(drivingMode->get_Sport())) {
            to=weightsList[i];
            break;
         }
         i++;
      }

      i=0;
      while (i < drivingMode->weightsList.size()) {
         if (drivingMode->weightsList[i]->isPortMatch(drivingMode->get_Sport())) {
            from=drivingMode->weightsList[i];
            break;
         }
         i++;
      }

      if (from && to) return to->calculateSij(from);
   }

   return complex<double> (DBL_MAX,DBL_MAX);
}

double Mode::get_maxReflection ()
{
   double reflection;
   double maxReflection=-DBL_MAX;

   long unsigned int i=0;
   while (i < weightsList.size()) {
      if (get_Sport() != weightsList[i]->get_drivenSport()) {
         reflection=20*log10(abs(1/weightsList[i]->calculateSii()));
         if (reflection > maxReflection) maxReflection=reflection;
      }
      i++;
   }

   return maxReflection;
}

double Mode::get_maxReflection (int drivingSport)
{
   double reflection;
   double maxReflection=-DBL_MAX;

   if (get_Sport() == drivingSport) {
      long unsigned int i=0;
      while (i < weightsList.size()) {
         if (get_Sport() != weightsList[i]->get_drivenSport()) {
            reflection=20*log10(abs(1/weightsList[i]->calculateSii()));
            if (reflection > maxReflection) maxReflection=reflection;
         }
         i++;
      }
   }

   return maxReflection;
}

void Mode::printPortReflections ()
{
   double reflection;
   PetscMPIInt rank;

   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   if (rank == 0) {
      if (weightsList.size() > 1) cout << "|      S-port " << get_Sport() << ":" << endl;
      long unsigned int i=0;
      while (i < weightsList.size()) {
         if (get_Sport() != weightsList[i]->get_drivenSport()) {
            reflection=20*log10(abs(1/weightsList[i]->calculateSii()));
            cout << "|         S-port " << weightsList[i]->get_drivenSport() << " reflection=" << reflection << " dB" << endl;
         }
         i++;
      }
   }
}

void Mode::transfer_2Dsolution_2Dgrids_to_3Dgrids()
{
   // E    

   ParTransferMap *port_to_full_ReEt=new ParTransferMap(*grid2DReEt,*grid3DReEt);
   port_to_full_ReEt->Transfer(*grid2DReEt,*grid3DReEt);
   delete port_to_full_ReEt; port_to_full_ReEt=nullptr;

   ParTransferMap *port_to_full_ImEt=new ParTransferMap(*grid2DImEt,*grid3DImEt);
   port_to_full_ImEt->Transfer(*grid2DImEt,*grid3DImEt);
   delete port_to_full_ImEt; port_to_full_ImEt=nullptr;

   ParTransferMap *port_to_full_ReEz=new ParTransferMap(*grid2DReEz,*grid3DReEz);
   port_to_full_ReEz->Transfer(*grid2DReEz,*grid3DReEz);
   delete port_to_full_ReEz; port_to_full_ReEz=nullptr;

   ParTransferMap *port_to_full_ImEz=new ParTransferMap(*grid2DImEz,*grid3DImEz);
   port_to_full_ImEz->Transfer(*grid2DImEz,*grid3DImEz);
   delete port_to_full_ImEz; port_to_full_ImEz=nullptr;

   // H    

   ParTransferMap *port_to_full_ReHt=new ParTransferMap(*grid2DReHt,*grid3DReHt);
   port_to_full_ReHt->Transfer(*grid2DReHt,*grid3DReHt);
   delete port_to_full_ReHt; port_to_full_ReHt=nullptr;

   ParTransferMap *port_to_full_ImHt=new ParTransferMap(*grid2DImHt,*grid3DImHt);
   port_to_full_ImHt->Transfer(*grid2DImHt,*grid3DImHt);
   delete port_to_full_ImHt; port_to_full_ImHt=nullptr;

   ParTransferMap *port_to_full_ReHz=new ParTransferMap(*grid2DReHz,*grid3DReHz);
   port_to_full_ReHz->Transfer(*grid2DReHz,*grid3DReHz);
   delete port_to_full_ReHz; port_to_full_ReHz=nullptr;

   ParTransferMap *port_to_full_ImHz=new ParTransferMap(*grid2DImHz,*grid3DImHz);
   port_to_full_ImHz->Transfer(*grid2DImHz,*grid3DImHz);
   delete port_to_full_ImHz; port_to_full_ImHz=nullptr;
}

// transfer from 3D back to 2D to capture the orientation operations applied by MFEM
void Mode::transfer_2Dsolution_3Dgrids_to_2Dgrids()
{
   // E

   ParTransferMap *port_to_full_ReEt=new ParTransferMap(*grid3DReEt,*grid2DmodalReEt);
   port_to_full_ReEt->Transfer(*grid3DReEt,*grid2DmodalReEt);
   delete port_to_full_ReEt; port_to_full_ReEt=nullptr;

   ParTransferMap *port_to_full_ImEt=new ParTransferMap(*grid3DImEt,*grid2DmodalImEt);
   port_to_full_ImEt->Transfer(*grid3DImEt,*grid2DmodalImEt);
   delete port_to_full_ImEt; port_to_full_ImEt=nullptr;

   ParTransferMap *port_to_full_ReEz=new ParTransferMap(*grid3DReEz,*grid2DmodalReEz);
   port_to_full_ReEz->Transfer(*grid3DReEz,*grid2DmodalReEz);
   delete port_to_full_ReEz; port_to_full_ReEz=nullptr;

   ParTransferMap *port_to_full_ImEz=new ParTransferMap(*grid3DImEz,*grid2DmodalImEz);
   port_to_full_ImEz->Transfer(*grid3DImEz,*grid2DmodalImEz);
   delete port_to_full_ImEz; port_to_full_ImEz=nullptr;

   // H

   ParTransferMap *port_to_full_ReHt=new ParTransferMap(*grid3DReHt,*grid2DmodalReHt);
   port_to_full_ReHt->Transfer(*grid3DReHt,*grid2DmodalReHt);
   delete port_to_full_ReHt; port_to_full_ReHt=nullptr;

   ParTransferMap *port_to_full_ImHt=new ParTransferMap(*grid3DImHt,*grid2DmodalImHt);
   port_to_full_ImHt->Transfer(*grid3DImHt,*grid2DmodalImHt);
   delete port_to_full_ImHt; port_to_full_ImHt=nullptr;

   ParTransferMap *port_to_full_ReHz=new ParTransferMap(*grid3DReHz,*grid2DmodalReHz);
   port_to_full_ReHz->Transfer(*grid3DReHz,*grid2DmodalReHz);
   delete port_to_full_ReHz; port_to_full_ReHz=nullptr;

   ParTransferMap *port_to_full_ImHz=new ParTransferMap(*grid3DImHz,*grid2DmodalImHz);
   port_to_full_ImHz->Transfer(*grid3DImHz,*grid2DmodalImHz);
   delete port_to_full_ImHz; port_to_full_ImHz=nullptr;
}

void Mode::save2DParaView(ParSubMesh *psubmesh2D, struct projectData *projData, double frequency, bool add_extension)
{
   if (!projData->debug_save_port_fields) return;

   stringstream ss;
   ss << projData->project_name << "_frequency_" << frequency << "_Sport_" << get_Sport();

   stringstream ssParaView;
   ssParaView << "ParaView_2D_port_" << projData->project_name;
   if (add_extension) ssParaView << "_test";

   ParaViewDataCollection *pd=nullptr;
   pd=new ParaViewDataCollection(ss.str(),psubmesh2D);
   pd->SetOwnData(false);
   pd->SetPrefixPath(ssParaView.str());
   pd->RegisterField("grid2DReEt",grid2DReEt);
   pd->RegisterField("grid2DReEz",grid2DReEz);
   pd->RegisterField("grid2DImEt",grid2DImEt);
   pd->RegisterField("grid2DImEz",grid2DImEz);
   pd->RegisterField("grid2DReHt",grid2DReHt);
   pd->RegisterField("grid2DReHz",grid2DReHz);
   pd->RegisterField("grid2DImHt",grid2DImHt);
   pd->RegisterField("grid2DImHz",grid2DImHz);
   pd->SetLevelsOfDetail(3);
   pd->SetDataFormat(VTKFormat::ASCII);
   pd->SetHighOrderOutput(true);
   pd->Save();
   delete pd;
}

void Mode::save3DParaView(ParMesh *pmesh, struct projectData *projData, double frequency, bool add_extension)
{
   if (!projData->debug_save_port_fields) return;

   stringstream ss;
   ss << projData->project_name << "_frequency_" << frequency << "_Sport_" << get_Sport();

   stringstream ssParaView;
   ssParaView << "ParaView_3D_port_" << projData->project_name;
   if (add_extension) ssParaView << "_test";

   ParaViewDataCollection *pd=nullptr;
   pd=new ParaViewDataCollection(ss.str(),pmesh);
   pd->SetOwnData(false);
   pd->SetPrefixPath(ssParaView.str());
   pd->RegisterField("grid3DReEt",grid3DReEt);
   pd->RegisterField("grid3DImEt",grid3DImEt);
   pd->RegisterField("grid3DReEz",grid3DReEz);
   pd->RegisterField("grid3DImEz",grid3DImEz);
   pd->RegisterField("grid3DReHt",grid3DReHt);
   pd->RegisterField("grid3DImHt",grid3DImHt);
   pd->RegisterField("grid3DReHz",grid3DReHz);
   pd->RegisterField("grid3DImHz",grid3DImHz);
   pd->SetLevelsOfDetail(3);
   pd->SetDataFormat(VTKFormat::ASCII);
   pd->SetHighOrderOutput(true);
   pd->Save();
   delete pd;
}

void Mode::save2DModalParaView(ParSubMesh *psubmesh2D, struct projectData *projData, double frequency, bool add_extension)
{
   if (!projData->debug_save_port_fields) return;

   stringstream ss;
   ss << projData->project_name << "_frequency_" << frequency << "_Sport_" << get_Sport();

   stringstream ssParaView;
   ssParaView << "ParaView_modal_2D_" << projData->project_name;
   if (add_extension) ssParaView << "_test";

   ParaViewDataCollection *pd=nullptr;
   pd=new ParaViewDataCollection(ss.str(),psubmesh2D);
   pd->SetOwnData(false);
   pd->SetPrefixPath(ssParaView.str());
   pd->RegisterField("grid2DmodalReEt",grid2DmodalReEt);
   pd->RegisterField("grid2DmodalImEt",grid2DmodalImEt);
   pd->RegisterField("grid2DmodalReEz",grid2DmodalReEz);
   pd->RegisterField("grid2DmodalImEz",grid2DmodalImEz);
   pd->RegisterField("grid2DmodalReHt",grid2DmodalReHt);
   pd->RegisterField("grid2DmodalImHt",grid2DmodalImHt);
   pd->RegisterField("grid2DmodalReHz",grid2DmodalReHz);
   pd->RegisterField("grid2DmodalImHz",grid2DmodalImHz);
   pd->SetLevelsOfDetail(3);
   pd->SetDataFormat(VTKFormat::ASCII);
   pd->SetHighOrderOutput(true);
   pd->Save();
   delete pd;
}

void Mode::resetElementNumbers ()
{
   long unsigned int i=0;
   while (i < pointsList.size()) {
      pointsList[i]->resetElementNumbers();
      i++;
   }
}

void Mode::reset()
{
   isSolutionLoaded=false;

   if (eVecReE) {free(eVecReE); eVecReE=nullptr;}
   if (eVecImE) {free(eVecImE); eVecImE=nullptr;}
   if (eVecReH) {free(eVecReH); eVecReH=nullptr;}
   if (eVecImH) {free(eVecImH); eVecImH=nullptr;}

   if (eigenVecReEt) {delete eigenVecReEt; eigenVecReEt=nullptr;}
   if (eigenVecImEt) {delete eigenVecImEt; eigenVecImEt=nullptr;}
   if (eigenVecReEz) {delete eigenVecReEz; eigenVecReEz=nullptr;}
   if (eigenVecImEz) {delete eigenVecImEz; eigenVecImEz=nullptr;}
   if (eigenVecReHt) {delete eigenVecReHt; eigenVecReHt=nullptr;}
   if (eigenVecImHt) {delete eigenVecImHt; eigenVecImHt=nullptr;}
   if (eigenVecReHz) {delete eigenVecReHz; eigenVecReHz=nullptr;}
   if (eigenVecImHz) {delete eigenVecImHz; eigenVecImHz=nullptr;}

   if (grid2DReEt) {delete grid2DReEt; grid2DReEt=nullptr;}
   if (grid2DImEt) {delete grid2DImEt; grid2DImEt=nullptr;}
   if (grid2DReEz) {delete grid2DReEz; grid2DReEz=nullptr;}
   if (grid2DImEz) {delete grid2DImEz; grid2DImEz=nullptr;}
   if (grid2DReHt) {delete grid2DReHt; grid2DReHt=nullptr;}
   if (grid2DImHt) {delete grid2DImHt; grid2DImHt=nullptr;}
   if (grid2DReHz) {delete grid2DReHz; grid2DReHz=nullptr;}
   if (grid2DImHz) {delete grid2DImHz; grid2DImHz=nullptr;}

   if (grid3DReEt) {delete grid3DReEt; grid3DReEt=nullptr;}
   if (grid3DImEt) {delete grid3DImEt; grid3DImEt=nullptr;}
   if (grid3DReEz) {delete grid3DReEz; grid3DReEz=nullptr;}
   if (grid3DImEz) {delete grid3DImEz; grid3DImEz=nullptr;}
   if (grid3DReHt) {delete grid3DReHt; grid3DReHt=nullptr;}
   if (grid3DImHt) {delete grid3DImHt; grid3DImHt=nullptr;}
   if (grid3DReHz) {delete grid3DReHz; grid3DReHz=nullptr;}
   if (grid3DImHz) {delete grid3DImHz; grid3DImHz=nullptr;}

   if (grid2DmodalReEt) {delete grid2DmodalReEt; grid2DmodalReEt=nullptr;}
   if (grid2DmodalImEt) {delete grid2DmodalImEt; grid2DmodalImEt=nullptr;}
   if (grid2DmodalReEz) {delete grid2DmodalReEz; grid2DmodalReEz=nullptr;}
   if (grid2DmodalImEz) {delete grid2DmodalImEz; grid2DmodalImEz=nullptr;}
   if (grid2DmodalReHt) {delete grid2DmodalReHt; grid2DmodalReHt=nullptr;}
   if (grid2DmodalImHt) {delete grid2DmodalImHt; grid2DmodalImHt=nullptr;}
   if (grid2DmodalReHz) {delete grid2DmodalReHz; grid2DmodalReHz=nullptr;}
   if (grid2DmodalImHz) {delete grid2DmodalImHz; grid2DmodalImHz=nullptr;}

   clearWeights();
}

Mode::~Mode()
{
   long unsigned int i=0;
   while (i < pathNameList.size()) {
      if (pathNameList[i]) delete pathNameList[i];
      i++;
   }

   i=0;
   while (i < pointsList.size()) {
      if (pointsList[i]) delete pointsList[i];
      i++;
   }

   i=0;
   while (i < weightsList.size()) {
      if (weightsList[i]) delete weightsList[i];
      i++;
   }

   if (eVecReE) {free(eVecReE); eVecReE=nullptr;}
   if (eVecImE) {free(eVecImE); eVecImE=nullptr;}
   if (eVecReH) {free(eVecReH); eVecReH=nullptr;}
   if (eVecImH) {free(eVecImH); eVecImH=nullptr;}

   if (eigenVecReEt) {delete eigenVecReEt; eigenVecReEt=nullptr;}
   if (eigenVecImEt) {delete eigenVecImEt; eigenVecImEt=nullptr;}
   if (eigenVecReEz) {delete eigenVecReEz; eigenVecReEz=nullptr;}
   if (eigenVecImEz) {delete eigenVecImEz; eigenVecImEz=nullptr;}
   if (eigenVecReHt) {delete eigenVecReHt; eigenVecReHt=nullptr;}
   if (eigenVecImHt) {delete eigenVecImHt; eigenVecImHt=nullptr;}
   if (eigenVecReHz) {delete eigenVecReHz; eigenVecReHz=nullptr;}
   if (eigenVecImHz) {delete eigenVecImHz; eigenVecImHz=nullptr;}

   if (grid2DReEt) {delete grid2DReEt; grid2DReEt=nullptr;}
   if (grid2DImEt) {delete grid2DImEt; grid2DImEt=nullptr;}
   if (grid2DReEz) {delete grid2DReEz; grid2DReEz=nullptr;}
   if (grid2DImEz) {delete grid2DImEz; grid2DImEz=nullptr;}
   if (grid2DReHt) {delete grid2DReHt; grid2DReHt=nullptr;}
   if (grid2DImHt) {delete grid2DImHt; grid2DImHt=nullptr;}
   if (grid2DReHz) {delete grid2DReHz; grid2DReHz=nullptr;}
   if (grid2DImHz) {delete grid2DImHz; grid2DImHz=nullptr;}

   if (grid3DReEt) {delete grid3DReEt; grid3DReEt=nullptr;}
   if (grid3DImEt) {delete grid3DImEt; grid3DImEt=nullptr;}
   if (grid3DReEz) {delete grid3DReEz; grid3DReEz=nullptr;}
   if (grid3DImEz) {delete grid3DImEz; grid3DImEz=nullptr;}
   if (grid3DReHt) {delete grid3DReHt; grid3DReHt=nullptr;}
   if (grid3DImHt) {delete grid3DImHt; grid3DImHt=nullptr;}
   if (grid3DReHz) {delete grid3DReHz; grid3DReHz=nullptr;}
   if (grid3DImHz) {delete grid3DImHz; grid3DImHz=nullptr;}

   if (grid2DmodalReEt) {delete grid2DmodalReEt; grid2DmodalReEt=nullptr;}
   if (grid2DmodalImEt) {delete grid2DmodalImEt; grid2DmodalImEt=nullptr;}
   if (grid2DmodalReEz) {delete grid2DmodalReEz; grid2DmodalReEz=nullptr;}
   if (grid2DmodalImEz) {delete grid2DmodalImEz; grid2DmodalImEz=nullptr;}
   if (grid2DmodalReHt) {delete grid2DmodalReHt; grid2DmodalReHt=nullptr;}
   if (grid2DmodalImHt) {delete grid2DmodalImHt; grid2DmodalImHt=nullptr;}
   if (grid2DmodalReHz) {delete grid2DmodalReHz; grid2DmodalReHz=nullptr;}
   if (grid2DmodalImHz) {delete grid2DmodalImHz; grid2DmodalImHz=nullptr;}

}

///////////////////////////////////////////////////////////////////////////////////////////
// Port
///////////////////////////////////////////////////////////////////////////////////////////

Port::Port(int startLine_, int endLine_)
{
   startLine=startLine_;
   endLine=endLine_;

   // port
   name.push_alias("name");
   name.set_loaded(false);
   name.set_positive_required(false);
   name.set_non_negative_required(false);
   name.set_lowerLimit(0);
   name.set_upperLimit(0);
   name.set_checkLimits(false);

   // impedance_definition
   impedance_definition.push_alias("impedance_definition");
   impedance_definition.set_loaded(false);
   impedance_definition.set_positive_required(false);
   impedance_definition.set_non_negative_required(false);
   impedance_definition.set_lowerLimit(0);
   impedance_definition.set_upperLimit(0);
   impedance_definition.set_checkLimits(false);

   // impedance_calculation
   impedance_calculation.push_alias("impedance_calculation");
   impedance_calculation.set_loaded(false);
   impedance_calculation.set_positive_required(false);
   impedance_calculation.set_non_negative_required(false);
   impedance_calculation.set_lowerLimit(0);
   impedance_calculation.set_upperLimit(0);
   impedance_calculation.set_checkLimits(false);

   normal.SetSize(3);
   normal=0.0;

   rotated_normal.SetSize(3);
   rotated_normal=0.0;
}

bool Port::load(string *indent, inputFile *inputs)
{
   bool fail=false;
   bool found_first_path=false;

   // Mode and Line
   if (findModeBlocks(inputs)) fail=true;
   if (findLineBlocks(inputs)) fail=true;

   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->load(indent, inputs)) fail=true;
      i++;
   }

   // now the keywords

   int lineNumber=inputs->get_next_lineNumber(startLine);
   int stopLineNumber=inputs->get_previous_lineNumber(endLine);
   while (lineNumber <= stopLineNumber) {

      if (!inModeBlocks(lineNumber)) {

         string token,value,line;
         line=inputs->get_line(lineNumber);
         get_token_pair(&line,&token,&value,&lineNumber,*indent);

         int recognized=0;

         if (name.match_alias(&token)) {
            if (name.is_loaded()) {
               PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3053: Duplicate entry at line %d for previous entry at line %d.\n",
                                            indent->c_str(),indent->c_str(),lineNumber,impedance_definition.get_lineNumber());
               fail=true;
            } else {
               name.set_keyword(token);
               name.set_value(value);
               name.set_lineNumber(lineNumber);
               name.set_loaded(true);
            }
            recognized++;
         }

         if (impedance_definition.match_alias(&token)) {
            if (impedance_definition.is_loaded()) {
               PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3054: Duplicate entry at line %d for previous entry at line %d.\n",
                                            indent->c_str(),indent->c_str(),lineNumber,impedance_definition.get_lineNumber());
               fail=true;
            } else {
               impedance_definition.set_keyword(token);
               impedance_definition.set_value(value);
               impedance_definition.set_lineNumber(lineNumber);
               impedance_definition.set_loaded(true);
            }
            recognized++;
         }

         if (impedance_calculation.match_alias(&token)) {
            if (impedance_calculation.is_loaded()) {
               PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3055: Duplicate entry at line %d for previous entry at line %d.\n",
                                            indent->c_str(),indent->c_str(),lineNumber,impedance_calculation.get_lineNumber());
               fail=true;
            } else {
               impedance_calculation.set_keyword(token);
               impedance_calculation.set_value(value);
               impedance_calculation.set_lineNumber(lineNumber);
               impedance_calculation.set_loaded(true);
            }
            recognized++;
         }

         if (token.compare("path") == 0) {
            if (found_first_path) {
               PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3057: Extraneous path= statement at line %d.\n",
                                            indent->c_str(),indent->c_str(),lineNumber);
               fail=true;
            } else {
               bool reverse=false;
               if (value.substr(0,1).compare("+") == 0) value=value.substr(1);
               else if (value.substr(0,1).compare("-") == 0) {value=value.substr(1); reverse=true;}

               keywordPair *path=new keywordPair();
               path->push_alias("path");
               path->set_keyword(token);
               path->set_value(value);
               path->set_lineNumber(lineNumber);
               path->set_positive_required(false);
               path->set_non_negative_required(false);
               path->set_lowerLimit(0);
               path->set_upperLimit(0);
               path->set_checkLimits(false);
               path->set_loaded(true);

               pathNameList.push_back(path);
               reverseList.push_back(reverse);
            }

            found_first_path=true;
            recognized++;
         }

         if (token.compare("path+") == 0) {
            if (found_first_path) {
               if (value.substr(0,1).compare("+")  == 0 || value.substr(0,1).compare("-") == 0) {
                  PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3058: Misformatted path at line %d.\n",
                                               indent->c_str(),indent->c_str(),lineNumber);
                  fail=true;
               } else {

                  keywordPair *path=new keywordPair();
                  path->push_alias("path");
                  path->set_keyword(token);
                  path->set_value(value);
                  path->set_lineNumber(lineNumber);
                  path->set_positive_required(false);
                  path->set_non_negative_required(false);
                  path->set_lowerLimit(0);
                  path->set_upperLimit(0);
                  path->set_checkLimits(false);
                  path->set_loaded(true);

                  pathNameList.push_back(path);
                  reverseList.push_back(false);
               }
            } else {
               PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3059: Missing path= statement before line %d.\n",
                                            indent->c_str(),indent->c_str(),lineNumber);
               fail=true;
            }
            recognized++;
         }

         if (token.compare("path-") == 0) {
            if (found_first_path) {
               if (value.substr(0,1).compare("+") == 0 || value.substr(0,1).compare("-") == 0) {
                  PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3060: Misformatted path at line %d.\n",
                                               indent->c_str(),indent->c_str(),lineNumber);
                  fail=true;
               } else {

                  keywordPair *path=new keywordPair();
                  path->push_alias("path");
                  path->set_keyword(token);
                  path->set_value(value);
                  path->set_lineNumber(lineNumber);
                  path->set_positive_required(false);
                  path->set_non_negative_required(false);
                  path->set_lowerLimit(0);
                  path->set_upperLimit(0);
                  path->set_checkLimits(false);
                  path->set_loaded(true);

                  pathNameList.push_back(path);
                  reverseList.push_back(true);
               }
            } else {
               PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3061: Missing path= statement before line %d.\n",
                                            indent->c_str(),indent->c_str(),lineNumber);
               fail=true;
            }
            recognized++;
         }

         // should recognize one keyword
         if (recognized != 1) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3062: Unrecognized keyword at line %d.\n",
                                         indent->c_str(),indent->c_str(),lineNumber);
            fail=true;
         }
      }
      lineNumber=inputs->get_next_lineNumber(lineNumber);
   }

   return fail;
}

bool Port::inBlock (int lineNumber)
{
   if (lineNumber >= startLine && lineNumber <= endLine) return true;
   return false;
}

bool Port::inModeBlocks (int lineNumber)
{
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->inModeBlock(lineNumber)) return true;
      i++;
   }
   return false;
}

bool Port::findModeBlocks(inputFile *inputs)
{
   bool fail=false;
   int start_lineNumber=startLine;
   int stop_lineNumber=endLine;
   int block_start,block_stop;

   while (start_lineNumber < stop_lineNumber) {
      if (inputs->findBlock(start_lineNumber,stop_lineNumber, &block_start, &block_stop,
                                "Mode", "EndMode", false)) {
         fail=true;
      } else {
         if (block_start >= 0 && block_stop >= 0) {
            Mode *newMode=new Mode(block_start,block_stop,"modal");
            modeList.push_back(newMode);
         }
      }
      start_lineNumber=inputs->get_next_lineNumber(block_stop);
   }
   return fail;
}

bool Port::findLineBlocks(inputFile *inputs)
{
   bool fail=false;
   int start_lineNumber=startLine;
   int stop_lineNumber=endLine;
   int block_start,block_stop;

   while (start_lineNumber < stop_lineNumber) {
      if (inputs->findBlock(start_lineNumber,stop_lineNumber, &block_start, &block_stop,
                                "Line", "EndLine", false)) {
         fail=true;
      } else {
         if (block_start >= 0 && block_stop >= 0) {
            Mode *newMode=new Mode(block_start,block_stop,"line");
            modeList.push_back(newMode);
         }
      }
      start_lineNumber=inputs->get_next_lineNumber(block_stop);
   }
   return fail;
}

bool Port::check(string *indent, vector<Path *> pathList, bool check_closed_loop)
{
   bool fail=false;

   if (!name.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3063: Port block at line %d must specify a name.\n",
                                   indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   // impedance_definition
   if (impedance_definition.is_loaded()) {
      if (!(impedance_definition.get_value().compare("VI") == 0 || 
            impedance_definition.get_value().compare("PV") == 0 ||
            impedance_definition.get_value().compare("PI") == 0)) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3064: Input at line %d must be \"VI\", \"PV\", or \"PI\".\n",
                                      indent->c_str(),indent->c_str(),impedance_definition.get_lineNumber());
         fail=true;      
      }
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3065: Port block at line %d must specify an impedance definition.\n",
                                   indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   // impedance_calculation
   if (impedance_calculation.is_loaded()) {
      if (!(impedance_calculation.get_value().compare("modal") == 0 ||
            impedance_calculation.get_value().compare("line") == 0)) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3066: Input at line %d must be \"modal\" or \"line\".\n",
                                      indent->c_str(),indent->c_str(),impedance_calculation.get_lineNumber());
         fail=true;
      }
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3067: Port block at line %d must specify an impedance calculation.\n",
                                   indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   // must have a path
   if (pathNameList.size() == 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3070: Port block at line %d must specify a path.\n",
                                   indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   // paths exist
   long unsigned int i=0;
   while (i < pathNameList.size()) {
      bool found=false;
      long unsigned int j=0;
      while (j < pathList.size()) {
         if (pathNameList[i]->get_value().compare(pathList[j]->get_name()) == 0) {found=true; break;}
         j++;
      }
      if (! found) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3071: Port block at line %d specifies a non-existent path.\n",
                                      indent->c_str(),indent->c_str(),startLine);
         fail=true;
      }
      i++;
   }

   // paths are not duplicated
   i=0;
   while (pathNameList.size() > 0 && i < pathNameList.size()-1) {
      long unsigned int j=i+1;
      while (j < pathNameList.size()) {
         if (pathNameList[i]->get_value().compare(pathNameList[j]->get_value()) == 0) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3072: Port block at line %d duplicates path \"%s\".\n",
                                         indent->c_str(),indent->c_str(),startLine,pathNameList[j]->get_value().c_str());
            fail=true;
         }
         j++;
      }
      i++;
   }

   // Mode

   // must have at least one mode
   if (modeList.size() == 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3073: Port block at line %d must specify at least one mode.\n",
                                   indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   // mode checks
   i=0;
   while (i < modeList.size()) {
      if (modeList[i]->check(indent,pathList)) fail=true;
      if (modeList[i]->check_current_paths(indent,&pathList,check_closed_loop)) fail=true;
      i++;
   }

   // each mode must have a voltage definition for S-parameter calculation even if the PI impedance definition is used
   if (!fail) {
      bool found=false;
      i=0;
      while (i < modeList.size()) {
         if (modeList[i]->is_voltage()) {found=true; break;}
         i++;
      }
      if (!found) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3074: Port block at line %d must specify a voltage line.\n",
                                      indent->c_str(),indent->c_str(),startLine);
         fail=true;
      }
   }

   // modes must be enclosed by the port boundary
   if (!fail && !is_modePathInside (indent,&pathList)) fail=true;

   return fail;
}

bool Port::is_overlapPath (Path *testPath)
{
   if (pathIndexList.size() != 1) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Port::is_overlapPath operation on a Port with an invalid path definition.\n");
      return false;
   }

   // any point interior constitutes an overlap
   long unsigned int i=0;
   while (i < testPath->get_points_size()) {
      if (rotated->is_point_interior(testPath->get_point_x(i),testPath->get_point_y(i),testPath->get_point_z(i))) return true;
      i++;
   }

   // all points inside constitutes an overlap
   bool outside=false;
   i=0;
   while (i < testPath->get_points_size()) {
      if (! rotated->is_point_inside(testPath->get_point_x(i),testPath->get_point_y(i),testPath->get_point_z(i))) {
         outside=true;
         break;
      }
      i++;
   }
   if (! outside) return true;

   // crossing lines constitutes an overlap
   if (rotated->is_path_overlap(testPath)) return true;

   return false;
}

bool Port::assignPathIndices (vector<Path *> *pathList)
{
   bool fail=false;

   long unsigned int i=0;
   while (i < pathNameList.size()) {
      bool found=false;
      long unsigned int j=0;
      while (j < (*pathList).size()) {
         if ((*pathList)[j]->get_name().compare(pathNameList[i]->get_value()) == 0) {
            pathIndexList.push_back(j);
            found=true;
            break;
         }
         j++;
      }
      if (! found) {
         // errors previously reported
         fail=true;
      }
      i++;
   }

   i=0;
   while (i < modeList.size()) {
      if (modeList[i]->assignPathIndices(pathList)) fail=true;
      i++;
   }

   return fail;
}

vector<int> Port::get_SportList ()
{
   vector<int> SportList;

   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
         SportList.push_back(modeList[i]->get_Sport());
      }
      i++;
   }

   return SportList;
}

bool Port::checkBoundingBox (Vector *lowerLeft, Vector *upperRight, string *indent, double tol, vector<Path *> *pathList)
{
   bool fail=false;

   long unsigned int i=0;
   while (i < pathIndexList.size()) {
      if ((*pathList)[pathIndexList[i]]->checkBoundingBox(lowerLeft, upperRight, indent, tol)) fail=true;
      i++;
   }

   i=0;
   while (i < modeList.size()) {
      if (modeList[i]->checkBoundingBox(lowerLeft,upperRight,indent,tol,pathList)) fail=true;
      i++;
   }

   return fail;
}

void Port::set2DModeNumbers()
{
   int modeNumber=1;
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
         modeList[i]->set_modeNumber2D(modeNumber);
         modeNumber++;
      }
      i++;
   }
}

void Port::fillUnused2DModeNumbers()
{
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
        long unsigned int j=0;
        while (j < modeList.size()) {
           if (!modeList[j]->get_isUsed() && modeList[i]->get_Sport() == modeList[j]->get_Sport()) {
              modeList[j]->set_modeNumber2D(modeList[i]->get_modeNumber2D());
           }
           j++;
         }
      }
      i++;
   }
}

void Port::print()
{
   PetscPrintf(PETSC_COMM_WORLD,"Port\n");
   PetscPrintf(PETSC_COMM_WORLD,"   name=%s\n",get_name().c_str());
   long unsigned int i=0;
   while (i < pathNameList.size()) {
      if (i == 0) {
         if (reverseList[i]) PetscPrintf(PETSC_COMM_WORLD,"   path=-%s\n",pathNameList[i]->get_value().c_str());
         else PetscPrintf(PETSC_COMM_WORLD,"   path=%s\n",pathNameList[i]->get_value().c_str());
      } else {
         if (reverseList[i]) PetscPrintf(PETSC_COMM_WORLD,"   path-=%s\n",pathNameList[i]->get_value().c_str());
         else PetscPrintf(PETSC_COMM_WORLD,"   path+=%s\n",pathNameList[i]->get_value().c_str());
      }
      i++;
   }

   PetscPrintf(PETSC_COMM_WORLD,"   impedance_definition=%s\n",impedance_definition.get_value().c_str());
   PetscPrintf(PETSC_COMM_WORLD,"   impedance_calculation=%s\n",impedance_calculation.get_value().c_str());

   i=0;
   while (i < modeList.size()) {
      modeList[i]->print("   ");
      i++;
   }

   PetscPrintf(PETSC_COMM_WORLD,"   attribute=%d\n",attribute);
   PetscPrintf(PETSC_COMM_WORLD,"   meshFilename=%s\n",meshFilename.c_str());
   PetscPrintf(PETSC_COMM_WORLD,"   modesFilename=%s\n",modesFilename.c_str());
   if (rotated) {
      PetscPrintf(PETSC_COMM_WORLD,"   rotated=%p:\n",rotated);
      rotated->print("      ");
   } else PetscPrintf(PETSC_COMM_WORLD,"   rotated=%p\n",rotated);
   PetscPrintf(PETSC_COMM_WORLD,"   outward normal=(%g,%g,%g)\n",normal(0),normal(1),normal(2));
   PetscPrintf(PETSC_COMM_WORLD,"   rotated outward normal=(%g,%g,%g)\n",rotated_normal(0),rotated_normal(1),rotated_normal(2));
   PetscPrintf(PETSC_COMM_WORLD,"EndPort\n");

   return;
}

void Port::printPaths(vector<Path *> *pathList)
{
   PetscPrintf(PETSC_COMM_WORLD,"Port=%s\n",get_name().c_str());
   PetscPrintf(PETSC_COMM_WORLD,"   Paths:\n");
   long unsigned int i=0;
   while (i < pathIndexList.size()) {
      (*pathList)[pathIndexList[i]]->print("      ");
      i++;
   }
   PetscPrintf(PETSC_COMM_WORLD,"   Rotated Path:\n");
   if (rotated) rotated->print("      ");
}

bool Port::createDirectory(string *tempDirectory)
{
   if (std::filesystem::exists(tempDirectory->c_str())) {
      stringstream PortDir;
      PortDir << *tempDirectory << "/" << "S" << get_name();
      if (std::filesystem::create_directory(PortDir.str().c_str())) return false;
      else PetscPrintf(PETSC_COMM_WORLD,"ERROR3075: Failed to create the port working directory \"%s\".\n",
                                        PortDir.str().c_str());
   }
   return true;
}

void Port::saveMesh(meshMaterialList *materials, string *directory, ParSubMesh *parSubMesh)
{
   PetscMPIInt rank,size;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   MPI_Comm_size(PETSC_COMM_WORLD, &size);

   stringstream parSaveMeshFilename,processed_parSaveMeshFilename;
   parSaveMeshFilename << *directory << "/S" << get_name() << "/par" << meshFilename;
   if (size > 1)  parSaveMeshFilename << "." << setw(6) << setfill('0') << rank;

   processed_parSaveMeshFilename << *directory << "/S" << get_name() << "/" << meshFilename;
   if (size > 1) processed_parSaveMeshFilename << "." << setw(6) << setfill('0') << rank;

   bool useSerial=false;
   if (useSerial) {
      ofstream serout(processed_parSaveMeshFilename.str().c_str());
      if (serout.is_open()) {
         serout.precision(15);

         RotatedMesh *mesh2D=new RotatedMesh();
         Mesh *mesh2Dbase;
         mesh2Dbase=mesh2D;
         *mesh2Dbase=parSubMesh->GetSerialMesh(0);
         mesh2D->rotate(rotated,spin180degrees);
         mesh2D->set_spaceDim(2);
         mesh2D->Print(serout);
         delete mesh2D;

         serout.close();
      } else {
         PetscPrintf(PETSC_COMM_WORLD,"ERROR3076: Failed to open file \"%s\" for writing.\n",
                                      parSaveMeshFilename.str().c_str());
      }
   } else {

      // save
      ofstream parout(parSaveMeshFilename.str().c_str());
      if (parout.is_open()) {
         parout.precision(15);
         parSubMesh->ParPrint(parout);
         parout.close();
      } else {
         PetscPrintf(PETSC_COMM_WORLD,"ERROR3077: Failed to open file \"%s\" for writing.\n",
                                      parSaveMeshFilename.str().c_str());
      }

      // post-process the mesh to make it fully 2D and to flip the direction the port faces, if needed
      if (postProcessMesh(parSaveMeshFilename.str(),processed_parSaveMeshFilename.str())) {
         PetscPrintf(PETSC_COMM_WORLD,"ERROR3078: Failed to post-process mesh file \"%s\".\n",
                                      parSaveMeshFilename.str().c_str());
      }
   }

   // save the materials associated with the mesh
   if (rank == 0) {
      stringstream saveRegionsFilename;
      saveRegionsFilename << *directory << "/S" << get_name() << "/materials_for_" << meshFilename;
      materials->saveRegionsFile(saveRegionsFilename.str().c_str());
   }
}

bool Port::postProcessMesh (string input_filename, string temp_output_filename)
{
   PetscMPIInt rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   ifstream INPUT;
   INPUT.open(input_filename.c_str(),ifstream::in);
   if (!INPUT.is_open()) {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3079: Failed to open file \"%s\" for reading.\n",input_filename.c_str());
      return true;
   }

   ofstream OUTPUT;
   OUTPUT.open(temp_output_filename.c_str(),ofstream::out);
   if (!OUTPUT.is_open()) {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3080: Failed to open file \"%s\" for writing.\n",temp_output_filename.c_str());
      INPUT.close();
      return true;
   }

   // process line-by-line
   bool inVertices=false;
   int dimension=-1;
   int vertexCount=-1;
   int vertexLineCount=0;
   string line;
   while (std::getline(INPUT,line)) {

      if (inVertices) {
        if (vertexLineCount == 0) {
           vertexCount=stoi(line);
           OUTPUT << vertexCount << endl;
        } else if (vertexLineCount == 1) {
           dimension=stoi(line);
           if (dimension == 3) dimension--;
           OUTPUT << dimension << endl;
        } else {

           // read the coordinate
           double x,y,z;
           int coordCount=0;
           stringstream ssLine(line);
           string value;
           while (std::getline(ssLine,value,' ')) {
              if (coordCount == 0) x=stod(value);
              else if (coordCount == 1) y=stod(value);
              else if (coordCount == 2) z=stod(value);
              coordCount++;
           }

           // rotate so the normal points in the +z direction
           rotated->rotatePoint(&x,&y,&z,spin180degrees);

           // skip z on output
           OUTPUT << setprecision(15) << x << " " << y << endl;

           ssLine.str("");
           ssLine.clear();
        }

        vertexLineCount++;
        if (vertexLineCount-2 == vertexCount) inVertices=false;

      } else {
         OUTPUT << line << endl;
      }

      if (line.compare("vertices") == 0) {
         inVertices=true;
      }
   }

   INPUT.close();
   OUTPUT.close();

   return false;
}

void Port::save2Dsetup(struct projectData *projData, string *directory, double frequency, Gamma *gamma)
{
   stringstream filename;
   filename << *directory << "/S" << get_name() << "/S" << get_name() << ".proj";

   ofstream out;
   out.open(filename.str().c_str(),ofstream::out);

   if (out.is_open()) {

      out << "#OpenParEM2Dproject 1.0" << endl << endl;
      out << "project.save.fields                      true" << endl << endl;            // must save fields, so override the proj file
      out << "mesh.file                                " << meshFilename << endl;
      out << "mesh.order                               " << projData->mesh_order << endl;
      out << "mesh.refinement.fraction                 " << projData->mesh_2D_refinement_fraction << endl;
      out << "mesh.uniform_refinement.count            " << "0" << endl;
      out << "mesh.enable.refine                       " << "0" << endl << endl; 
      out << "mode.definition.file                     " << modesFilename << endl;
      out << "materials.global.path                    " << "../../" << projData->materials_global_path << endl;
      out << "materials.global.name                    " << projData->materials_global_name << endl;
      out << "materials.local.path                     " << "../../" << projData->materials_local_path << endl;
      out << "materials.local.name                     " << projData->materials_local_name << endl << endl;
      // do not independently refine the ports
      // The refinement variables are included in case the user wants to use refinement when running OpenParEM2D manually.
      out << "refinement.frequency                     none" << endl;
      out << "refinement.variable                      |Zo|" << endl;
      out << "refinement.iteration.min                 " << projData->refinement_iteration_min << endl;
      out << "refinement.iteration.max                 " << projData->refinement_iteration_max << endl;
      out << "refinement.required.passes               " << projData->refinement_required_passes << endl;
      out << "refinement.tolerance                     " << projData->refinement_tolerance << endl << endl;

      out << setprecision(15) << "frequency.plan.point                     " << frequency << endl << endl;

      out << "solution.modes                           " << get_modeCount() << endl;
      out << "solution.temperature                     " << projData->solution_temperature << endl;
      out << "solution.tolerance                       " << projData->solution_2D_tolerance << endl;
      out << "solution.iteration.limit                 " << projData->solution_iteration_limit << endl;
      out << "solution.modes.buffer                    " << projData->solution_modes_buffer << endl;
      out << "solution.check.closed.loop               " << convertLogic(projData->solution_check_closed_loop) << endl;
      out << "solution.impedance.definition            " << impedance_definition.get_value() << endl;
      out << "solution.impedance.calculation           " << impedance_calculation.get_value() << endl;
      out << "solution.accurate.residual               " << convertLogic(projData->solution_accurate_residual) << endl;
      out << "solution.shift.invert                    " << convertLogic(projData->solution_shift_invert) << endl;
      out << "solution.use.initial.guess               " << convertLogic(projData->solution_use_initial_guess) << endl;
      if (gamma && gamma->get_alpha() > 0) {
         out << "solution.initial.alpha                   " << gamma->get_alpha() << endl;
      }
      if (gamma && gamma->get_beta() > 0) {
         // align with with OpenParEM2D
         out << "solution.initial.beta                    " << gamma->get_beta()*frequency/gamma->get_frequency()*2 << endl;
      }
      out << "solution.shift.factor                    " << projData->solution_shift_factor << endl << endl;

      out << "output.show.refining.mesh                " << convertLogic(projData->output_show_refining_mesh) << endl;
      out << "output.show.iterations                   " << convertLogic(projData->output_show_iterations) << endl;
      out << "output.show.license                      " << convertLogic(projData->output_show_license) << endl << endl;

      out << "test.create.cases                        " << convertLogic(projData->test_create_cases) << endl;
      out << "test.show.audit                          " << "false" << endl;
      out << "test.show.detailed.cases                 " << convertLogic(projData->test_show_detailed_cases) << endl << endl;

      out << "debug.show.memory                        " << convertLogic(projData->debug_show_memory) << endl;
      out << "debug.show.project                       " << convertLogic(projData->debug_show_project) << endl;
      out << "debug.show.frequency.plan                " << convertLogic(projData->debug_show_frequency_plan) << endl;
      out << "debug.show.materials                     " << convertLogic(projData->debug_show_materials) << endl;
      out << "debug.show.mode.definitions              " << convertLogic(projData->debug_show_port_definitions) << endl;
      out << "debug.show.impedance.details             " << convertLogic(projData->debug_show_impedance_details) << endl;
      out << "debug.skip.solve                         " << "false" << endl;
      out << "debug.tempfiles.keep                     " << "true" << endl << endl;

      out.close();
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3081: Failed to open file \"%s\" for writing.\n",filename.str().c_str());
   }
}

void Port::saveModeFile (struct projectData *projData, vector<Path *> *pathList, BoundaryDatabase *boundaryDatabase)
{
   int pathNumber=1;   
   double x1,y1,z1,x2,y2,z2;
   double xr1,yr1,zr1,xr2,yr2,zr2;

   if (pathIndexList.size() != 1) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Port::saveModeFile operation on a Port with an invalid path definition.\n");
   }

   if (!rotated) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Port::saveModeFile operation on a Port without a rotated path.\n");
   }

   Path *path=(*pathList)[pathIndexList[0]];

   if (path->get_points_size() == 0) return;

   // set up the save file
   stringstream modeFilename;
   modeFilename << boundaryDatabase->get_tempDirectory() << "/S" << get_name() << "/S" << get_name() << "_modes.txt";

   ofstream modeFile;
   modeFile.open(modeFilename.str().c_str(),ofstream::out);

   if (modeFile.is_open()) {

      modeFile << "#OpenParEMmodes 1.0" << endl << endl << "File" << endl << "   name=generated by OpenParEM3D" << endl << "EndFile" << endl << endl;

      // boundaries

      x1=path->get_point_x(0); y1=path->get_point_y(0); z1=path->get_point_z(0);

      long unsigned int i=0;
      while (i < path->get_points_size()) {
         if (!path->is_closed() && i == path->get_points_size()-1) break;

         if (i < path->get_points_size()-1) {x2=path->get_point_x(i+1); y2=path->get_point_y(i+1); z2=path->get_point_z(i+1);}
         else {x2=path->get_point_x(0); y2=path->get_point_y(0); z2=path->get_point_z(0);}

         Boundary* boundary=boundaryDatabase->get_matchBoundary (x1,y1,z1,x2,y2,z2);
         if (boundary) {

            xr1=x1; yr1=y1; zr1=z1;
            xr2=x2; yr2=y2; zr2=z2;

            rotated->rotatePoint(&xr1,&yr1,&zr1);
            rotated->rotatePoint(&xr2,&yr2,&zr2);

            modeFile << "Path" << endl;
            modeFile << "   name=path" << pathNumber << endl;
            modeFile << setprecision(16) << "   point=(" << xr1 << "," << yr1 << ")" << endl;
            modeFile << setprecision(16) << "   point=(" << xr2 << "," << yr2 << ")" << endl;
            modeFile << "   closed=false" << endl;
            modeFile << "EndPath" << endl << endl;

            modeFile << "Boundary" << endl;
            modeFile << "   name=boundary" << pathNumber << endl;
            if (boundary->is_perfect_electric_conductor()) modeFile << "   type=perfect_electric_conductor" << endl;
            if (boundary->is_perfect_magnetic_conductor()) modeFile << "   type=perfect_magnetic_conductor" << endl;
            if (boundary->is_surface_impedance()) {
               modeFile << "   type=surface_impedance" << endl;
               modeFile << "   material=" << boundary->get_material() << endl;
            }
            if (boundary->is_radiation()) {
               // A proper surface impedance boundary is not yet implemented in OpenParEM2D.
               // Assume PEC.
               // ToDo: update this when the surface impedance boundary is available.
               modeFile << "   type=perfect_electric_conductor" << endl;
            }
            modeFile << "   path=path" << pathNumber << endl;
            modeFile << "EndBoundary" << endl << endl;

            pathNumber++;
         }

         x1=x2; y1=y2; z1=z2;

         i++;
      }

      // modes
      i=0;
      while (i < modeList.size()) {
         modeList[i]->output(&modeFile,pathList,rotated,spin180degrees);
         i++;
      }
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3082: Failed to open file \"%s\" for writing.\n",modeFilename.str().c_str());
   }
}

bool Port::merge(vector<Path *> *pathList)
{
   Path *mergedPath=nullptr;

   // merge

   stringstream boundaryName;
   boundaryName << get_name();

   bool fail=mergePaths(pathList,&pathIndexList,&reverseList,"Port",boundaryName.str(),&mergedPath);
   if (fail) return fail;

   if (! mergedPath) return fail;

   stringstream pathName;
   pathName << "S" << get_name() << "_OpenParEM3D_generated";
   mergedPath->set_name(pathName.str());

   pathList->push_back(mergedPath);

   // update the tracking information

   pathNameList.clear();
   keywordPair *newName=new keywordPair();
   newName->set_value(pathName.str());
   pathNameList.push_back(newName);

   pathIndexList.clear();
   pathIndexList.push_back(pathList->size()-1);

   reverseList.clear();
   reverseList.push_back(false);

   return fail;
}

bool Port::createRotated (vector<Path *> *pathList, string indent)
{
   if (pathIndexList.size() != 1) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3180: Port at line %d is incorrectly formatted.\n",
                                    indent.c_str(),indent.c_str(),startLine);
      return true;
   }

   if (rotated != nullptr) delete rotated;
   rotated=(*pathList)[pathIndexList[0]]->rotateToXYplane();

   if (! rotated) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3179: Port at line %d does not form a closed polygon with nonzero area.\n",
                                    indent.c_str(),indent.c_str(),startLine);
      return true;
   }
   return false;
}

bool Port::is_point_inside (double x, double y, double z)
{
   if (! rotated) return false;  // not a closed boundary because it failed the rotation operation, so no point can be inside
   if (rotated->is_point_inside (x,y,z)) return true;
   return false;
}

bool Port::is_triangleInside (DenseMatrix *pointMat)
{
   if (is_point_inside (pointMat->Elem(0,0),pointMat->Elem(1,0),pointMat->Elem(2,0)) &&
       is_point_inside (pointMat->Elem(0,1),pointMat->Elem(1,1),pointMat->Elem(2,1)) &&
       is_point_inside (pointMat->Elem(0,2),pointMat->Elem(1,2),pointMat->Elem(2,2))) {
      return true;
   }
   return false;
}

bool Port::is_modePathInside (string *indent, vector<Path *> *pathList)
{
   bool fail=false;

   if (pathIndexList.size() != 1) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Port::is_modePathInside operation on a Port with an invalid path definition.\n");
      return false;
   }

   long unsigned i=0;
   while (i < modeList.size()) {
      if (! modeList[i]->is_enclosedByPath(pathList,rotated)) {
         if (modeList[i]->is_modal()) PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3083: Port %s Mode %d of type \"%s\" is not enclosed within the port.\n",
            indent->c_str(),indent->c_str(),get_name().c_str(),modeList[i]->get_Sport(),modeList[i]->get_type().c_str());
         if (modeList[i]->is_line()) PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3084: Port %s Line %d of type \"%s\" is not enclosed within the port.\n",
            indent->c_str(),indent->c_str(),get_name().c_str(),modeList[i]->get_Sport(),modeList[i]->get_type().c_str());
         fail=true;
      }
      i++;
   }

   if (fail) return false;
   return true;
}

void Port::create2Dmesh (int order, ParMesh *mesh3D, vector<ParSubMesh> *parSubMeshes, long unsigned int parSubMeshIndex, double tol)
{
   PetscMPIInt rank,size;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   MPI_Comm_size(PETSC_COMM_WORLD, &size);

   bool match;
   int foundNormal;
   int adjacentElemIndex,info;
   DenseMatrix pointMat2D(3,3),pointMat3D(3,3);
   Array<int> used3D(3),used2D(3);

   // find the normal
   foundNormal=0;
   int i=0;
   while (i < mesh3D->GetNBE()) {
      if (mesh3D->GetBdrAttribute(i) == attribute) {
         if (!foundNormal) {

            ElementTransformation *elemTr=mesh3D->GetBdrElementTransformation(i);
            elemTr->SetIntPoint(&Geometries.GetCenter(elemTr->GetGeometryType()));
            CalcOrtho(elemTr->Jacobian(), normal);

            double mag=sqrt(normal(0)*normal(0)+normal(1)*normal(1)+normal(2)*normal(2));
            normal(0)/=mag;
            normal(1)/=mag;
            normal(2)/=mag;

            foundNormal=1;
         }
      }
      i++;
   }

   // put the normal onto all ranks

   if (rank == 0) {
      double nx,ny,nz;

      i=1;
      while (i < size) {
         MPI_Recv(&nx,1,MPI_DOUBLE,i,0,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         MPI_Recv(&ny,1,MPI_DOUBLE,i,1,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         MPI_Recv(&nz,1,MPI_DOUBLE,i,2,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

         if (nx != 0 || ny != 0 || nz != 0) {
            normal(0)=nx;
            normal(1)=ny;
            normal(2)=nz;
            foundNormal=1;
         }

         i++;
      }

   } else {
      MPI_Send(&(normal(0)),1,MPI_DOUBLE,0,0,PETSC_COMM_WORLD);
      MPI_Send(&(normal(1)),1,MPI_DOUBLE,0,1,PETSC_COMM_WORLD);
      MPI_Send(&(normal(2)),1,MPI_DOUBLE,0,2,PETSC_COMM_WORLD);
   }

   MPI_Bcast(&(normal(0)),1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
   MPI_Bcast(&(normal(1)),1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
   MPI_Bcast(&(normal(2)),1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
   MPI_Bcast(&foundNormal,1,MPI_INT,0,PETSC_COMM_WORLD);

   if (!foundNormal) cout << "ASSERT: Port::create2Dmesh failed to find a normal." << endl;

   // reset the boundary attributes from indicating the port to indicating the material
   i=0;
   while (i < mesh3D->GetNBE()) {

      mesh3D->GetBdrPointMatrix(i,pointMat3D);

      // find the matching border element from the 2D mesh
      int j=0;
      while (j < (*parSubMeshes)[parSubMeshIndex].GetNE()) {

         (*parSubMeshes)[parSubMeshIndex].GetPointMatrix(j,pointMat2D);

         used3D=0; used2D=0;

         int m=0;
         while (m < 3) {
            int n=0;
            while (n < 3) {
               if (compare_xyz(pointMat3D.Elem(0,m),pointMat3D.Elem(1,m),pointMat3D.Elem(2,m),
                               pointMat2D.Elem(0,n),pointMat2D.Elem(1,n),pointMat2D.Elem(2,n),tol)) {
                  used3D[m]=1;
                  used2D[n]=1;
                  break;
               }
               n++;
            }
            m++;
         }

         match=true;
         m=0;
         while (m < 3) {
            if (used3D[m] == 0) {match=false; break;}
            if (used2D[m] == 0) {match=false; break;}
            m++;
         }

         if (match) {
            mesh3D->GetBdrElementAdjacentElement (i,adjacentElemIndex,info);
            (*parSubMeshes)[parSubMeshIndex].SetAttribute(j,mesh3D->GetAttribute(adjacentElemIndex));
            break;
         }
         j++;
      }
      i++;
   }

   // create finite element collections and spaces for this 2D port

   fec2D_ND=new ND_FECollection(order,(*parSubMeshes)[parSubMeshIndex].Dimension());
   fes2D_ND=new ParFiniteElementSpace(&((*parSubMeshes)[parSubMeshIndex]),fec2D_ND);

   fec2D_H1=new H1_FECollection(order,(*parSubMeshes)[parSubMeshIndex].Dimension());
   fes2D_H1=new ParFiniteElementSpace(&((*parSubMeshes)[parSubMeshIndex]),fec2D_H1);

   fec2D_L2=new L2_FECollection(order,(*parSubMeshes)[parSubMeshIndex].Dimension());
   fes2D_L2=new ParFiniteElementSpace(&((*parSubMeshes)[parSubMeshIndex]),fec2D_L2);

   // rotate the outward normal
   rotated_normal(0)=normal(0);
   rotated_normal(1)=normal(1);
   rotated_normal(2)=normal(2);

   rotated->rotatePoint(&rotated_normal(0),&rotated_normal(1),&rotated_normal(2));

   // spin 180 degrees about the y axis if needed to ensure that the outward normal points in the +z direction
   spin180degrees=false;
   if (rotated_normal(2) < 0) spin180degrees=true;

   set_filenames();
}

void Port::set_filenames()
{
   stringstream ssMeshFilename;
   ssMeshFilename << "S" << get_name() << "_mesh2D.msh";
   meshFilename=ssMeshFilename.str();

   stringstream ssModesFilename;
   ssModesFilename << "S" << get_name() << "_modes.txt";
   modesFilename=ssModesFilename.str();
}

void eh(MPI_Comm *comm, int *err, ...)
{
   PetscPrintf(PETSC_COMM_WORLD,"ERROR3085: Failed to spawn OpenParEM2D.  Manually remove the lock file.\n");
  
   // ToDo - Eliminate the MPI_Abort for a graceful exit.

   // Abort to avoid hanging ranks
   MPI_Abort (PETSC_COMM_WORLD,1);
}

bool Port::solve(string *directory, MPI_Comm *MPI_PORT_COMM)
{
   bool fail=false;
   PetscMPIInt size,rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   MPI_Comm_size(PETSC_COMM_WORLD, &size);

   // cd to the project directory
   stringstream projDirectory;
   projDirectory << *directory << "/S" << get_name();
   std::filesystem::current_path(projDirectory.str().c_str());

   // lock file name
   stringstream ssLock;
   ssLock << "." << "S" << get_name() << ".lock";

   // argv 
   char *argv[2];

   char project[64];
   sprintf(project,"S%s.proj",get_name().c_str());

   argv[0]=project;
   argv[1]=nullptr;

   // launch the jobs from rank 0
   MPI_Errhandler errorHandler;
   MPI_Comm_create_errhandler(eh,&errorHandler);
   MPI_Comm_set_errhandler(PETSC_COMM_WORLD,errorHandler);

   MPI_Comm_spawn ("OpenParEM2D",argv,size,MPI_INFO_NULL,0,PETSC_COMM_WORLD,MPI_PORT_COMM,MPI_ERRCODES_IGNORE);

   if (rank == 0) {

      // wait for the lock file to appear
      bool launch_success=true;
      int watchdog=0;
      while (!std::filesystem::exists(ssLock.str().c_str())) {
         sleep(0.001);
         watchdog+=1;
         if (watchdog == 100000) { // ~100 seconds
            launch_success=false;
            PetscPrintf(PETSC_COMM_WORLD,"ERROR3086: OpenParEM2D failed to launch for \"%s\"\n",project);
            break;
         }
      }

      if (launch_success) {
         // wait for the lock file to disappear - can hang
         while (std::filesystem::exists(ssLock.str().c_str())) {sleep(0.01);}
      }

      // send breakout tags
      int message=1;
      int i=1;
      while (i < size) {
         MPI_Send(&message,1,MPI_INT,i,10,PETSC_COMM_WORLD);
         i++;
      }

   } else {

      // wait for breakout tag
      int flag=0;
      while (!flag) {
         sleep(0.01);
         MPI_Iprobe(0,10,PETSC_COMM_WORLD,&flag,MPI_STATUS_IGNORE);
      }

      // get the tag to clear it
      int retVal;
      MPI_Recv(&retVal,1,MPI_INT,0,10,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
   }

   MPI_Barrier(*MPI_PORT_COMM);

   MPI_Comm_set_errhandler(PETSC_COMM_WORLD,MPI_ERRORS_RETURN);
   MPI_Errhandler_free(&errorHandler);

   // get the return values
   if (rank == 0) {
      int i=0;
      while (i < size) {
         int retVal;
         MPI_Recv(&retVal,1,MPI_INT,i,0,*MPI_PORT_COMM,MPI_STATUS_IGNORE);
         if (retVal > 0) fail=true;
         i++;
      }
      if (fail) PetscPrintf(PETSC_COMM_WORLD,"ERROR3087: OpenParEM2D error in execution for Port \"%s\".\n",get_name().c_str());
   }

   std::filesystem::current_path("../../");

   return fail;
}

bool Port::uses_current()
{
   if (impedance_definition.get_value().compare("VI") == 0) return true;
   if (impedance_definition.get_value().compare("PI") == 0) return true;
   return false;
}

bool Port::uses_voltage()
{
   if (impedance_definition.get_value().compare("VI") == 0) return true;
   if (impedance_definition.get_value().compare("PV") == 0) return true;
   return false;
}

void Port::markUsedModes ()
{
   long unsigned int i=0;
   while (i < modeList.size()) {
      if ((uses_current() && modeList[i]->is_current()) ||
          (uses_voltage() && modeList[i]->is_voltage())) {
         modeList[i]->set_isUsed();
      }
      i++;
   }
}

bool Port::loadSizes_tz (string *directory)
{
   char filename[128];

   PetscMPIInt size;
   MPI_Comm_size(PETSC_COMM_WORLD, &size);

   // cd to the project directory

   stringstream projDirectory;
   projDirectory << *directory << "/S" << get_name();

   try {
      std::filesystem::current_path(projDirectory.str().c_str());
   } catch (std::filesystem::filesystem_error const& ex) {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3088: Missing project directory for the 2D solution of Port %s.\n",get_name().c_str());
      return true;
   }

   // cd to the temp directory

   stringstream tempDirectory;
   tempDirectory << "temp_S" << get_name();

   try {
      std::filesystem::current_path(tempDirectory.str().c_str());
   } catch (std::filesystem::filesystem_error const& ex) {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3089: Missing temporary directory for the 2D solution of Port %s.\n",get_name().c_str());
      std::filesystem::current_path("../../");
      return true;
   }

   // t_size

   sprintf(filename,"St_mat.S%s.%05d",get_name().c_str(),size-1);
   ifstream ss_size_t;
   ss_size_t.open(filename,ifstream::in);
   if (ss_size_t.is_open()) {
      bool fail=false;
      string line;
      if (getline(ss_size_t,line)) {
         bool found=false;
         stringstream sstream(line);
         string entry;
         int count=0;
         while (std::getline(sstream,entry,' ')) {
            if (count == 3 && is_int(&entry)) {
               found=true;
               t_size=stoi(entry)+1;
            }
            count++;
         }
         if (! found) {
            PetscPrintf(PETSC_COMM_WORLD,"ERROR3090: Failed to parse data in file \"%s\".\n",filename);
            fail=true;
         }
      } else {
          PetscPrintf(PETSC_COMM_WORLD,"ERROR3091: Failed to read data in file \"%s\".\n",filename);
          fail=true;
      }
      ss_size_t.close();
      if (fail) return true;
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3092: Unable to open file \"%s\" for reading.\n",filename);
      std::filesystem::current_path("../../../");
      return true;
   }

   // z_size

   sprintf(filename,"Tz_eps_re_mat.S%s.%05d",get_name().c_str(),size-1);
   ifstream ss_size_z;
   ss_size_z.open(filename,ifstream::in);
   if (ss_size_z.is_open()) {
      bool fail=false;
      string line;
      if (getline(ss_size_z,line)) {
         bool found=false;
         stringstream sstream(line);
         string entry;
         int count=0;
         while (std::getline(sstream,entry,' ')) {
            if (count == 3 && is_int(&entry)) {
               found=true;
               z_size=stoi(entry)+1;
            }
            count++;
         }
         if (! found) {
            PetscPrintf(PETSC_COMM_WORLD,"ERROR3093: Failed to parse data in file \"%s\".\n",filename);
            fail=true;
         }
      } else {   
          PetscPrintf(PETSC_COMM_WORLD,"ERROR3094: Failed to read data in file \"%s\".\n",filename);
          fail=true;
      }
      ss_size_z.close();
      if (fail) return true;
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3095: Unable to open file \"%s\" for reading.\n",filename);
      std::filesystem::current_path("../../../");
      return true;
   }

   // cd back to the 3D project directory
   std::filesystem::current_path("../../../");

   return false;
}

bool Port::set_alphaBeta (double alpha, double beta, int mode)
{
   int assignCount=0;
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed() && modeList[i]->get_modeNumber2D() == mode) {
         modeList[i]->set_alpha(alpha);
         modeList[i]->set_beta(beta);
         assignCount++;
      }
      i++;
   }

   if (assignCount == 0) {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3096: Failed to assign propagation constant for Port \"%s\".\n",get_name().c_str());
      return true;
   }

   if (assignCount > 1) {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3097: Assigned propagation constant to more than one mode for Port \"%s\".\n",get_name().c_str());
      return true;
   }

   return false;
}

bool Port::load_gammaZ (string *directory, double frequency)
{
   char filename[128];
   double NpTodB=20*log10(exp(1));
   double eps0=8.8541878176e-12;
   double ko=2*M_PI*frequency*sqrt(4.0e-7*M_PI*eps0);
   int entry=0;

   PetscMPIInt size;
   MPI_Comm_size(PETSC_COMM_WORLD, &size);

   // cd to the project directory

   stringstream projDirectory;
   projDirectory << *directory << "/S" << get_name();

   try {
      std::filesystem::current_path(projDirectory.str().c_str());
   } catch (std::filesystem::filesystem_error const& ex) {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3098: Missing project directory for the 2D solution of Port %s.\n",get_name().c_str());
      return true;
   }

   // load the data for the given frequency

   vector<string> tokenList;
   sprintf(filename,"S%s_results.csv",get_name().c_str());
   ifstream CSV;
   CSV.open(filename,ifstream::in);
   if (CSV.is_open()) {
      string line;
      while (getline(CSV,line)) {
         stringstream ssLine(line);
         string value;
         bool load=false;
         bool is_first=true;
         while (std::getline(ssLine,value,',')) {
            if (load) tokenList.push_back(value);
            if (is_first) {
               if (is_double(&value) && double_compare(stod(value),frequency,1e-12)) load=true;
               if (! load) break;
               is_first=false;
            }
         }

         ssLine.str("");
         ssLine.clear();
         if (load) break;
      }
      CSV.close();
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3099: Unable to open file \"%s\" for reading.\n",filename);
      std::filesystem::current_path("../../");
      return true;
   }

   // pull out the needed data
   if (tokenList.size() > 0) {
      bool fail=false;
      int modal_impedance_calculation=0;
      int mode_count=0;
      double alpha=0,beta=0;

      // impedance type and mode count

      entry=6;

      if (tokenList.size() > 8) {
         if (is_int(&tokenList[entry])) modal_impedance_calculation=stoi(tokenList[entry]); else fail=true;
         entry++;
         if (is_int(&tokenList[entry])) mode_count=stoi(tokenList[entry]); else fail=true;
         entry++;
      } else fail=true;

      if (mode_count != get_modeCount()) {
         PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Port::load_gammaZ mismatched mode counts for Port \"%s\".\n",get_name().c_str());
         fail=true;
      }

      if ((is_modal() && modal_impedance_calculation == 0) ||
          (!is_modal() && modal_impedance_calculation == 1)) {
         PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Port::load_gammaZ mismatched impedance calculation for Port \"%s\".\n",get_name().c_str());
         fail=true;
      }

      // alpha and beta

      if ((int)tokenList.size() >= entry+mode_count*2) {
         int i=0;
         while (i < mode_count) {
            if (is_double(&tokenList[entry])) alpha=stod(tokenList[entry]); else fail=true;   // dB/m
            alpha/=NpTodB;   // MKS units
            entry++;

            if (is_double(&tokenList[entry])) beta=stod(tokenList[entry]); else fail=true;  // beta/ko
            beta*=ko;        // MKS units
            entry++;

            if (!fail) set_alphaBeta(alpha,beta,i+1);

            i++;
         }

      } else fail=true;

      // impedance count
      int impedanceEntries=0;
      if ((int)tokenList.size() >= entry && is_int(&tokenList[entry])) impedanceEntries=stoi(tokenList[entry]); else fail=true;
      entry++;

      if (modal_impedance_calculation && impedanceEntries != mode_count) {
         PetscPrintf (PETSC_COMM_WORLD,"ERROR3100: Port impedance results are inconsistent with the mode count.\n");
         fail=true;
      }

      if (!modal_impedance_calculation && impedanceEntries != mode_count*mode_count) {
         PetscPrintf (PETSC_COMM_WORLD,"ERROR3101: Port impedance results are inconsistent with the mode count.\n");
         fail=true;
      }

      if (fail) return fail;

      // impedance

      if ((int)tokenList.size() >= entry+impedanceEntries) {

         ReZ2D.SetSize(mode_count,mode_count);
         ReZ2D=0.;

         ImZ2D.SetSize(mode_count,mode_count);
         ImZ2D=0.;

         if (modal_impedance_calculation) {
            int i=0;
            while (i < mode_count) {
               if (is_double(&tokenList[entry])) ReZ2D(i,i)=stod(tokenList[entry]); else fail=true;
               entry++;

               if (is_double(&tokenList[entry])) ImZ2D(i,i)=stod(tokenList[entry]); else fail=true;
               entry++;

               i++;
            }
         } else {
            int i=0;
            while (i < mode_count) {
               int j=0;
               while (j < mode_count) {
                  if (is_double(&tokenList[entry])) ReZ2D(i,j)=stod(tokenList[entry]); else fail=true;
                  entry++;

                  if (is_double(&tokenList[entry])) ImZ2D(i,j)=stod(tokenList[entry]); else fail=true;
                  entry++;

                  j++;
               }
               i++;
            }
         }

      } else fail=true;

      // voltage count
      int voltageEntries=0;
      if ((int)tokenList.size() >= entry && is_int(&tokenList[entry])) voltageEntries=stoi(tokenList[entry]); else fail=true;
      entry++;

      if (modal_impedance_calculation && voltageEntries != mode_count) {
         PetscPrintf (PETSC_COMM_WORLD,"ERROR3102: Port voltage results are inconsistent with the mode count.\n");
         fail=true;
      }

      if (!modal_impedance_calculation && voltageEntries != mode_count*mode_count) {
         PetscPrintf (PETSC_COMM_WORLD,"ERROR3103: Port voltage results are inconsistent with the mode count.\n");
         fail=true;
      }

      if (fail) return fail;

      // voltage

      if ((int)tokenList.size() >= entry+voltageEntries) {

         ReV2D.SetSize(mode_count,mode_count);
         ReV2D=0.;

         ImV2D.SetSize(mode_count,mode_count);
         ImV2D=0.;

         if (modal_impedance_calculation) {
            int i=0;
            while (i < mode_count) {
               if (is_double(&tokenList[entry])) ReV2D(i,i)=stod(tokenList[entry]); else fail=true;
               entry++;

               if (is_double(&tokenList[entry])) ImV2D(i,i)=stod(tokenList[entry]); else fail=true;
               entry++;

               i++;
            }
         } else {
            int i=0;
            while (i < mode_count) {
               int j=0;
               while (j < mode_count) {
                  if (is_double(&tokenList[entry])) ReV2D(i,j)=stod(tokenList[entry]); else fail=true;
                  entry++;

                  if (is_double(&tokenList[entry])) ImV2D(i,j)=stod(tokenList[entry]); else fail=true;
                  entry++;

                  j++;
               }
               i++;
            }
         }

      } else fail=true;

      if (fail) {
         PetscPrintf(PETSC_COMM_WORLD,"ERROR3104: Failed to extract computed parameters for Port \"%s\".\n",get_name().c_str());
         std::filesystem::current_path("../../");
         return true;
      }
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3105: Missing computed results for Port \"%s\".\n",get_name().c_str());
   }

   // cd back to the 3D project directory
   std::filesystem::current_path("../../");

   return false;
}

bool Port::loadSolution (string *directory, double frequency)
{
   PetscMPIInt rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   loadSizes_tz(directory);
   load_gammaZ(directory,frequency);

   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->loadSolution(directory,get_name(),t_size,z_size)) return true;
      if (modeList[i]->scaleSolution()) return true;
      int index=modeList[i]->get_modeNumber2D()-1;
      modeList[i]->setSign(&(ReV2D(index,index)),&(ImV2D(index,index)));
      i++;
   }

   return false;
}

// for the first mode found
bool Port::getGamma (Gamma *gamma, double frequency)
{
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
         gamma->set(modeList[i]->get_alpha(),modeList[i]->get_beta(),frequency);
         return false;
      }
      i++;
   }
   return true;
}

void Port::build2Dgrids()
{
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
         modeList[i]->build2Dgrids(fes2D_ND,fes2D_H1);
      }
      i++;
   }
}

void Port::build3Dgrids(ParFiniteElementSpace *fes3D_ND, ParFiniteElementSpace *fes3D_H1)
{
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
         modeList[i]->build3Dgrids(fes3D_ND,fes3D_H1);
      }
      i++;
   }
}

void Port::build2DModalGrids()
{
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
         modeList[i]->build2DModalGrids(fes2D_ND,fes2D_H1);
      }
      i++;
   }
}

void Port::build2DSolutionGrids ()
{
   grid2DsolutionReEt=new ParGridFunction(fes2D_ND);
   *grid2DsolutionReEt=0.0;

   grid2DsolutionImEt=new ParGridFunction(fes2D_ND);
   *grid2DsolutionImEt=0.0;

   grid2DsolutionReEz=new ParGridFunction(fes2D_L2);
   *grid2DsolutionReEz=0.0;

   grid2DsolutionImEz=new ParGridFunction(fes2D_L2);
   *grid2DsolutionImEz=0.0;

   grid2DsolutionReHt=new ParGridFunction(fes2D_ND);
   *grid2DsolutionReHt=0.0;

   grid2DsolutionImHt=new ParGridFunction(fes2D_ND);
   *grid2DsolutionImHt=0.0;

   grid2DsolutionReHz=new ParGridFunction(fes2D_L2);
   *grid2DsolutionReHz=0.0;

   grid2DsolutionImHz=new ParGridFunction(fes2D_L2);
   *grid2DsolutionImHz=0.0;
}

int Port::get_modeCount()
{
   int count=0;
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) count++;
      i++;
   }
   return count;
}

int Port::get_SportCount()
{
   int count=0;
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
         if (modeList[i]->get_Sport() > count) count=modeList[i]->get_Sport();
      }
      i++;
   }
   return count;
}

void Port::printSolution (string indent)
{
   PetscMPIInt rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   if (rank == 0) {
      cout << indent << "Port " << get_name() << endl;
      cout << indent << indent << "t_size=" << t_size << endl;
      cout << indent << indent << "z_size=" << z_size << endl;

      if (is_modal()) {
         int i=0;
         while (i < get_modeCount()) {
            cout << indent << indent << "Z[" << i+1 << "]=(" << ReZ2D(i,i) << "," << ImZ2D(i,i) << ")" << endl;
            i++;
         }
      } else {
         int i=0;
         while (i < get_modeCount()) {
            cout << indent << indent;
            int j=0;
            while (j < get_modeCount()) {
               cout << "Z[" << i+1 << "," << j+1 << "]=(" << ReZ2D(i,j) << "," << ImZ2D(i,j) << "), ";
               j++;
            }
            cout << endl;
            i++;
         }
      }

      long unsigned int i=0;
      while (i < modeList.size()) {
         modeList[i]->printSolution();
         i++;
      }
   }
}

void Port::build_essTdofList (ParFiniteElementSpace *fespace, ParMesh *pmesh)
{
   Array<int> border_attributes;
   border_attributes.SetSize(pmesh->bdr_attributes.Max());
   border_attributes=0;               // default to nothing
   border_attributes[attribute-1]=1;  // enable this port

   if (ess_tdof_list) delete ess_tdof_list;
   ess_tdof_list=new Array<int>;
   fespace->GetEssentialTrueDofs(border_attributes, *ess_tdof_list);
   offset=fespace->GetTrueDofOffsets();
}

void Port::fillX (Vec *X, Vec *Xdofs, int drivingSport)
{
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
         if (modeList[i]->get_Sport() == drivingSport) {
            modeList[i]->fillX(X,Xdofs,ess_tdof_list,offset);
            break;
         }
      }
      i++;
   }
}

// ToDo: include Inv_mur
bool Port::addMassPortIntegrators (ParMesh *pmesh, ParMixedBilinearForm *pmblf, PWConstCoefficient *Inv_mur, vector<Array<int> *> &borderAttributesList,
                                   vector<ConstantCoefficient *> &alphaConstList, vector<ConstantCoefficient *> &betaConstList, bool isReal)
{
   bool fail=false;

   bool found=false;
   double alpha=0;
   double beta=0;
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
         alpha=modeList[i]->get_alpha();
         beta=modeList[i]->get_beta();
         found=true;
         break;
      }
      i++;
   }
   if (!found) {
      fail=true;
      cout << "ASSERT: Port::addMassPortIntegrators failed to find a mode to use." << endl;
   }

   Array<int> *border_attributes=new Array<int>;
   border_attributes->SetSize(pmesh->bdr_attributes.Max());
   (*border_attributes)=0;                // default to nothing
   (*border_attributes)[attribute-1]=1;   // enable this port

   ConstantCoefficient *alphaConst=new ConstantCoefficient(alpha);
   ConstantCoefficient *betaConst=new ConstantCoefficient(beta);

   if (isReal) pmblf->AddBoundaryIntegrator(new VectorFEMassIntegrator(*alphaConst),*border_attributes);
   else pmblf->AddBoundaryIntegrator(new VectorFEMassIntegrator(*betaConst),*border_attributes);

   // save for later deleting
   borderAttributesList.push_back(border_attributes);
   alphaConstList.push_back(alphaConst);
   betaConstList.push_back(betaConst);

   return fail;
}

void Port::flipModalHsign(long unsigned int mode)
{
   long unsigned int j=0;
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
         if (j == mode) {
            modeList[i]->flipModalHsign(!spin180degrees);
            break;
         }
         j++;
      }
      i++;
   }
}

void Port::extract2Dmesh (ParMesh *pmesh, vector<ParSubMesh> *parSubMeshes)
{
   Array<int> border_attributes(1);
   border_attributes[0]=attribute;

   ParSubMesh pmesh2D=ParSubMesh::CreateFromBoundary(*pmesh,border_attributes);
   parSubMeshes->push_back(pmesh2D);
}

void Port::set_V3Dsize (int size)
{
   ReV3D.SetSize(size,size);
   ReV3D=0.;

   ImV3D.SetSize(size,size);
   ImV3D=0.;
}

void Port::calculateWeights (Mode *drivingMode)
{
   int i,j;

   i=drivingMode->get_Sport()-1;

   long unsigned int k=0;
   while (k < modeList.size()) {
      if (modeList[k]->get_isUsed()) {

         // 2D
         int index2D=modeList[k]->get_modeNumber2D();

         // 3D
         j=modeList[k]->get_Sport()-1;
         modeList[k]->calculateWeights(fes2D_L2,grid2DsolutionReEt,grid2DsolutionImEt,grid2DsolutionReEz,grid2DsolutionImEz,grid2DsolutionReHt,
                                       grid2DsolutionImHt,grid2DsolutionReHz,grid2DsolutionImHz,normal,
                                       drivingMode->get_Sport(),
                                       complex<double>(ReV2D(index2D-1,index2D-1),ImV2D(index2D-1,index2D-1)),
                                       complex<double>(ReZ2D(index2D-1,index2D-1),ImZ2D(index2D-1,index2D-1)),
                                       complex<double>(ReV3D(i,j),ImV3D(i,j)));
      }
      k++;
   }
}

double Port::get_maxReflection ()
{
   double reflection;
   double maxReflection=-DBL_MAX;

   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
         reflection=modeList[i]->get_maxReflection();
         if (reflection > maxReflection) maxReflection=reflection;
      }
      i++;
   }

   return maxReflection;
}

double Port::get_maxReflection (int drivingSport)
{
   double reflection;
   double maxReflection=-DBL_MAX;

   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
         reflection=modeList[i]->get_maxReflection(drivingSport);
         if (reflection > maxReflection) maxReflection=reflection;
      }
      i++;
   }

   return maxReflection;
}

void Port::printPortReflections ()
{
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
         modeList[i]->printPortReflections();
      }
      i++;
   }
}

Mode* Port::getDrivingMode (int drivingSport)
{
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed() && modeList[i]->get_Sport() == drivingSport) {
         return modeList[i];
      }
      i++;
   }
   return nullptr;
}

void Port::fillIntegrationPoints (vector<Path *> *pathList)
{
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->is_voltage()) {
         modeList[i]->fillIntegrationPoints(pathList);
      }
      i++;
   }
}

void Port::calculateVoltages (ParMesh *pmesh, vector<Path *> *pathList, ParGridFunction *grid_re, ParGridFunction *grid_im, Mode* drivingMode)
{
   complex<double> voltage;

   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->is_voltage()) {
         voltage=-modeList[i]->calculateLineIntegral (pmesh,pathList,grid_re,grid_im);
         ReV3D(drivingMode->get_Sport()-1,modeList[i]->get_Sport()-1)=real(voltage);
         ImV3D(drivingMode->get_Sport()-1,modeList[i]->get_Sport()-1)=imag(voltage);
      }
      i++;
   }
}

void Port::calculateSparameter (Result *result, Mode *drivingMode)
{
   complex<double> S;

   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
         if (modeList[i]->get_Sport() == drivingMode->get_Sport()) {
            int index2D=modeList[i]->get_modeNumber2D();
            complex<double> Zo=complex<double>(ReZ2D(index2D-1,index2D-1),ImZ2D(index2D-1,index2D-1));
            result->set_Zo(Zo);
         }
         S=modeList[i]->calculateSparameter(drivingMode);
         result->push_S(S,modeList[i]->get_Sport(),drivingMode->get_Sport());
      }
      i++;
   }
}

void Port::transfer_2Dsolution_2Dgrids_to_3Dgrids()
{
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
         modeList[i]->transfer_2Dsolution_2Dgrids_to_3Dgrids();
      }
      i++;    
   }
}

void Port::transfer_2Dsolution_3Dgrids_to_2Dgrids()
{
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
         modeList[i]->transfer_2Dsolution_3Dgrids_to_2Dgrids();
      }
      i++;
   }
}

void Port::transfer_3Dsolution_3Dgrids_to_2Dgrids(fem3D *fem)
{
   // E

   ParTransferMap *full_to_port_ReEt=new ParTransferMap(*(fem->get_gridReE()),*grid2DsolutionReEt);
   full_to_port_ReEt->Transfer(*(fem->get_gridReE()),*grid2DsolutionReEt);
   delete full_to_port_ReEt; full_to_port_ReEt=nullptr;

   ParTransferMap *full_to_port_ImEt=new ParTransferMap(*(fem->get_gridImE()),*grid2DsolutionImEt);
   full_to_port_ImEt->Transfer(*(fem->get_gridImE()),*grid2DsolutionImEt);
   delete full_to_port_ImEt; full_to_port_ImEt=nullptr;

//   ParTransferMap *full_to_port_ReEz=new ParTransferMap(*(fem->get_gridReEz()),*grid2DsolutionReEz);
//   full_to_port_ReEz->Transfer(*(fem->get_gridReEz()),*grid2DsolutionReEz);
//   delete full_to_port_ReEz; full_to_port_ReEz=nullptr;

//   ParTransferMap *full_to_port_ImEz=new ParTransferMap(*(fem->get_gridImEz()),*grid2DsolutionImEz);
//   full_to_port_ImEz->Transfer(*(fem->get_gridImEz()),*grid2DsolutionImEz);
//   delete full_to_port_ImEz; full_to_port_ImEz=nullptr;

   // H

   ParTransferMap *full_to_port_ReHt=new ParTransferMap(*(fem->get_gridReH()),*grid2DsolutionReHt);
   full_to_port_ReHt->Transfer(*(fem->get_gridReH()),*grid2DsolutionReHt);
   delete full_to_port_ReHt; full_to_port_ReHt=nullptr;

   ParTransferMap *full_to_port_ImHt=new ParTransferMap(*(fem->get_gridImH()),*grid2DsolutionImHt);
   full_to_port_ImHt->Transfer(*(fem->get_gridImH()),*grid2DsolutionImHt);
   delete full_to_port_ImHt; full_to_port_ImHt=nullptr;

//   ParTransferMap *full_to_port_ReHz=new ParTransferMap(*(fem->get_gridReHz()),*grid2DsolutionReHz);
//   full_to_port_ReHz->Transfer(*(fem->get_gridReHz()),*grid2DsolutionReHz);
//   delete full_to_port_ReHz; full_to_port_ReHz=nullptr;

//   ParTransferMap *full_to_port_ImHz=new ParTransferMap(*(fem->get_gridImHz()),*grid2DsolutionImHz);
//   full_to_port_ImHz->Transfer(*(fem->get_gridImHz()),*grid2DsolutionImHz);
//   delete full_to_port_ImHz; full_to_port_ImHz=nullptr;
}

void Port::save2DParaView(ParSubMesh *psubmesh2D, struct projectData *projData, double frequency, bool add_extension)
{
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
         modeList[i]->save2DParaView(psubmesh2D,projData,frequency,add_extension);
      }
      i++;
   }
}

void Port::save3DParaView(ParMesh *pmesh, struct projectData *projData, double frequency, bool add_extension)
{
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
         modeList[i]->save3DParaView(pmesh,projData,frequency,add_extension);
      }
      i++;
   }
}

void Port::save2DSolutionParaView(ParSubMesh *psubmesh2D, struct projectData *projData, double frequency, int drivingSport, bool add_extension)
{
   if (!projData->debug_save_port_fields) return;

   stringstream ss;
   ss << projData->project_name << "_frequency_" << frequency << "_port_" << get_name() << "_driving_mode_" << drivingSport;

   stringstream ssParaView;
   ssParaView << "ParaView_solution_2D_" << projData->project_name;
   if (add_extension) ssParaView << "_test";

   ParaViewDataCollection *pd=nullptr;
   pd=new ParaViewDataCollection(ss.str(),psubmesh2D);
   pd->SetOwnData(false);
   pd->SetPrefixPath(ssParaView.str());
   pd->RegisterField("grid2DsolutionReEt",grid2DsolutionReEt);
   pd->RegisterField("grid2DsolutionImEt",grid2DsolutionImEt);
   pd->RegisterField("grid2DsolutionReEz",grid2DsolutionReEz);
   pd->RegisterField("grid2DsolutionImEz",grid2DsolutionImEz);
   pd->RegisterField("grid2DsolutionReHt",grid2DsolutionReHt);
   pd->RegisterField("grid2DsolutionImHt",grid2DsolutionImHt);
   pd->RegisterField("grid2DsolutionReHz",grid2DsolutionReHz);
   pd->RegisterField("grid2DsolutionImHz",grid2DsolutionImHz);
   pd->SetLevelsOfDetail(3);
   pd->SetDataFormat(VTKFormat::ASCII);
   pd->SetHighOrderOutput(true);
   pd->Save();
   delete pd;
}

void Port::save2DModalParaView(ParSubMesh *psubmesh2D, struct projectData *projData, double frequency, bool add_extension)
{
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
         modeList[i]->save2DModalParaView(psubmesh2D,projData,frequency,add_extension);
      }
      i++;
   }
}

void Port::resetElementNumbers ()
{
   long unsigned int i=0;
   while (i < modeList.size()) {
      if (modeList[i]->get_isUsed()) {
          modeList[i]->resetElementNumbers();
      }
      i++;
   }
}

void Port::reset()
{
   long unsigned int i=0;
   while (i < modeList.size()) {
      modeList[i]->reset();
      i++;
   }

   if (fec2D_ND) {delete fec2D_ND; fec2D_ND=nullptr;}
   if (fes2D_ND) {delete fes2D_ND; fes2D_ND=nullptr;}

   if (fec2D_H1) {delete fec2D_H1; fec2D_H1=nullptr;}
   if (fes2D_H1) {delete fes2D_H1; fes2D_H1=nullptr;}

   if (fec2D_L2) {delete fec2D_L2; fec2D_L2=nullptr;}
   if (fes2D_L2) {delete fes2D_L2; fes2D_L2=nullptr;}

   meshFilename="";
   modesFilename="";

   // do not delete rotated

   if (ess_tdof_list) {delete ess_tdof_list; ess_tdof_list=nullptr;}

   if (grid2DsolutionReEt) {delete grid2DsolutionReEt; grid2DsolutionReEt=nullptr;}
   if (grid2DsolutionImEt) {delete grid2DsolutionImEt; grid2DsolutionImEt=nullptr;}
   if (grid2DsolutionReEz) {delete grid2DsolutionReEz; grid2DsolutionReEz=nullptr;}
   if (grid2DsolutionImEz) {delete grid2DsolutionImEz; grid2DsolutionImEz=nullptr;}
   if (grid2DsolutionReHt) {delete grid2DsolutionReHt; grid2DsolutionReHt=nullptr;}
   if (grid2DsolutionImHt) {delete grid2DsolutionImHt; grid2DsolutionImHt=nullptr;}
   if (grid2DsolutionReHz) {delete grid2DsolutionReHz; grid2DsolutionReHz=nullptr;}
   if (grid2DsolutionImHz) {delete grid2DsolutionImHz; grid2DsolutionImHz=nullptr;}
}

Port::~Port()
{
   long unsigned int i=0;
   while (i < pathNameList.size()) {
      if (pathNameList[i]) delete pathNameList[i];
      i++;
   }

   i=0;
   while (i < modeList.size()) {
      if (modeList[i]) {
         delete modeList[i];
         modeList[i]=nullptr;
      }
      i++;
   }

   if (rotated) {delete rotated;}

   if (fec2D_ND) {delete fec2D_ND; fec2D_ND=nullptr;}
   if (fes2D_ND) {delete fes2D_ND; fes2D_ND=nullptr;}

   if (fec2D_H1) {delete fec2D_H1; fec2D_H1=nullptr;}
   if (fes2D_H1) {delete fes2D_H1; fes2D_H1=nullptr;}

   if (fec2D_L2) {delete fec2D_L2; fec2D_L2=nullptr;}
   if (fes2D_L2) {delete fes2D_L2; fes2D_L2=nullptr;}

   if (ess_tdof_list) {delete ess_tdof_list; ess_tdof_list=nullptr;}
}

///////////////////////////////////////////////////////////////////////////////////////////
// BoundaryDatabase
///////////////////////////////////////////////////////////////////////////////////////////

bool BoundaryDatabase::findSourceFileBlocks()
{
   bool fail=false;
   int start_lineNumber=inputs.get_first_lineNumber();
   int stop_lineNumber=inputs.get_last_lineNumber();
   int block_start,block_stop;

   // skip the version line, which must be the first line
   start_lineNumber=inputs.get_next_lineNumber(start_lineNumber);

   while (start_lineNumber < stop_lineNumber) {
      if (inputs.findBlock(start_lineNumber,stop_lineNumber, &block_start, &block_stop,
                                "File", "EndFile", false)) {
         fail=true;
      } else {
         if (block_start >= 0 && block_stop >= 0) {
            SourceFile *newSourceFile=new SourceFile(block_start,block_stop);
            sourceFileList.push_back(newSourceFile);
         }
      }
      start_lineNumber=inputs.get_next_lineNumber(block_stop);
   }
   return fail;
}

bool BoundaryDatabase::findPathBlocks()
{
   bool fail=false;
   int start_lineNumber=inputs.get_first_lineNumber();
   int stop_lineNumber=inputs.get_last_lineNumber();
   int block_start,block_stop;

   // skip the version line, which must be the first line
   start_lineNumber=inputs.get_next_lineNumber(start_lineNumber);

   while (start_lineNumber < stop_lineNumber) {
      if (inputs.findBlock(start_lineNumber,stop_lineNumber, &block_start, &block_stop,
                                "Path", "EndPath", false)) {
         fail=true;
      } else {
         if (block_start >= 0 && block_stop >= 0) {
            Path *newPath=new Path(block_start,block_stop);
            pathList.push_back(newPath);
         }
      }
      start_lineNumber=inputs.get_next_lineNumber(block_stop);
   }
   return fail;
}

bool BoundaryDatabase::findBoundaryBlocks()
{
   bool fail=false;
   int start_lineNumber=inputs.get_first_lineNumber();
   int stop_lineNumber=inputs.get_last_lineNumber();
   int block_start,block_stop;

   // skip the version line, which must be the first line
   start_lineNumber=inputs.get_next_lineNumber(start_lineNumber);

   while (start_lineNumber < stop_lineNumber) {
      if (inputs.findBlock(start_lineNumber,stop_lineNumber, &block_start, &block_stop,
                                "Boundary", "EndBoundary", false)) {
         fail=true;
      } else {
         if (block_start >= 0 && block_stop >= 0) {
            Boundary *newBoundary=new Boundary(block_start,block_stop);
            boundaryList.push_back(newBoundary);
         }
      }
      start_lineNumber=inputs.get_next_lineNumber(block_stop);
   }
   return fail;
}

bool BoundaryDatabase::findPortBlocks()
{
   bool fail=false;
   int start_lineNumber=inputs.get_first_lineNumber();
   int stop_lineNumber=inputs.get_last_lineNumber();
   int block_start,block_stop;

   // skip the version line, which must be the first line
   start_lineNumber=inputs.get_next_lineNumber(start_lineNumber);

   while (start_lineNumber < stop_lineNumber) {
      if (inputs.findBlock(start_lineNumber,stop_lineNumber, &block_start, &block_stop,
                                "Port", "EndPort", false)) {
         fail=true;
      } else {
         if (block_start >= 0 && block_stop >= 0) {
            Port *newPort=new Port(block_start,block_stop);
            portList.push_back(newPort);
         }
      }
      start_lineNumber=inputs.get_next_lineNumber(block_stop);
   }
   return fail;
}

bool BoundaryDatabase::load(const char *filename, bool check_closed_loop) {
   bool fail=false;
   int dim=3;

   if (strcmp(filename,"") == 0) return false;  // this file is optional

   PetscPrintf(PETSC_COMM_WORLD,"   loading port definition file \"%s\"\n",filename);

   if (inputs.load(filename)) return true;
   inputs.createCrossReference();

   if (inputs.checkVersion(version_name, version_value)) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3106: Version mismatch.  Expecting the first line to be: %s %s\n",
                                   indent.c_str(),indent.c_str(),version_name.c_str(),version_value.c_str());
      return true;
   }

   // Source
   if (findSourceFileBlocks()) {fail=true;}

   long unsigned int i=0;
   while (i < sourceFileList.size()) {
      if (sourceFileList[i]->load(&indent, &inputs)) fail=true;
      i++;
   }

   // Path
   if (findPathBlocks()) fail=true;

   i=0;
   while (i < pathList.size()) {
      if (pathList[i]->load(dim,&indent, &inputs)) fail=true;
      i++;
   }

   // Boundary

   if (findBoundaryBlocks()) fail=true;

   i=0;
   while (i < boundaryList.size()) {
      if (boundaryList[i]->load(&indent, &inputs)) fail=true;
      else {
         boundaryList[i]->assignPathIndices(&pathList);
         if (boundaryList[i]->merge(&pathList)) fail=true;
         else {
            if (boundaryList[i]->createRotated(&pathList,indent)) fail=true;
         }
      }
      i++;
   }

   // Port

   if (findPortBlocks()) fail=true;

   i=0;
   while (i < portList.size()) {
      if (portList[i]->load(&indent, &inputs)) fail=true;
      else {
         portList[i]->assignPathIndices(&pathList);
         if (portList[i]->merge(&pathList)) fail=true;
         else {
            if (portList[i]->createRotated(&pathList,indent)) fail=true;
         }
      }
      i++;
   }

   // database checks
   if (check(check_closed_loop)) fail=true;

   // subdivide the paths to eliminate partial overlaps
//xxx
// ToDo
// This is corrupting the path.  Need to fix.
//   if (!fail) subdivide_paths();


   if (fail) PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3107: Failed to load port definitions.\n",indent.c_str(),indent.c_str());

   return fail;
}

void BoundaryDatabase::print()
{
   PetscPrintf(PETSC_COMM_WORLD,"%s %s\n",version_name.c_str(),version_value.c_str());

   long unsigned int i=0;
   while (i < sourceFileList.size()) {
      sourceFileList[i]->print();
      i++;
   }

   i=0;
   while (i < pathList.size()) {
      pathList[i]->print("");
      i++;
   }

   i=0;
   while (i < boundaryList.size()) {
      boundaryList[i]->print();
      i++;
   }

   i=0;
   while (i < portList.size()) {
      portList[i]->print();
      i++;
   }
}

bool BoundaryDatabase::inBlocks(int lineNumber)
{
   long unsigned int i=0;
   while (i < pathList.size()) {
      if (pathList[i]->inBlock(lineNumber)) return true;
      i++;
   }

   i=0;
   while (i < boundaryList.size()) {
      if (boundaryList[i]->inBlock(lineNumber)) return true;
      i++;
   }

   i=0;
   while (i < sourceFileList.size()) {
      if (sourceFileList[i]->inBlock(lineNumber)) return true;
      i++;
   }

   i=0;
   while (i < portList.size()) {
      if (portList[i]->inBlock(lineNumber)) return true;
      i++;
   }

   return false;
}


bool BoundaryDatabase::check(bool check_closed_loop)
{
   bool fail=false;

   // Source
   if (sourceFileList.size() > 1) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3108: Only one File block is allowed.\n",indent.c_str(),indent.c_str());
      fail=true;
   }

   long unsigned int i=0;
   while (i < sourceFileList.size()) {
      if (sourceFileList[i]->check (&indent)) fail=true;
      i++;
   }

   // Path
   i=0;
   while (i < pathList.size()) {

      // individual block checks
      if (pathList[i]->check(&indent)) fail=true;

      // cross block checks

      // duplicated names
      long unsigned int j=i+1;
      while (pathList.size() > 0 && i < pathList.size()-1 && j < pathList.size()) {
         if (pathList[i]->name_is_loaded() && pathList[j]->name_is_loaded()) {
            if (pathList[i]->get_name().compare(pathList[j]->get_name()) == 0) {
               PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3109: name at line %d duplicates the name at line %d.\n",
                                            indent.c_str(),indent.c_str(),pathList[j]->get_name_lineNumber(),pathList[i]->get_name_lineNumber());
               fail=true;
            }
         }
         j++;
      }

      i++;
   }

   // Boundary
   i=0;
   while (i < boundaryList.size()) {

      // individual block checks
      if (boundaryList[i]->check(&indent,pathList)) fail=true;

      // cross block checks

      // duplicated names
      long unsigned int j=i+1;
      while (boundaryList.size() > 0 && i < boundaryList.size()-1 && j < boundaryList.size()) {
         if (boundaryList[i]->name_is_loaded() && boundaryList[j]->name_is_loaded()) {
            if (boundaryList[i]->get_name().compare(boundaryList[j]->get_name()) == 0) {
               PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3110: name at line %d duplicates the name at line %d.\n",
                                            indent.c_str(),indent.c_str(),boundaryList[j]->get_name_lineNumber(),boundaryList[i]->get_name_lineNumber());
               fail=true;
            }
         }
         j++;
      }

      i++;
   }

   // Port

   i=0;
   while (i < portList.size()) {
      if (portList[i]->check(&indent,pathList,check_closed_loop)) fail=true;
      i++;
   }

   // check for extraneous text
   i=1;  // skip the first line, which is the version information
   while (i < inputs.get_size()) {
      if (! inBlocks(inputs.get_lineNumber(i))) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3111: Invalid input at line %d.\n",indent.c_str(),indent.c_str(),inputs.get_lineNumber(i));
         fail=true;
      }
      i++;
   }

   return fail;
}

// check that the S-parameter ports are legally numbered across all ports and modes
// i.e. sequential starting at 1 with no duplicates or gaps
bool BoundaryDatabase::checkSportNumbering()
{
   bool fail=false;

   vector<int> subList;
   vector<int> fullList;

   long unsigned int i=0;
   while (i < portList.size()) {
      subList=portList[i]->get_SportList();

      long unsigned int j=0;
      while (j < subList.size()) {
         fullList.push_back(subList[j]);
         j++;
      }

      subList.clear();

      i++;
   }

   sort(fullList.begin(),fullList.end());

   if (fullList.size() == 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3112: No S-parameter ports are defined.\n",indent.c_str(),indent.c_str());
      fail=true;
   }

   if (!fail && fullList[0] != 1) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3113: S-parameter ports must start numbering with 1.\n",indent.c_str(),indent.c_str());
      fail=true;
   }

   if (fullList.size() > 0) {
      i=0;
      while (i < fullList.size()-1) {
         if (fullList[i] != fullList[i+1]-1) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3114: S-parameter ports are not numbered sequentially.\n",indent.c_str(),indent.c_str());
            fail=true;
         }
         i++;
      }
   }

   return fail;
}


// Make sure the paths fall within the bounding box
// to catch scaling errors.
bool BoundaryDatabase::check_scale (Mesh *mesh, int order)
{
   bool fail=false;

   Vector lowerLeft,upperRight;
   mesh->GetBoundingBox(lowerLeft,upperRight,max(order,1));

   long unsigned int i=0;
   while (i < boundaryList.size()) {
      if (boundaryList[i]->checkBoundingBox(&lowerLeft,&upperRight,&indent,tol,&pathList)) fail=true;
      i++;
   }

   i=0;
   while (i < portList.size()) {
      if (portList[i]->checkBoundingBox(&lowerLeft,&upperRight,&indent,tol,&pathList)) fail=true;
      i++;
   }

   if (fail) {
      PetscPrintf(PETSC_COMM_WORLD,"         Bounding box: (%g,%g,%g) to (%g,%g,%g)\n",
         lowerLeft[0],lowerLeft[1],lowerLeft[2],upperRight[0],upperRight[1],upperRight[2]);
   }

   return fail;
}

bool BoundaryDatabase::check_overlaps ()
{
   bool fail=false;

   // ports cannot overlap ports
   long unsigned int i=0;
   while (i < portList.size()) {
      long unsigned int j=0;
      while (j < portList.size()) {
         if (i != j && portList[i]->is_overlapPath(pathList[portList[j]->get_pathIndex(0)])) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3115: Port %s and Port %s overlap.\n",
               indent.c_str(),indent.c_str(),portList[i]->get_name().c_str(),portList[j]->get_name().c_str());
            fail=true;
         }
         j++;
      }
      i++;
   }

   // ports cannot overlap boundaries
   i=0;
   while (i < portList.size()) {
      long unsigned int j=0;
      while (j < boundaryList.size()) {
         if (boundaryList[j]->is_overlapPath(&pathList,pathList[portList[i]->get_pathIndex(0)])) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3116: Port %s and Boundary %s overlap.\n",
               indent.c_str(),indent.c_str(),portList[i]->get_name().c_str(),boundaryList[j]->get_name().c_str());
            fail=true;
         }
         j++;
      }
      i++;
   }

   // boundaries cannot overlap boundaries
   i=0;
   while (i < boundaryList.size()) {
      long unsigned int j=0;
      while (j < boundaryList.size()) {
         if (i != j && boundaryList[i]->is_overlapPath(&pathList,pathList[boundaryList[j]->get_pathIndex(0)])) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3117: Boundary %s and Boundary %s overlap.\n",
               indent.c_str(),indent.c_str(),boundaryList[i]->get_name().c_str(),boundaryList[j]->get_name().c_str());
            fail=true;
         }
         j++;
      }
      i++;
   }

   return fail;
}

// remove overlaps in paths
void BoundaryDatabase::subdivide_paths ()
{
   long unsigned int i=0;
   while (i < pathList.size()) {
      long unsigned int j=0;
      while (j < pathList.size()) {
         if (i != j) pathList[i]->subdivide3D(pathList[j]);
         j++;
      }
      i++;
   }
}

void BoundaryDatabase::create2Dmeshes(int order, ParMesh *mesh3D, vector<ParSubMesh> *parSubMeshes)
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->create2Dmesh (order,mesh3D,parSubMeshes,i,tol);
      i++;
   }
}

// assign unique attributes for the ports and boundaries
void BoundaryDatabase::assignAttributes ()
{
   // 0 is not allowed by MFEM
   // 1 is reserved for PEC as the default
   // so start with 2
   int attribute=2;

   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->set_attribute(attribute);
      attribute++;
      i++;
   }

   i=0;
   while (i < boundaryList.size()) {
      if (boundaryList[i]->is_perfect_electric_conductor()) {
         boundaryList[i]->set_attribute(1);  // PEC is always 1
      } else {
         boundaryList[i]->set_attribute(attribute);   // PMC gets an attribute, but nothing is done with it (the natural boundary condition)
         attribute++;
      }
      i++;
   }
}

Port* BoundaryDatabase::get_port (int attribute)
{
   long unsigned int i=0;
   while (i < portList.size()) {
      if (portList[i]->get_attribute() == attribute) return portList[i];
      i++;
   }
   return nullptr;
}

// mark the mesh boundaries with attributes into the ports and boundaries
void BoundaryDatabase::markMeshBoundaries (Mesh *mesh)
{
   DenseMatrix pointMat(3,3);

   // loop through the boundary elements
   int i=0;
   while (i < mesh->GetNBE()) {

      int attribute=-1;
      if (mesh->GetBdrElementType(i) == Element::TRIANGLE) {
         mesh->GetBdrPointMatrix(i,pointMat);

         // loop through the ports
         long unsigned int j=0;
         while (j < portList.size()) {
            if (portList[j]->is_triangleInside(&pointMat)) {
                attribute=portList[j]->get_attribute();
                break;
            }
            j++;
         }

         // loop through the boundaries
         j=0;
         while (j < boundaryList.size()) {
            if (boundaryList[j]->is_triangleInside(&pointMat)) { 
                attribute=boundaryList[j]->get_attribute();
                break;
            }
            j++;
         }
      }

      if (attribute > 0) {
         mesh->SetBdrAttribute(i,attribute);
      } else {
         // set to 1 as the PEC default
         mesh->SetBdrAttribute(i,1);          
      }

      i++;
   }
}

int BoundaryDatabase::getLastAttribute ()
{
   int attribute=-1;

   long unsigned int i=0;
   while (i < portList.size()) {
      if (portList[i]->get_attribute() > attribute) attribute=portList[i]->get_attribute();
      i++;
   }

   i=0;
   while (i < boundaryList.size()) {
      if (boundaryList[i]->get_attribute() > attribute) attribute=boundaryList[i]->get_attribute();
      i++;
   } 

   return attribute;
}

bool BoundaryDatabase::setDefaultMeshBoundary (struct projectData *projData, Mesh *mesh,
                                               MaterialDatabase *materialDatabase, BoundaryDatabase *boundaryDatabase)
{
   // nothing to do if PEC
   if (strcmp(projData->materials_default_boundary,"PEC") == 0) return false;

   // get the material
   string default_material_name=projData->materials_default_boundary;
   Material *default_material=materialDatabase->get(default_material_name);
   if (! default_material) {
      PetscPrintf(PETSC_COMM_WORLD,
         "%sERROR3007: \"%s\" specified by material.default.boundary in the project file does not exist in the materials database.\n",
         indent.c_str(),projData->materials_default_boundary);
      return true;
   }

   // create a Boundary to hold the default boundary - note that this includes a reduced set of information
   int new_attribute=getLastAttribute()+1;
   Boundary *default_boundary=new Boundary(0,0);
   default_boundary->set_attribute(new_attribute);
   default_boundary->set_name("default_boundary");
   default_boundary->set_type("surface_impedance");
   default_boundary->set_material(default_material_name);
   boundaryList.push_back(default_boundary);

   // sweep through the mesh and replace attributes of 1 with this attribute
   // converts from PEC to the new Boundary
   int i=0;
   while (i < mesh->GetNBE()) {
      int attribute=mesh->GetBdrAttribute(i);
      if (attribute == 1) mesh->SetBdrAttribute(i,new_attribute);
      i++;
   }

   return false;
}

void BoundaryDatabase::savePortMeshes (meshMaterialList *materials, vector<ParSubMesh> *parSubMeshes)
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->saveMesh(materials,&tempDirectory,&((*parSubMeshes)[i]));
      i++;
   }
}

void BoundaryDatabase::save2Dsetups (struct projectData *projData, double frequency, GammaDatabase *gammaDatabase)
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->save2Dsetup(projData,&tempDirectory,frequency,gammaDatabase->getGamma(i));
      i++;
   }
}

void BoundaryDatabase::saveModeFiles (struct projectData *projData)
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->saveModeFile (projData,&pathList,this);
      i++;
   }
}

// return the boundary that encloses the segment defined by (x1,y1,z1) to (x2,y2,z2) [unrotated]
Boundary* BoundaryDatabase::get_matchBoundary (double x1, double y1, double z1, double x2, double y2, double z2)
{
   long unsigned int i=0;
   while (i < boundaryList.size()) {
      if (boundaryList[i]->is_point_inside(x1,y1,z1) && boundaryList[i]->is_point_inside(x2,y2,z2)) return boundaryList[i];
      i++;
   }
   return nullptr;
}

bool BoundaryDatabase::createPortDirectories ()
{
   bool fail=false;
   long unsigned int i=0;
   while (i < portList.size()) {
      if (portList[i]->createDirectory(&tempDirectory)) fail=true;
      i++;
   }
   return fail;
}

void BoundaryDatabase::extract2Dmesh(ParMesh *pmesh, vector<ParSubMesh> *parSubMeshes)
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->extract2Dmesh(pmesh,parSubMeshes);
      i++;
   }
}

bool BoundaryDatabase::solvePorts (int mesh_order, ParMesh *pmesh, vector<ParSubMesh> *parSubMeshes, double frequency, meshMaterialList *meshMaterials,
                                   struct projectData *projData, GammaDatabase *gammaDatabase)
{
   bool fail=false;
   PetscMPIInt rank;
   chrono::duration<double> elapsed;
   chrono::system_clock::time_point start;
   chrono::system_clock::time_point current;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   MPI_Comm MPI_PORT_COMM;

   create2Dmeshes(mesh_order,pmesh,parSubMeshes); 

   savePortMeshes(meshMaterials,parSubMeshes);
   if (rank == 0) {
      save2Dsetups(projData,frequency,gammaDatabase);
      saveModeFiles (projData);
   }

   long unsigned int i=0;
   while (i < portList.size()) {
      PetscPrintf(PETSC_COMM_WORLD,"----------------------------------------------------------------------------------------------------------------------------------\n");
      PetscPrintf(PETSC_COMM_WORLD,"Port \"%s\"\n",portList[i]->get_name().c_str());
      PetscPrintf(PETSC_COMM_WORLD,"----------------------------------------------------------------------------------------------------------------------------------\n");

      if (portList[i]->solve(&tempDirectory,&MPI_PORT_COMM)) fail=true;

      // wait for the lock file to disappear
      stringstream ssLock;
      ssLock << tempDirectory << "/" << "S" << portList[i]->get_name() << "/" << "." << portList[i]->get_name() << ".lock";
      start=chrono::system_clock::now();
      while (std::filesystem::exists(ssLock.str().c_str())) {
         current=chrono::system_clock::now();
         elapsed=current-start;
         if (elapsed.count() > 60) {
            PetscPrintf(PETSC_COMM_WORLD,"ERROR3118: OpenParEM2D lock file is present, implying a failed 2D port simulation.\n");
            break;
         }
      }

      i++;
   }

   if (loadPortSolutions(frequency)) fail=true;
   if (getGamma(gammaDatabase,frequency)) fail=true;

   return fail;
}

void BoundaryDatabase::markPortModes ()
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->markUsedModes();
      i++;
   }
}

bool BoundaryDatabase::loadPortSolutions (double frequency)
{
   bool fail=false;
   long unsigned int i=0;
   while (i < portList.size()) {
      if (portList[i]->loadSolution(&tempDirectory,frequency)) fail=true;
      i++;
   }
   return fail;
}

bool BoundaryDatabase::getGamma (GammaDatabase *gammaDatabase, double frequency)
{
   bool fail=false;

   gammaDatabase->reset();

   long unsigned int i=0;
   while (i < portList.size()) {
      Gamma *gamma=new Gamma;
      if (portList[i]->getGamma(gamma,frequency)) fail=true;
      gammaDatabase->push(gamma);
      i++;
   }
   return fail;
}

void BoundaryDatabase::printPortSolutions()
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->printSolution("   ");
      i++;
   }
}

int BoundaryDatabase::get_totalModeCount()
{
   int count=0;
   long unsigned int i=0;
   while (i < portList.size()) {
      count+=portList[i]->get_modeCount();
      i++;
   }
   return count;
}

int BoundaryDatabase::get_SportCount()
{
   int count=0;
   long unsigned int i=0;
   while (i < portList.size()) {
      if (portList[i]->get_SportCount() > count) count=portList[i]->get_SportCount();
      i++;
   }
   return count;
}

void BoundaryDatabase::set2DModeNumbers()
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->set2DModeNumbers();
      i++;
   }
}

void BoundaryDatabase::fillUnused2DModeNumbers()
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->fillUnused2DModeNumbers();
      i++;
   }
}

void BoundaryDatabase::build_portEssTdofLists (ParFiniteElementSpace *fespace, ParMesh *pmesh)
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->build_essTdofList(fespace,pmesh);
      i++;
   }
}

void BoundaryDatabase::fillX (Vec *X, Vec *Xdofs, int drivingSport)
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->fillX(X,Xdofs,drivingSport);
      i++;
   }
}

bool BoundaryDatabase::addMassPortIntegrators (ParMesh *pmesh, ParMixedBilinearForm *pmblf, PWConstCoefficient *Inv_mur, vector<Array<int> *> &borderAttributesList,
                                               vector<ConstantCoefficient *> &alphaConstList, vector<ConstantCoefficient *> &betaConstList, bool isReal)
{
   bool fail=false;
   long unsigned int i=0;
   while (i < portList.size()) {
      if (portList[i]->addMassPortIntegrators(pmesh,pmblf,Inv_mur,borderAttributesList,alphaConstList,betaConstList,isReal)) fail=true;
      i++;
   }

   return fail;
}

void BoundaryDatabase::addMassImpedanceIntegrators (double frequency, double temperature, ParMesh *pmesh, ParMixedBilinearForm *pmblf, MaterialDatabase *materialDatabase,
                                                    vector<Array<int> *> &borderAttributesList, vector<ConstantCoefficient *> &ZconstList, bool isReal)
{
   long unsigned int i=0;
   while (i < boundaryList.size()) {
      boundaryList[i]->addMassImpedanceIntegrator(frequency,temperature,pmesh,pmblf,materialDatabase,borderAttributesList,ZconstList,isReal);
      i++;
   }
}


void BoundaryDatabase::set_V3Dsize ()
{
   int size=get_totalModeCount();

   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->set_V3Dsize(size);
      i++;
   }
}

void BoundaryDatabase::calculateWeights(Mode *drivingMode)
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->calculateWeights(drivingMode);
      i++;
   }
}

double BoundaryDatabase::get_maxReflection ()
{
   double reflection;
   double maxReflection=-DBL_MAX;

   long unsigned int i=0;
   while (i < portList.size()) {
      reflection=portList[i]->get_maxReflection();
      if (reflection > maxReflection) maxReflection=reflection;
      i++;
   }

   return maxReflection;
}

double BoundaryDatabase::get_maxReflection (int drivingSport)
{
   double reflection;
   double maxReflection=-DBL_MAX;

   long unsigned int i=0;
   while (i < portList.size()) {
      reflection=portList[i]->get_maxReflection(drivingSport);
      if (reflection > maxReflection) maxReflection=reflection;
      i++;
   }

   return maxReflection;
}

void BoundaryDatabase::printPortReflections ()
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->printPortReflections();
      i++;
   }
}

Mode* BoundaryDatabase::getDrivingMode (int drivingSport)
{
   Mode *mode=nullptr;

   long unsigned int i=0;
   while (i < portList.size()) {
      mode=portList[i]->getDrivingMode(drivingSport);
      if (mode) break;
      i++;
   }

   return mode;
}

void BoundaryDatabase::fillIntegrationPoints ()
{
   long unsigned int i=0;
   while (i < portList.size()) {
         portList[i]->fillIntegrationPoints(&pathList);
      i++;
   }
}

void BoundaryDatabase::calculateVoltages(ParMesh *pmesh, ParGridFunction *grid_re, ParGridFunction *grid_im, Mode *drivingMode)
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->calculateVoltages(pmesh,&pathList,grid_re,grid_im,drivingMode);
      i++;
   }
}

void BoundaryDatabase::calculateSparameters(Result *result, Mode *drivingMode)
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->calculateSparameter(result,drivingMode);
      i++;
   }
}

void BoundaryDatabase::build2Dgrids()
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->build2Dgrids();
      i++;
   }
}

void BoundaryDatabase::build2DModalGrids()
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->build2DModalGrids();
      i++;
   }
}

void BoundaryDatabase::build3Dgrids(ParFiniteElementSpace *fes3D_ND, ParFiniteElementSpace *fes3D_H1)
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->build3Dgrids(fes3D_ND,fes3D_H1);
      i++;
   }
}

void BoundaryDatabase::build2DSolutionGrids()
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->build2DSolutionGrids();
      i++;
   }
}

void BoundaryDatabase::transfer_2Dsolution_2Dgrids_to_3Dgrids()
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->transfer_2Dsolution_2Dgrids_to_3Dgrids();
      i++;
   }
}

void BoundaryDatabase::transfer_2Dsolution_3Dgrids_to_2Dgrids()
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->transfer_2Dsolution_3Dgrids_to_2Dgrids();
      i++;
   }
}

void BoundaryDatabase::transfer_3Dsolution_3Dgrids_to_2Dgrids(fem3D *fem)
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->transfer_3Dsolution_3Dgrids_to_2Dgrids(fem);
      i++;
   }
}

void BoundaryDatabase::save2DParaView(vector<ParSubMesh> *parSubMeshes, struct projectData *projData, double frequency, bool add_extension)
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->save2DParaView(&((*parSubMeshes)[i]),projData,frequency,add_extension);
      i++;
   }
}

void BoundaryDatabase::save3DParaView(ParMesh *pmesh, struct projectData *projData, double frequency, bool add_extension)
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->save3DParaView(pmesh,projData,frequency,add_extension);
      i++;
   }
}

void BoundaryDatabase::save2DSolutionParaView(vector<ParSubMesh> *parSubMeshes, struct projectData *projData, double frequency, int drivingSport, bool add_extension)
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->save2DSolutionParaView(&((*parSubMeshes)[i]),projData,frequency,drivingSport,add_extension);
      i++;
   }
}

void BoundaryDatabase::save2DModalParaView(vector<ParSubMesh> *parSubMeshes, struct projectData *projData, double frequency, bool add_extension)
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->save2DModalParaView(&((*parSubMeshes)[i]),projData,frequency,add_extension);
      i++;
   }
}

void BoundaryDatabase::buildGrids (fem3D *fem)
{
   build3Dgrids(fem->get_fespace_ND(),fem->get_fespace_H1());
   transfer_2Dsolution_2Dgrids_to_3Dgrids();
   build2DModalGrids();
   transfer_2Dsolution_3Dgrids_to_2Dgrids();
   build2DSolutionGrids();
}

bool BoundaryDatabase::solve2Dports (ParMesh *pmesh, vector<ParSubMesh> *parSubMeshes, struct projectData *projData, double frequency, meshMaterialList *meshMaterials, GammaDatabase *gammaDatabase)
{
   bool fail=false;
   PetscMPIInt rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   // extract 2D meshes for the ports from the 3D mesh
   parSubMeshes->clear();
   extract2Dmesh(pmesh,parSubMeshes);

   // solve 2D ports

   stringstream ssTemp;
   ssTemp << "temp_" << projData->project_name;
   set_tempDirectory(ssTemp.str());

   if (rank == 0) {
      if (std::filesystem::exists(get_tempDirectory().c_str())) {
        std::filesystem::remove_all(get_tempDirectory().c_str());
      }

      // create the temp directory
      if (! std::filesystem::create_directory(get_tempDirectory().c_str())) {
         PetscPrintf(PETSC_COMM_WORLD,"ERROR3119: Failed to create results directory \"%s\".\n",get_tempDirectory().c_str());
         fail=true;
      }
   }
   if (fail) return fail;

   if (rank == 0 && createPortDirectories()) fail=true;
   if (fail) return fail;

   if (solvePorts(projData->mesh_order,pmesh,parSubMeshes,frequency,meshMaterials,projData,gammaDatabase)) fail=true;
   if (fail) return fail;

   build2Dgrids();

   save2DParaView(parSubMeshes,projData,frequency,false);

   return fail;
}

void BoundaryDatabase::resetElementNumbers ()
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->resetElementNumbers();
      i++;
   }
}

void BoundaryDatabase::reset ()
{
   long unsigned int i=0;
   while (i < portList.size()) {
      portList[i]->reset();
      i++;
   }

   tempDirectory="";
}

BoundaryDatabase::~BoundaryDatabase ()
{
   long unsigned int i=0;
   while (i < sourceFileList.size()) {
      if (sourceFileList[i]) {
         delete sourceFileList[i];
         sourceFileList[i]=nullptr;
      }
      i++;
   }

   i=0;
   while (i < portList.size()) {
      if (portList[i]) {
         delete portList[i];
         portList[i]=nullptr;
      }
      i++;
   }

   i=0;
   while (i < boundaryList.size()) {
      if (boundaryList[i]) {
         delete boundaryList[i];
         boundaryList[i]=nullptr;
      }
      i++;
   }

   i=0;
   while (i < pathList.size()) {
      if (pathList[i]) {
         delete pathList[i];
         pathList[i]=nullptr;
      }
      i++;    
   }
}

