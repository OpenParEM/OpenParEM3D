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

#include "results.hpp"
#include "fem3D.hpp"

bool isClose (double a, double b)
{
   double tolerance=1e-12;

   if (a == b) return true;
   if (a == 0 && abs(b) < tolerance) return true;
   if (b == 0 && abs(a) < tolerance) return true;
   if (abs((b-a)/a) < tolerance) return true;
   return false;
}

//---------------------------------------------------------------------------------------------------------------------------------
// Sparm
//---------------------------------------------------------------------------------------------------------------------------------

Sparm::Sparm (complex<double> S_, int portOut_, int portIn_)
{
   S=S_;
   portOut=portOut_;
   portIn=portIn_;
}

bool Sparm::is_match (int portOut_, int portIn_)
{
   if (portOut_ == portOut && portIn_ == portIn) return true;
   return false;
}

void Sparm::print ()
{
   cout << "      S(" << portOut << "," << portIn << ")=" << S << endl;
}

void Sparm::save_as_test (ofstream *out, const char *casename, int frequency_index, int *casenumber, double frequency, 
                          double magLimitdB, double equalMagLimit, double argLimitdeg, double equalArgLimit)
{

   // mag

   *out << casename << "_" << frequency_index+1 << "_" << (*casenumber)++ << "_result" << ",";
   *out << setprecision(15) << frequency << ",";
   *out << "magS(dB),";
   *out << portOut << ",";
   *out << portIn << ",";

   if (20*log10(abs(S)) > magLimitdB) {
      *out << "equal" << ",";
      *out << setprecision(15) << 20*log10(abs(S)) << ",";
      *out << equalMagLimit << endl;
   } else {
      *out << "lessthan" << ",";
      *out << magLimitdB << endl;
   }

   // phase

   if (20*log10(abs(S)) > magLimitdB) {

      *out << casename << "_" << frequency_index+1 << "_" << (*casenumber)++ << "_result" << ",";
      *out << setprecision(15) << frequency << ",";
      *out << "argS(deg),";
      *out << portOut << ",";
      *out << portIn << ",";

      if (abs(arg(S)*180/M_PI) > argLimitdeg) {
         *out << "equal" << ",";
         *out << setprecision(15) << arg(S)*180/M_PI << ",";
         *out << equalArgLimit << endl;
      } else {
         *out << "lessthan" << ",";
         *out << argLimitdeg << endl;
      }

   }
}

//---------------------------------------------------------------------------------------------------------------------------------
// SparmList
//---------------------------------------------------------------------------------------------------------------------------------

void SparmList::push (complex<double> S, int portOut, int portIn)
{
   Sparm *newSparm=new Sparm(S,portOut,portIn);
   parmList.push_back(newSparm);
}

complex<double> SparmList::get_S (int portOut, int portIn)
{
   complex<double> retval=complex<double>(-DBL_MAX,-DBL_MAX);
   long unsigned int i=0;
   while (i < parmList.size()) {
      if (parmList[i]->is_match(portOut,portIn)) {
         retval=parmList[i]->get_S();
         break;
      }
      i++;
   }

   return retval;
}

void SparmList::set_S (complex<double> Sparam, int portOut, int portIn)
{
   bool found=false;
   long unsigned int i=0;
   while (i < parmList.size()) {
      if (parmList[i]->is_match(portOut,portIn)) {
         parmList[i]->set_S(Sparam);
         found=true;
         break;
      }
      i++;
   }
   if (!found) cout << "ASSERT: SparmList::set_S failed to set a value for S(" << portOut << "," << portIn << ")." << endl;
}

void SparmList::print ()
{
   cout << "   SparmList: " << this << endl;
   long unsigned int i=0;
   while (i < parmList.size()) {
      parmList[i]->print();
      i++;
   }
}

void SparmList::save_as_test (ofstream *out, const char *casename, int frequency_index, int *casenumber, double frequency,
                              double magLimitdB, double equalMagLimit, double argLimitdeg, double equalArgLimit)
{
   long unsigned int i=0;
   while (i < parmList.size()) {
      parmList[i]->save_as_test (out,casename,frequency_index,casenumber,frequency,
                                 magLimitdB,equalMagLimit,argLimitdeg,equalArgLimit);
      i++;
   }
}

int commaCount (string a)
{
   int count=0;
   long unsigned int i=0;
   while (i < a.length()) {
      if (a[i] == ',') count ++;
      i++;
   }
   return count;
}

SparmList::~SparmList ()
{
   long unsigned int i=0;
   while (i < parmList.size()) {
      delete parmList[i];
      i++;
   }
}

//---------------------------------------------------------------------------------------------------------------------------------
// Result
//---------------------------------------------------------------------------------------------------------------------------------
Result::Result ()
{
   isNormalized=false;
}

void Result::set_S (complex<double> Sparam, int portOut, int PortIn)
{
   S.set_S(Sparam,portOut,PortIn);
}

void Result::save (ostream *out, struct projectData *projData, int SportCount)
{
   int rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   if (rank == 0) {
      *out << "[Result]" << endl;
      *out << "   drivingSport=" << drivingSport << endl;
      *out << "   active=" << active << endl;
      *out << "   iteration=" << iteration << endl;
      *out << setprecision(15) << "   frequency=" << frequency << endl;
      *out << "   shortestPerWavelength=" << shortestPerWavelength << endl;
      *out << "   longestPerWavelength=" << longestPerWavelength << endl;
      *out << "   maxReflection=" << maxReflection << endl;
      *out << "   Zo=(" << real(Zo) << "," << imag(Zo) << ")" << endl;
      *out << "   normalizeZo=" << normalizeZo << endl;
      *out << "   isNormalized=" << isNormalized << endl;
      if (strcmp(projData->touchstone_format,"RI") == 0) *out << "   S(R,I)=";
      if (strcmp(projData->touchstone_format,"MA") == 0) *out << "   S(mag,deg)=";
      if (strcmp(projData->touchstone_format,"DB") == 0) *out << "   S(dB,deg)=";

      *out << setprecision(15);
      int i=0;
      while (i < SportCount) {
         complex<double> Sparam=S.get_S(i+1,drivingSport);
         if (strcmp(projData->touchstone_format,"RI") == 0) *out << "(" << real(Sparam) << "," << imag(Sparam) << "),";
         if (strcmp(projData->touchstone_format,"MA") == 0) *out << "(" << abs(Sparam) << "," << arg(Sparam)*180/M_PI << "),";
         if (strcmp(projData->touchstone_format,"DB") == 0) *out << "(" << 20*log10(abs(Sparam)) << "," << arg(Sparam)*180/M_PI << "),";
         i++;
      }
      *out << endl;
      *out << setprecision(6);

      *out << "   mesh_size=" << meshSize << endl;
      *out << "   matrix_size(DOF_count)=" << matrixSize << endl;
      *out << "   sparse_width=" << sparseWidth << endl;
      *out << "   EfieldError=" << EfieldError << endl;
      *out << "   HfieldError=" << HfieldError << endl;
      *out << "   EfieldConverged=" << EfieldConverged << endl;
      *out << "   HfieldConverged=" << HfieldConverged << endl;
      if (iteration > 1) *out << "   maxRelativeError=" << maxRelativeError << endl;
      if (iteration > 1) *out << "   hasIterations=" << hasIterations << endl;
      if (iteration > 1) *out << "   isConverged=" << isConverged << endl;
      *out << "   fem_setup_time=" << get_fem_setup_time() << endl;
      *out << "   mesh_refine_time" << get_mesh_error_time() << endl;
      *out << "   solve_time=" << get_solve_time() << endl;
      *out << "   refine_time=" << get_refine_time() << endl;

      *out << "[EndResult]" << endl;
   }
}

void Result::saveFormatted (ostream *out, double solveElapsedTime, double meshErrorTime, double femSetupTime, double refineTime, double *priorFrequency)
{
   int rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   if (rank == 0) {
      if (isClose(frequency,*priorFrequency)) *out << setw(15) << "";
      else *out << setw(15) << frequency;

      *out << setw(12) << iteration;
      *out << setw(12) << meshSize;
      *out << setw(12) << matrixSize;   // DOF count
      *out << setw(12) << setprecision(4) << femSetupTime;
      *out << setw(12) << setprecision(4) << solveElapsedTime;
      *out << setw(12) << setprecision(4) << meshErrorTime;
      *out << setw(12) << setprecision(4) << refineTime;
      *out << setw(12) << setprecision(4) << solveElapsedTime+femSetupTime+meshErrorTime+refineTime;

      if (isClose(frequency,*priorFrequency)) *out << setw(17) << maxRelativeError;
      else  *out << setw(17) << "";

      *priorFrequency=frequency;

      *out << endl;
   }
}

void Result::saveCSV (ostream *out, struct projectData *projData, double scale, int SportCount)
{
   int rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   if (rank == 0) {

      *out << setprecision(15);
      if (drivingSport == 1) *out << frequency*scale;

      int i=0;
      while (i < SportCount) {
         complex<double> Sparam=S.get_S(i+1,drivingSport);
         if (strcmp(projData->touchstone_format,"RI") == 0) *out << "," << real(Sparam) << "," << imag(Sparam);
         if (strcmp(projData->touchstone_format,"MA") == 0) *out << "," << abs(Sparam) << "," << arg(Sparam)*180/M_PI;
         if (strcmp(projData->touchstone_format,"DB") == 0) *out << "," << 20*log10(abs(Sparam)) << "," << arg(Sparam)*180/M_PI;
         i++;
      }
      *out << setprecision(6);
   }
}

void Result::set (int drivingSport_, double frequency_, double shortestPerWavelength_, double longestPerWavelength_, int iteration_, double normalizeZo_, fem3D *fem, ParMesh *pmesh)
{
   drivingSport=drivingSport_;
   frequency=frequency_;
   shortestPerWavelength=shortestPerWavelength_;
   longestPerWavelength=longestPerWavelength_;
   iteration=iteration_+1;
   active=true;
   normalizeZo=normalizeZo_;
   EfieldError=fem->get_EfieldError();
   EfieldConverged=fem->get_EfieldConverged();
   HfieldError=fem->get_HfieldError();
   HfieldConverged=fem->get_HfieldConverged();
   meshSize=getGlobalNE(pmesh);
   matrixSize=fem->get_matrixSize();
   sparseWidth=fem->get_sparseWidth();
}

void Result::print ()
{
   cout << "Result: " << this << endl;
   cout << "   iteration=" << iteration << endl;
   cout << "   frequency=" << frequency << endl;
   cout << "   shortestPerWavelength=" << shortestPerWavelength << endl;
   cout << "   longestPerWavelength=" << longestPerWavelength << endl;
   cout << "   maxReflection=" << maxReflection << endl;
   cout << "   drivingSport=" << drivingSport << endl;
   cout << "   Zo=" << Zo << endl;
   cout << "   normalizeZo=" << normalizeZo << endl;
   cout << "   isNormalized=" << isNormalized << endl;
   cout << "   active=" << active << endl;
   cout << "   magLimitdB=" << magLimitdB << endl;
   cout << "   equalMagLimit=" << equalMagLimit << endl;
   cout << "   argLimitdeg=" << argLimitdeg << endl;
   cout << "   equalArgLimit=" << equalArgLimit << endl;
   cout << "   EfieldError=" << EfieldError << endl;
   cout << "   HfieldError=" << HfieldError << endl;
   cout << "   EfieldConverged=" << EfieldConverged << endl;
   cout << "   HfieldConverged=" << HfieldConverged << endl;
   cout << "   meshSize=" << meshSize << endl;
   cout << "   matrixSize(DOF_count)=" << matrixSize << endl;
   cout << "   sparseWidth=" << sparseWidth << endl;
   cout << "   maxRelativeError=" << maxRelativeError << endl;
   cout << "   hasIterations=" << hasIterations << endl;
   cout << "   isConverged=" << isConverged << endl;
   cout << "   solve_time=" << solve_time << endl;
   cout << "   fem_setup_time=" << fem_setup_time << endl;
   cout << "   refine_time=" << refine_time << endl;
   S.print();
}

void Result::save_as_test (ofstream *out, const char *casename, int frequency_index, int *casenumber)
{
   S.save_as_test (out,casename,frequency_index,casenumber,frequency,magLimitdB,equalMagLimit,argLimitdeg,equalArgLimit); 
}

bool Result::extractS (string line, string frequency_unit, int number_of_ports)
{
   if (line.length() == 0) return true;
   if (line[0] == '#') return true;

   bool foundFrequency=false;
   int portIn=1;
   int portOut=1;
   double data1=0;
   double data2=0;
   bool loaded_data1=false;

   stringstream ssLine(line);
   string value;
   while (std::getline(ssLine,value,',')) {
      if (foundFrequency) {
         if (loaded_data1) {
            data2=stod(value);
            push_S(complex<double>(data1,data2),portOut,portIn);

            portOut++;
            if (portOut > number_of_ports) {portIn++; portOut=1;}

            loaded_data1=false;
         } else {
            data1=stod(value);
            loaded_data1=true;
         }
      } else {
         double frequency=stod(value);
         if (frequency_unit.compare("Hz") == 0) frequency*=1;
         if (frequency_unit.compare("kHz") == 0) frequency*=1e3;
         if (frequency_unit.compare("MHz") == 0) frequency*=1e6;
         if (frequency_unit.compare("GHz") == 0) frequency*=1e9;
         set_frequency(frequency);
         foundFrequency=true;
      }
   }
   return false;
}

Result::~Result ()
{
}

//---------------------------------------------------------------------------------------------------------------------------------
// ResultDatabase
//---------------------------------------------------------------------------------------------------------------------------------

// result given iteration
Result* ResultDatabase::get_Result (double frequency, int drivingSport, int iteration)
{
   long unsigned int i=0;
   while (i < results.size()) {
      if (isClose(results[i]->get_frequency(),frequency) && 
          results[i]->get_drivingSport() == drivingSport &&
          results[i]->get_iteration() == iteration) return results[i];
      i++;
   }
   return nullptr;
}

// active result
Result* ResultDatabase::get_Result (double frequency, int drivingSport)
{
   long unsigned int i=0;
   while (i < results.size()) {
      if (results[i]->is_active() &&
          isClose(results[i]->get_frequency(),frequency) &&
          results[i]->get_drivingSport() == drivingSport) return results[i];
      i++;
   }
   return nullptr;
}

// for use with process3D.cpp
Result* ResultDatabase::get_Result (double frequency)
{
   long unsigned int i=0;
   while (i < results.size()) {
      if (isClose(results[i]->get_frequency(),frequency)) return results[i];
      i++;
   }
   return nullptr;
}

void ResultDatabase::push (Result *result)
{
   Result *test=get_Result(result->get_frequency(),result->get_drivingSport(),result->get_iteration()-1);
   if (test) {
      test->set_inactive();
   }

   result->set_active();
   results.push_back(result);

   // keep a sorted list of unique frequencies

   // unique?
   bool found=false;
   long unsigned int i=0;
   while (i < unique_frequencies.size()) {
      if (fabs(unique_frequencies[i]-result->get_frequency())/unique_frequencies[i] <= tol) {
         found=true;
         break;
      }
      i++;
   }

   if (! found) {

      // save
      unique_frequencies.push_back(result->get_frequency());

      // ascending sort
      found=true;
      while (found) {
         found=false;
         i=0;
         while (i < unique_frequencies.size()-1) {
            long unsigned int j=i+1;
            while (j < unique_frequencies.size()) {
               if (unique_frequencies[j] < unique_frequencies[i]) {
                  found=true;
                  double temp=unique_frequencies[i];
                  unique_frequencies[i]=unique_frequencies[j];
                  unique_frequencies[j]=temp;
               }
               j++;
            }
            i++;
         }
      }
   }

}

double ResultDatabase::calculate_maxRelativeError (double frequency, int iteration)
{
   PetscErrorCode ierr=0;
   double maxRelativeError=-1;

   if (iteration > 1) {

      Mat S;
      ierr=get_S(iteration,frequency,&S,false,true); if (ierr) return maxRelativeError;

      Mat difference;
      ierr=get_S(iteration-1,frequency,&difference,false,true); if (ierr) return maxRelativeError;

      ierr=MatAXPY(difference,-1,S,SAME_NONZERO_PATTERN);

      ierr=MatInvert(&S,0); if (ierr) return maxRelativeError;
  
      Mat error;  
      ierr=MatMatMult(S,difference,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&error); if (ierr) return maxRelativeError;
      MatDestroy(&S);
      MatDestroy(&difference);

      PetscReal *norms=(PetscReal *) malloc(SportCount*sizeof(PetscReal));
      if (norms) {
         ierr=MatGetColumnNorms(error,NORM_2,norms); if (ierr) return maxRelativeError;

         int i=0;
         while (i < SportCount) {
            if (norms[i] > maxRelativeError) maxRelativeError=norms[i];
            i++;
         }

         free(norms);
      }

      MatDestroy(&error);
   }

   return maxRelativeError;
}

void ResultDatabase::set_refine_time (double elapsed, double frequency, int iteration)
{
   int i=1;
   while (i <= SportCount) {
      Result *result=get_Result(frequency,i,iteration);
      if (result) result->set_refine_time(elapsed);
      i++;
   }
}

void ResultDatabase::set_maxRelativeError (double maxRelativeError, double frequency, int iteration)
{
   if (iteration > 1) {
      int i=1;
      while (i <= SportCount) {
         Result *result=get_Result(frequency,i,iteration);
         if (result) result->set_maxRelativeError(maxRelativeError);
         i++;
      }
   }
}

void ResultDatabase::set_hasIterations (bool hasIterations, double frequency, int iteration)
{
   if (iteration > 1) {
      int i=1;
      while (i <= SportCount) {
         Result *result=get_Result(frequency,i,iteration);
         if (result) result->set_hasIterations(hasIterations);
         i++;
      }
   }
}

void ResultDatabase::set_isConverged (bool isConverged, double frequency, int iteration)
{
   if (iteration > 1) {
      int i=1;
      while (i <= SportCount) {
         Result *result=get_Result(frequency,i,iteration);
         if (result) result->set_isConverged(isConverged);
         i++;
      }
   }
}

bool ResultDatabase::hasIterations ()
{
   long unsigned int i=0;
   while (i < results.size()) {
      if (results[i]->get_hasIterations()) {
          return true;
      }
      i++;
   }
   return false;
}

bool ResultDatabase::isAllConverged ()
{
   long unsigned int i=0;
   while (i < results.size()) {
      if (results[i]->is_active() && 
          results[i]->get_hasIterations() && 
         !results[i]->get_isConverged()) {
         return false;
      }
      i++;
   }
   return true;
}

bool ResultDatabase::isSequentialConverged (struct projectData *projData, double frequency)
{
   Result *result=get_Result(frequency,1);
   int iteration=result->get_iteration();

   int i=0;
   while (i < projData->refinement_required_passes) {
      Result *iteration_result=get_Result(frequency,1,iteration-i);
      if (iteration_result) {
         if (!iteration_result->get_isConverged()) return false;
      } else return false;
      i++;
   }

   return true;
}

// across all S-ports at a single frequency
void ResultDatabase::get_totalSolveTime (int index, double *solve_time)
{
   *solve_time=0;

   long unsigned int i=0;
   while (i < results.size()) {
      if (isClose(results[index]->get_frequency(),results[i]->get_frequency()) &&
          results[index]->get_iteration() == results[i]->get_iteration()) {
         (*solve_time)+=results[i]->get_solve_time();
      }
      i++;
   }
}

// across all S-ports at a single frequency
void ResultDatabase::get_totalMeshErrorTime (int index, double *mesh_error_time)
{
   *mesh_error_time=0;

   long unsigned int i=0;
   while (i < results.size()) {
      if (isClose(results[index]->get_frequency(),results[i]->get_frequency()) &&
          results[index]->get_iteration() == results[i]->get_iteration()) {
         (*mesh_error_time)+=results[i]->get_mesh_error_time();
      }
      i++;
   }
}

// across all S-ports at a single frequency
void ResultDatabase::get_totalFEMSetupTime (int index, double *fem_setup_time)
{
   *fem_setup_time=results[index]->get_fem_setup_time();
}

// across all S-ports at a single frequency
void ResultDatabase::get_totalRefineTime (int index, double *refine_time)
{
   *refine_time=0;

   long unsigned int i=0;
   while (i < results.size()) {
      if (isClose(results[index]->get_frequency(),results[i]->get_frequency()) &&
          results[index]->get_iteration() == results[i]->get_iteration()) {
         (*refine_time)+=results[i]->get_refine_time();
      }
      i++;
   }
}

//ToDo: add error checking to ensure that the file wrote out completely
bool ResultDatabase::save (struct projectData *projData)
{
   bool fail=false;
   int rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   if (rank == 0) {
      stringstream ss;
      ss << projData->project_name << "_results.txt";

      ofstream out;
      out.open(ss.str().c_str(),ofstream::out);

      if (out.is_open()) {
         save(&out,projData);
         out << "[Time]" << endl;
         out << "   job_time=" << solve_time << endl;
         out << "[EndTime]" << endl;
         out.close();
      } else {
         PetscPrintf(PETSC_COMM_WORLD,"ERROR3127: Failed to open file \"%s\" for writing.\n",ss.str().c_str());      
         fail=true;
      }
   }

   MPI_Barrier(PETSC_COMM_WORLD);

   return fail;
}

void ResultDatabase::save (ostream *out, struct projectData *projData)
{
   long unsigned int i=0;
   while (i < results.size()) {
      results[i]->save(out,projData,SportCount);
      i++;
   }
}

//ToDo: add error checking to ensure that the file wrote out completely
bool ResultDatabase::saveFormatted (struct projectData *projData)
{
   bool fail=false;
   int rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   if (rank == 0) {
      stringstream ss;
      ss << projData->project_name << "_iterations.txt";

      ofstream out;
      out.open(ss.str().c_str(),ofstream::out);

      if (out.is_open()) {
         saveFormatted(&out);
         out.close();
      } else {
         PetscPrintf(PETSC_COMM_WORLD,"ERROR3128: Failed to open file \"%s\" for writing.\n",ss.str().c_str());
         fail=true;
      }
   }

   MPI_Barrier(PETSC_COMM_WORLD);

   return fail;
}

void ResultDatabase::saveCSV (ostream *out, struct projectData *projData, int portCount, bool allIterations)
{
   double scale=1;

   // header

   *out << "#Touchstone format," << projData->touchstone_format << endl;
   *out << "#frequency unit," << projData->touchstone_frequency_unit << endl;
   *out << "#number of frequencies," << unique_frequencies.size() << endl;
   *out << "#number of ports," << portCount << endl;

   if (allIterations) *out << "#prior iterations pre-pended by #" << endl;

   if (strcmp(projData->touchstone_frequency_unit,"Hz") == 0) {*out << "#Frequency(Hz)"; scale=1;}
   if (strcmp(projData->touchstone_frequency_unit,"kHz") == 0) {*out << "#Frequency(kHz)"; scale=1e-3;}
   if (strcmp(projData->touchstone_frequency_unit,"MHz") == 0) {*out << "#Frequency(MHz)"; scale=1e-6;}
   if (strcmp(projData->touchstone_frequency_unit,"GHz") == 0) {*out << "#Frequency(GHz)"; scale=1e-9;}

   // S-parameters are stored by columns
   int col=0;
   while (col < portCount) {
      int row=0;
      while (row < portCount) {
         if (strcmp(projData->touchstone_format,"RI") == 0) *out << ",Re(S(" << row+1 << ";" << col+1 << ")),Im(S(" << row+1 << ";" << col+1 << "))";
         if (strcmp(projData->touchstone_format,"MA") == 0) *out << ",mag(S(" << row+1 << ";" << col+1 << ")),deg(S(" << row+1 << ";" << col+1 << "))";
         if (strcmp(projData->touchstone_format,"DB") == 0) *out << ",dB(S(" << row+1 << ";" << col+1 << ")),deg(S(" << row+1 << ";" << col+1 << "))";
         row++;
      }
      col++;
   }
   *out << endl;

   // data

   long unsigned int i=0;
   while (i < unique_frequencies.size()) {
      int lastIteration=get_lastIteration(unique_frequencies[i]);

      int k=lastIteration;
      if (allIterations) k=1;

      while (k <= lastIteration) {
         if (k != lastIteration) *out << "#";
         int drivingPort=1;
         while (drivingPort <= portCount) {
            Result *result=get_Result(unique_frequencies[i],drivingPort,k);
            if (result) result->saveCSV(out,projData,scale,SportCount);
            drivingPort++;
         }
         *out << endl;
         k++;
      }

      i++;
   }

}

bool ResultDatabase::saveTouchstone (struct projectData *projData, int portCount)
{
   double scale=1;
   int rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   if (rank == 0) {
      if (strcmp(projData->touchstone_version,"1.1") == 0) {
         stringstream ss;
         ss << projData->project_name << ".s" << portCount << "p";

         ofstream out;
         out.open(ss.str().c_str(),ofstream::out);

         if (!out.is_open()) {
            PetscPrintf(PETSC_COMM_WORLD,"ERROR3175: Failed to open file \"%s\" for writing.\n",ss.str().c_str());
            return true;
         }

         // header
         out << "! " << portCount << "-port S-parameter data" << endl;
         out << "# " << projData->touchstone_frequency_unit << " S " << projData->touchstone_format << " R " << projData->reference_impedance << endl;

         if (strcmp(projData->touchstone_frequency_unit,"Hz") == 0) scale=1;
         if (strcmp(projData->touchstone_frequency_unit,"kHz") == 0) scale=1e-3;
         if (strcmp(projData->touchstone_frequency_unit,"MHz") == 0) scale=1e-6;
         if (strcmp(projData->touchstone_frequency_unit,"GHz") == 0) scale=1e-9;

         // comment line for 1 and 2 ports

         if (portCount < 3) {
            out << "!freq";

            int col=0;
            while (col < portCount) {
               int row=0;
               while (row < portCount) {
                  if (strcmp(projData->touchstone_format,"RI") == 0) out << " ReS"  << row+1 << col+1 << " ImS"  << row+1 << col+1;
                  if (strcmp(projData->touchstone_format,"MA") == 0) out << " magS" << row+1 << col+1 << " angS" << row+1 << col+1;
                  if (strcmp(projData->touchstone_format,"DB") == 0) out << " dBS"  << row+1 << col+1 << " angS" << row+1 << col+1;
                  row++;
               }
               col++;
            }
            out << endl;
         }

         // data

         long unsigned int k=0;
         while (k < unique_frequencies.size()) {
            int lastIteration=get_lastIteration(unique_frequencies[k]);

            out << setprecision(15);
            out << unique_frequencies[k]*scale;

            int i=1;
            while (i <= portCount) {
               int count=0;
               bool row_printed=false;

               int j=1;
               while (j <= portCount) {
                  complex<double> Sparam;
                  if (portCount < 3) {
                     Result *result=get_Result(unique_frequencies[k],i,lastIteration);
                     Sparam=result->get_S(j,i);
                  } else {
                     Result *result=get_Result(unique_frequencies[k],j,lastIteration);
                     Sparam=result->get_S(i,j);
                  }

                  if (strcmp(projData->touchstone_format,"RI") == 0) out << " " << real(Sparam) << " " << imag(Sparam);
                  if (strcmp(projData->touchstone_format,"MA") == 0) out << " " << abs(Sparam) << " " << arg(Sparam)*180/M_PI;
                  if (strcmp(projData->touchstone_format,"DB") == 0) out << " " << 20*log10(abs(Sparam)) << " " << arg(Sparam)*180/M_PI;

                  count++;
                  if (portCount > 2 && (portCount == count || count == 4)) {
                     if (!row_printed) {
                        out << " !row " << i;
                        row_printed=true;
                     }
                     out << endl;
                     count=0;
                  }

                  j++;
               }

               if (portCount > 2 && !row_printed) out << " !row " << i;
               if (portCount > 2 && count > 0) out << endl;

               i++;
            }

            if (portCount <= 2) out << endl;

            k++;
         }

         out.close();
      } else {
         PetscPrintf(PETSC_COMM_WORLD,"ERROR3176: Touchstone version \"%s\" is not supported.\n",projData->touchstone_version);
      }
   }

   MPI_Barrier(PETSC_COMM_WORLD);

   return 0;
}

bool ResultDatabase::saveCSV (struct projectData *projData, int portCount, bool allIterations)
{
   bool fail=false;
   int rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   if (rank == 0) {
      stringstream ss;
      ss << projData->project_name << "_results.csv";

      ofstream out;
      out.open(ss.str().c_str(),ofstream::out);

      if (out.is_open()) {
         saveCSV(&out,projData,portCount,allIterations);
         out.close();
      } else {
         PetscPrintf(PETSC_COMM_WORLD,"ERROR3129: Failed to open file \"%s\" for writing.\n",ss.str().c_str());
         fail=true;
      }
   }

   MPI_Barrier(PETSC_COMM_WORLD);

   return fail;
}

void ResultDatabase::saveFormatted (ostream *out)
{
   *out << setw(15) << "---------------"
        << setw(12) << "------------"
        << setw(12) << "------------"
        << setw(12) << "------------"
        << setw(12) << "------------"
        << setw(12) << "------------"
        << setw(12) << "------------"
        << setw(12) << "------------"
        << setw(12) << "------------"
        << setw(17) << "-----------------"
        << endl;

   *out << setw(15) << "frequency"
        << setw(12) << "iteration"
        << setw(12) << "mesh size"
        << setw(12) << "DOF count"
        << setw(12) << "setup,s"
        << setw(12) << "solve,s"
        << setw(12) << "mesh err,s"
        << setw(12) << "refine,s"
        << setw(12) << "total,s"
        << setw(17) << "relative error"
        << endl;

   *out << setw(15) << "---------------"
        << setw(12) << "------------"
        << setw(12) << "------------"
        << setw(12) << "------------"
        << setw(12) << "------------"
        << setw(12) << "------------"
        << setw(12) << "------------"
        << setw(12) << "------------"
        << setw(12) << "------------"
        << setw(17) << "-----------------"
        << endl;

   double priorFrequency=-DBL_MAX;
   double solveElapsed,meshError,femSetup,refine;
   long unsigned int i=0;
   while (i < results.size()) {
      if (results[i]->get_drivingSport() == 1) {
         get_totalSolveTime(i,&solveElapsed);
         get_totalMeshErrorTime(i,&meshError);
         get_totalFEMSetupTime(i,&femSetup);
         get_totalRefineTime(i,&refine);
         results[i]->saveFormatted(out,solveElapsed,meshError,femSetup,refine,&priorFrequency);
      }
      i++;
   }
}

int ResultDatabase::get_lastIteration (double frequency)
{
   int iteration=0;

   long unsigned int i=0;
   while (i < results.size()) {
      if (isClose(results[i]->get_frequency(),frequency) &&
          results[i]->get_iteration() > iteration) {
         iteration=results[i]->get_iteration();
      }
      i++;
   }
   return iteration;
}

int ResultDatabase::get_SportCount (int iteration, double frequency)
{
   int Sport=0;

   long unsigned int i=0;
   while (i < results.size()) {
      if (results[i]->get_iteration() == iteration &&
          isClose(results[i]->get_frequency(),frequency) &&
          results[i]->get_drivingSport() > Sport) {
         Sport=results[i]->get_drivingSport();
      }
      i++;
   }
   return Sport; 
}

Result* ResultDatabase::get_col (int iteration, double frequency, int drivingSport)
{
   long unsigned int i=0;
   while (i < results.size()) {
      if (results[i]->get_iteration() == iteration &&
         isClose(results[i]->get_frequency(),frequency) &&
         results[i]->get_drivingSport() == drivingSport) {
         return results[i];
      }
      i++;
   }
   return nullptr;
}

complex<double> ResultDatabase::get_Zo (int iteration, double frequency, int drivingSport)
{
   long unsigned int i=0;
   while (i < results.size()) {
      if (results[i]->get_iteration() == iteration &&
         isClose(results[i]->get_frequency(),frequency) &&
         results[i]->get_drivingSport() == drivingSport) {
         return results[i]->get_Zo();
      }
      i++;
   }
   return complex<double>(-DBL_MAX,-DBL_MAX);
}

complex<double> ResultDatabase::get_normalizeZo (int iteration, double frequency, int drivingSport)
{
   long unsigned int i=0;
   while (i < results.size()) {
      if (results[i]->get_iteration() == iteration &&
         isClose(results[i]->get_frequency(),frequency) &&
         results[i]->get_drivingSport() == drivingSport) {
         return results[i]->get_normalizeZo();
      }
      i++;
   }
   return complex<double>(-DBL_MAX,-DBL_MAX);
}

// S must be destroyed elsewhere
PetscErrorCode  ResultDatabase::get_S (int iteration, double frequency, Mat *S, bool flipSign, bool assemble)
{
   PetscErrorCode ierr=0;
   long unsigned int j;
   PetscScalar value;
   PetscInt low,high;

   if (SportCount == 0) return 1;

   ierr=MatCreate(PETSC_COMM_WORLD,S); if (ierr) return 3;
   ierr=MatSetType(*S,MATAIJ); if (ierr) return 4;
   ierr=MatSetSizes(*S,PETSC_DECIDE,PETSC_DECIDE,SportCount,SportCount); if (ierr) return 5;
   ierr=MatSeqAIJSetPreallocation(*S,SportCount,NULL); if (ierr) return 6;
   ierr=MatMPIAIJSetPreallocation(*S,SportCount,NULL,SportCount,NULL); if (ierr) return 7;
   ierr=MatZeroEntries(*S); if (ierr) return 8;

   MatGetOwnershipRange(*S,&low,&high);

   int i=low+1;  // portIn
   while (i < high+1) {
      Result *row=get_col(iteration,frequency,i);
      if (row == nullptr) {
         cout << "ASSERT: ResultDatabase::get_S missing row." << endl;
      } else {
         j=1;  // portOut
         while (j <= (long unsigned int)SportCount) {
            complex<double> Svalue=row->get_S(j,i);
            if (Svalue != complex<double>(-DBL_MAX,-DBL_MAX)) {
               value=real(Svalue)+PETSC_i*imag(Svalue);
               if (flipSign) value=-value;
               ierr=MatSetValue(*S,i-1,j-1,value,ADD_VALUES); if (ierr) return ierr;
            }
            j++;
         }
      }
      i++;
   }

   if (assemble) {
      ierr=MatAssemblyBegin(*S,MAT_FINAL_ASSEMBLY); if (ierr) return ierr;
      ierr=MatAssemblyEnd(*S,MAT_FINAL_ASSEMBLY); if (ierr) return ierr;
   }

   return ierr;
}

PetscErrorCode ResultDatabase::set_S (int iteration, double frequency, Mat *S)
{
   int rank,size;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   MPI_Comm_size(PETSC_COMM_WORLD, &size);

   PetscErrorCode ierr=0;
   int isValid,irow,icol;
   PetscInt m,n,low,high;
   complex<double> Svalue;
   double ReValue,ImValue;

   int SportCount=get_SportCount(iteration,frequency);
   if (SportCount == 0) return 1;

   ierr=MatGetSize(*S,&m,&n);
   if (m != n) return 2;
   if (m != SportCount) return 3;

   // transfer the local rows 
   MatGetOwnershipRange(*S,&low,&high);
   int i=low+1;
   while (i < high+1) {
      Result *row=get_col(iteration,frequency,i);
      if (row == nullptr) {
         cout << "ASSERT: ResultDatabase::set_S missing row." << endl;
      } else {
         int j=1;
         while (j <= SportCount) {
            ierr=MatGetValue(*S,i-1,j-1,&Svalue); if (ierr) return 9;
            Svalue=complex<double>(PetscRealPart(Svalue),PetscImaginaryPart(Svalue));
            row->set_S(Svalue,j,i);
            j++;
         }
         row->set_isNormalized();
      }
      i++;
   }

   // collect data at rank 0
   // Could do some of this above, but this code is simpler and speed/memory is not a factor for
   // the relatively small S-parameter matrices.
   i=1;
   while (i <= SportCount) {
      Result *row=get_col(iteration,frequency,i);
      if (rank == 0) {
         int k=1;
         while (k < size) {
            int j=1;
            while (j <= SportCount) {
               MPI_Recv(&isValid,1,MPI_INT,k,1000,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
               MPI_Recv(&irow,1,MPI_INT,k,1001,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
               MPI_Recv(&icol,1,MPI_INT,k,1002,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
               MPI_Recv(&ReValue,1,MPI_DOUBLE,k,1003,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
               MPI_Recv(&ImValue,1,MPI_DOUBLE,k,1004,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

               Svalue=complex<double>(ReValue,ImValue);
               if (isValid) row->set_S(Svalue,icol,irow);

               j++;
            }
            k++;
         }
      } else {
         int j=1;
         while (j <= SportCount) {
            Svalue=row->get_S(j,i);
            ReValue=PetscRealPart(Svalue);
            ImValue=PetscImaginaryPart(Svalue);

            isValid=0;
            if (i-1 >= low && i-1 < high) isValid=1;

            MPI_Send(&isValid,1,MPI_INT,0,1000,PETSC_COMM_WORLD);
            MPI_Send(&i,1,MPI_INT,0,1001,PETSC_COMM_WORLD);
            MPI_Send(&j,1,MPI_INT,0,1002,PETSC_COMM_WORLD);
            MPI_Send(&ReValue,1,MPI_DOUBLE,0,1003,PETSC_COMM_WORLD);
            MPI_Send(&ImValue,1,MPI_DOUBLE,0,1004,PETSC_COMM_WORLD);

            j++;
         }
      }
      i++;
   }

   // copy the data out to all the other ranks
   if (rank == 0) {
      int i=1;
      while (i <= SportCount) {
         Result *row=get_col(iteration,frequency,i);
         int j=1;
         while (j <= SportCount) {
            Svalue=row->get_S(j,i);
            ReValue=PetscRealPart(Svalue);
            ImValue=PetscImaginaryPart(Svalue);

            int k=1;
            while (k < size) {
               MPI_Send(&i,1,MPI_INT,k,1000,PETSC_COMM_WORLD);
               MPI_Send(&j,1,MPI_INT,k,1001,PETSC_COMM_WORLD);
               MPI_Send(&ReValue,1,MPI_DOUBLE,k,1002,PETSC_COMM_WORLD);
               MPI_Send(&ImValue,1,MPI_DOUBLE,k,1003,PETSC_COMM_WORLD);
               k++;
            }
 
            j++;
         }
         i++;
      }
   } else {
      int i=1;
      while (i <= SportCount) {
         Result *row=get_col(iteration,frequency,i);
         int j=1;
         while (j <= SportCount) {
            MPI_Recv(&irow,1,MPI_INT,0,1000,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(&icol,1,MPI_INT,0,1001,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(&ReValue,1,MPI_DOUBLE,0,1002,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(&ImValue,1,MPI_DOUBLE,0,1003,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
            Svalue=complex<double>(ReValue,ImValue);
            row->set_S(Svalue,icol,irow);
            j++;
         }
         i++;
      }
   }

   return ierr;
}

// k must be destroyed elsewhere
PetscErrorCode  ResultDatabase::get_k (int iteration, double frequency, bool normalize, Mat *k)
{
   PetscErrorCode ierr=0;
   complex<double> Zo;
   PetscScalar value;
   PetscInt low,high;

   int SportCount=get_SportCount(iteration,frequency);
   if (SportCount == 0) return 1;

   ierr=MatCreate(PETSC_COMM_WORLD,k); if (ierr) return ierr;
   ierr=MatSetType(*k,MATAIJ); if (ierr) return ierr;
   ierr=MatSetSizes(*k,PETSC_DECIDE,PETSC_DECIDE,SportCount,SportCount); if (ierr) return ierr;
   ierr=MatSeqAIJSetPreallocation(*k,1,NULL); if (ierr) return ierr;
   ierr=MatMPIAIJSetPreallocation(*k,1,NULL,1,NULL); if (ierr) return ierr;
   ierr=MatZeroEntries(*k); if (ierr) return ierr;

   MatGetOwnershipRange(*k,&low,&high);

   int i=low+1;
   while (i < high+1) {
      if (normalize) Zo=get_normalizeZo(iteration,frequency,i);
      else Zo=get_Zo(iteration,frequency,i);
      ierr=MatSetValue(*k,i-1,i-1,sqrt(Zo),INSERT_VALUES); if (ierr) return ierr;
      i++;
   }

   ierr=MatAssemblyBegin(*k,MAT_FINAL_ASSEMBLY); if (ierr) return ierr;
   ierr=MatAssemblyEnd(*k,MAT_FINAL_ASSEMBLY); if (ierr) return ierr;

   return ierr;
}

PetscErrorCode ResultDatabase::S2Z (int iteration, double frequency)
{
   PetscErrorCode ierr=0;
   bool renormalize=false;
   PetscInt low,high;

   // to become I+S
   Mat IpS;
   ierr=get_S(iteration,frequency,&IpS,false,false); if (ierr) return ierr;

   // I+S
   MatGetOwnershipRange(IpS,&low,&high);
   int i=low;
   while (i < high) {
      ierr=MatSetValue(IpS,i,i,1.0,ADD_VALUES); if (ierr) return ierr;
      i++;
   }

   ierr=MatAssemblyBegin(IpS,MAT_FINAL_ASSEMBLY); if (ierr) return ierr;
   ierr=MatAssemblyEnd(IpS,MAT_FINAL_ASSEMBLY); if (ierr) return ierr;

   // to become I-S
   Mat ImS;
   ierr=get_S(iteration,frequency,&ImS,true,false); if (ierr) return ierr;

   // I-S 
   MatGetOwnershipRange(ImS,&low,&high);
   i=low;
   while (i < high) {
      ierr=MatSetValue(ImS,i,i,1.0,ADD_VALUES); if (ierr) return ierr;
      i++;
   }

   ierr=MatAssemblyBegin(ImS,MAT_FINAL_ASSEMBLY); if (ierr) return ierr;
   ierr=MatAssemblyEnd(ImS,MAT_FINAL_ASSEMBLY); if (ierr) return ierr;

   // k 
   Mat k;
   ierr=get_k(iteration,frequency,renormalize,&k); if (ierr) return ierr;

   // (I-S)^-1
   ierr=MatInvert(&ImS,0); if (ierr) return ierr;

   // calculate Z
   Mat C,D;
   ierr=MatMatMult(k,IpS,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&C); if (ierr) return ierr;
   ierr=MatMatMult(C,ImS,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&D); if (ierr) return ierr;
   ierr=MatDestroy(&C); if (ierr) return ierr;
   ierr=MatMatMult(D,k,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&C); if (ierr) return ierr;

   // save it in the S data location
   ierr=set_S(iteration,frequency,&C); if (ierr) return ierr;

   // cleanup
   ierr=MatDestroy(&IpS); if (ierr) return ierr;
   ierr=MatDestroy(&ImS); if (ierr) return ierr;
   ierr=MatDestroy(&C); if (ierr) return ierr;
   ierr=MatDestroy(&D); if (ierr) return ierr;
   ierr=MatDestroy(&k); if (ierr) return ierr;

   return ierr;
}

PetscErrorCode ResultDatabase::Z2S (int iteration, double frequency, bool renormalize)
{
   PetscErrorCode ierr=0;

   // k
   Mat k;
   ierr=get_k(iteration,frequency,renormalize,&k); if (ierr) return ierr;

   // k^(-1)
   Mat ki;
   ierr=get_k(iteration,frequency,renormalize,&ki); if (ierr) return ierr;
   ierr=MatInvert(&ki,1); if (ierr) return ierr;

   // Z
   Mat Z;
   ierr=get_S(iteration,frequency,&Z,false,true); if (ierr) return ierr;

   Mat Zkipk;
   ierr=MatMatMult(Z,ki,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Zkipk); if (ierr) return ierr;
   ierr=MatAXPY(Zkipk,1,k,DIFFERENT_NONZERO_PATTERN); if (ierr) return ierr;
   ierr=MatInvert(&Zkipk,0); if (ierr) return ierr;

   Mat Zkimk;
   ierr=MatMatMult(Z,ki,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Zkimk); if (ierr) return ierr;
   ierr=MatAXPY(Zkimk,-1,k,DIFFERENT_NONZERO_PATTERN); if (ierr) return ierr;

   // calculate S
   Mat C;
   ierr=MatMatMult(Zkipk,Zkimk,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&C); if (ierr) return ierr;

   // save it in the S data location
   ierr=set_S(iteration,frequency,&C); if (ierr) return ierr;

   // cleanup
   ierr=MatDestroy(&Z); if (ierr) return ierr;
   ierr=MatDestroy(&Zkipk); if (ierr) return ierr;
   ierr=MatDestroy(&Zkimk); if (ierr) return ierr;
   ierr=MatDestroy(&C); if (ierr) return ierr;
   ierr=MatDestroy(&k); if (ierr) return ierr;
   ierr=MatDestroy(&ki); if (ierr) return ierr;

   return ierr;
}

PetscErrorCode ResultDatabase::renormalize (int iteration, double frequency)
{
   PetscErrorCode ierr=0;
   ierr=S2Z (iteration,frequency); if (ierr) return ierr;
   ierr=Z2S (iteration,frequency,true); if (ierr) return ierr;

   return ierr;
}

void ResultDatabase::print ()
{
   long unsigned int i=0;
   while (i < results.size()) {
      results[i]->print();
      i++;
   }
}

void ResultDatabase::save_as_test (struct projectData *projData)
{
   ofstream out;

   stringstream ssTests;
   ssTests << projData->project_name << "_prototype_test_cases.csv";

   out.open(ssTests.str().c_str(),ofstream::out);
   int casenumber=0;

   if (out.is_open()) {
      out << "# ResultDatabase::save_as_test" << endl;

      int drivingSport=1;
      while (drivingSport <= SportCount) {
         long unsigned int i=0;
         while (i < unique_frequencies.size()) {
            Result *result=this->get_Result(unique_frequencies[i],drivingSport);  // gets the active result
            if (result) result->save_as_test (&out, projData->project_name, i, &casenumber);
            else PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Failed to find a result at frequency %g.\n",unique_frequencies[i]);
            i++;
         }
         drivingSport++;
      }

      out.close();
   } else {
       PetscPrintf(PETSC_COMM_WORLD,"ERROR3130: Failed to open file \"%s\" for writing.\n",ssTests.str().c_str());
   }
}

bool ResultDatabase::loadCSV (const char *filename)
{
   bool fail=false;
   string touchstone_format="";
   string frequency_unit="";
   //int number_of_frequencies=-1;
   int number_of_ports=-1;

   bool load_touchstone_format=false;
   bool load_frequency_unit=false;
   //bool load_number_of_frequencies=false;
   bool load_number_of_ports=false;

   // load header information

   ifstream CSV;
   CSV.open(filename,ifstream::in);
   if (CSV.is_open()) {
      string line;
      while (getline(CSV,line)) {
         stringstream ssLine(line);
         string value;
         while (std::getline(ssLine,value,',')) {
            if (load_touchstone_format) {touchstone_format=value; load_touchstone_format=false;}
            if (load_frequency_unit) {frequency_unit=value; load_frequency_unit=false;}
            //if (load_number_of_frequencies) {number_of_frequencies=stoi(value); load_number_of_frequencies=false;}
            if (load_number_of_ports) {number_of_ports=stoi(value); load_number_of_ports=false;}

            if (value.compare("#Touchstone format") == 0) load_touchstone_format=true;
            if (value.compare("#frequency unit") == 0) load_frequency_unit=true;
            //if (value.compare("#number of frequencies") == 0) load_number_of_frequencies=true;
            if (value.compare("#number of ports") == 0) load_number_of_ports=true;
         }

         ssLine.str("");
         ssLine.clear();
      }
      CSV.close();
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3131: Unable to open file \"%s\" for reading.\n",filename);
      fail=true;
   }

   // load S-parameters

   CSV.open(filename,ifstream::in);
   if (CSV.is_open()) {
      string line;
      while (getline(CSV,line)) {
         Result *newResult=new Result();
         if (newResult->extractS (line,frequency_unit,number_of_ports)) {
            delete newResult;
         } else {
            results.push_back(newResult);
         }
      }
      CSV.close();
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR3132: Unable to open file \"%s\" for reading.\n",filename);
      fail=true;
   }

   SportCount=number_of_ports;

   return fail;
}

ResultDatabase::~ResultDatabase ()
{
   long unsigned int i=0;
   while (i < results.size()) {
      delete results[i];
      i++;
   }
}



