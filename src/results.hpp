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


#ifndef RESULTS_H
#define RESULTS_H

#include "mfem.hpp"
#include "petscsys.h"
#include <petsc.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <chrono>
#include "fem3D.hpp"

extern "C" PetscErrorCode MatInvert (Mat *, int);

using namespace std;
using namespace mfem;

class fem3D;

class Sparm
{
   private:
      complex<double> S;
      int portIn;            // column
      int portOut;           // row
   public:
      Sparm (complex<double>, int, int);
      complex<double> get_S () {return S;}
      bool is_match (int, int);
      void set_S (complex<double> S_) {S=S_;}
      void print ();
      void save_as_test (ofstream *, const char *, int, int *, double, double, double, double, double);
};

class SparmList
{
   private:
      vector<Sparm *> parmList;
   public:
      ~SparmList ();
      void push (complex<double>, int, int);
      void set_S (complex<double>, int, int);
      complex<double> get_S (int, int);
      void print ();
      void save_as_test (ofstream *, const char *, int, int *, double, double, double, double, double);
};

class Result
{
   private:
      int iteration;                       // starts at 1
      double frequency;                    // Hz
      int drivingSport;                    // S-paramter port number for the driving port
      SparmList S;                         // S-parameter column for the driving port
      complex<double> Zo;                  // impedance for the driving port
      double normalizeZo;                  // impedance for normalization of the S-parameter matrix
      bool isNormalized;
      bool active;                         // false means that a frequency was re-calculated and a newer result is available
      PetscReal EfieldError,HfieldError;
      PetscInt EfieldConverged,HfieldConverged;
      PetscInt meshSize;
      PetscInt matrixSize;                 // global true dof size of the finite element space (keep the matrixSize variable name to align with fem3D.hpp)
      PetscInt sparseWidth;                // allocated width of the main compuational sparse matrices; memory usage is then matrixSize*sparseWidth*16
      double shortestPerWavelength;        // min edge length / wavelength in across all materials
      double longestPerWavelength;         // max edge length / wavelength in across all materials
      double maxRelativeError;             // iterative solution: error from prior iteration across all ports (multiple results)
      bool hasIterations=false;            //                     indicates that this is part of an iterative solution
      bool isConverged=false;              //                     indicates if the iterative solution converged
      double maxReflection;                // from the non-driven ports
      double solve_time=0;
      double fem_setup_time=0;
      double mesh_error_time=0;
      double refine_time=0;
      double magLimitdB=-150;              // breakover point for "equal" vs. "lessthan" comparisons for results
      double equalMagLimit=1e-12;
      double argLimitdeg=0.1;
      double equalArgLimit=1e-12;
   public:
      Result();
      ~Result();
      void set (int, double, double, double, int, double, fem3D *, ParMesh *);
      void set_iteration (int iteration_) {iteration=iteration_;}
      void set_frequency (double frequency_) {frequency=frequency_;}
      void set_drivingSport (int drivingSport_) {drivingSport=drivingSport_;}
      void push_S (complex<double> s_, int portOut_, int portIn_) {S.push(s_,portOut_,portIn_);}  
      void set_S (complex<double>, int, int);
      void set_Zo (complex<double> Zo_) {Zo=Zo_;}
      void set_active() {active=true;}
      void set_inactive() {active=false;}
      void set_maxRelativeError (double maxRelativeError_) {maxRelativeError=maxRelativeError_;}
      void set_hasIterations (bool hasIterations_) {hasIterations=hasIterations_;}
      void set_isConverged (bool isConverged_) {isConverged=isConverged_;}
      void set_maxReflection (double maxReflection_) {maxReflection=maxReflection_;}
      int get_iteration() {return iteration;}
      double get_frequency() {return frequency;}
      int get_drivingSport() {return drivingSport;}
      bool get_isNormalized() {return isNormalized;}
      complex<double> get_Zo () {return Zo;}
      complex<double> get_normalizeZo () {return normalizeZo;}
      complex<double> get_S (int portOut, int portIn) {return S.get_S(portOut,portIn);}
      void set_isNormalized() {isNormalized=true;}
      bool is_active() {return active;}
      void save (ostream *, struct projectData *, int);
      void saveFormatted (ostream *, double, double, double, double, double *);
      void saveCSV (ostream *, struct projectData *, double, int);
      void save_result_component (ofstream *, const char *, int, int *, const char *, long unsigned int);
      bool get_hasIterations () {return hasIterations;}
      bool get_isConverged () {return isConverged;}
      double get_solve_time () {return solve_time;}
      double get_fem_setup_time () {return fem_setup_time;}
      double get_refine_time () {return refine_time;}
      double get_mesh_error_time () {return mesh_error_time;}
      void set_solve_time (double t) {solve_time=t;}
      void set_fem_setup_time (double t) {fem_setup_time=t;}
      void set_refine_time (double t) {refine_time=t;}
      void set_mesh_error_time (double t) {mesh_error_time=t;}
      void print ();
      void save_as_test (ofstream *, const char *, int, int *);
      bool extractS (string, string, int);
};

class ResultDatabase
{
   private:
      vector<Result *> results;
      vector<double> unique_frequencies;  // sorted in increasing order
      int SportCount;
      double solve_time=0;
      double tol=1e-12;
   public:
      ResultDatabase(){}
      ~ResultDatabase();
      void push(Result *);
      Result* get_Result (double, int, int);
      Result* get_Result (double, int);
      Result* get_Result (double);
      double calculate_maxRelativeError (double, int);
      void set_maxRelativeError (double, double, int);
      void set_hasIterations (bool, double, int);
      void set_isConverged (bool, double, int);
      void set_SportCount (int SportCount_) {SportCount=SportCount_;}
      void save (ostream *, struct projectData *);
      bool save (struct projectData *);
      void saveFormatted(ostream *);
      bool saveFormatted(struct projectData *);
      void saveCSV (ostream *, struct projectData *, int, bool);
      bool saveCSV (struct projectData *, int, bool);
      bool saveTouchstone (struct projectData *, int);
      void set_refine_time (double, double, int);
      bool hasIterations ();
      bool isAllConverged ();
      bool isSequentialConverged (struct projectData *, double);
      void get_totalSolveTime (int, double *);
      void get_totalMeshErrorTime (int, double *);
      void get_totalFEMSetupTime (int, double *);
      void get_totalRefineTime (int, double *);
      void set_solve_time (double t) {solve_time=t;}
      double get_solve_time () {return solve_time;}
      int get_lastIteration (double);
      int get_SportCount (int, double);
      Result* get_col (int, double, int);
      complex<double> get_Zo (int, double, int);
      complex<double> get_normalizeZo (int, double, int);
      PetscErrorCode get_S (int, double, Mat *, bool, bool);
      PetscErrorCode get_k (int, double, bool, Mat *);
      PetscErrorCode set_S (int, double, Mat *);
      PetscErrorCode S2Z (int, double);
      PetscErrorCode Z2S (int, double, bool);
      PetscErrorCode renormalize (int, double);
      void print();
      void save_as_test(struct projectData *);
      bool loadCSV (const char *);
};

#endif

