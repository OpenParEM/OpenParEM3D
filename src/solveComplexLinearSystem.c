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

#include "solveComplexLinearSystem.h"

double current_residual;

void showReason (KSPConvergedReason reason) {
   if (reason < 0) PetscPrintf(PETSC_COMM_WORLD," NOT CONVERGED: ");
   else PetscPrintf(PETSC_COMM_WORLD," Converged: ");

   if (reason == KSP_CONVERGED_ITERATING) PetscPrintf(PETSC_COMM_WORLD,"CONVERGED_ITERATING");
   if (reason == KSP_CONVERGED_RTOL_NORMAL) PetscPrintf(PETSC_COMM_WORLD,"RTOL_NORMAL");
   if (reason == KSP_CONVERGED_ATOL_NORMAL) PetscPrintf(PETSC_COMM_WORLD,"ATOL_NORMAL");
   if (reason == KSP_CONVERGED_RTOL) PetscPrintf(PETSC_COMM_WORLD,"RTOL");
   if (reason == KSP_CONVERGED_ATOL) PetscPrintf(PETSC_COMM_WORLD,"ATOL");
   if (reason == KSP_CONVERGED_ITS) PetscPrintf(PETSC_COMM_WORLD,"ITS");
   if (reason == KSP_CONVERGED_NEG_CURVE) PetscPrintf(PETSC_COMM_WORLD,"NEG_CURVE");
   if (reason == KSP_CONVERGED_STEP_LENGTH) PetscPrintf(PETSC_COMM_WORLD,"STEP_LENGTH");
   if (reason == KSP_CONVERGED_HAPPY_BREAKDOWN) PetscPrintf(PETSC_COMM_WORLD,"HAPPY_BREAKDOWN");
   if (reason == KSP_DIVERGED_NULL) PetscPrintf(PETSC_COMM_WORLD,"NULL");
   if (reason == KSP_DIVERGED_ITS) PetscPrintf(PETSC_COMM_WORLD,"ITS - Hit the iteration limit.");
   if (reason == KSP_DIVERGED_DTOL) PetscPrintf(PETSC_COMM_WORLD,"DTOL");
   if (reason == KSP_DIVERGED_BREAKDOWN) PetscPrintf(PETSC_COMM_WORLD,"BREAKDOWN - Generic breakdown during solution.");
   if (reason == KSP_DIVERGED_BREAKDOWN_BICG) PetscPrintf(PETSC_COMM_WORLD,"BREAKDOWN_BICG");
   if (reason == KSP_DIVERGED_NONSYMMETRIC) PetscPrintf(PETSC_COMM_WORLD,"NONSYMMETRIC");
   if (reason == KSP_DIVERGED_INDEFINITE_PC) PetscPrintf(PETSC_COMM_WORLD,"INDEFINITE_PC");
   if (reason == KSP_DIVERGED_NANORINF) PetscPrintf(PETSC_COMM_WORLD,"NANORINF - Detected a nan or inf.");
   if (reason == KSP_DIVERGED_INDEFINITE_MAT) PetscPrintf(PETSC_COMM_WORLD,"INDEFINITE_MAT");
   if (reason == KSP_DIVERGED_PC_FAILED) PetscPrintf(PETSC_COMM_WORLD,"PC_FAILED - Could not build or use the requested preconditioner.");

   PetscPrintf(PETSC_COMM_WORLD,"\n");
}

PetscErrorCode monitorEM3D(KSP ksp, PetscInt its, PetscReal residual, void *vf)
{
   if (residual < current_residual/10) {
      PetscPrintf (PETSC_COMM_WORLD,"|               iteration=%ld   residual=%g\n",its,residual);
      current_residual=residual;
   }
   PetscFunctionReturn(0);
}

PetscErrorCode eliminatePEC (Mat *A, PetscInt nPEC, PetscInt *PEC) {
   PetscErrorCode ierr=0;
   ierr=MatZeroRowsColumns(*A,nPEC,PEC,1.0,NULL,NULL);
   return ierr;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// solveComplexLinearSystem
//
// on input:
//    PortDof is list of dof indices for the driven port 
//    x is the dof strengths over the 3D space including just the dof strengths for the driven port
//
// on output;
//    x is the dof strengths for all dofs in the 3D space 
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

PetscErrorCode solveComplexLinearSystem (const char *directory, struct projectData *projData, PetscMPIInt rank,
    PetscInt nPortDof, PetscInt *PortDof, Mat *A, Vec *x, PetscInt *matrixSize, PetscReal *rnorm, PetscInt *converged)
{
   PetscErrorCode ierr=0;
   PetscInt i,maxits,m,n,*ivals;;
   PetscReal rtol,atol,dtol;
   PetscScalar *vals;
   KSP ksp;
   PC pc;
   KSPConvergedReason reason;

   MatGetSize(*A,&m,&n);
   *matrixSize=m;

   // calculate the right-hand-side using the given x, which has the known dofs from the driven port
   Vec b;
   MatCreateVecs(*A,NULL,&b);
   ierr=MatMult(*A,*x,b); if (ierr) return 1;
   ierr=VecScale(b,-1); if (ierr) return 2;

   // set the known dof values from x in b

   ierr=PetscMalloc(nPortDof*sizeof(PetscInt),&ivals); if (ierr) return 3;
   i=0;
   while (i < nPortDof) {
      ivals[i]=PortDof[i];
      i++;
   }

   ierr=PetscMalloc(nPortDof*sizeof(PetscScalar),&vals); if (ierr) return 4;

   ierr=VecGetValues(*x,nPortDof,ivals,vals); if (ierr) return 5;
   ierr=VecSetValues(b,nPortDof,ivals,vals,INSERT_VALUES); if (ierr) return 6;

   ierr=VecAssemblyBegin(b); if (ierr) return 7;
   ierr=VecAssemblyEnd(b); if (ierr) return 8;

   ierr=PetscFree(ivals); if (ierr) return 9;
   ierr=PetscFree(vals); if (ierr) return 10;

   // eliminate the known port dofs from A
   ierr=MatZeroRowsColumns(*A,nPortDof,PortDof,1.0,NULL,NULL); if (ierr) return 11;

   // solve Ax=b

   i=0;
   while (i < 5) {

      ierr=KSPCreate(PETSC_COMM_WORLD, &ksp); if (ierr) return 12;
      ierr=KSPSetType(ksp,KSPGMRES); if (ierr) return 13;
      ierr=KSPSetOperators(ksp,*A,*A); if (ierr) return 14;
      ierr=KSPGetPC(ksp,&pc); if (ierr) return 15;
      ierr=PCFactorSetShiftType(pc,MAT_SHIFT_NONZERO); if (ierr) return 16;
      ierr=KSPSetFromOptions(ksp); if (ierr) return 17;

      maxits=projData->solution_iteration_limit;       // PETSc default = 1e4
      rtol=projData->solution_3D_tolerance;            // PETSc default = 1e-5
      //atol=projData->solution_3D_tolerance;          // PETSC default = 1e-50
      atol=1e-50;                                      // PETSC default = 1e-50
      dtol=1e5;                                        // PETSC default = 1e5
      ierr=KSPSetTolerances(ksp,rtol,atol,dtol,maxits); if (ierr) return 18;

      current_residual=DBL_MAX;
      if (projData->output_show_iterations) KSPMonitorSet(ksp,monitorEM3D,NULL,NULL);

      if (i == 0) {
         ierr=PCSetType(pc,PCCHOLESKY); if (ierr) return 19;
         PetscPrintf(PETSC_COMM_WORLD,"|            using Cholesky preconditioner ...\n");
      } else if (i == 1) {
         ierr=PCSetType(pc,PCCHOLESKY); if (ierr) return 20;
         PetscPrintf(PETSC_COMM_WORLD,"|            re-trying ...\n");
      } else if (i == 2) {
         ierr=PCSetType(pc,PCLU); if (ierr) return 21;
         PetscPrintf(PETSC_COMM_WORLD,"|            re-trying using LU preconditioner ...\n");
      } else if (i == 3) {
         ierr=PCSetType(pc,PCBJACOBI); if (ierr) return 22;
         PetscPrintf(PETSC_COMM_WORLD,"|            re-trying using Jacobi preconditioner ...\n");
      } else {
         ierr=PCSetType(pc,PCNONE); if (ierr) return 23;
         PetscPrintf(PETSC_COMM_WORLD,"|            re-trying using no preconditioner ...\n");
      }

      // switch x to get the needed local distribution for the solve
      VecDestroy(x);
      MatCreateVecs(*A,NULL,x);

      // solve
      ierr=KSPSolve(ksp,b,*x); if (ierr) return 24;
      // get stats
      ierr=KSPGetResidualNorm(ksp,rnorm); if (ierr) return 25;
      ierr=KSPGetConvergedReason(ksp,&reason); if (ierr) return 26;
      ierr=KSPDestroy(&ksp); if (ierr) return 27;

      // exit criteria

      if (reason >= 0) {      // converged
         PetscPrintf(PETSC_COMM_WORLD,"|              ");
         showReason(reason);
         break;
      } else if (reason == KSP_DIVERGED_PC_FAILED) {
         if (i >= 2) {
            PetscPrintf(PETSC_COMM_WORLD,"|              ");
            showReason(reason);
         }
      } else {
         PetscPrintf(PETSC_COMM_WORLD,"|              ");
         showReason(reason);
      }

      if (reason == KSP_DIVERGED_ITS) break;   // break out if a preconditioner ran

      i++;
   }

   ierr=VecDestroy(&b); if (ierr) return 28;

   *converged=1;
   if (reason < 0) *converged=0;

   return 0;
}

PetscErrorCode buildHb (Mat *Q, Vec *x, Vec *b)
{
   PetscErrorCode ierr=0;

   ierr=VecDuplicate(*x,b); if (ierr) return 1;

   // rhs
   ierr=MatMult(*Q,*x,*b); if (ierr) return 2;

   return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// solveHfield
//
// on input:
//    x is the dofs solution from solveComplexLinearSystem
//
// on output;
//    hdof is the dof solution for the H field
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

PetscErrorCode solveHfield (const char *directory, struct projectData *projData, PetscMPIInt rank, 
                            Mat *P, Vec *b, Vec* hdofs, PetscReal *rnorm, PetscInt *converged)
{
   PetscErrorCode ierr=0;
   PetscInt i,maxits;
   PetscReal rtol,atol,dtol;
   KSP ksp;
   PC pc;
   KSPConvergedReason reason;

   i=0;
   while (i < 5) {

      ierr=KSPCreate(PETSC_COMM_WORLD, &ksp); if (ierr) return 1;
      ierr=KSPSetType(ksp,KSPGMRES); if (ierr) return 2;
      ierr=KSPSetOperators(ksp,*P,*P); if (ierr) return 3;
      ierr=KSPGetPC(ksp,&pc); if (ierr) return 4;
      ierr=PCFactorSetShiftType(pc,MAT_SHIFT_NONZERO); if (ierr) return 5;
      ierr=KSPSetFromOptions(ksp); if (ierr) return 6;

      maxits=projData->solution_iteration_limit;       // PETSc default = 1e4
      rtol=projData->solution_3D_tolerance;            // PETSc default = 1e-5
      //atol=projData->solution_3D_tolerance;          // PETSC default = 1e-50
      atol=1e-50;                                      // PETSC default = 1e-50
      dtol=1e5;                                        // PETSC default = 1e5
      ierr=KSPSetTolerances(ksp,rtol,atol,dtol,maxits); if (ierr) return 7;

      current_residual=DBL_MAX;
      if (projData->output_show_iterations) KSPMonitorSet(ksp,monitorEM3D,NULL,NULL);

      if (i == 0) {
         ierr=PCSetType(pc,PCCHOLESKY); if (ierr) return 8;
         PetscPrintf(PETSC_COMM_WORLD,"|            using Cholesky preconditioner ...\n");
      } else if (i == 1) {
         ierr=PCSetType(pc,PCCHOLESKY); if (ierr) return 9;
         PetscPrintf(PETSC_COMM_WORLD,"|            re-trying ...\n");
      } else if (i == 2) {
         ierr=PCSetType(pc,PCLU); if (ierr) return 10;
         PetscPrintf(PETSC_COMM_WORLD,"|            re-trying using LU preconditioner ...\n");
      } else if (i == 3) {
         ierr=PCSetType(pc,PCBJACOBI); if (ierr) return 11;
         PetscPrintf(PETSC_COMM_WORLD,"|            re-trying using Jacobi preconditioner ...\n");
      } else {
         ierr=PCSetType(pc,PCNONE); if (ierr) return 12;
         PetscPrintf(PETSC_COMM_WORLD,"|            re-trying using no preconditioner ...\n");
      }

      // solve
      ierr=KSPSolve(ksp,*b,*hdofs); if (ierr) return 13;

      // get stats
      ierr=KSPGetConvergedReason(ksp,&reason); if (ierr) return 14;
      ierr=KSPGetResidualNorm(ksp,rnorm); if (ierr) return 15;
      ierr=KSPDestroy(&ksp); if (ierr) return 16;

      // exit criteria

      if (reason >= 0) {      // converged
         PetscPrintf(PETSC_COMM_WORLD,"|              ");
         showReason(reason);
         break;
      } else if (reason == KSP_DIVERGED_PC_FAILED) {
         if (i >= 2) {
            PetscPrintf(PETSC_COMM_WORLD,"|              ");
            showReason(reason);
         }
      } else {
         PetscPrintf(PETSC_COMM_WORLD,"|              ");
         showReason(reason);
      }

      if (reason == KSP_DIVERGED_ITS) break;   // break out if a preconditioner ran

      i++;
   }

   *converged=1;
   if (reason < 0) *converged=0;

   return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// miscellaneous routines
//
///////////////////////////////////////////////////////////////////////////////////////////////////////


PetscErrorCode MatInvert (Mat *A, int isDiagonal)
{
   PetscErrorCode ierr=0;
   PetscInt m,n,i,j,high,low;
   //PetscScalar det,A11,A12,A21,A22;
   PetscScalar value,*idx;
   KSP ksp;
   PC pc;
   Mat Ap,Ai;
   Vec b,x;

   if (A == NULL) return 1;

   ierr=MatGetSize(*A,&m,&n); if (ierr) return 2;
   if (m != n) return 3;

   if (isDiagonal || m == 1) {

      ierr=MatGetOwnershipRange(*A,&low,&high); if (ierr) return 4;

      ierr=PetscMalloc((high-low)*sizeof(PetscScalar),&idx);  if (ierr) return 5;
      i=low;
      while (i < high) {
         ierr=MatGetValue(*A,i,i,&value); if (ierr) return 6;
         idx[i-low]=value;
         i++;
      }

      i=low;
      while (i < high) {
         ierr=MatSetValue(*A,i,i,1/idx[i-low],INSERT_VALUES); if (ierr) return 7;
         i++;
      }

      ierr=PetscFree(idx); if (ierr) return 8;

/* ToDo: parallelize the m=2 case
   } else if (m == 2) {

      // analytic inversion

      ierr=MatGetValue(*A,0,0,&A11); if (ierr) return 9;
      ierr=MatGetValue(*A,0,1,&A12); if (ierr) return 10;
      ierr=MatGetValue(*A,1,0,&A21); if (ierr) return 11;
      ierr=MatGetValue(*A,1,1,&A22); if (ierr) return 12;

      det=A11*A22-A21*A12;
      if (det == 0) return 13;

      ierr=MatSetValue(*A,0,0,A22/det,INSERT_VALUES); if (ierr) return 14;
      ierr=MatSetValue(*A,0,1,-A12/det,INSERT_VALUES); if (ierr) return 15;
      ierr=MatSetValue(*A,1,0,-A21/det,INSERT_VALUES); if (ierr) return 16;
      ierr=MatSetValue(*A,1,1,A11/det,INSERT_VALUES); if (ierr) return 17;
*/
   } else {

      // iterative inversion: inefficient, but ok for usage here with small S-parameter matrices

      ierr=KSPCreate(PETSC_COMM_WORLD, &ksp); if (ierr) return 18;
      ierr=KSPSetType(ksp,KSPGMRES); if (ierr) return 19;
      ierr=KSPGetPC(ksp,&pc); if (ierr) return 20;
      ierr=PCSetType(pc,PCLU); if (ierr) return 21;                    // works fine, but a little slower

      ierr=MatDuplicate(*A,MAT_DO_NOT_COPY_VALUES,&Ap); if (ierr) return 23;
      ierr=MatDuplicate(*A,MAT_DO_NOT_COPY_VALUES,&Ai); if (ierr) return 24;
      ierr=MatCreateVecs(*A,NULL,&x); if (ierr) return 25;
      ierr=MatCreateVecs(*A,NULL,&b); if (ierr) return 26;

      i=0;
      while (i < m) {
         ierr=MatCopy(*A,Ap,SAME_NONZERO_PATTERN); if (ierr) return 27;
         ierr=KSPSetOperators(ksp,Ap,Ap); if (ierr) return 28;

         ierr=VecZeroEntries(b); if (ierr) return 29;

         ierr=VecGetOwnershipRange(b,&low,&high); if (ierr) return 30;
         if (i >= low && i < high) {
            ierr=VecSetValue(b,i,1,INSERT_VALUES); if (ierr) return 31;
         }

         ierr=VecAssemblyBegin(b); if (ierr) return 32;
         ierr=VecAssemblyEnd(b); if (ierr) return 33;

         ierr=VecAssemblyBegin(x); if (ierr) return 34;
         ierr=VecAssemblyEnd(x); if (ierr) return 35;

         ierr=KSPSolve(ksp,b,x); if (ierr) return 36;

         j=0;
         while (j < m) {
            if (j >= low && j < high) {
               ierr=VecGetValues(x,1,&j,&value); if (ierr) return 37;
               ierr=MatSetValue(Ai,j,i,value,INSERT_VALUES); if (ierr) return 38;
            }
            j++;
         }

         i++;
      }

      ierr=KSPDestroy(&ksp); if (ierr) return 39;

      ierr=MatAssemblyBegin(Ai,MAT_FINAL_ASSEMBLY); if (ierr) return 40;
      ierr=MatAssemblyEnd(Ai,MAT_FINAL_ASSEMBLY); if (ierr) return 41;
   
      ierr=MatCopy(Ai,*A,SAME_NONZERO_PATTERN); if (ierr) return 42;

      ierr=MatDestroy(&Ap); if (ierr) return 43;
      ierr=MatDestroy(&Ai); if (ierr) return 44;
      ierr=VecDestroy(&x); if (ierr) return 45;
      ierr=VecDestroy(&b); if (ierr) return 46;
   }

   ierr=MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY); if (ierr) return 47;
   ierr=MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY); if (ierr) return 48;

   return ierr;
}
