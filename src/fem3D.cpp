//////////////////////////////////////////////////////////////////////////
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

#include "fem3D.hpp"

extern "C" PetscErrorCode eliminatePEC (Mat *, PetscInt, PetscInt *);
extern "C" PetscErrorCode solveComplexLinearSystem (const char *, struct projectData *, PetscMPIInt,
    PetscInt, PetscInt *, Mat *, Vec *, PetscInt *, PetscReal *, PetscInt *);
extern "C" PetscErrorCode buildHb (Mat *, Vec *, Vec *);
extern "C" PetscErrorCode solveHfield (const char *, struct projectData *, PetscMPIInt, Mat *, Vec *, Vec*, PetscReal *, PetscInt *);
extern "C" PetscErrorCode hypre_ParCSRMatrixToMat(hypre_ParCSRMatrix *, Mat *, PetscInt, int, int, int);
extern "C" PetscErrorCode hypre_getSparseWidth (hypre_ParCSRMatrix *, PetscInt *);

int getGlobalNE (ParMesh *pmesh)
{
   int localNE=pmesh->GetNE();
   int globalNE=0;

   MPI_Allreduce (&localNE,&globalNE,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);

   return globalNE;
}

int getGlobalNBE (ParMesh *pmesh)
{
   int localNBE=pmesh->GetNBE();
   int globalNBE=0;

   MPI_Allreduce (&localNBE,&globalNBE,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);

   return globalNBE;
}

//--------------------------------------------------------------------------------------------------------------------
// fem3D
//--------------------------------------------------------------------------------------------------------------------

void fem3D::set_data (ParMesh **pmesh_, struct projectData *projData_, double frequency_)
{
   pmesh=pmesh_;
   projData=projData_;
   frequency=frequency_;
}

void fem3D::build_fe_spaces ()
{
   fec_ND=new ND_FECollection(projData->mesh_order,(*pmesh)->Dimension());
   fespace_ND=new ParFiniteElementSpace(*pmesh,fec_ND);

   fec_RT=new RT_FECollection(projData->mesh_order,(*pmesh)->Dimension());
   fespace_RT=new ParFiniteElementSpace(*pmesh,fec_RT);

   fec_H1=new H1_FECollection(projData->mesh_order,(*pmesh)->Dimension());
   fespace_H1=new ParFiniteElementSpace(*pmesh,fec_H1);

   fec_L2=new L2_FECollection(projData->mesh_order,(*pmesh)->Dimension(),BasisType::GaussLobatto);
   fespace_L2=new ParFiniteElementSpace(*pmesh,fec_L2);
}

void fem3D::build_PEC_dofs ()
{
   Array<int> ess_bdr_PEC;
   ess_bdr_PEC.SetSize((*pmesh)->bdr_attributes.Max());
   ess_bdr_PEC=0;       // disable all
   ess_bdr_PEC[0]=1;    // enable PEC

   Array<int> ess_dofs_ND;
   fespace_ND->GetEssentialTrueDofs(ess_bdr_PEC, ess_dofs_ND);
   HYPRE_BigInt *offset=fespace_ND->GetTrueDofOffsets();

   nPEC=ess_dofs_ND.Size();
   PetscMalloc(nPEC*sizeof(PetscInt),&PEC);

   int i=0;
   while (i < ess_dofs_ND.Size()) {
      PEC[i]=ess_dofs_ND[i]+offset[0];
      i++;
   }
}

// A is a complex matrix that is matrixSize x sparseWidth.
// Building A costs 50% more memory ReA or ReB must be allocated at the same time as A.
bool fem3D::build_A (BoundaryDatabase *boundaryDatabase, MaterialDatabase *materialDatabase, double temperature,
                     PWConstCoefficient *neg_ko2_Re_er, PWConstCoefficient *neg_ko2_Im_er,PWConstCoefficient *Inv_mur)
{
   vector<Array<int> *> borderAttributeList;
   vector<ConstantCoefficient *> alphaConstList,betaConstList,ZconstList;

   PetscMPIInt rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   // real

   ParMixedBilinearForm *pmblfReA=new ParMixedBilinearForm(fespace_ND,fespace_ND);
   pmblfReA->AddDomainIntegrator(new CurlCurlIntegrator (*Inv_mur));
   pmblfReA->AddDomainIntegrator(new VectorFEMassIntegrator (*neg_ko2_Re_er));
   boundaryDatabase->addMassPortIntegrators(*pmesh,pmblfReA,Inv_mur,borderAttributeList,alphaConstList,betaConstList,true);
   pmblfReA->Assemble();
   pmblfReA->Finalize();
   HypreParMatrix *ReA=pmblfReA->ParallelAssemble();

   long unsigned int i=0;
   while (i < borderAttributeList.size()) {
      if (borderAttributeList[i]) delete borderAttributeList[i];
      if (alphaConstList[i]) delete alphaConstList[i];
      if (betaConstList[i]) delete betaConstList[i];
      i++;
   }
   borderAttributeList.clear();
   alphaConstList.clear();
   betaConstList.clear();

   delete pmblfReA;
   hypre_getSparseWidth((hypre_ParCSRMatrix *) *ReA,&sparseWidth);
   if (hypre_ParCSRMatrixToMat((hypre_ParCSRMatrix *) *ReA, &A, sparseWidth, 1, 0, 0)) return true;
   delete ReA;

   // imag

   ParMixedBilinearForm *pmblfImA=new ParMixedBilinearForm(fespace_ND,fespace_ND);
   pmblfImA->AddDomainIntegrator(new VectorFEMassIntegrator (*neg_ko2_Im_er));
   boundaryDatabase->addMassPortIntegrators(*pmesh,pmblfImA,Inv_mur,borderAttributeList,alphaConstList,betaConstList,false);
   boundaryDatabase->addMassImpedanceIntegrators(frequency,temperature,*pmesh,pmblfImA,materialDatabase,
                                                 borderAttributeList,ZconstList,false);
   pmblfImA->Assemble(); 
   pmblfImA->Finalize();
   HypreParMatrix *ImA=pmblfImA->ParallelAssemble();

   i=0;
   while (i < borderAttributeList.size()) {
      if (borderAttributeList[i]) delete borderAttributeList[i];
      i++;
   }

   i=0;
   while (i < alphaConstList.size()) {
      if (alphaConstList[i]) delete alphaConstList[i];
      if (betaConstList[i]) delete betaConstList[i];
      i++;
   }

   i=0;
   while (i < ZconstList.size()) {
      if (ZconstList[i]) delete ZconstList[i];
      i++;
   }

   delete pmblfImA;
   if (hypre_ParCSRMatrixToMat((hypre_ParCSRMatrix *) *ImA, &A, sparseWidth, 0, 1, 1)) return true;
   delete ImA;

   return false;
}

void fem3D::build_P ()
{
   ConstantCoefficient one(1.0);
   ParBilinearForm *pblfReP=new ParBilinearForm(fespace_ND);
   pblfReP->AddDomainIntegrator(new VectorFEMassIntegrator(one));
   pblfReP->Assemble();
   pblfReP->Finalize();
   HypreParMatrix *ReP=pblfReP->ParallelAssemble();
   delete pblfReP;

   if (hypre_ParCSRMatrixToMat((hypre_ParCSRMatrix *) *ReP, &P, sparseWidth, 1, 0, 1)) {
      cout << "ASSERT: fem3D::build_P failed to build Mat P." << endl;
   }

   delete ReP;
}

void fem3D::build_Q (PWConstCoefficient *Inv_w_mu)
{
   ParMixedBilinearForm *pmblfReQ=new ParMixedBilinearForm(fespace_ND,fespace_ND);
   pmblfReQ->AddDomainIntegrator(new MixedVectorCurlIntegrator (*Inv_w_mu));
   pmblfReQ->Assemble();
   pmblfReQ->Finalize();

   HypreParMatrix *ReQ=pmblfReQ->ParallelAssemble();
   delete pmblfReQ;

   // anticipate solving PH=jQE and put Q in the imaginary position
   if (hypre_ParCSRMatrixToMat((hypre_ParCSRMatrix *) *ReQ, &Q, sparseWidth, 1, 1, 1)) {
      cout << "ASSERT: fem3D::build_Q failed to build Mat Q." << endl;
   }

   delete ReQ;
}

PetscErrorCode fem3D::build_X()
{
   PetscErrorCode ierr=0;
   ierr=MatCreateVecs(A,&X,NULL);
   return ierr;
}

PetscErrorCode fem3D::build_x()
{
   PetscErrorCode ierr=0;
   ierr=MatCreateVecs(A,&x,NULL);
   return ierr;
}

void fem3D::build_grids()
{
   if (gridReE) delete gridReE;
   if (gridImE) delete gridImE;
   if (gridReH) delete gridReH;
   if (gridImH) delete gridImH;

   if (gridReExH) delete gridReExH;
   if (gridImExH) delete gridImExH;

//   if (gridReEz) delete gridReEz;
//   if (gridImEz) delete gridImEz;
//   if (gridReHz) delete gridReHz;
//   if (gridImHz) delete gridImHz;

   gridReE=new ParGridFunction(fespace_ND);
   gridImE=new ParGridFunction(fespace_ND);
   gridReH=new ParGridFunction(fespace_ND);
   gridImH=new ParGridFunction(fespace_ND);

   gridReExH=new ParGridFunction(fespace_RT);
   gridImExH=new ParGridFunction(fespace_RT);

//   gridReEz=new ParGridFunction(fespace_L2);
//   gridImEz=new ParGridFunction(fespace_L2);
//   gridReHz=new ParGridFunction(fespace_L2);
//   gridImHz=new ParGridFunction(fespace_L2);
}

// fill out the dofs on the ports
void fem3D::build_portEssTdofLists(BoundaryDatabase *boundaryDatabase)
{
   boundaryDatabase->build_portEssTdofLists(fespace_ND,*pmesh);
}

PetscErrorCode fem3D::build_Xdofs()
{
   PetscErrorCode ierr=0;
   ierr=MatCreateVecs(A,&Xdofs,NULL);
   return ierr;
}

// x uses a local distribution on A (which is destroyed before here)
// e_re and e_im must use a local distribution on fespace_ND
// The two local spaces do not exactly line up, so they must be realigned from A to fespace_ND for e_re and e_im.
void fem3D::build_e_re_e_im ()
{
   PetscMPIInt rank,size;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   MPI_Comm_size(PETSC_COMM_WORLD, &size);

   // local distribution of fespace_ND - for e_re and e_im
   HYPRE_BigInt *fespace_ND_offsets=fespace_ND->GetTrueDofOffsets();
   int local_ND_size=fespace_ND_offsets[1]-fespace_ND_offsets[0];

   if (e_re) delete e_re;
   if (e_im) delete e_im;

   e_re=new Vector(local_ND_size);
   e_im=new Vector(local_ND_size);

   // local distribution of A - for x

   PetscInt i,low,high;
   VecGetOwnershipRange(x,&low,&high);
   int local_x_size=high-low;

   PetscInt *ixvals;
   PetscMalloc(local_x_size*sizeof(PetscInt),&ixvals);

   i=low;
   while (i < high) {
      ixvals[i-low]=i;
      i++;
   }

   PetscScalar *local_x_vals;
   PetscMalloc(local_x_size*sizeof(PetscScalar),&local_x_vals);

   VecGetValues(x,local_x_size,ixvals,local_x_vals);

   // global x

   int global_x_size=0;
   MPI_Allreduce (&local_x_size,&global_x_size,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);

   PetscScalar *global_x_vals;
   PetscMalloc(global_x_size*sizeof(PetscScalar),&global_x_vals);

   // collect x at rank 0

   if (rank == 0) {

      // local
      int i=0;
      while (i < local_x_size) {
         global_x_vals[low+i]=local_x_vals[i];
         i++;
      }

      // collected
      i=1;
      while (i < size) {
         int transfer_size=0;
         MPI_Recv(&transfer_size,1,MPI_INT,i,10000,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

         int k=0;
         while (k < transfer_size) {
            int location=0;
            double transfer_real=0;
            double transfer_imag=0;
            MPI_Recv(&location,1,MPI_INT,i,10001,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(&transfer_real,1,MPI_DOUBLE,i,10002,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(&transfer_imag,1,MPI_DOUBLE,i,10003,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

            global_x_vals[location]=transfer_real+PETSC_i*transfer_imag;

            k++;
         }
         i++;
      }
   } else {
      MPI_Send(&local_x_size,1,MPI_INT,0,10000,PETSC_COMM_WORLD);

      int i=0;
      while (i < local_x_size) {
         int location=low+i;
         double transfer_real=real(local_x_vals[i]);
         double transfer_imag=imag(local_x_vals[i]);
         MPI_Send(&location,1,MPI_INT,0,10001,PETSC_COMM_WORLD);
         MPI_Send(&transfer_real,1,MPI_DOUBLE,0,10002,PETSC_COMM_WORLD);
         MPI_Send(&transfer_imag,1,MPI_DOUBLE,0,10003,PETSC_COMM_WORLD);
         i++;
      }
   }

   // send global x to all ranks

   if (rank == 0) {
      int i=1;
      while (i < size) {
         int k=0;
         while (k < global_x_size) {
            double transfer_real=real(global_x_vals[k]);
            double transfer_imag=imag(global_x_vals[k]);
            MPI_Send(&transfer_real,1,MPI_DOUBLE,i,20001,PETSC_COMM_WORLD);
            MPI_Send(&transfer_imag,1,MPI_DOUBLE,i,20002,PETSC_COMM_WORLD);
            k++;
         }
         i++;
      }
   } else {
      int k=0;
      while (k < global_x_size) {
         double transfer_real=0;
         double transfer_imag=0;
         MPI_Recv(&transfer_real,1,MPI_DOUBLE,0,20001,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         MPI_Recv(&transfer_imag,1,MPI_DOUBLE,0,20002,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

         global_x_vals[k]=transfer_real+PETSC_i*transfer_imag;

         k++;
      }
   }

   // transfer data to e_re and e_im - completes the re-alignment

   int index=0;
   i=0;
   while (i < global_x_size) {
      if (i >= fespace_ND_offsets[0] && i < fespace_ND_offsets[1]) {
         e_re->Elem(index)=real(global_x_vals[i]);
         e_im->Elem(index)=imag(global_x_vals[i]);
         index++;
      }
      i++;
   }

   // clean up
   if (ixvals) {free(ixvals); ixvals=nullptr;}
   if (local_x_vals) {free(local_x_vals); local_x_vals=nullptr;}
   if (global_x_vals) {free(global_x_vals); global_x_vals=nullptr;}
}

void fem3D::buildEgrids (BoundaryDatabase *boundaryDatabase)
{
   gridReE->Distribute(*e_re);
   gridImE->Distribute(*e_im);
}

// see notes for build_e_re_e_im
void fem3D::build_h_re_h_im (Vec *hdofs)
{
   PetscMPIInt rank,size;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   MPI_Comm_size(PETSC_COMM_WORLD, &size);

   // local distribution of fespace_ND - for h_re and h_im
   HYPRE_BigInt *fespace_ND_offsets=fespace_ND->GetTrueDofOffsets();
   int local_ND_size=fespace_ND_offsets[1]-fespace_ND_offsets[0];

   if (h_re) delete h_re;
   if (h_im) delete h_im;

   h_re=new Vector(local_ND_size);
   h_im=new Vector(local_ND_size);

   // local distribution of A - for hdofs

   PetscInt i,low,high;
   VecGetOwnershipRange(*hdofs,&low,&high);
   int local_hdof_size=high-low;

   PetscInt *ixvals;
   PetscMalloc(local_hdof_size*sizeof(PetscInt),&ixvals);

   i=low;
   while (i < high) {
      ixvals[i-low]=i;
      i++;
   }

   PetscScalar *local_hdof_vals;
   PetscMalloc(local_hdof_size*sizeof(PetscScalar),&local_hdof_vals);

   VecGetValues(*hdofs,local_hdof_size,ixvals,local_hdof_vals);

   // global hdofs

   int global_hdof_size=0;
   MPI_Allreduce (&local_hdof_size,&global_hdof_size,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);

   PetscScalar *global_hdof_vals;
   PetscMalloc(global_hdof_size*sizeof(PetscScalar),&global_hdof_vals);

   // collect hdofs at rank 0

   if (rank == 0) {

      // local
      int i=0;
      while (i < local_hdof_size) {
         global_hdof_vals[low+i]=local_hdof_vals[i];
         i++;
      }

      // collected
      i=1;
      while (i < size) {
         int transfer_size=0;
         MPI_Recv(&transfer_size,1,MPI_INT,i,10000,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

         int k=0;
         while (k < transfer_size) {
            int location=0;
            double transfer_real=0;
            double transfer_imag=0;
            MPI_Recv(&location,1,MPI_INT,i,10001,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(&transfer_real,1,MPI_DOUBLE,i,10002,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(&transfer_imag,1,MPI_DOUBLE,i,10003,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

            global_hdof_vals[location]=transfer_real+PETSC_i*transfer_imag;

            k++;
         }
         i++;
      }
   } else {
      MPI_Send(&local_hdof_size,1,MPI_INT,0,10000,PETSC_COMM_WORLD);

      int i=0;
      while (i < local_hdof_size) {
         int location=low+i;
         double transfer_real=real(local_hdof_vals[i]);
         double transfer_imag=imag(local_hdof_vals[i]);
         MPI_Send(&location,1,MPI_INT,0,10001,PETSC_COMM_WORLD);
         MPI_Send(&transfer_real,1,MPI_DOUBLE,0,10002,PETSC_COMM_WORLD);
         MPI_Send(&transfer_imag,1,MPI_DOUBLE,0,10003,PETSC_COMM_WORLD);
         i++;
      }
   }

   // send global hdofs to all ranks

   if (rank == 0) {
      int i=1;
      while (i < size) {
         int k=0;
         while (k < global_hdof_size) {
            double transfer_real=real(global_hdof_vals[k]);
            double transfer_imag=imag(global_hdof_vals[k]);
            MPI_Send(&transfer_real,1,MPI_DOUBLE,i,20001,PETSC_COMM_WORLD);
            MPI_Send(&transfer_imag,1,MPI_DOUBLE,i,20002,PETSC_COMM_WORLD);
            k++;
         }
         i++;
      }
   } else {
      int k=0;
      while (k < global_hdof_size) {
         double transfer_real=0;
         double transfer_imag=0;
         MPI_Recv(&transfer_real,1,MPI_DOUBLE,0,20001,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         MPI_Recv(&transfer_imag,1,MPI_DOUBLE,0,20002,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

         global_hdof_vals[k]=transfer_real+PETSC_i*transfer_imag;

         k++;
      }
   }

   // transfer data to h_re and h_im - completes the re-alignment

   int index=0;
   i=0;
   while (i < global_hdof_size) {
      if (i >= fespace_ND_offsets[0] && i < fespace_ND_offsets[1]) {
         h_re->Elem(index)=real(global_hdof_vals[i]);
         h_im->Elem(index)=imag(global_hdof_vals[i]);
         index++;
      }
      i++;
   }

   // clean up
   if (ixvals) {free(ixvals); ixvals=nullptr;}
   if (local_hdof_vals) {free(local_hdof_vals); local_hdof_vals=nullptr;}
   if (global_hdof_vals) {free(global_hdof_vals); global_hdof_vals=nullptr;}
}

void fem3D::buildHgrids(BoundaryDatabase *boundaryDatabase, PWConstCoefficient *Inv_w_mu)
{
   int ierr=0;
   PetscMPIInt rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   build_Q(Inv_w_mu);

   // hdofs
   Vec hdofs;
   VecDuplicate(x,&hdofs);

   // rhs b
   Vec b;
   buildHb(&Q,&x,&b);

   MatDestroy(&Q);

   build_P();

   // solve for hofs
   ierr=solveHfield (boundaryDatabase->get_tempDirectory().c_str(),projData,rank,&P,&b,&hdofs,&HfieldError,&HfieldConverged);
   if (ierr) {
      cout << "fem3D::buildHgrids: ASSERT: solveHfield returned error code " <<  ierr << endl;
   }

   MatDestroy(&P);

   // Vec hdofs to HypreParVector h_re,h_im
   build_h_re_h_im (&hdofs);

   // fill grids
   gridReH->Distribute(*h_re);
   gridImH->Distribute(*h_im);

   // clean up
   VecDestroy(&hdofs);
   VecDestroy(&b);
}

// Re(E x H*)
void fem3D::build_ReExHgrid()
{
   VectorGridFunctionCoefficient ReEcoef(gridReE);
   VectorGridFunctionCoefficient ImEcoef(gridImE);
   VectorGridFunctionCoefficient ReHcoef(gridReH);
   VectorGridFunctionCoefficient ImHcoef(gridImH);

   // ReE x ReH
   ParDiscreteLinearOperator ReExReH (fespace_ND,fespace_RT);
   ReExReH.AddDomainInterpolator(new VectorCrossProductInterpolator(ReEcoef));
   ReExReH.Assemble();
   ReExReH.Finalize();
   ParGridFunction gridReExReH=ParGridFunction(fespace_RT);
   ReExReH.Mult(*gridReH,gridReExReH);
   VectorGridFunctionCoefficient ReExReHcoef(&gridReExReH);

   // ImE x ImH
   ParDiscreteLinearOperator ImExImH (fespace_ND,fespace_RT);
   ImExImH.AddDomainInterpolator(new VectorCrossProductInterpolator(ImEcoef));
   ImExImH.Assemble();
   ImExImH.Finalize();
   ParGridFunction gridImExImH=ParGridFunction(fespace_RT);
   ImExImH.Mult(*gridImH,gridImExImH);
   VectorGridFunctionCoefficient ImExImHcoef(&gridImExImH);

   VectorSumCoefficient ReExH=VectorSumCoefficient(ReExReHcoef,ImExImHcoef,0.5,0.5);

   gridReExH->ProjectCoefficient(ReExH);
}

// Im(E x H*)
void fem3D::build_ImExHgrid()
{
   VectorGridFunctionCoefficient ReEcoef(gridReE);
   VectorGridFunctionCoefficient ImEcoef(gridImE);
   VectorGridFunctionCoefficient ReHcoef(gridReH);
   VectorGridFunctionCoefficient ImHcoef(gridImH);

   // ReE x ImH
   ParDiscreteLinearOperator ReExImH (fespace_ND,fespace_RT);
   ReExImH.AddDomainInterpolator(new VectorCrossProductInterpolator(ReEcoef));
   ReExImH.Assemble();
   ReExImH.Finalize();
   ParGridFunction gridReExImH=ParGridFunction(fespace_RT);
   ReExImH.Mult(*gridImH,gridReExImH);
   VectorGridFunctionCoefficient ReExImHcoef(&gridReExImH);

   // ImE x ReH
   ParDiscreteLinearOperator ImExReH (fespace_ND,fespace_RT);
   ImExReH.AddDomainInterpolator(new VectorCrossProductInterpolator(ImEcoef));
   ImExReH.Assemble();
   ImExReH.Finalize();
   ParGridFunction gridImExReH=ParGridFunction(fespace_RT);
   ImExReH.Mult(*gridReH,gridImExReH);
   VectorGridFunctionCoefficient ImExReHcoef(&gridImExReH);

   VectorSumCoefficient ImExH=VectorSumCoefficient(ReExImHcoef,ImExReHcoef,-0.5,0.5);

   gridImExH->ProjectCoefficient(ImExH);
}

bool fem3D::solve(BoundaryDatabase *boundaryDatabase, MaterialDatabase *materialDatabase, double temperature,
                  PWConstCoefficient *neg_ko2_Re_er, PWConstCoefficient *neg_ko2_Im_er,
                  PWConstCoefficient *Inv_mur, PWConstCoefficient *Inv_w_mu, int drivingSport_)
{
   PetscInt i,low,high;
   PetscScalar value;
   PetscMPIInt rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   drivingSport=drivingSport_;

   if (build_A(boundaryDatabase,materialDatabase,temperature,neg_ko2_Re_er,neg_ko2_Im_er,Inv_mur)) return true;;
   build_PEC_dofs();
   eliminatePEC(&A,nPEC,PEC);
   build_X();
   build_x();
   build_grids();
   build_Xdofs();
   boundaryDatabase->fillX(&X,&Xdofs,drivingSport);
   VecCopy(X,x);

   // port dofs

   PetscInt nPortDof=0;
   VecGetOwnershipRange(Xdofs,&low,&high);
   i=low;
   while (i < high) {
      VecGetValues(Xdofs,1,&i,&value);
      if (value != 0) nPortDof++;
      i++;
   }

   PetscInt *PortDof;
   PetscMalloc(nPortDof*sizeof(PetscInt),&PortDof); 

   i=0;
   while (i < nPortDof) {
      PortDof[i]=0;
      i++;
   }

   PetscInt index=0;
   i=low;
   while (i < high) {
      VecGetValues(Xdofs,1,&i,&value);
      if (value != 0) {
         PortDof[index]=i; index++;
      }
      i++;
   }

   solveComplexLinearSystem (boundaryDatabase->get_tempDirectory().c_str(),projData,rank,
                                  nPortDof,PortDof,&A,&x,&matrixSize,&EfieldError,&EfieldConverged);

   MatDestroy(&A);
   VecDestroy(&X);
   VecDestroy(&Xdofs);
   PetscFree(PortDof); PortDof=nullptr;
   PetscFree(PEC); PEC=nullptr;
   build_e_re_e_im();
   buildEgrids(boundaryDatabase);

   PetscPrintf(PETSC_COMM_WORLD,"|         solving H field in 3D volume ...\n");
   buildHgrids (boundaryDatabase,Inv_w_mu);

   VecDestroy(&x);
//   buildZgrids(boundaryDatabase);

   if (projData->project_calculate_poynting) {
      PetscPrintf(PETSC_COMM_WORLD,"|         calculating Poynting vector field in 3D volume ...\n");
      build_ReExHgrid();
      build_ImExHgrid();
   }

   saveParaView();

   return false;
}

void fem3D::calculateMeshErrors (struct projectData *projData, BoundaryDatabase *boundaryDatabase, PWConstCoefficient *Inv_mur)
{
   PetscMPIInt rank,size;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   MPI_Comm_size(PETSC_COMM_WORLD, &size);

   CurlCurlIntegrator flux_integrator(*Inv_mur);

   RT_FECollection flux_fec(projData->mesh_order-1, (*pmesh)->SpaceDimension());
   ParFiniteElementSpace flux_fes(*pmesh, &flux_fec);

   ND_FECollection smooth_flux_fec(projData->mesh_order,(*pmesh)->Dimension());
   ParFiniteElementSpace smooth_flux_fes(*pmesh, &smooth_flux_fec);

   Vector ReLocalErrors,ImLocalErrors;
   L2ZZErrorEstimator(flux_integrator,*gridReE,smooth_flux_fes,flux_fes,ReLocalErrors,1);
   L2ZZErrorEstimator(flux_integrator,*gridImE,smooth_flux_fes,flux_fes,ImLocalErrors,1);

   // combine the real and imaginary errors into the real Vector
   int i=0; 
   while (i < ReLocalErrors.Size()) {
      ReLocalErrors[i]=sqrt(ReLocalErrors[i]*ReLocalErrors[i]+ImLocalErrors[i]*ImLocalErrors[i]);
      i++;
   }

   // merge in the errors from the prior pass (i.e. different driven port) - keep the larger error
   i=0;
   while (i < (int)errors.size()) {
      if (errors[i] > ReLocalErrors[elements[i]]) ReLocalErrors[elements[i]]=errors[i];
      i++;
   }

   // transfer the data
   errors.clear();
   elements.clear();
   ranks.clear();
   i=0;
   while (i < ReLocalErrors.Size()) {
      errors.push_back(ReLocalErrors[i]);
      elements.push_back(i);
      ranks.push_back(rank);
      i++;
   }

   // collect the errors at rank 0

   if (rank == 0) {
      int i=1;
      while (i < size) {
         int transfer_size;
         MPI_Recv(&transfer_size,1,MPI_INT,i,1000,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

         int j=0;
         while (j < transfer_size) {

            double transfer_error=0;
            MPI_Recv(&transfer_error,1,MPI_DOUBLE,i,1001,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
            errors.push_back(transfer_error);

            int transfer_element=0;
            MPI_Recv(&transfer_element,1,MPI_INT,i,1002,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
            elements.push_back(transfer_element);

            int transfer_rank=0;
            MPI_Recv(&transfer_rank,1,MPI_INT,i,1003,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
            ranks.push_back(transfer_rank);
            
            j++;
         }
         i++;
      }
   } else {
      int local_size=errors.size();
      MPI_Send(&local_size,1,MPI_INT,0,1000,PETSC_COMM_WORLD);

      int i=0;
      while (i < local_size) {
         MPI_Send(&(errors[i]),1,MPI_DOUBLE,0,1001,PETSC_COMM_WORLD);
         MPI_Send(&(elements[i]),1,MPI_INT,0,1002,PETSC_COMM_WORLD);
         MPI_Send(&(ranks[i]),1,MPI_INT,0,1003,PETSC_COMM_WORLD);
         i++;
      }
   }

   // sort and find the top errors

   refinementCount=getGlobalNE(*pmesh)*projData->mesh_3D_refinement_fraction;
   if (refinementCount == 0) refinementCount=1;

   if (rank == 0) {
      int i=0;
      while (i < refinementCount) {
         int j=i+1;
         while (j < (int)errors.size()) {
            if (errors[i] < errors[j]) {
               double temp_error=errors[i];
               errors[i]=errors[j];
               errors[j]=temp_error;

               int temp_element=elements[i];
               elements[i]=elements[j];
               elements[j]=temp_element;

               int temp_rank=ranks[i];
               ranks[i]=ranks[j];
               ranks[j]=temp_rank;
            }
            j++;
         }
         i++;
      }
   }

   // send back to the ranks

   if (rank == 0) {

      int i=1;
      while (i < size) {

         // count the elements to send
         int count=0;
         long unsigned int j=0;
         while (j < (long unsigned int)refinementCount) {
            if (ranks[j] == i) count++;
            j++;
         }

         MPI_Send(&count,1,MPI_INT,i,2000,PETSC_COMM_WORLD);

         j=0;
         while (j < (long unsigned int)refinementCount) {
            if (ranks[j] == i) {
               MPI_Send(&(errors[j]),1,MPI_DOUBLE,i,2001,PETSC_COMM_WORLD);
               MPI_Send(&(elements[j]),1,MPI_INT,i,2002,PETSC_COMM_WORLD);
               MPI_Send(&(ranks[j]),1,MPI_INT,i,2003,PETSC_COMM_WORLD);
            }
            j++;
         }

         i++;
      }

      // reduce the rank 0 elements to just the local elements
      vector<double> errors_copy=errors; errors.clear();
      vector<int> elements_copy=elements; elements.clear();
      vector<int> ranks_copy=ranks; ranks.clear();
      int j=0;
      while (j < (int)refinementCount) {
         if (ranks[j] == 0) {
            errors.push_back(errors_copy[j]);
            elements.push_back(elements_copy[j]);
            ranks.push_back(ranks_copy[j]);
         }
         j++;
      }

   } else {
      errors.clear();
      elements.clear();
      ranks.clear();

      int transfer_count=0;
      MPI_Recv(&transfer_count,1,MPI_INT,0,2000,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

      int i=0;
      while (i < transfer_count) {
         double local_error=0;
         MPI_Recv(&local_error,1,MPI_DOUBLE,0,2001,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         errors.push_back(local_error);

         int local_element=0;
         MPI_Recv(&local_element,1,MPI_INT,0,2002,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         elements.push_back(local_element);

         int local_rank=0;
         MPI_Recv(&local_rank,1,MPI_INT,0,2003,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         ranks.push_back(local_rank);

         i++;
      }
   }
}

void fem3D::refineMesh ()
{
   PetscMPIInt rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   // transfer the data for use with GeneralRefinement
   Array<int> localRefineList(elements.size());
   int i=0;
   while (i < (int)elements.size()) {
      localRefineList[i]=elements[i];
      i++;
   }

   // refine the mesh
   int previousMeshSize=getGlobalNE(*pmesh);
   (*pmesh)->GeneralRefinement(localRefineList);
   int newMeshSize=getGlobalNE(*pmesh);

   if (projData->output_show_refining_mesh) PetscPrintf(PETSC_COMM_WORLD,"|             new mesh size: %d\n",newMeshSize);
   if (newMeshSize > 3*previousMeshSize) PetscPrintf(PETSC_COMM_WORLD,
      "|             Warning: The mesh size jumped from %d to %d.\n",previousMeshSize,newMeshSize);
}

void fem3D::saveParaView()
{
   if (!projData->project_save_fields) return;

   stringstream ss;
   ss << projData->project_name << "_frequency_" << frequency << "_Sport_" << drivingSport;

   stringstream ssParaView;
   ssParaView << "ParaView_" << projData->project_name;

   ParaViewDataCollection *pd=nullptr;
   pd=new ParaViewDataCollection(ss.str(),*pmesh);
   pd->SetOwnData(false);
   pd->SetPrefixPath(ssParaView.str());
   pd->RegisterField("gridReE",gridReE);
   pd->RegisterField("gridImE",gridImE);
   pd->RegisterField("gridReH",gridReH);
   pd->RegisterField("gridImH",gridImH);
   if (projData->project_calculate_poynting) {
      pd->RegisterField("gridReExH",gridReExH);
      pd->RegisterField("gridImExH",gridImExH);
   }
   pd->SetLevelsOfDetail(3);
   pd->SetDataFormat(VTKFormat::ASCII);
   pd->SetHighOrderOutput(true);
   pd->Save();
   delete pd;
}

void fem3D::saveFieldValuesHeader (struct projectData *projData)
{
   PetscMPIInt rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   // nothing to do
   if (projData->field_points_count == 0) return;

   if (rank == 0) {
      stringstream ss;
      ss << projData->project_name << "_fields.csv";

      ofstream out;
      out.open(ss.str().c_str(),ofstream::app);

      if (out.is_open()) {
         out << "#frequency,iteration,driving Sport,"
             << "real(Ex),real(Ey),real(Ez),"
             << "imag(Ex),imag(Ey),imag(Ez),"
             << "real(Hx),real(Hy),real(Hz),"
             << "imag(Hx),imag(Hy),imag(Hz)"
             << endl;

         out.close();
      } else {
         PetscPrintf(PETSC_COMM_WORLD,"ERROR3009: Failed to open file \"%s\" for writing.\n",ss.str().c_str());
      }
   }
}

void fem3D::saveFieldValues (struct projectData *projData, ParMesh *pmesh, int iteration, int drivingSport)
{
   PetscMPIInt rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   // nothing to do
   if (projData->field_points_count == 0) return;

   // transfer the data into a useful format - make a copy for E and a copy for H
   OPEMIntegrationPointList Elist,Hlist;
   int i=0;
   while (i < projData->field_points_count) {

      OPEMIntegrationPoint *point=new OPEMIntegrationPoint(i,
         projData->field_points_x[i],projData->field_points_y[i],projData->field_points_z[i]);
      Elist.push(point);

      point=new OPEMIntegrationPoint(i,
         projData->field_points_x[i],projData->field_points_y[i],projData->field_points_z[i]);
      Hlist.push(point);

      i++;
   }

   // get the values

   Elist.update(pmesh);
   Elist.get_fieldValues(gridReE,gridImE);
   Elist.assemble();

   Hlist.update(pmesh);
   Hlist.get_fieldValues(gridReH,gridImH);
   Hlist.assemble();

   if (rank == 0) {

      stringstream ss;
      ss << projData->project_name << "_fields.csv";

      ofstream out;
      out.open(ss.str().c_str(),ofstream::app);

      if (out.is_open()) {
         long unsigned int i=0;
         while (i < Elist.get_size()) {
            double x1,y1,z1,x2,y2,z2;
            Elist.get_point(i)->get_location(&x1,&y1,&z1);
            Hlist.get_point(i)->get_location(&x2,&y2,&z2);
            if (!double_compare(x1,x2,1e-12) || !double_compare(y1,y2,1e-12) || !double_compare(z1,z2,1e-12)) {
               cout << "ASSERT: fem3D::saveFieldValues: Mismatched coordinates." << endl;
            }

            complex<double> Ex,Ey,Ez,Hx,Hy,Hz;
            Elist.get_point(i)->get_fields(&Ex,&Ey,&Ez);
            Hlist.get_point(i)->get_fields(&Hx,&Hy,&Hz);

            // match ParaView's "probe" output format
            // Note that the results do not exactly match Paraview.
            out << frequency << ","
                << iteration << ","
                << drivingSport << ", "
                << x1 << "," << y1 << "," << z1 << ", "
                << real(Ex) << "," << real(Ey) << "," << real(Ez) << ", "
                << imag(Ex) << "," << imag(Ey) << "," << imag(Ez) << ", "
                << real(Hx) << "," << real(Hy) << "," << real(Hz) << ", "
                << imag(Hx) << "," << imag(Hy) << "," << imag(Hz) 
                << endl;

            i++;
         }

         out.close();
      } else {
         PetscPrintf(PETSC_COMM_WORLD,"ERROR3004: Failed to open file \"%s\" for writing.\n",ss.str().c_str());
      }
   }
}

/*
// alternative serial calculation for debugging
void fem3D::saveFieldValues2 (struct projectData *projData, ParMesh *pmesh, int iteration, int drivingSport)
{
   int i=0;
   while (i < projData->field_points_count) {

      Vector point;
      point.SetSize(3);
      point(0)=projData->field_points_x[i];
      point(1)=projData->field_points_y[i];
      point(2)=projData->field_points_z[i];

      IntegrationPoint integrationPoint;

      int j=0;
      while (j < pmesh->GetNE()) {

         ElementTransformation *eltransf=pmesh->GetElementTransformation(j);

         InverseElementTransformation *inv_tr=new InverseElementTransformation;
         inv_tr->SetPrintLevel(-1);         // -1 - never print
         inv_tr->SetInitialGuessType(InverseElementTransformation::Center);
         inv_tr->SetSolverType(InverseElementTransformation::Newton);
         inv_tr->SetTransformation(*eltransf);
         int res=inv_tr->Transform(point,integrationPoint);
         delete inv_tr;

         if (res == InverseElementTransformation::Inside) break;

         j++;
      }

      if (j < pmesh->GetNE()) {
         Vector vectorReE,vectorImE,vectorReH,vectorImH;
         gridReE->GetVectorValue(j,integrationPoint,vectorReE);
         gridImE->GetVectorValue(j,integrationPoint,vectorImE);
         gridReH->GetVectorValue(j,integrationPoint,vectorReH);
         gridImH->GetVectorValue(j,integrationPoint,vectorImH);
         cout << frequency << ","
              << iteration << ","
              << drivingSport << ","
              << projData->field_points_x[i] << "," << projData->field_points_y[i] << "," << projData->field_points_z[i] << ","
              << vectorReE[0] << "," << vectorReE[1] << "," << vectorReE[2] << ","
              << vectorImE[0] << "," << vectorImE[1] << "," << vectorImE[2] << ","
              << vectorReH[0] << "," << vectorReH[1] << "," << vectorReH[2] << ","
              << vectorImH[0] << "," << vectorImH[1] << "," << vectorImH[2] << ","
              << endl;
      } else {
         PetscPrintf (PETSC_COMM_WORLD,"ERROR3006: Failed to locate a test point within an element.\n");
      }

      i++;
   }
}
*/

fem3D::~fem3D()
{
   if (fespace_ND) {delete fespace_ND; fespace_ND=nullptr;}
   if (fespace_RT) {delete fespace_RT; fespace_RT=nullptr;}
   if (fespace_H1) {delete fespace_H1; fespace_H1=nullptr;}
   if (fespace_L2) {delete fespace_L2; fespace_L2=nullptr;}
   if (fec_ND) {delete fec_ND; fec_ND=nullptr;}
   if (fec_RT) {delete fec_RT; fec_RT=nullptr;}
   if (fec_H1) {delete fec_H1; fec_H1=nullptr;}
   if (fec_L2) {delete fec_L2; fec_L2=nullptr;}

   if (e_re) {delete e_re; e_re=nullptr;}
   if (e_im) {delete e_im; e_im=nullptr;}
   if (h_re) {delete h_re; h_re=nullptr;}
   if (h_im) {delete h_im; h_im=nullptr;}

   if (gridReE) {delete gridReE; gridReE=nullptr;}
   if (gridImE) {delete gridImE; gridImE=nullptr;}
   if (gridReH) {delete gridReH; gridReH=nullptr;}
   if (gridImH) {delete gridImH; gridImH=nullptr;}
   if (gridReExH) {delete gridReExH; gridReExH=nullptr;}
   if (gridImExH) {delete gridImExH; gridImExH=nullptr;}

//   if (gridReEz) {delete gridReEz; gridReEz=nullptr;}
//   if (gridImEz) {delete gridImEz; gridImEz=nullptr;}
//   if (gridReHz) {delete gridReHz; gridReHz=nullptr;}
//   if (gridImHz) {delete gridImHz; gridImHz=nullptr;}

}

