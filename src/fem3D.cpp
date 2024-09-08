////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//    OpenParEM3D - A fullwave 3D electromagnetic simulator.                  //
//    Copyright (C) 2024 Brian Young                                          //
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
                     PWConstCoefficient *Inv_mur, PWConstCoefficient *neg_ko2_Re_er, PWConstCoefficient *neg_ko2_Im_er,
                     bool solution_check_homogeneous, string indent)
{
   vector<Array<int> *> borderAttributeList;
   vector<ConstantCoefficient *> ReInvGammaConstList,ImInvGammaConstList,RekConstList,ImkConstList,ZconstList;

   PetscMPIInt rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   // real

   ParBilinearForm *pmblfReA=new ParBilinearForm(fespace_ND);
   pmblfReA->AddDomainIntegrator(new CurlCurlIntegrator (*Inv_mur));
   pmblfReA->AddDomainIntegrator(new VectorFEMassIntegrator (*neg_ko2_Re_er));
   if (boundaryDatabase->addPortIntegrators(*pmesh,pmblfReA,Inv_mur,neg_ko2_Re_er,neg_ko2_Im_er,borderAttributeList,
                                            ReInvGammaConstList,ImInvGammaConstList,RekConstList,ImkConstList,
                                            true,drivingSet,solution_check_homogeneous,indent)) return true;
   pmblfReA->Assemble();
   pmblfReA->Finalize();
   HypreParMatrix *ReA=pmblfReA->ParallelAssemble();

   long unsigned int i=0;
   while (i < borderAttributeList.size()) {
      if (borderAttributeList[i]) delete borderAttributeList[i];
      i++;
   }
   borderAttributeList.clear();

   i=0;
   while (i < ReInvGammaConstList.size()) {
      if (ReInvGammaConstList[i]) delete ReInvGammaConstList[i];
      if (ImInvGammaConstList[i]) delete ImInvGammaConstList[i];
      if (RekConstList[i]) delete RekConstList[i];
      if (ImkConstList[i]) delete ImkConstList[i];
      i++;
   }
   ReInvGammaConstList.clear();
   ImInvGammaConstList.clear();
   RekConstList.clear();
   ImkConstList.clear();


   delete pmblfReA;
   hypre_getSparseWidth((hypre_ParCSRMatrix *) *ReA,&sparseWidth);
   if (hypre_ParCSRMatrixToMat((hypre_ParCSRMatrix *) *ReA, &A, sparseWidth, 1, 0, 0)) return true;
   delete ReA;

   // imag

   ParBilinearForm *pmblfImA=new ParBilinearForm(fespace_ND);
   pmblfImA->AddDomainIntegrator(new VectorFEMassIntegrator (*neg_ko2_Im_er));
   if (boundaryDatabase->addPortIntegrators(*pmesh,pmblfImA,Inv_mur,neg_ko2_Re_er,neg_ko2_Im_er,borderAttributeList,
                                            ReInvGammaConstList,ImInvGammaConstList,RekConstList,ImkConstList,
                                            false,drivingSet,solution_check_homogeneous,indent)) return true;
   boundaryDatabase->addImpedanceIntegrators(frequency,temperature,*pmesh,pmblfImA,materialDatabase,
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
   while (i < ReInvGammaConstList.size()) {
      if (ReInvGammaConstList[i]) delete ReInvGammaConstList[i];
      if (ImInvGammaConstList[i]) delete ImInvGammaConstList[i];
      if (RekConstList[i]) delete RekConstList[i];
      if (ImkConstList[i]) delete ImkConstList[i];
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
   ParBilinearForm *pmblfReQ=new ParBilinearForm(fespace_ND);
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

void fem3D::build_R (PWConstCoefficient *ReInv_w_er, PWConstCoefficient *ImInv_w_er)
{
   // real

   ParBilinearForm *pmblfReR=new ParBilinearForm(fespace_ND);
   pmblfReR->AddDomainIntegrator(new MixedVectorCurlIntegrator (*ReInv_w_er));
   pmblfReR->Assemble();
   pmblfReR->Finalize();
   HypreParMatrix *ReR=pmblfReR->ParallelAssemble();
   delete pmblfReR;

   // anticipate solving PE=-jRH and put ReR in the imaginary position as a negative
   (*ReR)*=-1.0;
   if (hypre_ParCSRMatrixToMat((hypre_ParCSRMatrix *) *ReR, &R, sparseWidth, 1, 1, 0)) {
      cout << "ASSERT: fem3D::build_R failed to build Mat ReR." << endl;
   }

   delete ReR;

   // imag

   ParBilinearForm *pmblfImR=new ParBilinearForm(fespace_ND);
   pmblfImR->AddDomainIntegrator(new MixedVectorCurlIntegrator (*ImInv_w_er));
   pmblfImR->Assemble();
   pmblfImR->Finalize();
   HypreParMatrix *ImR=pmblfImR->ParallelAssemble();
   delete pmblfImR;

   // anticipate solving PE=-jRH and put ImR in the real position
   if (hypre_ParCSRMatrixToMat((hypre_ParCSRMatrix *) *ImR, &R, sparseWidth, 0, 0, 1)) {
      cout << "ASSERT: fem3D::build_R failed to build Mat ImR." << endl;
   }

   delete ImR;
}

PetscErrorCode fem3D::build_X()
{
   PetscErrorCode ierr=0;
   ierr=MatCreateVecs(A,&X,NULL); if (ierr) return ierr;
   ierr=VecZeroEntries(X);
   return ierr;
}

PetscErrorCode fem3D::build_x()
{
   PetscErrorCode ierr=0;
   ierr=MatCreateVecs(A,&x,NULL); if (ierr) return ierr;
   ierr=VecZeroEntries(x);
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
   ierr=MatCreateVecs(A,&Xdofs,NULL);  if (ierr) return ierr;
   ierr=VecZeroEntries(Xdofs);
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
   build_h_re_h_im(&hdofs);

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
                  PWConstCoefficient *Inv_mur, PWConstCoefficient *Inv_w_mu, int drivingSet_, bool calculateH,
                  bool solution_check_homogeneous, string indent)
{
   PetscInt i,low,high;
   PetscScalar value;
   PetscMPIInt rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   drivingSet=drivingSet_;

   if (build_A(boundaryDatabase,materialDatabase,temperature,Inv_mur,neg_ko2_Re_er,neg_ko2_Im_er,solution_check_homogeneous,indent)) return true;

   build_PEC_dofs();
   eliminatePEC(&A,nPEC,PEC);
   build_X();
   build_x();
   build_grids();
   build_Xdofs();
   boundaryDatabase->fillX(&X,&Xdofs,drivingSet);
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

   if (calculateH) {
      prefix(); PetscPrintf(PETSC_COMM_WORLD,"         solving H field in 3D volume ...\n");
      buildHgrids (boundaryDatabase,Inv_w_mu);
   }

   VecDestroy(&x);

   //buildZgrids(boundaryDatabase);

   if (projData->project_calculate_poynting) {
      if (calculateH) {
         prefix(); PetscPrintf(PETSC_COMM_WORLD,"         calculating Poynting vector field in 3D volume ...\n");
         build_ReExHgrid();
         build_ImExHgrid();
      } else {
         prefix(); PetscPrintf(PETSC_COMM_WORLD,"         skipping Poynting vector field in 3D volume due to omitted H-field calculation ...\n");
      }
   }

   saveParaView(calculateH,projData->project_calculate_poynting,boundaryDatabase->get_drivingSetName());

   return false;
}

void fem3D::printElementVertices ()
{
   PetscMPIInt rank,size;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   MPI_Comm_size(PETSC_COMM_WORLD, &size);

   Array<int> vertices;
   real_t *coord;
   int i=0;
   while (i < (*pmesh)->GetNE()) {
      (*pmesh)->GetElementVertices(i,vertices);

      cout << rank << ": element " << i << ":" << endl;
      int j=0;
      while (j < vertices.Size()) {
         coord=(*pmesh)->GetVertex(j);
         cout << rank << ":   vertex j x=" << coord[0] << " y=" << coord[1] << " z=" << coord[2] << endl;
         j++;
      }
      i++;
   }
}

// calculate the mesh errors using a Zienkiewicz-Zhu estimator
//
// Calculates the mesh errors on the H fields.  Since these are calculated from E, any error in E is magnified, and adaptive refinement
// does a better job of refining the mesh than if the errors are calculated on E.  To switch to calculating the error on E, change
// gridReH and gridImH to gridReE and gridImE in the two calls to OPEM_L2ZZErrorEstimator.
bool fem3D::calculateMeshErrors (struct projectData *projData, BoundaryDatabase *boundaryDatabase, PWConstCoefficient *Inv_mur, double *maxMeshError_)
{
   PetscMPIInt rank,size;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   MPI_Comm_size(PETSC_COMM_WORLD, &size);

   CurlCurlIntegrator flux_integrator(*Inv_mur);

   RT_FECollection flux_fec(projData->mesh_order-1, (*pmesh)->SpaceDimension());
   ParFiniteElementSpace flux_fes(*pmesh, &flux_fec);

   ND_FECollection smooth_flux_fec(projData->mesh_order,(*pmesh)->Dimension());
   ParFiniteElementSpace smooth_flux_fes(*pmesh, &smooth_flux_fec);

   double solution_tolerance=1e-12;
   double solution_tolerance_message_limit=1e-3;
   int iteration_limit=2000;
   double ReFinalResidualNorm=0;
   double ImFinalResidualNorm=0;
   bool printError=false;
   Vector ReLocalErrors,ImLocalErrors;
   if (OPEM_L2ZZErrorEstimator(flux_integrator,*gridReH,smooth_flux_fes,flux_fes,ReLocalErrors,1,solution_tolerance,iteration_limit,ReFinalResidualNorm)) printError=true;
   if (OPEM_L2ZZErrorEstimator(flux_integrator,*gridImH,smooth_flux_fes,flux_fes,ImLocalErrors,1,solution_tolerance,iteration_limit,ImFinalResidualNorm)) printError=true;
   if (printError) {
      if (ImFinalResidualNorm > ReFinalResidualNorm) ReFinalResidualNorm=ImFinalResidualNorm;
      if (ReFinalResidualNorm > solution_tolerance_message_limit) {
         prefix(); PetscPrintf (PETSC_COMM_WORLD,"            INFO: Calculated mesh errors are approximate with the final error of %g > the target of %g.\n",ReFinalResidualNorm,solution_tolerance_message_limit);
      }
   }

   int i=0; 
   while (i < ReLocalErrors.Size()) {
      ReLocalErrors[i]=sqrt(ReLocalErrors[i]*ReLocalErrors[i]+ImLocalErrors[i]*ImLocalErrors[i]);
      i++;
   }

   /* Does not work as well as with no scaling
   // scale the local errors by the longest element edge
   Array<int> vertices;
   real_t *coord1,*coord2;
   i=0;
   while (i < (*pmesh)->GetNE()) {
      (*pmesh)->GetElementVertices(i,vertices);
      real_t maxLength=0;
      int j=0;
      while (j < vertices.Size()) {
         coord1=(*pmesh)->GetVertex(j);
         int k=0;
         while (k < vertices.Size()) {
            if (k != j) {
               coord2=(*pmesh)->GetVertex(k);
               real_t length=sqrt(pow(coord1[0]-coord2[0],2)+pow(coord1[1]-coord2[1],2)+pow(coord1[2]-coord2[2],2));
               if (length > maxLength) maxLength=length;
            }
            k++;
         }
         j++;
      }
      //ReLocalErrors[i]*=pow(maxLength,3);
      ReLocalErrors[i]*=maxLength;
      i++;
   }
   */

   /* Does not work as well as with no scaling
   // scale the local errors by the element volume
   i=0;
   while (i < (*pmesh)->GetNE()) {
      real_t volume=(*pmesh)->GetElementVolume(i);
      ReLocalErrors[i]*=volume;
      //ReLocalErrors[i]*=pow(volume,0.3333);
      i++;
   }
   */

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

   double maxMeshError=-1;
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
      maxMeshError=errors[0];
   }

   MPI_Bcast(&maxMeshError,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
   *maxMeshError_=maxMeshError;

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

   return false;
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

   if (projData->output_show_refining_mesh) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"             new mesh size: %d\n",newMeshSize);}
   if (newMeshSize > 3*previousMeshSize) {prefix(); PetscPrintf(PETSC_COMM_WORLD,
      "             Warning: The mesh size jumped from %d to %d.\n",previousMeshSize,newMeshSize);}
}

void fem3D::saveParaView (bool plot_H, bool plot_poynting, string set_name)
{
   if (!projData->project_save_fields) return;

   stringstream ss;
   ss << projData->project_name << "_frequency_" << frequency << "_" << set_name << "_" << drivingSet;

   stringstream ssParaView;
   ssParaView << "ParaView_" << projData->project_name;

   ParaViewDataCollection *pd=nullptr;
   pd=new ParaViewDataCollection(ss.str(),*pmesh);
   pd->SetOwnData(false);
   pd->SetPrefixPath(ssParaView.str());
   pd->RegisterField("gridReE",gridReE);
   pd->RegisterField("gridImE",gridImE);
   if (plot_H) pd->RegisterField("gridReH",gridReH);
   if (plot_H) pd->RegisterField("gridImH",gridImH);
   if (plot_poynting) pd->RegisterField("gridReExH",gridReExH);
   if (plot_poynting) pd->RegisterField("gridImExH",gridImExH);
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
         prefix(); PetscPrintf(PETSC_COMM_WORLD,"ERROR3009: Failed to open file \"%s\" for writing.\n",ss.str().c_str());
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
         prefix(); PetscPrintf(PETSC_COMM_WORLD,"ERROR3004: Failed to open file \"%s\" for writing.\n",ss.str().c_str());
      }
   }
}

// conversion = 0 => conversion from modal to single-ended
// conversion = 1 => conversion from single-ended to mixed-mode
PetscErrorCode fem3D::build_Mc_Ms (double referenceZ, BoundaryDatabase *boundaryDatabase, vector<DifferentialPair *> *differentialPairList, int conversion)
{
   PetscErrorCode ierr=0;
   bool SINGLEENDED=false,MIXEDMODE=false;

   if (conversion == 0) SINGLEENDED=true;
   if (conversion == 1) MIXEDMODE=true;

   // for MIXEDMODE
   double Zoe=referenceZ/2;
   double Zoo=referenceZ*2;

   int SportCount=boundaryDatabase->get_SportCount();

   // make a list of the modes across all ports - ordered in terms of increasing Sport to align with the S-parameters
   vector<Mode *> modeList;
   if (boundaryDatabase->buildAggregateModeList(&modeList)) return true;
   if (SportCount != (int)modeList.size()) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: fem3D::build_Mc_Ms found inconsistent number of Sports.\n");
      return true;
   }

   // for convenience, make lists of the ports for each mode
   vector<Port *> portList;
   long unsigned int i=0;
   while (i < modeList.size()) {
      Port *port;
      long unsigned int index;
      if (boundaryDatabase->get_port_from_mode(modeList[i],&port,&index)) {
         PetscPrintf(PETSC_COMM_WORLD,"ASSERT: fem3D::build_Mc_Ms failed to find a port.\n");
         return true;
      }
      portList.push_back(port);
      i++;
   }

   // allocate space
   Mat Kv;
   ierr=MatCreate(PETSC_COMM_WORLD,&Kv); if (ierr) return ierr;
   ierr=MatSetType(Kv,MATDENSE); if (ierr) return ierr;
   ierr=MatSetSizes(Kv,PETSC_DECIDE,PETSC_DECIDE,SportCount,SportCount); if (ierr) return ierr;
   ierr=MatZeroEntries(Kv); if (ierr) return ierr;
   PetscInt Kvlow,Kvhigh;
   ierr=MatGetOwnershipRange(Kv,&Kvlow,&Kvhigh); if (ierr) return ierr;

   Mat Ki;
   ierr=MatCreate(PETSC_COMM_WORLD,&Ki); if (ierr) return ierr;
   ierr=MatSetType(Ki,MATDENSE); if (ierr) return ierr;
   ierr=MatSetSizes(Ki,PETSC_DECIDE,PETSC_DECIDE,SportCount,SportCount); if (ierr) return ierr;
   ierr=MatZeroEntries(Ki); if (ierr) return ierr;
   PetscInt Kilow,Kihigh;
   ierr=MatGetOwnershipRange(Ki,&Kilow,&Kihigh); if (ierr) return ierr;

   Mat rootZa;
   ierr=MatCreate(PETSC_COMM_WORLD,&rootZa); if (ierr) return ierr;
   ierr=MatSetType(rootZa,MATDENSE); if (ierr) return ierr;
   ierr=MatSetSizes(rootZa,PETSC_DECIDE,PETSC_DECIDE,SportCount,SportCount); if (ierr) return ierr;
   ierr=MatZeroEntries(rootZa); if (ierr) return ierr;
   PetscInt rootZalow,rootZahigh;
   ierr=MatGetOwnershipRange(rootZa,&rootZalow,&rootZahigh); if (ierr) return ierr;

   Mat rootZb;
   ierr=MatCreate(PETSC_COMM_WORLD,&rootZb); if (ierr) return ierr;
   ierr=MatSetType(rootZb,MATDENSE); if (ierr) return ierr;
   ierr=MatSetSizes(rootZb,PETSC_DECIDE,PETSC_DECIDE,SportCount,SportCount); if (ierr) return ierr;
   ierr=MatZeroEntries(rootZb); if (ierr) return ierr;
   PetscInt rootZblow,rootZbhigh;
   ierr=MatGetOwnershipRange(rootZb,&rootZblow,&rootZbhigh); if (ierr) return ierr;

   Mat invRootZa;
   ierr=MatCreate(PETSC_COMM_WORLD,&invRootZa); if (ierr) return ierr;
   ierr=MatSetType(invRootZa,MATDENSE); if (ierr) return ierr;
   ierr=MatSetSizes(invRootZa,PETSC_DECIDE,PETSC_DECIDE,SportCount,SportCount); if (ierr) return ierr;
   ierr=MatZeroEntries(invRootZa); if (ierr) return ierr;
   PetscInt invRootZalow,invRootZahigh;
   ierr=MatGetOwnershipRange(invRootZa,&invRootZalow,&invRootZahigh); if (ierr) return ierr;

   Mat invRootZb;
   ierr=MatCreate(PETSC_COMM_WORLD,&invRootZb); if (ierr) return ierr;
   ierr=MatSetType(invRootZb,MATDENSE); if (ierr) return ierr;
   ierr=MatSetSizes(invRootZb,PETSC_DECIDE,PETSC_DECIDE,SportCount,SportCount); if (ierr) return ierr;
   ierr=MatZeroEntries(invRootZb); if (ierr) return ierr;
   PetscInt invRootZblow,invRootZbhigh;
   ierr=MatGetOwnershipRange(invRootZb,&invRootZblow,&invRootZbhigh); if (ierr) return ierr;


   // fill rootZa and rootZb

   PetscScalar value;

   if (SINGLEENDED) {
      PetscInt i=0;
      while (i < (PetscInt)modeList.size()) {
         value=real(sqrt(modeList[i]->get_impedance()))+PETSC_i*imag(sqrt(modeList[i]->get_impedance()));
         if (i >= rootZalow && i < rootZahigh) {ierr=MatSetValue(rootZa,i,i,value,INSERT_VALUES); if (ierr) return ierr;}
         SportZoList.push_back(referenceZ);
         value=sqrt(referenceZ);
         if (i >= rootZblow && i < rootZbhigh) {ierr=MatSetValue(rootZb,i,i,value,INSERT_VALUES); if (ierr) return ierr;}
         i++;
      }
   }

   if (MIXEDMODE) {
      PetscInt i=0;
      while (i < (PetscInt)modeList.size()) {
         value=sqrt(referenceZ);
         if (i >= rootZalow && i < rootZahigh) {ierr=MatSetValue(rootZa,i,i,value,INSERT_VALUES); if (ierr) return ierr;}
         SportZoList.push_back(referenceZ);
         value=sqrt(referenceZ);
         if (i >= rootZblow && i < rootZbhigh) {ierr=MatSetValue(rootZb,i,i,value,INSERT_VALUES); if (ierr) return ierr;} // default to single-ended
         i++;
      }

      i=0;
      while (i < (PetscInt)differentialPairList->size()) {
         PetscInt p=(*differentialPairList)[i]->get_Sport_P()-1;
         PetscInt n=(*differentialPairList)[i]->get_Sport_N()-1;

         // set to mixed mode
         SportZoList[p]=Zoe;
         value=sqrt(Zoe);
         if (p >= rootZblow && p < rootZbhigh) {ierr=MatSetValue(rootZb,p,p,value,INSERT_VALUES); if (ierr) return ierr;}
         SportZoList[n]=Zoo;
         value=sqrt(Zoo);
         if (n >= rootZblow && n < rootZbhigh) {ierr=MatSetValue(rootZb,n,n,value,INSERT_VALUES); if (ierr) return ierr;}
         i++;
      }
   }

   ierr=MatAssemblyBegin(rootZa,MAT_FINAL_ASSEMBLY); if (ierr) return ierr;
   ierr=MatAssemblyEnd(rootZa,MAT_FINAL_ASSEMBLY); if (ierr) return ierr;

   ierr=MatAssemblyBegin(rootZb,MAT_FINAL_ASSEMBLY); if (ierr) return ierr;
   ierr=MatAssemblyEnd(rootZb,MAT_FINAL_ASSEMBLY); if (ierr) return ierr;

   // fill invRootZa and invRootZb

   if (SINGLEENDED) {
      PetscInt i=0;
      while (i < (PetscInt)modeList.size()) {
         value=real(1/sqrt(modeList[i]->get_impedance()))+PETSC_i*imag(1/sqrt(modeList[i]->get_impedance()));
         if (i >= invRootZalow && i < invRootZahigh) {ierr=MatSetValue(invRootZa,i,i,value,INSERT_VALUES); if (ierr) return ierr;}
         value=1/sqrt(referenceZ);
         if (i >= invRootZblow && i < invRootZbhigh) {ierr=MatSetValue(invRootZb,i,i,value,INSERT_VALUES); if (ierr) return ierr;}
         i++;
      }
   }

   if (MIXEDMODE) {
      PetscInt i=0;
      while (i < (PetscInt)modeList.size()) {
         value=1/sqrt(referenceZ);
         if (i >= invRootZalow && i < invRootZahigh) {ierr=MatSetValue(invRootZa,i,i,value,INSERT_VALUES); if (ierr) return ierr;}
         value=1/sqrt(referenceZ);
         if (i >= invRootZblow && i < invRootZbhigh) {ierr=MatSetValue(invRootZb,i,i,value,INSERT_VALUES); if (ierr) return ierr;} // default to single-ended
         i++;
      }

      i=0;
      while (i < (PetscInt)differentialPairList->size()) {
         PetscInt p=(*differentialPairList)[i]->get_Sport_P()-1;
         PetscInt n=(*differentialPairList)[i]->get_Sport_N()-1;

         // set to mixed mode
         value=1/sqrt(Zoe);
         if (p >= invRootZblow && p < invRootZbhigh) {ierr=MatSetValue(invRootZb,p,p,value,INSERT_VALUES); if (ierr) return ierr;} // common mode
         value=1/sqrt(Zoo);
         if (n >= invRootZblow && n < invRootZbhigh) {ierr=MatSetValue(invRootZb,n,n,value,INSERT_VALUES); if (ierr) return ierr;} // differential mode
         i++;
      }
   }

   ierr=MatAssemblyBegin(invRootZa,MAT_FINAL_ASSEMBLY); if (ierr) return ierr;
   ierr=MatAssemblyEnd(invRootZa,MAT_FINAL_ASSEMBLY); if (ierr) return ierr;

   ierr=MatAssemblyBegin(invRootZb,MAT_FINAL_ASSEMBLY); if (ierr) return ierr;
   ierr=MatAssemblyEnd(invRootZb,MAT_FINAL_ASSEMBLY); if (ierr) return ierr;

   // Kv and Ki

   // Note: From OpenParEM2D for line calculations, Ti/Tv goes from single-ended voltages to modal voltages.

   if (SINGLEENDED) {
      PetscInt i=0;
      while (i < (PetscInt)portList.size()) {

         // select the row from Ti/Tv to insert
         // ports may be in sequential order
         int row=-1;
         PetscInt j=0;
         while (j <= i) {
            if (portList[j] == portList[i]) row++;
            j++;
         }

         // fill Ki/Kv from Ti/Tv

         int TiTvSize=portList[i]->get_TiTvSize();
         j=0;
         while (j < (PetscInt) TiTvSize) {

            // get the mode number for this TiTv value
            int count=-1;
            PetscInt k=0;
            while (k < (PetscInt)portList.size()) {
               if (portList[k] == portList[i]) count++;
               if (count == (int)j) break;
               k++;
            }

            value=portList[i]->get_ReTi(row,j)+PETSC_i*portList[i]->get_ImTi(row,j);
            if (i >= Kilow && i < Kihigh) {ierr=MatSetValue(Ki,i,k,value,INSERT_VALUES); if (ierr) return ierr;} 

            value=portList[i]->get_ReTv(row,j)+PETSC_i*portList[i]->get_ImTv(row,j);
            if (i >= Kvlow && i < Kvhigh) {ierr=MatSetValue(Kv,i,k,value,INSERT_VALUES); if (ierr) return ierr;}

            j++;
         }
         i++;
      }
   }

   if (MIXEDMODE) {

      // default for single ended
      PetscInt i=0;
      while (i < (PetscInt)modeList.size()) {
         value=1;
         if (i >= Kvlow && i < Kvhigh) {ierr=MatSetValue(Kv,i,i,value,INSERT_VALUES); if (ierr) return ierr;}
         value=1;
         if (i >= Kilow && i < Kihigh) {ierr=MatSetValue(Ki,i,i,value,INSERT_VALUES); if (ierr) return ierr;}
         i++;
      }

      // mixed mode
      i=0;
      while (i < (PetscInt)differentialPairList->size()) {
         PetscInt p=(*differentialPairList)[i]->get_Sport_P()-1;
         PetscInt n=(*differentialPairList)[i]->get_Sport_N()-1;

         // common-mode voltage
         value=0.5;
         if (p >= Kvlow && p < Kvhigh) {ierr=MatSetValue(Kv,p,p,value,INSERT_VALUES); if (ierr) return ierr;}
         value=0.5;
         if (p >= Kvlow && p < Kvhigh) {ierr=MatSetValue(Kv,p,n,value,INSERT_VALUES); if (ierr) return ierr;}

         // differential voltage
         value=+1;
         if (n >= Kvlow && n < Kvhigh) {ierr=MatSetValue(Kv,n,p,value,INSERT_VALUES); if (ierr) return ierr;}
         value=-1;
         if (n >= Kvlow && n < Kvhigh) {ierr=MatSetValue(Kv,n,n,value,INSERT_VALUES); if (ierr) return ierr;}

         // common-mode current
         value=1;
         if (p >= Kilow && p < Kihigh) {ierr=MatSetValue(Ki,p,p,value,INSERT_VALUES); if (ierr) return ierr;}
         value=1;
         if (p >= Kilow && p < Kihigh) {ierr=MatSetValue(Ki,p,n,value,INSERT_VALUES); if (ierr) return ierr;}

         // differential current
         value=+0.5;
         if (n >= Kilow && n < Kihigh) {ierr=MatSetValue(Ki,n,p,value,INSERT_VALUES); if (ierr) return ierr;}
         value=-0.5;
         if (n >= Kilow && n < Kihigh) {ierr=MatSetValue(Ki,n,n,value,INSERT_VALUES); if (ierr) return ierr;}

         i++;
      }
   }

   ierr=MatAssemblyBegin(Kv,MAT_FINAL_ASSEMBLY); if (ierr) return ierr;
   ierr=MatAssemblyEnd(Kv,MAT_FINAL_ASSEMBLY); if (ierr) return ierr;

   ierr=MatAssemblyBegin(Ki,MAT_FINAL_ASSEMBLY); if (ierr) return ierr;
   ierr=MatAssemblyEnd(Ki,MAT_FINAL_ASSEMBLY); if (ierr) return ierr;

   // The calculation is going from modal to single-ended, so invert since
   // Ti/Tv goes from single-ended to modal.
   if (SINGLEENDED) {
      ierr=MatInvert(&Kv,0); if (ierr) return ierr;
      ierr=MatInvert(&Ki,0); if (ierr) return ierr;
   }

   // intermediate terms

   Mat term1a,term1;
   ierr=MatMatMult(invRootZb,Kv,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&term1a); if (ierr) return ierr;
   ierr=MatMatMult(term1a,rootZa,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&term1); if (ierr) return ierr;
   ierr=MatDestroy(&term1a);  if (ierr) return ierr;

   Mat term2a,term2;
   ierr=MatMatMult(rootZb,Ki,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&term2a); if (ierr) return ierr;
   ierr=MatMatMult(term2a,invRootZa,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&term2); if (ierr) return ierr;
   ierr=MatDestroy(&term2a);  if (ierr) return ierr;

   // Ms
   ierr=MatConvert(term1,MATSAME,MAT_INITIAL_MATRIX,&Ms); if (ierr) return ierr;
   ierr=MatAXPY(Ms,1,term2,SAME_NONZERO_PATTERN); if (ierr) return ierr;
   ierr=MatScale(Ms,0.5); if (ierr) return ierr;

   // Mc
   ierr=MatConvert(term1,MATSAME,MAT_INITIAL_MATRIX,&Mc); if (ierr) return ierr;
   ierr=MatAXPY(Mc,-1,term2,SAME_NONZERO_PATTERN); if (ierr) return ierr;
   ierr=MatScale(Mc,0.5); if (ierr) return ierr;

   // clean up
   ierr=MatDestroy(&Ki); if (ierr) return ierr;
   ierr=MatDestroy(&Kv); if (ierr) return ierr;
   ierr=MatDestroy(&rootZa); if (ierr) return ierr;
   ierr=MatDestroy(&rootZb); if (ierr) return ierr;
   ierr=MatDestroy(&invRootZa); if (ierr) return ierr;
   ierr=MatDestroy(&invRootZb); if (ierr) return ierr;
   ierr=MatDestroy(&term1); if (ierr) return ierr;
   ierr=MatDestroy(&term2); if (ierr) return ierr;

   return ierr;
}


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

   MatDestroy(&Mc);
   MatDestroy(&Ms);
}

