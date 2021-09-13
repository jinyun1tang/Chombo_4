#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _chombo_ebpetscsolver_H_
#define _chombo_ebpetscsolver_H_

#ifdef CH_USE_PETSC
#include "PetscSolver.H"
#include "BaseDomainBC.H"
#include "petsc.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscviewer.h"

namespace Chombo4
{

///
/**
   EBCellFAB version of the PetscLinearSolverPoisson, iterative solver which only uses the operator's preconditioner. Form a matrix in Petsc mat format which includes only domain cells, EB or non-EB, at the finest level (so that there is no multi-valued cells). Boundary conditions which are coded in a_phi boundary cells are then translated into the RHS of the matrix.
   Handles two types of Boundary conditions: (1) Domain BC (2) EBBC
   (1) Domain BC, (a) for Dirichlet, the ghost cell has been evaluated outside in applyDomainFlux. Since we do u_B = u_(i+1/2) = (u_i+u_i+1)/2, there are two parts that need to be bring into the matrix, the diagonal term and the RHS. They are done, respectively, in addBCdiagValue and addBCrhsValue; (b) for Neumann, u_i == u_i+1, only diag needs to be changed (addBCdiagValue)
   (2) EBBC are completely handled in opStencil. Even at the domain boundary, BC is written in flux form and encoded in opStencil. Therefore nothing needs to be done for the RHS.
   Therefore, when we have a regular domain boundary we set bccode=true, when it is irregular bccode=false
   diag terms are treated in formMatrix; RHS needs to be done each time step in solveprivate;

*/

////////////////////////////////////////////////////////////////////////
// an instantialtion of PetscLinearSolver with Poisson
//     solve: alpha u - div( beta grad u) = f
//         - alpha and beta are constants
////////////////////////////////////////////////////////////////////////
  class EBPetscSolver : public PetscSolver< EBLevelBoxData<CELL, 1> >
  {
  public:
    EBPetscSolver();
    virtual void define(LinearOp<EBLevelBoxData<CELL, 1> >* a_operator, bool a_homogeneous = false);
    BaseFab<Real>& getRegFab(EBLevelBoxData<CELL, 1>& a_fab, const DataIndex& a_datInd)
    {
      BaseFab<Real>& fabref = a_fab[a_datInd].getSingleValuedFAB();
      return fabref;
    }
    const BaseFab<Real>& getRegFab(const EBLevelBoxData<CELL, 1>& a_fab, const DataIndex& a_datInd) const
    {
      const BaseFab<Real>& fabref = a_fab[a_datInd].getSingleValuedFAB();
      return fabref;
    }
    const BaseFab<Real>& getRegFab(const EBLevelBoxData<CELL, 1>& a_fab, const DataIndex& a_datInd, Box& a_box) const
    {
      const EBCellFAB& ebfab = a_fab[a_datInd];
      const BaseFab<Real>& fabref = ebfab.getSingleValuedFAB();
      a_box = ebfab.getRegion();
      return fabref;
    }
    virtual void defineData(EBLevelBoxData<CELL, 1>& a_data, const EBLevelBoxData<CELL, 1>& a_phi)
    {
      IntVect idghosts = a_phi.ghostVect();
      const DisjointBoxLayout &dbl = a_phi.disjointBoxLayout();
      EBCellFactory ebcellfact(m_finestEBISL);
      a_data.define(dbl, 1, idghosts, ebcellfact);
    }
    Real addBCrhsValue(const IntVect& a_iv, const EBLevelBoxData<CELL, 1>& a_phi, const DataIndex& a_datInd, const Real& coeff = 1);
    Real addBCdiagValue(const IntVect a_iv, const IntVect a_jv, const EBLevelBoxData<CELL, 1>& a_rhs, const DataIndex& a_datInd, const Real coeff);
  protected:
    ///
/**
 // Limitations/Assumptions:
 // Only for Pressure (MAC or CC Projector)
 // Assumes Neumann BCs have derivative == 0
 */
    virtual int formMatrix( Mat a_Mat, const EBLevelBoxData<CELL, 1> *a_rhs=0,
                            PetscInt n=0, PetscInt z=0, 
                            PetscInt *d_nnz=0, PetscInt *o_nnz=0 );
  public:
    Real m_alpha;
    Real m_beta;
    EBISLayout m_finestEBISL;
  private:
    //making them pointers since they only link to the items from EBAMRPoissonOp
    //these pointers are made const* since they shall not modify whoever they are pointing to, but can change their targets.
    LayoutData<BaseIVFAB<VoFStencil> > const* m_vofStencil;
    LayoutData<BaseIVFAB<Real> >       const* m_alphaDiagWeight;
    EBLevelGrid m_eblg;
    EBAMRPoissonOp* m_op;
    //m_rhsAdd store the values that need to be added to RHS due to BCs
    RefCountedPtr<BaseDomainBC> m_domainBC;
  public:
    /**
       get an estimate of the number of nnz/row for matrix allocation
    */
    virtual int getNNZPerRow() const
    {
      if (CH_SPACEDIM==2) return 50;
      else return 500;
    }
    /**
       do I support having formMatrix precompute exact NNZ/row
    */
    virtual bool supportNNZExact() const
    {
      return true;
    }
  };
}
#endif
#endif /*_PETSCSOLVER_H_*/
