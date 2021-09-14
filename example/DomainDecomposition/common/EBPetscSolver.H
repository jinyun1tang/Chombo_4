#ifndef _chombo_ebpetscsolver_H_
#define _chombo_ebpetscsolver_H_

#ifdef CH_USE_PETSC
#include "PetscSolver.H"
#include "BaseDomainBC.H"
#include "petsc.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscviewer.h"
#include "Chombo_EBChombo.H"
#include "Chombo_EBLevelBoxData.H"

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
  template <int order>
  class EBPetscSolver : public PetscSolver< EBLevelBoxData<CELL, 1> >
  {
  public:
    typedef GraphConstructorFactory< EBBoxData<CELL, int, 1> >  devifactint_t;
    typedef GraphConstructorFactory<EBHostData<CELL, int, 1> >  hostfactint_t;
    typedef GraphConstructorFactory<EBHostData<CELL,Real, 1> > hostfactreal_t;
    typedef EBDictionary<order, Real, CELL, CELL>                dictionary_t;

    EBPetscSolver(const shared_ptr<GeometryService<order> >  & a_geoserv,
                  const shared_ptr<dictionary_t           >  & a_ebdictionary,
                  const shared_ptr<LevelData<EBGraph> >      & a_graphs,
                  const Chombo4::DisjointBoxLayout                    & a_grids,
                  const Chombo4::Box                                  & a_domain,
                  string a_stencilName,
                  string a_domainBCName[2*DIM],
                  string a_ebbcName,
                  string a_prestring,
                  Real a_dx,   Real a_alpha, Real a_beta, Point a_ghost)
    {
      m_geoserv      =  a_geoserv;
      m_ebdictionary =  a_ebdictionary;
      m_graphs       =  a_graphs;
      m_grids        =  a_grids;
      m_domain       =  a_domain;
      m_stencilName  =  a_stencilName;
      m_dx           =  a_dx;
      m_alpha        =  a_alpha;
      m_beta         =  a_beta;
      m_prestring    =  a_prestring;
      for(int iface = 0; iface < 2*DIM; iface++)
      {
        m_domainBCName[iface] = a_domainBCName[iface];
      }
      m_ebbcName = a_ebbcName;
    
      m_ivghost      =    ProtoCh::getIntVect(a_ghost);
      m_ptghost      =    a_ghost;
      //create map of locations in space to matrix row.
      defineGIDS();
      //create space for a matrix and the necessary vectors.
      createMatrixAndVectors();
      //put actual values into the matrix
      formMatrix();
      setupSolver();
    }

    virtual void define(LinearOp<EBLevelBoxData<CELL, 1> >* a_operator, bool a_homogeneous = false)
    {
      PetscSolver<EBLevelBoxData<CELL, 1> >::define(a_operator, a_homogeneous);
      m_op = static_cast<EBAMRPoissonOp*>(a_operator);
      //We need to be sure it is EBAMRPoissonOp that is being passed in
      // this is defined in AMRMultiGrid::define->a_factory.AMRnewOp with
      //How to do a check?

      //getting a pointer to the vofStencil
      m_op->getVoFStencil(m_vofStencil);
      m_op->getAlphaDiagWeight(m_alphaDiagWeight);
      m_eblg = (*m_op).getEBLG();
      m_finestEBISL= m_eblg.getEBISL();
      m_domainBC = m_op->getDomainBC();

      // m_op->getAlphaBeta(m_alpha,m_beta);//too early
    }

    Real addBCrhsValue(const IntVect& a_iv, const EBLevelBoxData<CELL, 1>& a_phi, const DataIndex& a_datInd, const Real& coeff = 1);
    Real addBCdiagValue(const IntVect a_iv, const IntVect a_jv, const EBLevelBoxData<CELL, 1>& a_rhs, const DataIndex& a_datInd, const Real coeff);
    
  private:    
    EBPetscSolver();
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
    
    shared_ptr<GeometryService<order> >  m_geoserv;
    shared_ptr< dictionary_t>            m_ebdictionary;
    shared_ptr<LevelData<EBGraph> >      m_graphs;
    Chombo4::DisjointBoxLayout                    m_grids;
    Chombo4::Box                                  m_domain;
    Real                                 m_dx;
  
    string m_stencilName;
    string m_domainBCName[2*DIM];
    string m_ebbcName;

    Vector<int> m_numPtsAllProc;
    Point   m_ptghost;
    IntVect m_ivghost;
    Mat m_mat;
    void *m_ctx; // pointer for nonlnear solver call backs

    Vec m_xx, m_rr, m_bb;
    SNES m_snes;
    KSP m_ksp;
    PetscInt m_defined;
    Real m_alpha, m_beta;
    
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
#include "EBPetscSolverI.H"
#endif
