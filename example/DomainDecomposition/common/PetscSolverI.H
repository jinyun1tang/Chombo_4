
#ifndef _PETSCSOLVERI_H_
#define _PETSCSOLVERI_H_
#ifdef CH_USE_PETSC

#include "Chombo_LevelData.H"
#include "Chombo_BoxIterator.H"

namespace Chombo4
{

// *******************************************************
// PetscSolver Implementation
// *******************************************************
  template <class T>
  PetscSolver<T>::PetscSolver()
    :m_homogeneous(false),
     m_mat(0), // m_xx, m_rr, m_bb;
     m_snes(0),
     m_ksp(0),
     m_defined(0),
     m_function(0),
     m_jacobian(0),
     m_null(false),
     m_nz_init_guess(false),
     m_gid0(0)
  {
    m_dx = 1.;
    m_prestring[0] = '\0';
  }
// *******************************************************
  template <class T>
  void PetscSolver<T>::destroy()
  {
    if ( m_defined )
    {
      if (m_mat)
      {
        MatDestroy(&m_mat);
        m_mat = 0;
      }
      VecDestroy(&m_bb);
      VecDestroy(&m_xx);
      VecDestroy(&m_rr);
      m_defined = 0;
    }
    if ( m_ksp )
    {
      KSPDestroy( &m_ksp );
      m_ksp = 0;
    }
    if ( m_snes )
    {
      SNESDestroy( &m_snes );
      m_snes = 0;
    }
  }
// *******************************************************
  template <class T>
  PetscSolver<T>::~PetscSolver()
  {
    destroy();
  }

// *******************************************************
  template <class T>
  void PetscSolver<T>::define( Real a_dx,
                               bool a_homogeneous )
  {
    m_homogeneous = a_homogeneous; // not used!!!
    m_dx = a_dx;
    CH_assert(m_dx!=0.);
  }

// *******************************************************
  template <class T>
  void PetscSolver<T>::define( LinearOp<T> *a_op,
                               bool a_homogeneous )
  {
    define( a_op->dx(), a_homogeneous);
  }

// *******************************************************
  template <class T>
  void PetscSolver<T>::setNull( bool n /*= true*/ )
  {
    m_null = n; CH_assert(0);
  }

// *******************************************************
  template <class T>
  void PetscSolver<T>::solve( T & a_phi,
                              const T & a_rhs )
  {
    T& phi = a_phi;
    const T& rhs = a_rhs;
    solve_private( phi, rhs );
  }


// *******************************************************
//   create_mat_vec
//     Create 'm_mat', 'm_xx', ....  Constructs 'm_gids'.
//
  template <class T>
  int PetscSolver<T>::create_mat_vec( const T& a_phi )
  {
    CH_TIME("PetscSolver::create_mat_vec");
    const DisjointBoxLayout &dbl = a_phi.disjointBoxLayout();
    DataIterator dit( dbl );
    PetscErrorCode ierr;
    const PetscInt nc = a_phi.nComp();
#ifdef CH_MPI
    MPI_Comm wcomm = Chombo_MPI::comm;
#else
    MPI_Comm wcomm = PETSC_COMM_SELF;
#endif

    if ( !m_mat )
    {
      // print_memory_line("Before AMG set up");
      m_defined = 2;
 
      IntVect idghosts = a_phi.ghostVect();
      if (idghosts == IntVect::Zero)
      {
        MayDay::Error("PetscSolver<T>::create_mat_vec: No ghost cells in input LevelData<>.");
      }
      m_bccode.define(dbl, 1, idghosts);

      // global ids with ghost cells
      m_gids.define(dbl, 1, idghosts);

      // get first (zero based) id on this processor      
      PetscInt data = 0;
      for ( dit = a_phi.dataIterator() ; dit.ok() ; ++dit )
      {
        const Box &box = dbl.get(dit()); // not ghosted
        data += box.numPts();
        BaseFab<PetscInt> &gidsFab = this->m_gids[dit]; 
        gidsFab.setVal(-1); // flag for BC eventually
      }
      const PetscInt NN = nc*data;
#ifdef CH_MPI
      PetscInt result;
      MPI_Datatype mtype;
      PetscDataTypeToMPIDataType(PETSC_INT,&mtype);
      MPI_Scan( &data, &result, 1, mtype, MPI_SUM, wcomm );
      m_gid0 = result - data;
#else
      m_gid0 = 0;
#endif
      PetscInt gid = m_gid0;
      for ( dit = a_phi.dataIterator() ; dit.ok() ; ++dit )
      {
        BaseFab<PetscInt> &gidsFab = this->m_gids[dit]; 
        const Box& box = dbl.get(dit());
        BoxIterator bit(box);
        for (bit.begin(); bit.ok(); bit.next(), gid++ )
        {
          IntVect iv = bit();
          gidsFab(iv,0) = gid;
        }
      }
      m_gids.exchange();

      // create matrix
      PetscInt nnzrow = getNNZPerRow();
      PetscInt *d_nnz=PETSC_NULL, *o_nnz=PETSC_NULL;
      
      ierr = MatCreate(wcomm,&m_mat);CHKERRQ(ierr);
      ierr = MatSetOptionsPrefix(m_mat,"");CHKERRQ(ierr);

      ierr = MatSetSizes(m_mat,NN,NN,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
      ierr = MatSetBlockSize(m_mat,nc);CHKERRQ(ierr);
      ierr = MatSetType(m_mat,MATAIJ);CHKERRQ(ierr);
      // ierr = MatSetOption(m_mat,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE) ;CHKERRQ(ierr);
      ierr = MatSetFromOptions( m_mat ); CHKERRQ(ierr);
      ierr = MatSeqAIJSetPreallocation(m_mat,nnzrow, d_nnz);CHKERRQ(ierr);
      ierr = MatMPIAIJSetPreallocation(m_mat,nnzrow, d_nnz, nnzrow/2, o_nnz);CHKERRQ(ierr);
      ierr = MatSetOption(m_mat,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE) ;CHKERRQ(ierr);
      
      if ( d_nnz )
      {
        ierr = PetscFree( d_nnz );  CHKERRQ(ierr);
        ierr = PetscFree( o_nnz );  CHKERRQ(ierr);
      }

      // create vectors
      ierr = MatCreateVecs(m_mat,&m_bb,&m_xx);CHKERRQ(ierr);
      //ierr = VecCreate( wcomm, &m_bb ); CHKERRQ(ierr);
      //ierr = VecSetFromOptions( m_bb ); CHKERRQ(ierr);
      //ierr = VecSetSizes( m_bb, NN, PETSC_DECIDE ); CHKERRQ(ierr);
      ierr = VecDuplicate( m_bb, &m_rr ); CHKERRQ(ierr);
      //ierr = VecDuplicate( m_bb, &m_xx ); CHKERRQ(ierr);
    }

    return 0;
  }

// *******************************************************
//    setup_solver
//      - creates solver if needed.  forms matrix, sets up KSP
//
  template <class T>
  int PetscSolver<T>::setup_solver( const T& a_phi )
  {
    CH_TIMERS("PetscSolver::setup_solver");
    CH_TIMER("solve-setup-1st", t1);
    CH_TIMER("solve-setup-rest", t2);

    const DisjointBoxLayout &dbl = a_phi.disjointBoxLayout();
    DataIterator dit( dbl );
    PetscErrorCode ierr;
    KSP ksp;
#ifdef CH_MPI
    MPI_Comm wcomm = Chombo_MPI::comm;
#else
    MPI_Comm wcomm = PETSC_COMM_SELF;
#endif

    if ( m_defined == 0 )
    {
      m_defined = 2;
      // print_memory_line("Before AMG set up");

      ierr = create_mat_vec( a_phi ); CHKERRQ(ierr);

      CH_START(t1);

      // Add values to A
      ierr = formMatrix( m_mat, &a_phi ); CHKERRQ(ierr);
      {
        char str[256];
        strcpy (str,"-");
        strcat (str,m_prestring);
#if PETSC_VERSION_GE(3,6,0)
        strcat (str,"pc_gamg_square_graph 20");
#else
        strcat (str,"pc_gamg_square_graph true");
#endif
#if PETSC_VERSION_GE(3,7,0)
        ierr = PetscOptionsInsertString(PETSC_NULL,str);CHKERRQ(ierr);
#else
        ierr = PetscOptionsInsertString(str);CHKERRQ(ierr);
#endif
      }

      // create solvers
      PetscBool ism = PETSC_FALSE;
#if PETSC_VERSION_GE(3,7,0)
      PetscOptionsGetBool(PETSC_NULL,m_prestring,"-ksp_monitor",&ism,PETSC_NULL);
#else 
      PetscOptionsGetBool(m_prestring,"-ksp_monitor",&ism,PETSC_NULL);
#endif
      if ( m_function && !m_snes )
      {
        ierr = SNESCreate( wcomm, &m_snes );                                       CHKERRQ(ierr);
        // these casts get rid of the 'const'. These are used in SNES as temp vector!  
        ierr = SNESSetFunction( m_snes, m_rr, m_function, (void*)this );         CHKERRQ(ierr);
        ierr = SNESSetJacobian( m_snes, m_mat, m_mat, m_jacobian, (void*)this ); CHKERRQ(ierr);
        //ierr = SNESSetApplicationContext( m_snes, (void*)this );                   CHKERRQ(ierr);
        ierr = SNESSetFromOptions( m_snes );                                       CHKERRQ(ierr);
        ierr = SNESGetKSP( m_snes, &ksp );                                         CHKERRQ(ierr);
        if (ism)
        {
          ierr = KSPMonitorSet(ksp,ksp_monitor_pout,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
        }
      }
      else if ( !m_function && !m_ksp )
      {
        // create the KSP so that we can set KSP parameters
        KSPCreate( wcomm, &m_ksp );
        if ( strlen(m_prestring) > 0 )
        {
          ierr = KSPSetOptionsPrefix( m_ksp, m_prestring );    CHKERRQ(ierr);
        }
        ierr = KSPSetFromOptions(m_ksp);CHKERRQ(ierr);
        if (ism)
        {
          ierr = KSPMonitorSet(m_ksp,ksp_monitor_pout,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
        }
#if PETSC_VERSION_GE(3,5,0)
        ierr = KSPSetOperators(m_ksp,m_mat,m_mat);CHKERRQ(ierr);
#else
        ierr = KSPSetOperators(m_ksp,m_mat,m_mat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
#endif
        ierr = KSPSetInitialGuessNonzero(m_ksp, m_nz_init_guess ? PETSC_TRUE : PETSC_FALSE );CHKERRQ(ierr);
        ksp = m_ksp;
      }
      else CH_assert(0);

      { // coordinates
        PC pc; 
        PetscInt sz,ind,bs,n,m;
#if PETSC_VERSION_LT(3,4,0) & PETSC_VERSION_RELEASE
        const PCType type;
#else
        PCType type;
#endif  
        ierr = KSPGetPC( ksp, &pc );     CHKERRQ(ierr);
        ierr = PCGetType( pc, &type );    CHKERRQ(ierr);
        ierr = MatGetBlockSize( m_mat, &bs );               CHKERRQ( ierr );
        if ( strcmp(type,PCGAMG) == 0 && bs > 1 )
        {
          PetscReal    *coords;
          DataIterator dit(a_phi.disjointBoxLayout());
        
          ierr = MatGetLocalSize( m_mat, &m, &n );  CHKERRQ(ierr);
          sz = CH_SPACEDIM*(m/bs);
          ierr = PetscMalloc( (sz+1)*sizeof(PetscReal), &coords ); CHKERRQ(ierr);
          for ( dit = a_phi.dataIterator(), ind = 0 ; dit.ok() ; ++dit )
          {
            const Box &box = a_phi.getBoxes()[dit];
            BoxIterator bit(box);
            for (bit.begin(); bit.ok(); bit.next())
            {
              IntVect iv = bit(); // coordinate in any scaled, shifted, rotated frame.
              for (PetscInt k=0; k<CH_SPACEDIM; k++) coords[ind++] = (PetscReal)iv[k];
            }
          }
          CH_assert(ind==sz);
          ierr = PCSetCoordinates( pc, CH_SPACEDIM, sz/CH_SPACEDIM, coords ); CHKERRQ(ierr);
          ierr = PetscFree( coords );  CHKERRQ(ierr);
        }
      }

      CH_STOP(t1);
      // print_memory_line("After AMG set up");
    }
    else if ( m_defined == 1 )
    {
      m_defined = 2;
      // form A -- m_mat
      CH_START(t2);
      ierr = MatZeroEntries( m_mat );        CHKERRQ(ierr);
      ierr = formMatrix( m_mat, &a_phi );     CHKERRQ(ierr);
      if ( m_ksp )
      {
        ksp = m_ksp;
      }
      else
      {
        ierr = SNESGetKSP( m_snes, &ksp ); CHKERRQ(ierr);
      }
#if PETSC_VERSION_GE(3,5,0)
      ierr = KSPSetOperators(ksp,m_mat,m_mat); CHKERRQ(ierr);
#else
      ierr = KSPSetOperators(ksp,m_mat,m_mat,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
#endif
      CH_STOP(t2);
    }

    return 0;
  }

#endif
