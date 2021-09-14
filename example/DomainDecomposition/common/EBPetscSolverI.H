
#ifdef CH_USE_PETSC

namespace Chombo4
{

// *******************************************************

///
  template <unsigned int order>
  int
  EBPetscSolver<order>::
  formMatrix( Mat a_mat, const EBLevelBoxData<CELL, 1> *a_phi,
              PetscInt my0, PetscInt nloc,
              PetscInt *d_nnz, PetscInt *o_nnz )
  {

    CH_TIME("EBPetscSolver::formMatrix");

    const DisjointBoxLayout& dbl = this->m_gids.disjointBoxLayout();
    CH_assert(this->m_dx!=0.0);
  
    PetscErrorCode ierr;
    PetscInt nc=1,nreg=0;
    m_op->getAlphaBeta(this->m_alpha,this->m_beta);
    Real idx2 = 1.e0/(this->m_dx*this->m_dx) * this->m_beta;
    Real vo = 1.e0 * idx2;
    Real va = -vo;
    Real vd = -CH_SPACEDIM * 2.e0 * idx2 + this->m_alpha;
    //Real nil=  1e-16 * (abs(vd)+abs(vo));
    //GHM sign change
    Real nil=  1e-16 * (abs(vd)+abs(vo)) * ( vd < 0 ? -1 : 1 );
    // Real zero=0.;
    Real diagWeight =0.;
 

    DataIterator dit = this->m_gids.dataIterator();
    int nbox=dit.size();
    for (int mybox=0; mybox<nbox; mybox++)
    {
      const DataIndex& datInd = dit[mybox];
      const BaseFab<PetscInt> &gidsFab = this->m_gids[datInd]; //getRegFab(this->m_gids,dit);
      
      //for EB
      const EBCellFAB& ebfab = (*a_phi)[datInd]; // we could use m_gids for this and get rid of a_phi
      const EBISBox& ebbox   = ebfab.getEBISBox();
      const Box& box = dbl.get(datInd);
      const IntVectSet ivsIrreg = ebbox.getIrregIVS(box);
      const BaseIVFAB<VoFStencil>& vofSfab = (*m_vofStencil)[datInd];
      const BaseIVFAB<Real>& alphaDiagWeight = (*m_alphaDiagWeight)[datInd];
      BaseFab<bool> &bccode = m_bccode[datInd];
      // EBCellFAB & rhsAdd = m_rhsAdd[datInd];

      BoxIterator bit(box);
      for(bit.begin(); bit.ok(); bit.next())
      { // this loop only does the regular cells
        IntVect iv = bit();
        PetscInt i = nc*gidsFab(iv,0);
        if( d_nnz ) d_nnz[i-my0]++; // diag
        if (ebbox.isRegular(iv))
        {
          nreg++;
          // diag
          if( !d_nnz )
          {
            ierr = MatSetValues(a_mat,1,&i,1,&i,&vd,ADD_VALUES); CHKERRQ(ierr); 
          }	      

          for( int dim_dir = 0 ; dim_dir < CH_SPACEDIM ; dim_dir++ )
          {
            for( int i_dir = -1 ; i_dir <= 1 ; i_dir += 2 )
            {
              IntVect shiftiv = IntVect::Zero ;
              shiftiv[dim_dir] = i_dir;
              IntVect jv = iv + shiftiv;
              PetscInt j = nc*gidsFab(jv,0);
              if( j >= 0 )
              {
                if( !d_nnz )
                {
                  ierr = MatSetValues(a_mat,1,&i,1,&j,&vo,ADD_VALUES); CHKERRQ(ierr);
                  //ierr = MatSetValues(a_mat,1,&j,1,&i,&zero,ADD_VALUES); CHKERRQ(ierr);
                }
                else
                {
                  if( j < my0 || j >= my0+nloc ) o_nnz[i-my0]++;
                  else if (i != j)
                  {
                    d_nnz[i-my0]++;
                  }
                }
              }
              else if( !d_nnz )
              {
                va = addBCdiagValue(iv,jv, *a_phi,datInd,vo);
                ierr = MatSetValues(a_mat,1,&i,1,&i,&va,ADD_VALUES); CHKERRQ(ierr);
                bccode(iv,0) = true;
                // vo *= -1.;
              }
            }
          }
        }
        // else if (ebbox.isCovered(iv))
        //   { 
        //     // otherwise, for some reason it breaks with more CPUs
        //     ierr = MatSetValues(this->m_mat,1,&i,1,&i,&nil,ADD_VALUES); CHKERRQ(ierr); 
        //   }
        else if( !d_nnz ) // !isRegular -- add something everywhere.
        {
          ierr = MatSetValues(a_mat,1,&i,1,&i,&nil,ADD_VALUES); CHKERRQ(ierr); // 
        }            
      }
      
      //Now do the irregular cells
      IVSIterator ivsit (ivsIrreg);
      for (ivsit.begin(); ivsit.ok(); ++ivsit)
      {
        IntVect iv = ivsit();
        VolIndex vofc = VolIndex(iv,0);
        PetscInt i = nc*gidsFab(iv,0);
        const VoFStencil& vofSten = vofSfab(vofc,0);
        for (PetscInt k = 0; k < vofSten.size(); k++)
        {
          VolIndex vof = vofSten.vof(k);
          IntVect jv  = vof.gridIndex();
          Real w = vofSten.weight(k) * this->m_beta;
          PetscInt j = nc*gidsFab(jv,0);
          if (j >= 0)
          {
            PetscInt j = nc*gidsFab(jv,0);
            if( !d_nnz )
            {
              ierr = MatSetValues(a_mat,1,&i,1,&j,&w,ADD_VALUES); CHKERRQ(ierr);
              //ierr = MatSetValues(a_mat,1,&j,1,&i,&zero,ADD_VALUES); CHKERRQ(ierr);
            }
            else
            {
              if( j < my0 || j >= my0+nloc ) o_nnz[i-my0]++; 
              else if (i != j)
              {
                d_nnz[i-my0]++;
                //d_nnz[j-my0]++; // this is to be safe - stencil is not symm.
              }
            }
          }
          else
          { // this is a physical domain ghost cell
            // now we are dealing with a domain boundary 
            MayDay::Error("EB accessing domain ghost cells is not anticipated");
          }
        }
        if( !d_nnz )
        {
          diagWeight = this->m_alpha* alphaDiagWeight(vofc,0);
          ierr = MatSetValues(a_mat,1,&i,1,&i,&diagWeight,ADD_VALUES); CHKERRQ(ierr);
          bccode(iv,0)=false;
        }
      }
    }
  
    //applyInhomDomBC();
    if( !d_nnz )
    {	  
      ierr = MatAssemblyBegin(a_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); 
      ierr = MatAssemblyEnd(a_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); 
    }
    else if(true)
    {
      const PetscInt num_b = 20;
      PetscInt ob[num_b],ib[num_b],m,n;
      for(m=0;m<num_b-1;m++) ob[m]=0; // bins
      for(PetscInt k=0;k<nloc;k++) {
	for(m=0,n=d_nnz[k]+o_nnz[k];m<num_b-1;m++) {
	  n = n/2;
	  if(n==0) { ob[m]++;  break; }
	}
	if(m==num_b-1) pout() << "HIT limit: " << (d_nnz[k]+o_nnz[k]) << endl;
      }
      ob[num_b-1] = nreg*nc; // real cells
      ob[num_b-2] = nloc;    // real cells
#ifdef CH_MPI
      MPI_Datatype mtype;
      PetscDataTypeToMPIDataType(PETSC_INT,&mtype);
      MPI_Allreduce(ob,ib,num_b,mtype,MPI_SUM,PETSC_COMM_WORLD);
#else
      for(m=0;m<num_b;m++) ib[m] = ob[m]; 
#endif
      pout() << "NNZ bins: (regular cells = " << ib[num_b-1]/nc << ", " << (Real)ib[num_b-1]/(Real)(ib[num_b-2]/nc)*100. << " % regular)" << endl;
      pout() << "NNZ bin counts ";
      for(m=0;m<num_b-3;m++) pout() << ib[m]<< ", " ;
      pout()  << ib[num_b-3] << endl;
      for(m=0;m<num_b-3;m++) pout() << ob[m]<< ", " ;
      pout()  << ob[num_b-3] << endl;
    }

    return 0;
  }

////////////////////
  template <unsigned int order>
  Real
  EBPetscSolver<order>::
  addBCdiagValue(const IntVect a_iv, const IntVect a_jv, const EBLevelBoxData<CELL, 1>& a_phi, const DataIndex& a_datInd, const Real coeff)
  //Bring BC into diag
  //rhsAdd may get modified in someone else's addBCdiagValue, but not here
  // because rhsAdd value in our settings depends on a_phi from time step to time step
  // it has to be done 
  {
    const EBCellFAB& rhsfab = a_phi[a_datInd];
    const VolIndex ivof(a_iv,0);
    const VolIndex jvof(a_jv,0);
    const bool isDirichlet = (*m_domainBC).isDirichletDom(ivof,jvof,rhsfab);
    if (isDirichlet)
    {
      
      return -coeff;
    }
    else
    {
      return coeff;
    }
  }
 
//
// EBPetscSolver::addBCrhsValue
//
  template <unsigned int order>
  Real 
  EBPetscSolver<order>::
  addBCrhsValue(const IntVect& a_iv, const EBLevelBoxData<CELL, 1>& a_phi, const DataIndex& a_datInd, const Real& coeff)
  //Bring BC into RHS of matrix
  //one cell can have multiple types of BCs if they are at the domain corner
  {
    Real rhsAdd= 0.0;
    int nc = 1;
    const IntVect& iv = a_iv;
    const EBCellFAB &phiFab = a_phi[a_datInd];
    const EBISBox& ebbox   = phiFab.getEBISBox();
    const BaseFab<Real>& regPhiFAB = phiFab.getSingleValuedFAB();
    const BaseFab<PetscInt> &gidsFab = this->m_gids[a_datInd]; 
    for( int dim_dir = 0 ; dim_dir < CH_SPACEDIM ; dim_dir++ )
    {
      for( int i_dir = -1 ; i_dir <= 1 ; i_dir += 2 )
      {
        IntVect shiftiv = IntVect::Zero ;
        shiftiv[dim_dir] = i_dir;
        IntVect jv = iv + shiftiv;
        PetscInt j = nc*gidsFab(jv,0);
        if( j < 0 )
        {
          const VolIndex ivof(iv,0);
          const VolIndex jvof(jv,0);
          const bool isDirichlet = (*m_domainBC).isDirichletDom(ivof,jvof,phiFab);
          if (isDirichlet)
          {
            //1. reverse-engneering to get BC value
            //this is because whatever the BC is, it has been coded in the ghost cell, for dirichlet, it has been figured out to be phiB = (a_phiFab(ivof,0)+a_phiFab(jvof,0))/2
            Real phiB = regPhiFAB(iv, 0);
            phiB = (regPhiFAB(iv, 0) + regPhiFAB(jv, 0))/2.0;
            //  Real phiB    = (phiFab(ivof,0)+phiFab(jvof,0))/2.0;
            //2. throw it to the RHS
            //rhsAdd -= 2.0 * coeff * m_beta * phiB; 
            // if(phiB != 0.)
            //   {
            //     pout() << "non-zero add rhs" << endl;
            //   }
            rhsAdd -= 2.0 * coeff * m_beta * phiB; 
            //m_beta is due to helmholtz
            // in projection step, this is automatically 1
            // in PetscLinearSolverI.H, coeff doesn't know m_beta
          }
          else
          {
            // nothing need to be done to RHS for Neumann (and it always ==0)
          }
        }
      }
    }
    if (ebbox.isCovered(iv))
    {
      rhsAdd = BaseFabRealSetVal; //std::numeric_limits<double>::quiet_NaN();
    }
    return rhsAdd;
    // EB Cells are not supposed to access outside cells so there is nothing that needs to be added to RHS, it has already been coded into opStencil
  }
}
#endif //CH_USE_PETSC
