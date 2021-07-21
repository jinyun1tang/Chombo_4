#ifndef __hoeb_utilities__
#define __hoeb_utilities__
#include "Proto.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_GeometryService.H"
#include "Chombo_EBDictionary.H"
#include "Hoeb_LAPACKMatrix.H"
#include "Hoeb_Neighborhood.H"
#include "Hoeb_Data.H"


namespace hoeb
{
  ///
  inline RealVect
  getVoFLocation(const EBIndex<CELL>& a_vof,
                 Real                 a_dx)
  {
    RealVect retval;
    auto pt = a_vof.m_pt; 
    for(int idir = 0; idir < DIM; idir++)
    {
      retval[idir] = a_dx*(0.5 + Real(pt[idir]));
    }
    return retval;
  }
  ///
  template <CENTERING cent>
  inline RealVect
  getFaceLocation(const EBIndex<cent>& a_face,
                  Real                 a_dx,
                  int                  a_facedir)
  {
    RealVect retval;
    auto pt = a_face.m_pt; //high size of the face
    //this will work for cell centers or boundary faces as well
    for(int idir = 0; idir < DIM; idir++)
    {
      if(idir == a_facedir)
      {
        retval[idir] = a_dx*Real(pt[idir]);
      }
      else
      {
        retval[idir] = a_dx*(0.5 + Real(pt[idir]));
      }
    }
    return retval;
  }
  
  ///Annoying translation between data structures I do too much
  inline 
  IndexTM<Real, DIM> getIndexDistance(const RealVect& a_vofloc,
                                      const RealVect& a_xbar)
  {
    IndexTM<Real, DIM> distance;
    for(int idir = 0; idir < DIM; idir++)
    {
      distance[idir] = a_vofloc[idir] - a_xbar[idir];
    }
    return distance;
  }


  ////
  template <unsigned int OPERATOR_ORDER, unsigned int GEOMETRY_ORDER>
  inline void
  shiftMomentAndFillRow(LAPACKMatrix                             & a_mat,
                        IndexedMoments<     DIM, GEOMETRY_ORDER> & a_bcmom,
                        const IndexTM<Real, DIM>                 & a_xbardistance,
                        const unsigned int                       & a_currentRow)
  {
                          
    //shift so it has the same xbar
    a_bcmom.shift(a_xbardistance);
    for(MomentIterator<DIM, OPERATOR_ORDER> momit; momit.ok(); ++momit)
    {
      auto momind = momit();
      unsigned int currentCol = IndexedMoments<DIM, OPERATOR_ORDER>::indexOf(momind);
      a_mat(a_currentRow, currentCol) = a_bcmom[momind];
    }
  }
  /////
  template <unsigned int OPERATOR_ORDER, unsigned int GEOMETRY_ORDER>
  inline void
  fillMatrixRowForEBBC(LAPACKMatrix                                        & a_mat,
                       const shared_ptr< GeometryService<GEOMETRY_ORDER> > & a_geoserv,
                       const Chombo4::Box                                  & a_srcDomain,
                       const EBIndex<BOUNDARY>                             & a_face,
                       const EBGraph                                       & a_graph,
                       const IndexTM<Real, DIM>                            & a_xbardistance,
                       const Chombo4::DataIndex                            & a_datind,
                       const unsigned int                                  & a_currentRow,
                       int a_facedir, Real a_dx, string a_bcname)
  {
    IndexedMoments<    DIM, GEOMETRY_ORDER> bcmom;
    if(a_bcname == string("Dirichlet"))
    {
      bcmom =  hoeb_data::getEBIndexedMoment<GEOMETRY_ORDER>
        (a_face, a_geoserv, a_srcDomain,a_datind, a_graph, a_dx);
    }
    else if(a_bcname == string("Neumann"))
    {
      bcmom.setVal(0.);
      for(unsigned int normDir = 0; normDir < DIM; normDir++)
      {
        for(MomentIterator<DIM, OPERATOR_ORDER> momit; momit.ok(); ++momit)
        {
          auto momind = momit();
          Real bcval = 0;
          
          IndexedMoments<DIM, GEOMETRY_ORDER> normMom=
            hoeb_data::getEBNormalMoment<GEOMETRY_ORDER>
            (a_face, a_geoserv, a_srcDomain,a_datind, a_graph, normDir, a_dx);

          if(momind[normDir] > 0)
          {
            int pd = momind[normDir];
            IndexTM<int, DIM> derind =momind - IndexTM<int, DIM>::BASISV(normDir);
            Real derval = normMom[derind];
            bcmom[momind] += pd*derval;
          }
        }
      }
    }
    shiftMomentAndFillRow<OPERATOR_ORDER, GEOMETRY_ORDER>
      (a_mat, bcmom, a_xbardistance, a_currentRow);
  }

  ///
  template <unsigned int D>
  inline IndexTM<int, D-1>
  getFaceIndexFromVol(const  IndexTM<int, D>& a_volind, int a_facedir)
  {
    IndexTM<int, D-1> retval;
    int ifac = 0;
    for(int ivol = 0; ivol < D; ivol++)
    {
      if(ivol != a_facedir)
      {
        retval[ifac] = a_volind[ivol];
        ifac++;
      }
    }
    return retval;
  }
  
  ///
  template <unsigned int D, unsigned int GEOMETRY_ORDER>
  inline IndexedMoments<D, GEOMETRY_ORDER>
  getVolMomFromFacMom(const IndexedMoments<D-1, GEOMETRY_ORDER>& a_facmom,
                      int a_facedir)
  {
    IndexedMoments<D, GEOMETRY_ORDER> retval;
    for(MomentIterator< D, GEOMETRY_ORDER> momit; momit.ok(); ++momit)
    {
      IndexTM<int, D> volpow = momit();
      Real momval;
      if(volpow[a_facedir] > 0) //no derivs in this direction
      {
        momval = 0;
      }
      else
      {
        IndexTM<int, D-1> facpow = getFaceIndexFromVol<DIM>(volpow, a_facedir);
        momval = a_facmom[facpow];
      }
      retval[volpow] = momval;
    }
    return retval;
  }
  ///
  template <CENTERING cent, unsigned int OPERATOR_ORDER, unsigned int GEOMETRY_ORDER>
  inline void
  fillMatrixRowForDomBC
  (LAPACKMatrix                                                         & a_mat,
   const HostIrregData<cent, IndexedMoments<DIM-1, GEOMETRY_ORDER>, 1 > & a_facdatfab,
   const EBIndex<cent>                                                  & a_face,
   const EBGraph                                                        & a_graph,
   const IndexTM<Real, DIM>                                             & a_xbardistance,
   const unsigned int                                                   & a_currentRow,
   int a_facedir, Real  a_dx, string a_bcname)
  {
                       
    IndexedMoments<DIM-1, GEOMETRY_ORDER> facmom =
      hoeb_data::getIndexedMoment<DIM-1, cent, GEOMETRY_ORDER>(a_face, a_facdatfab, a_graph, a_dx);
    //need to translate this into higher dim because that is what the matrix sees
    int facedir = 0;
    IndexedMoments<DIM, GEOMETRY_ORDER> volmom = getVolMomFromFacMom<DIM, GEOMETRY_ORDER>(facmom, facedir);
    IndexedMoments<    DIM, GEOMETRY_ORDER> bcmom;
    if(a_bcname == string("Dirichlet"))
    {
      bcmom = volmom;
    }
    else if(a_bcname == string("Neumann"))
    {
      for(MomentIterator<DIM, OPERATOR_ORDER> momit; momit.ok(); ++momit)
      {
        auto momind = momit();
        Real bcval = 0;
        if(momind[a_facedir] > 0)
        {
          int pd = momind[a_facedir];
          IndexTM<int, DIM> derind =momind - IndexTM<int, DIM>::BASISV(a_facedir);
          Real derval = volmom[derind];
          bcval = pd*derval;
        }
        bcmom[momind] = bcval;
      }
    }
    else
    {
      PROTO_ASSERT(false, "fillmatrixrow: bogus bc name");
    }
    shiftMomentAndFillRow<OPERATOR_ORDER, GEOMETRY_ORDER>
      (a_mat, bcmom, a_xbardistance, a_currentRow);
  }


  ///
  /******/  
  template <CENTERING cent, unsigned int OPERATOR_ORDER, unsigned int GEOMETRY_ORDER>
  inline LAPACKMatrix
  getMomentMatrix(LAPACKMatrix                                        & a_weights,
                  const Neighborhood<cent>                            & a_block,
                  const shared_ptr< GeometryService<GEOMETRY_ORDER> > & a_geoserv,
                  const Chombo4::Box                                  & a_srcDomain,
                  const RealVect                                      & a_xbar,
                  unsigned int                                          a_ibox,
                  Real                                                  a_dx)
  {

    unsigned int Pcols = IndexedMoments<DIM  , OPERATOR_ORDER>::size();
    unsigned int Nrows = a_block.size();
    LAPACKMatrix retval(Nrows, Pcols);
    a_weights.define(Nrows, Nrows);
    a_weights.setVal(0.);
    Real cellvol = 1;
    for(int idir = 0; idir < DIM; idir++)
    {
      cellvol *= a_dx;
    }
    auto dit = a_geoserv->getDBL(a_srcDomain).dataIterator();
    
    const auto & graphsldptr =  a_geoserv->getGraphs(  a_srcDomain);
    const auto & volmomldptr =  a_geoserv->getVoluData(a_srcDomain);
    
    const auto & graph     = (*graphsldptr)[dit[a_ibox]];
    const auto & volmomfab = (*volmomldptr)[dit[a_ibox]];

    unsigned int currentRow = 0;
    {
      //volumes
      for(unsigned int ivol = 0; ivol < a_block.m_volumes.size(); ivol++)
      {
        EBIndex<cent> vof    = a_block.m_volumes[ivol];
        RealVect vofloc = hoeb::getVoFLocation(vof, a_dx);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        IndexedMoments<DIM  , GEOMETRY_ORDER> volmom = hoeb_data::
          getIndexedMoment<DIM, cent, GEOMETRY_ORDER>(vof, volmomfab, graph, a_dx);
        
        Real volmag = volmom[IndexTM<int, DIM>::Zero];
        if(volmag > 1.0e-6*cellvol)
        {
          volmom.shift(distance);
          for(MomentIterator<DIM, OPERATOR_ORDER> momit; momit.ok(); ++momit)
          {
            auto momind = momit();
            unsigned int currentCol = IndexedMoments<DIM, OPERATOR_ORDER>::indexOf(momind);
            retval(currentRow, currentCol) = volmom[momind]/volmag;
          }
        }
        else
        {
          //for tiny cells, set moments to their limits (1 0 0 0...)
          retval.setSmallCellRow(currentRow);
        }
        a_weights(currentRow, currentRow) = a_block.m_volumesWgt[currentRow];
        currentRow++;
      }
    }
    {
      //now for the  domain x faces
      const auto & xfacldptr  =  a_geoserv->getXFaceData(a_srcDomain);
      const auto & xfacdatfab = (*xfacldptr)[dit[a_ibox]];
      for(unsigned int ivol = 0; ivol < a_block.m_xfacesDiri.size(); ivol++)
      {
        const auto & face = a_block.m_xfacesDiri[ivol];
        RealVect vofloc = getFaceLocation<XFACE>(face, a_dx, 0);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<XFACE, OPERATOR_ORDER, GEOMETRY_ORDER>
          (retval, xfacdatfab, face, graph, distance,
           currentRow, 0, a_dx, string("Dirichlet"));

        a_weights(currentRow,currentRow) = a_block.m_xfacesDiriWgt[ivol];
        currentRow++;
      }
      for(unsigned int ivol = 0; ivol < a_block.m_xfacesNeum.size(); ivol++)
      {
        const auto & face = a_block.m_xfacesNeum[ivol];
        RealVect vofloc = getFaceLocation<XFACE>(face, a_dx, 0);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<XFACE, OPERATOR_ORDER, GEOMETRY_ORDER>
          (retval, xfacdatfab, face, graph, distance,
           currentRow, 0, a_dx, string("Neumann"));

        a_weights(currentRow,currentRow) = a_block.m_xfacesNeumWgt[ivol];
        currentRow++;
      }
    }
    {
      //now for the  domain y faces
      const auto & yfacldptr  =  a_geoserv->getYFaceData(a_srcDomain);
      const auto & yfacdatfab = (*yfacldptr)[dit[a_ibox]];
      for(unsigned int ivol = 0; ivol < a_block.m_yfacesDiri.size(); ivol++)
      {
        const auto & face = a_block.m_yfacesDiri[ivol];
        RealVect vofloc = getFaceLocation<YFACE>(face, a_dx, 1);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<YFACE, OPERATOR_ORDER, GEOMETRY_ORDER>
          (retval, yfacdatfab, face, graph, distance,
           currentRow, 1, a_dx, string("Dirichlet"));
      
        a_weights(currentRow, currentRow) = a_block.m_yfacesDiriWgt[ivol];
        currentRow++;
      }
      for(unsigned int ivol = 0; ivol < a_block.m_yfacesNeum.size(); ivol++)
      {
        const auto & face = a_block.m_yfacesNeum[ivol];
        RealVect vofloc = getFaceLocation<YFACE>(face, a_dx, 1);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<YFACE, OPERATOR_ORDER, GEOMETRY_ORDER>
          (retval, yfacdatfab, face, graph, distance,
           currentRow, 1, a_dx, string("Neumann"));
      
        a_weights(currentRow, currentRow) = a_block.m_yfacesNeumWgt[ivol];
        currentRow++;
      }
    }
#if DIM==3
    {
      //now for the  domain z faces
      const auto & zfacldptr  =  a_geoserv->getYFaceData(a_srcDomain);
      const auto & zfacdatfab = (*yfacldptr)[dit[a_ibox]];
      for(unsigned int ivol = 0; ivol < a_block.m_zfacesDiri.size(); ivol++)
      {
        const auto & face = a_block.m_zFacesDiri[ivol];
        RealVect vofloc = getFaceLocation<ZFACE>(face, a_dx, 2);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<ZFACE, OPERATOR_ORDER, GEOMETRY_ORDER>
          (retval, zfacdatfab, face, graph, distance,
           currentRow, 2, a_dx, string("Dirichlet"));
      
        a_weights(currentRow,currentRow) = a_block.m_zfacesDiriWgt[ivol];
        currentRow++;
      }
      for(unsigned int ivol = 0; ivol < a_block.m_zfacesNeum.size(); ivol++)
      {
        const auto & face = a_block.m_zFacesNeum[ivol];
        RealVect vofloc = getFaceLocation<ZFACE>(face, a_dx, 2);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<ZFACE, OPERATOR_ORDER, GEOMETRY_ORDER>
          (retval, zfacdatfab, face, graph, distance,
           currentRow, 2, a_dx, string("Neumann"));
      
        a_weights(currentRow,currentRow) = a_block.m_zfacesNeumWgt[ivol];
        currentRow++;
      }
    }
#endif
    {
      //Cut faces
      for(unsigned int ivol = 0; ivol < a_block.m_cutfacesDiri.size(); ivol++)
      {
        const auto & face = a_block.m_cutfacesDiri[ivol];
        //-1 gets the cell centered location
        RealVect vofloc = getFaceLocation<BOUNDARY>(face, a_dx, -1);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForEBBC< OPERATOR_ORDER, GEOMETRY_ORDER >
          (retval, a_geoserv, a_srcDomain, face, graph, distance,
           dit[a_ibox], currentRow, 2, a_dx, 
           string("Dirichlet"));
      
        a_weights(currentRow,currentRow) = a_block.m_cutfacesDiriWgt[ivol];
        currentRow++;
      }
      for(unsigned int ivol = 0; ivol < a_block.m_cutfacesNeum.size(); ivol++)
      {
        const auto & face = a_block.m_cutfacesNeum[ivol];
        //-1 gets the cell centered location
        RealVect vofloc = getFaceLocation<BOUNDARY>(face, a_dx, -1);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForEBBC< OPERATOR_ORDER, GEOMETRY_ORDER >
          (retval, a_geoserv, a_srcDomain, face, graph, distance,
           dit[a_ibox], currentRow, 2, a_dx, 
           string("Neumann"));
      
        a_weights(currentRow,currentRow) = a_block.m_cutfacesNeumWgt[ivol];
        currentRow++;
      }
    }

    return retval;
  }  


  //Moore-Penrose inverse of the moment matrix is calulcated using weighted least squares
  inline LAPACKMatrix
  getMoorePenroseInverse(const LAPACKMatrix & a_M)
  {
    LAPACKMatrix MT = a_M;
    MT.transpose();

    LAPACKMatrix retval;
    multiply(retval, MT, a_M);
    int maxiter = 100;
    Real tol = 1.0e-16;
    retval.invertUsingSVD(maxiter, tol);
    return retval;
  }

  ///
  LAPACKMatrix
  getStencilMatrix(const LAPACKMatrix & a_Qmat,
                   const LAPACKMatrix & a_Gmat,
                   const LAPACKMatrix & a_GTGinvMat)
  {
    LAPACKMatrix GTmat = a_Gmat;
    GTmat.transpose();
    //thi9s is (gtg)-1 GT
    LAPACKMatrix gtginvgtMat;
    multiply(gtginvgtMat, a_GTGinvMat, GTmat);
    
    LAPACKMatrix Smat;
    multiply(Smat, a_Qmat, gtginvgtMat);
    return Smat;
  }

   
  ///
  template <CENTERING cent, unsigned int OPERATOR_ORDER, unsigned int GEOMETRY_ORDER>
  inline LAPACKMatrix
  getDharshiQMatrix(const EBIndex<cent>                                 & a_face,
                    const Neighborhood<CELL>                            & a_block,
                    const shared_ptr< GeometryService<GEOMETRY_ORDER> > & a_geoserv,
                    const shared_ptr<LevelData<HostIrregData<cent,
                    IndexedMoments<DIM-1, GEOMETRY_ORDER>, 1> > >       & a_irregdata,
                    const Chombo4::Box                                  & a_srcDomain,
                    const RealVect                                      & a_xbar,
                    unsigned int                                          a_ibox,
                    Real                                                  a_dx,
                    int                                                   a_facedir,
                    const Chombo4::DataIndex                            & a_datind)
  {
    const auto& irregfab = (*a_irregdata)[a_datind];
    const auto & graphsldptr = a_geoserv->getGraphs(    a_srcDomain);
    const auto & graph = (*graphsldptr)[a_datind];
    
    IndexedMoments<DIM-1, GEOMETRY_ORDER>  facmom = hoeb_data::
      getIndexedMoment<DIM-1, cent, GEOMETRY_ORDER>(a_face, irregfab, graph, a_dx);
    IndexedMoments<DIM  , GEOMETRY_ORDER> volmom
      = getVolMomFromFacMom<DIM, GEOMETRY_ORDER>(facmom, a_facedir);

    RealVect facloc = hoeb::getFaceLocation<cent>(a_face, a_dx, a_facedir);
    IndexTM<Real, DIM> distance = getIndexDistance(facloc, a_xbar);
    volmom.shift(distance);

    
    LAPACKMatrix Qmat(1, IndexedMoments<DIM, OPERATOR_ORDER>::size());

    for(MomentIterator<DIM, OPERATOR_ORDER> momit; momit.ok(); ++momit)
    {
      auto momind = momit();
      Real matval = 0;
      if(momind[a_facedir] > 0)
      {
        int pd = momind[a_facedir];
        IndexTM<int, DIM> derind =momind - IndexTM<int, DIM>::BASISV(a_facedir);
        Real derval = volmom[derind];
        matval = pd*derval;
      }
      unsigned int currentCol = IndexedMoments<DIM, OPERATOR_ORDER>::indexOf(momind);
      Qmat(0,currentCol) = matval;
    }
    
    return Qmat;
  }

  ///template specialization to get the eb flux q matrix
  template <unsigned int OPERATOR_ORDER, unsigned int GEOMETRY_ORDER>
  inline LAPACKMatrix
  getDharshiQMatrix(const EBIndex<BOUNDARY>                             & a_face,
                    const Neighborhood<CELL>                            & a_block,
                    const shared_ptr< GeometryService<GEOMETRY_ORDER> > & a_geoserv,
                    const shared_ptr<LevelData<HostIrregData<BOUNDARY,
                    IndexedMoments<DIM-1, GEOMETRY_ORDER>, 1> > >       & a_irregdata,
                    const Chombo4::Box                                  & a_srcDomain,
                    const RealVect                                      & a_xbar,
                    unsigned int                                          a_ibox,
                    Real                                                  a_dx,
                    int                                                   a_facedir,
                    const Chombo4::DataIndex                            & a_datind)
  {
    LAPACKMatrix Qmat(1, IndexedMoments<DIM, OPERATOR_ORDER>::size());
    Qmat.setVal(0.);
    const auto & graphsldptr = a_geoserv->getGraphs(    a_srcDomain);
    const auto & graph = (*graphsldptr)[a_datind];
    for(unsigned int normDir = 0; normDir < DIM; normDir++)
    {
      IndexedMoments<DIM, GEOMETRY_ORDER> normMom=
        getEBNormalMoment(a_face, a_geoserv, a_srcDomain,
                          a_datind, graph, normDir, a_dx);
        
      for(MomentIterator<DIM, OPERATOR_ORDER> momit; momit.ok(); ++momit)
      {
        auto momind = momit();
        if(momind[normDir] > 0)
        {
          int pd = momind[normDir];
          IndexTM<int, DIM> derind =momind - IndexTM<int, DIM>::BASISV(normDir);
          Real derval = normMom[derind];
          unsigned int currentCol = IndexedMoments<DIM, OPERATOR_ORDER>::indexOf(momind);
          Qmat(0,currentCol) += pd*derval;
        }
      }
    }
    return Qmat;
  }
  
                       
  ///
  inline LocalStencil<CELL, Real>
  getStencilFromMatrix(const LAPACKMatrix      & a_stenmat,
                       const Neighborhood<CELL>& a_block)
  {
    LocalStencil<CELL, Real> retval;
    for(unsigned int ivol = 0; ivol < a_block.m_volumes.size(); ivol++)
    {
      const auto& vof = a_block.m_volumes[ivol];
      //this depends upon the volumes being the first bits
      Real wgt = a_stenmat(0, ivol);
      retval.add(vof, wgt);
    }

    return retval;
  }
  
  ///
  /**
     Return integral F dA  (do NOT divide out the area).
     For EB Fluxes, it is Integral F-n dA with an outward-facing normal.
  */
  template <CENTERING cent, unsigned int OPERATOR_ORDER, unsigned int GEOMETRY_ORDER>
  inline LocalStencil<CELL, Real> 
  getDharshiIntFluxDAStencil(const EBIndex<cent>                                 & a_face,
                             const EBIndex<CELL>                                 & a_startingVoF,
                             const std::string                                     a_dombcname[2*DIM],
                             const std::string                                   & a_ebbcname,
                             const shared_ptr< GeometryService<GEOMETRY_ORDER> > & a_geoserv,
                             const shared_ptr<LevelData<HostIrregData<cent,
                             IndexedMoments<DIM-1, GEOMETRY_ORDER>, 1> > >       & a_irregdata,
                             const Chombo4::Box                                  & a_srcDomain,
                             unsigned int a_ibox,                              
                             Real a_alpha, Real a_beta, Real a_dx, int a_facedir,
                             Side::LoHiSide a_vofside)
  {
    //facedir == -1 means the EB
    const auto & graphsldptr = a_geoserv->getGraphs(    a_srcDomain);
    auto dit = a_geoserv->getDBL(a_srcDomain).dataIterator();
    const auto & datind = dit[a_ibox];
    const auto & graph = (*graphsldptr)[datind];
    LocalStencil<CELL, Real>  retval;
    //when I know the answer, short circuit all the other stuff and return it
    if(a_facedir==-1)
    {
      if(a_ebbcname  == string("Neumann"))
      {
        return retval;
      }
    }
    else if ((a_facedir >= 0) && a_face.m_isBoundary)
    {
      int iindex = ebp_index(a_facedir, a_vofside);
      if(a_dombcname[iindex] == string("Neumann"))
      {
        return retval;
      }
    }
    else
    {
      //get the neighborhood near the target face
      hoeb::Neighborhood<CELL>  block = getCellNeighbors<cent, OPERATOR_ORDER>
        (a_face, graph, a_dombcname, a_ebbcname, a_facedir);
      const EBIndex<cent>& face = a_face;
      Real                 dx = a_dx;
      int                  facedir = a_facedir;
      RealVect xbar =   getFaceLocation(face, dx, facedir);

      //this matrix contains all the moments and boundary condition constraints
      LAPACKMatrix weight;
      LAPACKMatrix Mmat =  getMomentMatrix<CELL, OPERATOR_ORDER, GEOMETRY_ORDER>
        (weight, block, a_geoserv, a_srcDomain, xbar, a_ibox, a_dx);

      LAPACKMatrix Gmat; //G = WM
      multiply(Gmat, weight, Mmat);
      
      //Moore-Penrose inverse of the moment matrix
      //is calulcated using weighted least squares
      //because G is weighted version of M
      LAPACKMatrix GTGinvMat = getMoorePenroseInverse(Gmat);

      //this is the matrix that defines the operator
      LAPACKMatrix Qmat =
        getDharshiQMatrix<cent, OPERATOR_ORDER, GEOMETRY_ORDER>
        (a_face, block, a_geoserv, a_irregdata,
         a_srcDomain, xbar, a_ibox, a_dx, a_facedir, datind);

      LAPACKMatrix Smat = getStencilMatrix(Qmat, Gmat, GTGinvMat);
      
      retval = hoeb::getStencilFromMatrix(Smat, block);
    }
    return retval;
  }

  
  /******/  
  template<unsigned int OPERATOR_ORDER, unsigned int GEOMETRY_ORDER>
  inline LocalStencil<CELL, Real> 
  getFullDharshiStencil(const EBIndex<CELL>                                 & a_vof,
                        const std::string                                     a_dombcname[2*DIM],
                        const std::string                                   & a_ebbcname,
                        const shared_ptr< GeometryService<GEOMETRY_ORDER> > & a_geoserv,
                        Proto::Box                                          & a_srcDomain,
                        unsigned int a_ibox, Real a_alpha, Real a_beta, Real a_dx) 
  {
    LocalStencil<CELL, Real> vofsten;
    //use the age-old trick of building up the stencil by construction
    auto dbl = a_geoserv->getDBL(a_srcDomain);
    auto dit = dbl.dataIterator();
    const auto & graphsldptr = a_geoserv->getGraphs(    a_srcDomain);
    const auto & graph = (*graphsldptr)[dit[a_ibox]];
    for(SideIterator sit; sit.ok(); ++sit)
    {
      int isign = sign(sit());
      {                                                // 
        auto xfaces = graph.getXFaces(a_vof, sit());
        for(unsigned int iface = 0; iface < xfaces.size(); iface++)
        {
          auto face = xfaces[iface];
          LocalStencil<CELL, Real> fluxsten =
            getDharshiIntFluxDAStencil<XFACE, OPERATOR_ORDER, GEOMETRY_ORDER>
            (face,
             a_vof,a_dombcname,a_ebbcname,
             a_geoserv, a_geoserv->getXFaceData(a_srcDomain),
             a_srcDomain, a_ibox,
             a_alpha, a_beta, a_dx, 0, sit());
          fluxsten *= Real(isign);
          vofsten += fluxsten;
        }
      }
      {                                                // 
        auto yfaces = graph.getYFaces(a_vof, sit());
        for(unsigned int iface = 0; iface < yfaces.size(); iface++)
        {
          auto face = yfaces[iface];
          LocalStencil<CELL, Real> fluxsten =
            getDharshiIntFluxDAStencil<YFACE, OPERATOR_ORDER, GEOMETRY_ORDER>
            (face,
             a_vof,a_dombcname,a_ebbcname,
             a_geoserv, a_geoserv->getYFaceData(a_srcDomain),
             a_srcDomain, a_ibox,
             a_alpha, a_beta, a_dx, 1, sit());
          fluxsten *= Real(isign);
          vofsten += fluxsten;
        }
      }
#if DIM==3      
      {                                                // 
        auto zfaces = graph.getZFaces(a_vof, sit());
        for(unsigned int iface = 0; iface < zfaces.size(); iface++)
        {
          auto face = zfaces[iface];
          LocalStencil<CELL, Real> fluxsten =
            getDharshiIntFluxDAStencil<ZFACE, OPERATOR_ORDER, GEOMETRY_ORDER>
            (face,
             a_vof, a_dombcname, a_ebbcname,
             a_geoserv, , a_geoserv->getYFaceData(a_srcDomain),
             a_srcDomain, a_ibox,
             a_alpha, a_beta, a_dx, 2, sit());
          fluxsten *= Real(isign);
          vofsten += fluxsten;
        }
      }
#endif
    }
    {
      EBIndex<BOUNDARY> face = a_vof.getCutFace();
      //dummy not used
      shared_ptr<LevelData<HostIrregData<BOUNDARY, IndexedMoments<DIM-1, GEOMETRY_ORDER>, 1> > >  dummy;
      LocalStencil<CELL, Real>
        fluxsten =
        getDharshiIntFluxDAStencil<BOUNDARY, OPERATOR_ORDER, GEOMETRY_ORDER>
        (face,
         a_vof,a_dombcname,a_ebbcname,
         a_geoserv, dummy,
         a_srcDomain, a_ibox,
         a_alpha, a_beta, a_dx, -1, Side::Invalid);
      vofsten += fluxsten;
    }
    //need to divide by dx^d to get kappa*lapl(phi)
    Real cellvolinv = 1;
    for(int idir = 0; idir < DIM; idir++)
    {
      cellvolinv /= a_dx;
    }
    vofsten *= cellvolinv;

    // set diagonal weight of the stencil so we have the choice of using multigrid
    Real diagweight;
    bool foundVoF = vofsten.getValueAtVoF(diagweight, a_vof);
    PR_assert(foundVoF); //if this fails, the stencil has no diagonal value
    vofsten.setBetaWeight(diagweight);
    
    return vofsten;
  }

  /******/
  template<unsigned int OPERATOR_ORDER, unsigned int GEOMETRY_ORDER>
  inline void
  dharshiLaplStencil(string                                              & a_stencilName,
                     string                                              & a_ebbcName,
                     vector<EBIndex<CELL> >                              & a_dstVoFs,
                     vector<LocalStencil<CELL, Real> >                   & a_stencil,
                     Proto::Box                                          & a_srcValid,
                     Proto::Box                                          & a_dstValid,
                     Proto::Box                                          & a_srcDomain,
                     Proto::Box                                          & a_dstDomain,
                     Proto::Point                                        & a_srcGhost,
                     Proto::Point                                        & a_dstGhost,
                     bool                                                & a_needDiagonalWeights,
                     const shared_ptr< GeometryService<GEOMETRY_ORDER> > & a_geoserv,
                     const Chombo4::DisjointBoxLayout                    & a_grids,
                     const Chombo4::Box                                  & a_domain,
                     const Real                                          & a_dx,
                     unsigned int                                          a_ibox)

  {
    /* just in case a fit of madness causes me to use geometric multigrid */
    a_needDiagonalWeights = true;
    a_srcDomain = ProtoCh::getProtoBox(a_domain);
    a_dstDomain = ProtoCh::getProtoBox(a_domain);
    Chombo4::ParmParse pp;
    
    Real alpha = 1.0;
    Real beta = -0.001;
    string dombcname;
    int nghost;
    pp.get("num_ghost_cells", nghost);
    pp.get("domainBC"  , dombcname);
    pp.get("EBBC"      , a_ebbcName);
    pp.get("alpha"     , alpha);
    pp.get("beta"      , beta);
    static bool printed = false;
    if(!printed)
    {
      printed = true;
      pout() << "domainBC"  << " = " <<  dombcname      << endl;
      pout() << "EBBC"      << " = " <<  a_ebbcName     << endl;
    }
    
    a_stencilName = string("Dharshi_Laplacian");

    const auto & graphsldptr = a_geoserv->getGraphs(    a_srcDomain);
    const auto & voldatldptr = a_geoserv->getVoluData(  a_srcDomain);
    const auto & ebfdatldptr = a_geoserv->getEBFaceData(a_srcDomain);
    const auto & xfadatldptr = a_geoserv->getXFaceData( a_srcDomain);
    const auto & yfadatldptr = a_geoserv->getYFaceData( a_srcDomain);
    const auto & zfadatldptr = a_geoserv->getZFaceData( a_srcDomain);
    const auto & dbl         = a_geoserv->getDBL(a_srcDomain);
    auto dit = dbl.dataIterator();

    const auto & graph = (*graphsldptr)[dit[a_ibox]];
    a_srcValid  = graph.validBox();
    a_dstValid  = graph.validBox();
    a_srcGhost  = Point::Ones(nghost);
    a_dstGhost  = Point::Ones(nghost);
    a_needDiagonalWeights = true;
    
    string dombcarray[2*DIM];
    for(int ivec = 0; ivec < 2*DIM; ivec++)
    {
      dombcarray[ivec] = dombcname;
    }
    for(auto bit = a_dstValid.begin(); bit != a_dstValid.end(); ++bit)
    {
      auto vofs = graph.getVoFs(*bit);
      for(unsigned int ivof = 0; ivof < vofs.size(); ivof++)
      {
        LocalStencil<CELL, Real> vofsten =
          getFullDharshiStencil<OPERATOR_ORDER, GEOMETRY_ORDER>(vofs[ivof],
                                                                dombcarray, a_ebbcName,
                                                                a_geoserv, a_srcDomain, a_ibox,
                                                                alpha, beta, a_dx);
        a_dstVoFs.push_back(vofs[ivof]);
        a_stencil.push_back(vofsten);
      }
    }
  }

  inline shared_ptr<BaseIF> getImplicitFunction()
  {
    ParmParse pp;
    string whichgeom;
    pp.get("which_geom", whichgeom);
    
    shared_ptr<BaseIF> retval;
    if(whichgeom == string("sphere"))
    {
      RealVect center;
      Real radius = 0.1;
      bool inside = false;
      pp.get("radius", radius);
      pp.get("inside", inside);
      pp.getRV("center", center);

      pout() << "Implicit function is a sphere with radius = " << radius
             << ", and center = " << center << endl;
      SimpleSphereIF* sphereptr = new SimpleSphereIF(center, radius, inside);
      retval = shared_ptr<BaseIF>(static_cast<BaseIF*>(sphereptr));
    }
    else if(whichgeom == string("plane"))
    {
      RealVect startpt, normal;
      pp.getRV("plane_startpt", startpt);
      pp.getRV("plane_normal", normal);
    }
    return retval;
  }
  ///
  /**
     Fill in scalar from exactsolutions class
  */
  template< unsigned int OPERATOR_ORDER, unsigned int GEOMETRY_ORDER>
  inline void
  fillPhi(EBLevelBoxData<CELL, 1>                                &  a_phi,
          const shared_ptr<LevelData<EBGraph> >                  &  a_graphs,
          const Chombo4::DisjointBoxLayout                       &  a_grids,
          const Chombo4::Box                                     &  a_domain,
          const Real                                             &  a_dx,
          const shared_ptr< GeometryService<GEOMETRY_ORDER> >    &  a_geoserv)
  {
    typedef IndexedMoments<DIM  , GEOMETRY_ORDER> IndMomDIM;
    typedef HostIrregData<CELL,      IndMomDIM , 1>  VoluData;
    typedef GraphConstructorFactory<EBHostData<CELL, Real, 1> > hostfactory_t;
    
    LevelData<EBHostData<CELL, Real, 1> > hostdata
      (a_grids, 1,  a_phi.ghostVect(), hostfactory_t(a_graphs));
    
    RealVect center = 0.5*RealVect::Unit();
    Real radius = 0.1;
    bool inside = false;
    ParmParse pp;
    string which_phi;
    pp.get("which_phi", which_phi);
    shared_ptr<BaseExactSolution<OPERATOR_ORDER, GEOMETRY_ORDER> > phigenptr;
    if(which_phi  == string("SineSphere"))
    {
      pp.get("radius", radius);
      pp.get("inside", inside);
      pp.getRV("center", center);

      phigenptr =shared_ptr<BaseExactSolution<OPERATOR_ORDER, GEOMETRY_ORDER> >
        (new  hoeb::SineSphereEF<OPERATOR_ORDER, GEOMETRY_ORDER>(radius, center));
    }
    else if(which_phi  == string("Polynomial"))
    {
      int num_terms;
      pp.get("polynomial_num_terms", num_terms);

      std::vector<std::pair<Point, Real> > entries(num_terms);
      for(int iterm = 0; iterm < num_terms; iterm++)
      {
        Real coeff;
        std::vector<int>  powervec(DIM);
        string powerstr = string("polynomial_powers_")      + std::to_string(iterm);
        string coeffstr = string("polynomial_coefficient_") + std::to_string(iterm);
        pp.getarr(powerstr.c_str(), powervec, 0, DIM);
        pp.get   (coeffstr.c_str(), coeff);
        Point power(powervec);
        entries[iterm].first  = power;
        entries[iterm].second = coeff;
      }
      phigenptr = shared_ptr<BaseExactSolution<OPERATOR_ORDER, GEOMETRY_ORDER> >
        (new hoeb::PolynomialEF<OPERATOR_ORDER, GEOMETRY_ORDER>(entries));
    }

    auto& phigen = *phigenptr;
    auto voldatldptr = a_geoserv->getVoluData(a_domain);
    const auto& voldatld = *voldatldptr;
    Chombo4::DataIterator dit = a_grids.dataIterator();
    for(unsigned int ibox = 0; ibox < dit.size(); ibox++)
    {
      const auto& graph  = (*a_graphs)[dit[ibox]];
      const auto& voldat = voldatld[dit[ibox]];
      Bx validbx = graph.validBox();
      for(auto bit = validbx.begin(); bit != validbx.end(); ++bit)
      {
        auto vofs= graph.getVoFs(*bit);
        for(int ivof = 0; ivof < vofs.size(); ++ivof)
        {
          auto& vof = vofs[ivof];
          Real phival = phigen(graph, a_dx, voldat, vof);
          hostdata[dit[ibox]](vof, 0) = phival;
        }
      }
    }

    EBLevelBoxData<CELL, 1>::copyToDevice(a_phi, hostdata);
  }

  ///
  /**
   **/
  template<unsigned int GEOMETRY_ORDER>
  inline void
  restrictKappaLphi(EBLevelBoxData<CELL, 1>                                           &  a_klpFToC,
                    const EBLevelBoxData<CELL, 1>                                     &  a_klpFine,
                    const shared_ptr<LevelData<EBGraph> >                             &  a_graphsFine,
                    const Chombo4::DisjointBoxLayout                                  &  a_gridsFine,
                    const Chombo4::Box                                                &  a_domFine,
                    const Real                                                        &  a_dxFine,
                    const shared_ptr<LevelData<EBGraph> >                             &  a_graphsCoar,
                    const Chombo4::DisjointBoxLayout                                  &  a_gridsCoar,
                    const Chombo4::Box                                                &  a_domCoar,
                    const Real                                                        &  a_dxCoar,
                    const shared_ptr<EBDictionary<GEOMETRY_ORDER, Real, CELL, CELL> > &  a_dictionary,
                    const shared_ptr< GeometryService<GEOMETRY_ORDER> >               &  a_geoserv)
  {
    std::string  nobcname        = string("no_bcs");
    std::string  restrictionName = string("Multigrid_Restriction");
    string dombc[2*DIM];
    for(unsigned int idom = 0; idom < 2*DIM;  idom++)
    {
      dombc[idom] = nobcname;
    }
    a_dictionary->registerStencil(restrictionName, dombc, nobcname, a_domFine, a_domCoar, false);
    Chombo4::DataIterator dit = a_gridsCoar.dataIterator();
    for(unsigned int ibox = 0; ibox < dit.size(); ++ibox)
    {
      auto      & coarfab = a_klpFToC[dit[ibox]];
      const auto& finefab = a_klpFine[dit[ibox]];
      auto stencil = a_dictionary->getEBStencil(restrictionName, nobcname, a_domFine, a_domCoar, ibox);
      //set resc = Ave(resf) (true is initToZero)
      stencil->apply(coarfab, finefab,  true, 1.0);
    }
  }

  ///
  template< unsigned int GEOMETRY_ORDER >
  inline void  
  fillKappa(  EBLevelBoxData<CELL, 1>                              & a_kappa,
              const Chombo4::DisjointBoxLayout                     & a_grids,
              const Chombo4::Box                                   & a_domain,
              const shared_ptr<LevelData<EBGraph> >                & a_graphs,
              const shared_ptr<GeometryService<GEOMETRY_ORDER> >   & a_geoserv)
  {
    Chombo4::DataIterator dit = a_grids.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ++ibox)
    {
      auto& kappdat = a_kappa[dit[ibox]];
      auto grid =a_grids[dit[ibox]];
      Bx  grbx = kappdat.inputBox();
      const EBGraph  & graph = (*a_graphs)[dit[ibox]];
      EBHostData<CELL, Real, 1> hostdat(grbx, graph);
      //fill kappa on the host then copy to the device
      a_geoserv->fillKappa(hostdat, grid, dit[ibox], a_domain);
      // now copy to the device
      EBLevelBoxData<CELL, 1>::copyToDevice(kappdat, hostdat);
    }
  }
  ///
  /**
   **/
  template<unsigned int GEOMETRY_ORDER>
  inline void
  restrictPhi(EBLevelBoxData<CELL, 1>                                           &  a_klpFToC,
              const EBLevelBoxData<CELL, 1>                                     &  a_klpFine,
              const shared_ptr<LevelData<EBGraph> >                             &  a_graphsFine,
              const Chombo4::DisjointBoxLayout                                  &  a_gridsFine,
              const Chombo4::Box                                                &  a_domFine,
              const Real                                                        &  a_dxFine,
              const shared_ptr<LevelData<EBGraph> >                             &  a_graphsCoar,
              const Chombo4::DisjointBoxLayout                                  &  a_gridsCoar,
              const Chombo4::Box                                                &  a_domCoar,
              const Real                                                        &  a_dxCoar,
              const shared_ptr<EBDictionary<GEOMETRY_ORDER, Real, CELL, CELL> > &  a_dictionary,
              const shared_ptr< GeometryService<GEOMETRY_ORDER> >               &  a_geoserv)
  {
    std::string  nobcname        = string("no_bcs");
    std::string  restrictionName = string("Volume_Weighted_Restriction");
    string dombc[2*DIM];
    for(unsigned int idom = 0; idom < 2*DIM;  idom++)
    {
      dombc[idom] = nobcname;
    }
    a_dictionary->registerStencil(restrictionName, dombc, nobcname, a_domFine, a_domCoar, false);
    Chombo4::DataIterator dit = a_gridsCoar.dataIterator();
    for(unsigned int ibox = 0; ibox < dit.size(); ++ibox)
    {
      auto      & coarfab = a_klpFToC[dit[ibox]];
      const auto& finefab = a_klpFine[dit[ibox]];
      auto stencil = a_dictionary->getEBStencil(restrictionName, nobcname, a_domFine, a_domCoar, ibox);
      //set resc = Ave(resf) (true is initToZero)
      stencil->apply(coarfab, finefab,  true, 1.0);
    }
  }
}
///

#endif
