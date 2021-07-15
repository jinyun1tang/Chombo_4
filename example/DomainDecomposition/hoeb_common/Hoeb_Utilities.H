#ifndef __hoeb_utilities__
#define __hoeb_utilities__
#include "Proto.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_GeometryService.H"
#include "Chombo_EBDictionary.H"
#include "Hoeb_LAPACKMatrix.H"
#include "Hoeb_Neighborhood.H"


namespace hoeb
{
  /// extract data if it is there, use regular if not
  template <unsigned int D, CENTERING cent>
  inline IndexedMoments<D, HOEB_MAX_ORDER>   
  getIndexedMoment(const EBIndex<cent>                                              &   a_index,
                   const HostIrregData<cent, IndexedMoments<D, HOEB_MAX_ORDER>, 1 > &   a_irregfab,
                   const EBGraph                                                    &   a_graph,
                   const Real                                                       &   a_dx);
  
  /// get the dreaded M from the neighborhood
  template <CENTERING cent>
  LAPACKMatrix
  getMomentMatrix(const Neighborhood<cent>                            & a_block,
                  const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                  const Chombo4::Box                                  & a_srcDomain,
                  const RealVect                                      & a_xbar,
                  unsigned int                                          a_ibox,
                  Real                                                  a_dx);

  //Moore-Penrose inverse of the moment matrix is calulcated using weighted least squares
  inline LAPACKMatrix
  getMoorePenroseInverse(const LAPACKMatrix & a_M);
  ///do matrix algebra to get stencil from operator
  inline LAPACKMatrix
  getStencilMatrix(const LAPACKMatrix & a_Qmat,
                   const LAPACKMatrix & a_Gmat,
                   const LAPACKMatrix & a_GTGinvMat);
  
  ///template specialization to get the eb flux q matrix
  template <unsigned int OPERATOR_ORDER, unsigned int GEOMETRY_ORDER>
  inline LAPACKMatrix
  getDharshiQMatrix(const EBIndex<BOUNDARY>                             & a_face,
                    const Neighborhood<CELL>                            & a_block,
                    const shared_ptr< GeometryService<GEOMETRY_ORDER> > & a_geoserv,
                    const Chombo4::Box                                  & a_srcDomain,
                    const RealVect                                      & a_xbar,
                    unsigned int                                          a_ibox,
                    Real                                                  a_dx,
                    int                                                   a_facedir,
                    const Chombo4::DataIndex                            & a_datind);
  
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
                             const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                             const Chombo4::Box                                  & a_srcDomain,
                             unsigned int a_ibox,                              
                             Real a_alpha, Real a_beta, Real a_dx, int a_facedir,
                             Side::LoHiSide a_vofside);

  ///get location of face center
  template <CENTERING cent>
  inline RealVect
  getFaceLocation(const EBIndex<cent>& a_face,
                  Real                 a_dx,
                  int                  a_facedir);
  

  inline 
  LocalStencil<CELL, Real>
  getStencilFromMatrix(const LAPACKMatrix      & a_stenmat,
                       const Neighborhood<CELL>& a_block);
  
  ///dharshi stencil at  single vof
  template<unsigned int OPERATOR_ORDER, unsigned int GEOMETRY_ORDER>
  inline LocalStencil<CELL, Real> 
  getFullDharshiStencil(const EBIndex<CELL>                                 & a_vof,
                        const std::string                                     a_dombcname[2*DIM],
                        const std::string                                   & a_ebbcname,
                        const shared_ptr< GeometryService<GEOMETRY_ORDER> > & a_geoserv,
                        Proto::Box                                          & a_srcDomain,
                        unsigned int a_ibox, Real a_alpha, Real a_beta, Real a_dx) ;
  

  ///get the stencil for the laplacian described in Devendran, et al. for whole box
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
                     unsigned int                                          a_ibox);

  ///return an implicit function based on stuff in the input file
  inline shared_ptr<BaseIF> getImplicitFunction();

  ///
  /**
     Fill in scalar from exactsolutions class
  */
  template< unsigned int GEOMETRY_ORDER>
  inline void
  fillPhi(EBLevelBoxData<CELL, 1>                                &  a_phi,
          const shared_ptr<LevelData<EBGraph> >                  &  a_graphs,
          const Chombo4::DisjointBoxLayout                       &  a_grids,
          const Chombo4::Box                                     &  a_domFine,
          const Real                                             &  a_dx,
          const shared_ptr< GeometryService<HOEB_MAX_ORDER> >    &  a_geoserv);

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
                    const shared_ptr< GeometryService<GEOMETRY_ORDER> >               &  a_geoserv);

  ///
  template <unsigned int D>
  inline 
  IndexTM<int, D-1>
  getFaceIndexFromVol(const  IndexTM<int, D>& a_volind, int a_facedir);

  ///
  template <unsigned int D, unsigned int GEOMETRY_ORDER>
  inline 
  IndexedMoments<D, GEOMETRY_ORDER>
  getVolMomFromFacMom(const IndexedMoments<D-1, GEOMETRY_ORDER>& a_facmom,
                      int a_facedir);
}
#include "Hoeb_UtilitiesImplem.H"
#endif
