#include <cmath>
#include <cstdio>
#include <iostream>


#include "EBProto.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_LevelData.H"
#include "Chombo_BaseFab.H"

#include "Chombo_ParmParse.H"
#include "Chombo_LoadBalance.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_BRMeshRefine.H"
#include "Chombo_GeometryService.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "EBMultigrid.H"
#include "Proto_DebugHooks.H"
#include "DebugFunctions.H"
#include "Hoeb_ExactSolutions.H"
//this one defines HOEB_MAX_ORDER
#include "Hoeb_Utilities.H"
#include <iomanip>

//geometry_order has to be >= operator_order
#define GEOMETRY_ORDER 6
#define OPERATOR_ORDER 4

/****/
//phi not const because of exchange
void
getKappaLphi(EBLevelBoxData<CELL, 1>                                          & a_klp,
             EBLevelBoxData<CELL, 1>                                          & a_phi,
             const shared_ptr<LevelData<EBGraph> >                            & a_graphs,
             const Chombo4::DisjointBoxLayout                                 & a_grids,
             const Chombo4::Box                                               & a_domain,
             const Real                                                       & a_dx,
             const shared_ptr<EBDictionary<GEOMETRY_ORDER, Real, CELL, CELL> >& a_dictionary,
             const shared_ptr< GeometryService<GEOMETRY_ORDER> >              & a_geoserv)
{
  PROTO_ASSERT((GEOMETRY_ORDER >= OPERATOR_ORDER), "Need at least operator accuracy in geometric moments");
  string stencilName;
  string ebbcName;
  Chombo4::Copier exchangeCopier;
  exchangeCopier.exchangeDefine(a_grids, a_phi.ghostVect());
  a_phi.exchange(exchangeCopier);

  //register it for every box
  Chombo4::DataIterator dit = a_grids.dataIterator();
  for(unsigned int ibox = 0; ibox < dit.size(); ++ibox)
  {
    vector<     EBIndex<CELL>  >          dstVoFs;
    vector<LocalStencil<CELL, Real> >     stencil;
    Proto::Box                            srcValid;
    Proto::Box                            dstValid;
    Proto::Box                            srcDomain;
    Proto::Box                            dstDomain;
    Point                                 srcGhost;
    Point                                 dstGhost;
    bool                                  needDiagonalWeights;
    
    hoeb::
      dharshiLaplStencil<OPERATOR_ORDER, GEOMETRY_ORDER>
      (stencilName,        
       ebbcName,           
       dstVoFs,            
       stencil,            
       srcValid,           
       dstValid,           
       srcDomain,          
       dstDomain,          
       srcGhost,           
       dstGhost,
       needDiagonalWeights,
       a_geoserv,
       a_grids,
       a_domain,
       a_dx,
       ibox);
  
    ///registering stencil
    a_dictionary->registerStencil(stencilName,        
                                  ebbcName,           
                                  dstVoFs,            
                                  stencil,            
                                  srcValid,           
                                  dstValid,           
                                  srcDomain,          
                                  dstDomain,          
                                  srcGhost,           
                                  dstGhost,
                                  needDiagonalWeights,
                                  ibox);
  }
  //now apply it
  for(unsigned int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto      & lphfab = a_klp[dit[ibox]];
    const auto& phifab = a_phi[dit[ibox]];
    auto stencil = a_dictionary->getEBStencil(stencilName, ebbcName, a_domain, a_domain, ibox);
    //set resc = Ave(resf) (true is initToZero)
    stencil->apply(lphfab, phifab,  true, 1.0);
  }
}
/*******/ 
PROTO_KERNEL_START 
void subtractionTractionF(Var<Real, 1>    a_error,
                          Var<Real, 1>    a_klpFtoC,
                          Var<Real, 1>    a_klpCoar)
{
  a_error(0) = a_klpFtoC(0) - a_klpCoar(0);
}
PROTO_KERNEL_END(subtractionTractionF, subtractionTraction)
/****/
void
getKLPhiError(EBLevelBoxData<CELL,   1>                                           &  a_errCoar, 
              const shared_ptr<LevelData<EBGraph> >                               &  a_graphsFine,
              const Chombo4::DisjointBoxLayout                                    &  a_gridsFine,
              const Chombo4::Box                                                  &  a_domFine,
              const Real                                                          &  a_dxFine,
              const shared_ptr<LevelData<EBGraph> >                               &  a_graphsCoar,
              const Chombo4::DisjointBoxLayout                                    &  a_gridsCoar,
              const Chombo4::Box                                                  &  a_domCoar,
              const Real                                                          &  a_dxCoar,
              const shared_ptr<EBDictionary<GEOMETRY_ORDER, Real, CELL, CELL> >   &  a_dictionary,
              const shared_ptr< GeometryService<GEOMETRY_ORDER> >                 &  a_geoserv,
              int a_nxlev)
{
  ParmParse pp;
  int nghost;
  pp.get("num_ghost_cells", nghost);
  IntVect dataGhostIV =   nghost*IntVect::Unit;
  EBLevelBoxData<CELL,   1>  phiFine(a_gridsFine, dataGhostIV, a_graphsFine);
  EBLevelBoxData<CELL,   1>  phiCoar(a_gridsCoar, dataGhostIV, a_graphsCoar);
  EBLevelBoxData<CELL,   1>  klpFine(a_gridsFine, dataGhostIV, a_graphsFine);
  EBLevelBoxData<CELL,   1>  klpCoar(a_gridsCoar, dataGhostIV, a_graphsCoar);
  EBLevelBoxData<CELL,   1>  klpFtoC(a_gridsCoar, dataGhostIV, a_graphsCoar);
  
  hoeb::fillPhi<OPERATOR_ORDER,GEOMETRY_ORDER>
    (phiFine, a_graphsFine, a_gridsFine, a_domFine, a_dxFine, a_geoserv);
  hoeb::fillPhi<OPERATOR_ORDER,GEOMETRY_ORDER>
    (phiCoar, a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar, a_geoserv);

  getKappaLphi(klpFine, phiFine, a_graphsFine, a_gridsFine, a_domFine, a_dxFine, a_dictionary, a_geoserv);
  getKappaLphi(klpCoar, phiCoar, a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar, a_dictionary, a_geoserv);

  hoeb::restrictKappaLphi<GEOMETRY_ORDER>
    (klpFtoC, klpFine,
     a_graphsFine, a_gridsFine, a_domFine, a_dxFine,                    
     a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar,
     a_dictionary, a_geoserv);

  //error = Ave(klphifine) - klphicoar 
  Chombo4::DataIterator dit = a_gridsCoar.dataIterator();
  for(unsigned int ibox = 0; ibox < dit.size(); ibox++)
  {
    auto& ftocfab =   klpFtoC[dit[ibox]];
    auto& coarfab =   klpCoar[dit[ibox]];
    auto& errfab  = a_errCoar[dit[ibox]];
    auto  inputbx = ftocfab.inputBox();
    auto  validbx = (*a_graphsCoar)[dit[ibox]].validBox();
    ebforall(inputbx, subtractionTraction, validbx, errfab, ftocfab, coarfab);
  }

  string phiFineName = string("phiFine_") + to_string(a_nxlev) + string(".hdf5");
  string phiCoarName = string("phiCoar_") + to_string(a_nxlev) + string(".hdf5");
  string klpFineName = string("klpFine_") + to_string(a_nxlev) + string(".hdf5");
  string klpCoarName = string("klpCoar_") + to_string(a_nxlev) + string(".hdf5");
  string klpFtoCName = string("klpFtoC_") + to_string(a_nxlev) + string(".hdf5");
  string errCoarName = string("errCoar_") + to_string(a_nxlev) + string(".hdf5");
  Real coveredVal = 0;
  phiFine.writeToFileHDF5(  phiFineName, coveredVal);
  phiCoar.writeToFileHDF5(  phiCoarName, coveredVal);
  klpFine.writeToFileHDF5(  klpFineName, coveredVal);
  klpCoar.writeToFileHDF5(  klpCoarName, coveredVal);
  klpFtoC.writeToFileHDF5(  klpFtoCName, coveredVal);
  a_errCoar.writeToFileHDF5(errCoarName, coveredVal);
}


void getTestMetaData(Real                                          &  coveredval,
                     int                                           &  nx,
                     int                                           &  maxGrid,
                     int                                           &  nghost,
                     Chombo4::ProblemDomain                        &  domain,
                     Chombo4::ProblemDomain                        &  domFine,
                     Chombo4::ProblemDomain                        &  domMedi,
                     Chombo4::ProblemDomain                        &  domCoar,
                     Chombo4::DisjointBoxLayout                    & gridsFine,
                     Chombo4::DisjointBoxLayout                    & gridsMedi,
                     Chombo4::DisjointBoxLayout                    & gridsCoar,
                     shared_ptr<LevelData<EBGraph> >               & graphsFine,
                     shared_ptr<LevelData<EBGraph> >               & graphsMedi,
                     shared_ptr<LevelData<EBGraph> >               & graphsCoar,
                     shared_ptr< GeometryService<GEOMETRY_ORDER> > & geomptr,
                     Vector<Chombo4::DisjointBoxLayout>            & vecgrids,
                     vector<Chombo4::Box>                          & vecdomain,
                     vector<Real>                                  & vecdx,
                     Real & dxFine, Real & dxMedi, Real & dxCoar)
{
  ParmParse pp;

  pp.get("nx"        , nx);
  pp.get("num_ghost_cells"        , nghost);
  pp.get("maxGrid"  , maxGrid);
  pp.get("coveredval", coveredval);         


  static bool printedStuff = false;
  if(!printedStuff)
  {
    printedStuff = true;
    pout() << "nx"        << " = " <<  nx         << endl;
    pout() << "max_grid"  << " = " <<  maxGrid    << endl;
    pout() << "coveredval"<< " = " <<  coveredval << endl;
  }
  

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

  domain = Chombo4::ProblemDomain(domLo, domHi);

  GeometryService<GEOMETRY_ORDER>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

  gridsFine = vecgrids[0];
  gridsMedi = vecgrids[1];
  gridsCoar = vecgrids[2];
  int geomGhost = 6;
  RealVect origin = RealVect::Zero();
  
  shared_ptr<BaseIF>    impfunc = hoeb::getImplicitFunction();
  pout() << "defining geometry" << endl;
  Real dx = 1.0/(Real(nx));
  vecdomain.resize(vecgrids.size(), domain.domainBox());
  vecdx    .resize(vecgrids.size(), dx);
  for(int ilev = 1; ilev < vecgrids.size(); ilev++)
  {
    vecdomain[ilev] = coarsen(vecdomain[ilev-1], 2);
    vecdx    [ilev] =           2*vecdx[ilev-1];
  }
  
  dxFine = vecdx[0];
  dxMedi = vecdx[1];
  dxCoar = vecdx[2];
  
  geomptr = shared_ptr< GeometryService<GEOMETRY_ORDER> >
    (new GeometryService<GEOMETRY_ORDER>
     (impfunc, origin, dxFine, domain.domainBox(), vecgrids, geomGhost));
  
  domFine = vecdomain[0];
  domMedi = vecdomain[1];
  domCoar = vecdomain[2];

  graphsFine = geomptr->getGraphs(domFine.domainBox());
  graphsMedi = geomptr->getGraphs(domMedi.domainBox());
  graphsCoar = geomptr->getGraphs(domCoar.domainBox());
}
/****/
int
runDharhsiTruncationErrorTest()
{
  Real                                            coveredval;
  int                                             nx;
  int                                             maxGrid;
  int                                             nghost;
  Chombo4::ProblemDomain                          domain;
  Chombo4::ProblemDomain                          domFine;
  Chombo4::ProblemDomain                          domMedi;
  Chombo4::ProblemDomain                          domCoar;
  Chombo4::DisjointBoxLayout                      gridsFine;
  Chombo4::DisjointBoxLayout                      gridsMedi;
  Chombo4::DisjointBoxLayout                      gridsCoar;
  shared_ptr<LevelData<EBGraph> >                 graphsFine;
  shared_ptr<LevelData<EBGraph> >                 graphsMedi;
  shared_ptr<LevelData<EBGraph> >                 graphsCoar;
  shared_ptr< GeometryService<GEOMETRY_ORDER> >   geomptr;
  Real  dxFine; Real  dxMedi; Real  dxCoar;
  Vector<Chombo4::DisjointBoxLayout> vecgrids;
  Vector<Chombo4::Box>               vecdomain;
  vector<Real>                       vecdx;
  getTestMetaData( coveredval, nx, maxGrid, nghost, domain,
                   domFine  ,    domMedi,    domCoar,
                   gridsFine,  gridsMedi,  gridsCoar,
                   graphsFine,graphsMedi, graphsCoar,
                   geomptr, vecgrids, vecdomain,vecdx,
                   dxFine, dxMedi, dxCoar);
    
  shared_ptr<EBDictionary<GEOMETRY_ORDER, Real, CELL, CELL> >  dictionary
    (new     EBDictionary<GEOMETRY_ORDER, Real, CELL, CELL>
     (geomptr, vecgrids, vecdomain, vecdx, nghost*Point::Ones(1)));

  EBLevelBoxData<CELL,   1>  errMedi(gridsMedi, nghost*IntVect::Unit, graphsMedi);
  EBLevelBoxData<CELL,   1>  errCoar(gridsCoar, nghost*IntVect::Unit, graphsCoar);
  

  getKLPhiError(errMedi, 
                graphsFine, gridsFine, domFine.domainBox(), dxFine,
                graphsMedi, gridsMedi, domMedi.domainBox(), dxMedi,
                dictionary, geomptr, nx);

  getKLPhiError(errCoar, 
                graphsMedi, gridsMedi, domMedi.domainBox(), dxMedi,
                graphsCoar, gridsCoar, domCoar.domainBox(), dxCoar,
                dictionary, geomptr, nx/2);


  //Norm!
  Real normMedi = errMedi.maxNorm(0);
  Real normCoar = errCoar.maxNorm(0);

  Real tol = 1.0e-16;
  Real order = 0;
  if((normCoar > tol) && (normMedi > tol))
  {
    order = log(normCoar/normMedi)/log(2.0);
  }
  pout() << "||klphi errMedi||_max = " << normMedi << std::endl;
  pout() << "||klphi errCoar||_max = " << normCoar << std::endl;
  pout() << "Richardson truncation error order for kappa(L(phi))= " << order << std::endl;
  return 0;
}
/****/
int
runInitialPhiConvergenceTest()
{
  Real                                            coveredval;
  int                                             nx;
  int                                             maxGrid;
  int                                             nghost;
  Chombo4::ProblemDomain                          domain;
  Chombo4::ProblemDomain                          domFine;
  Chombo4::ProblemDomain                          domMedi;
  Chombo4::ProblemDomain                          domCoar;
  Chombo4::DisjointBoxLayout                      gridsFine;
  Chombo4::DisjointBoxLayout                      gridsMedi;
  Chombo4::DisjointBoxLayout                      gridsCoar;
  shared_ptr<LevelData<EBGraph> >                 graphsFine;
  shared_ptr<LevelData<EBGraph> >                 graphsMedi;
  shared_ptr<LevelData<EBGraph> >                 graphsCoar;
  shared_ptr< GeometryService<GEOMETRY_ORDER> >   geomptr;
  Real  dxFine; Real  dxMedi; Real  dxCoar;
  Vector<Chombo4::DisjointBoxLayout> vecgrids;
  Vector<Chombo4::Box>               vecdomain;
  vector<Real>                       vecdx;

  getTestMetaData( coveredval, nx, maxGrid, nghost, domain,
                   domFine  ,    domMedi,    domCoar,
                   gridsFine,  gridsMedi,  gridsCoar,
                   graphsFine,graphsMedi, graphsCoar,
                   geomptr, vecgrids, vecdomain,vecdx,
                   dxFine, dxMedi, dxCoar);

  return 0;
}

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI  
  MPI_Init(&a_argc, &a_argv);
  pout() << "MPI INIT called" << std::endl;
#endif

  //needs to be called after MPI_Init
  CH_TIMER_SETFILE("trunc.time.table");
  {
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
    char* in_file = a_argv[1];
    ParmParse  pp(a_argc-2,a_argv+2,NULL,in_file);
    runInitialPhiConvergenceTest();
    runDharhsiTruncationErrorTest();
  }

  pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI  
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif  
  return 0;
}
