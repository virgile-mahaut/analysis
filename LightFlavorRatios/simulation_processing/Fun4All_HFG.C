#ifndef MACRO_FUN4ALLG4SPHENIX_C
#define MACRO_FUN4ALLG4SPHENIX_C

#include "G4_Input.C"
#include <G4_Global.C>
#include <G4Setup_sPHENIX.C>

#include <Trkr_RecoInit.C>
#include <Trkr_Clustering.C>
#include <Trkr_TruthTables.C>
#include <Trkr_Reco.C>

#include <phpythia8/PHPy8ParticleTrigger.h>

#include <decayfinder/DecayFinder.h>
#include <hftrackefficiency/HFTrackEfficiency.h>
#include <kfparticle_sphenix/KFParticle_sPHENIX.h>

#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>
#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <simqa_modules/QAG4SimulationTracking.h>
#include <qautils/QAHistManagerDef.h>

#include "HF_selections.C"

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libdecayfinder.so)
R__LOAD_LIBRARY(libhftrackefficiency.so)
R__LOAD_LIBRARY(libsimqa_modules.so)

int Fun4All_HFG(std::string processID = "0", std::string channel = "Kshort2pipi")
{
  int nEvents = 0;
  std::string outDir = "./" + channel + "_20260324/";

  string makeDirectory = "mkdir -p " + outDir + "hfEff";
  system(makeDirectory.c_str());

  //F4A setup
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);

  Fun4AllInputManager *infile = new Fun4AllDstInputManager("DSTin");
  string fileName = "lists/g4hits_" + processID + ".list";
  infile->AddListFile(fileName.c_str());
  se->registerInputManager(infile);

  PHRandomSeed::Verbosity(1);
  recoConsts *rc = recoConsts::instance();

  int particleID = 421;
  if (channel == "Kshort2pipi")
  {
    run_pipi_reco = true;
    particleID = 310;
  }
  else if (channel == "Lambda2ppi")
  {
    run_ppi_reco = true;
    particleID = 3122;
  }
  else if (channel == "minBias")
  {
    std::cout << "Min bias simulations" << std::endl; 
  }
  else
  {
    std::cout << "Your decay channel " << channel << " is not known" << std::endl;
    exit(1); 
  }
  //Input::BEAM_CONFIGURATION = Input::pp_COLLISION;

  //InputInit();

  //CDB flags and such

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG","ProdA_2024");
  rc->set_uint64Flag("TIMESTAMP",1);
  rc->set_IntFlag("RUNNUMBER",1);
/*
  Enable::MVTX_APPLYMISALIGNMENT = true;
  ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;

  SyncReco *sync = new SyncReco();
  se->registerSubsystem(sync);

  HeadReco *head = new HeadReco();
  se->registerSubsystem(head);

  FlagHandler *flag = new FlagHandler();
  se->registerSubsystem(flag);

  //Simulation setup
  Enable::MBDFAKE = true;
  Enable::PIPE = true;
  Enable::PIPE_ABSORBER = true;
  Enable::MVTX = true;
  Enable::INTT = true;
  Enable::TPC = true;
  Enable::MICROMEGAS = true;

  //Tracking setup

  G4Init();
  MagnetInit();
  MagnetFieldInit();

  G4Setup();
*/
  //Tagging stuff
  float abs_eta = 1.2;
  DecayFinder *myFinder = new DecayFinder("myFinder");
  myFinder->Verbosity(INT_MAX);
  if (channel == "Kshort2pipi") myFinder->setDecayDescriptor("K_S0 -> pi^- pi^+");
  else myFinder->setDecayDescriptor("[Lambda0 -> proton^+ pi^-]cc");
  myFinder->saveDST(1);
  myFinder->allowPi0(1);
  myFinder->allowPhotons(1);
  myFinder->triggerOnDecay(1);
  myFinder->setPTmin(0.); //Note: sPHENIX min pT is 0.2 GeV for tracking
  myFinder->setEtaRange(-1*abs_eta, abs_eta); //Note: sPHENIX acceptance is |eta| <= 1.1
  myFinder->useDecaySpecificEtaRange(false);
  if (channel != "minBias") se->registerSubsystem(myFinder);  

  Mbd_Reco();
  Mvtx_Cells();
  Intt_Cells();
  TPC_Cells();
  Micromegas_Cells();

  TrackingInit();

  Mvtx_Clustering();
  Intt_Clustering();
  TPC_Clustering();
  Micromegas_Clustering();

  Tracking_Reco();

  auto vtxfinder = new PHSimpleVertexFinder;
  vtxfinder->Verbosity(0);
  vtxfinder->setDcaCut(0.5);
  vtxfinder->setTrackPtCut(-99999.);
  vtxfinder->setBeamLineCut(1);
  vtxfinder->setTrackQualityCut(1000000000);
  vtxfinder->setNmvtxRequired(2);
  vtxfinder->setOutlierPairCut(0.1);  
  se->registerSubsystem(vtxfinder);

  Global_Reco();

  build_truthreco_tables();

  HFTrackEfficiency *myTrackEff = new HFTrackEfficiency("myTrackEff");
  myTrackEff->Verbosity(INT_MAX);
  myTrackEff->setDFNodeName("myFinder");
  myTrackEff->triggerOnDecay(1);
  myTrackEff->writeSelectedTrackMap(true);
  myTrackEff->writeOutputFile(true);
  std::string outputHFEffFile = outDir + "/hfEff/outputHFTrackEff_" + channel + "_" + processID + ".root";
  myTrackEff->setOutputFileName(outputHFEffFile);
  if (channel != "minBias") se->registerSubsystem(myTrackEff);

  output_dir = outDir;

  if (run_pipi_reco) create_hf_directories(pipi_reconstruction_name, pipi_output_dir, pipi_output_reco_file);
  if (run_ppi_reco) create_hf_directories(ppi_reconstruction_name, ppi_output_dir, ppi_output_reco_file);


  if (run_pipi_reco) reconstruct_pipi_mass();
  if (run_ppi_reco) reconstruct_ppi_mass();

  se->run(nEvents);

  se->End();

  if (run_pipi_reco) end_kfparticle(pipi_output_reco_file, pipi_output_dir);
  if (run_ppi_reco) end_kfparticle(ppi_output_reco_file, ppi_output_dir);

  gSystem->Exit(0);

  return 0;
}

#endif
