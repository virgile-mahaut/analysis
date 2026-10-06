#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

#include <TSystem.h>

#include <upcmeson/UPCMeson.h>
#include <trackingdiagnostics/BeamCrossingAnalysis.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4dst.so)
R__LOAD_LIBRARY(libUPCMeson.so)
R__LOAD_LIBRARY(libTrackingDiagnostics.so)

//! Simple macro to analyze simDSTs or real DSTs for UPC analysis
void Fun4All_UPCMeson(const int nevnt = 0,
    const std::string& inputfile = "G4sPHENIX.root",
    std::string outputfile = "")
{
  Fun4AllServer *se = Fun4AllServer::instance();

  if ( outputfile.size() == 0 )
  {
    outputfile = "upcmeson_" + inputfile;
  }

  UPCMeson *upcmeson = new UPCMeson("upcmeson", outputfile);
  //upcmeson->Verbosity(10);
  upcmeson->analyzeTracks(true);
  upcmeson->analyzeTruth(false);
  upcmeson->SetGuessMass(0.13957039);
  se->registerSubsystem(upcmeson);

  BeamCrossingAnalysis *bcross = new BeamCrossingAnalysis("BEAMCROSS");
  bcross->set_output_file("beamcross.root");
  se->registerSubsystem(bcross);

  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTin");
  //Fun4AllInputManager *in2 = new Fun4AllDstInputManager("DSTin2");

  if ( inputfile.substr(inputfile.size()-5,inputfile.size()-1) == ".root" )
  {
    //in->fileopen(inputfile);
    in->AddFile(inputfile);
    //TString seedsfile = inputfile.c_str();
    //seedsfile.ReplaceAll("TRACKS","SEED");
    //std::cout << "SEEDS " << seedsfile << std::endl;
    //in2->AddFile(seedsfile.Data());
  }
  else if ( inputfile.substr(inputfile.size()-5,inputfile.size()-1) == ".list" )
  {
    in->AddListFile(inputfile);
  }
  se->registerInputManager(in);

  se->run(nevnt);
  se->End();
  delete se;

  cout <<"All done. Exiting..."<<endl;

  gSystem->Exit(0);
}

