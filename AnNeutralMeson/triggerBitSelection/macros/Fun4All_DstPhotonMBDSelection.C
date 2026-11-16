#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/SubsysReco.h>

#include <phool/recoConsts.h>

#include <ffamodules/CDBInterface.h>
#include <GlobalVariables.C>

#include <ctime>
#include <iostream>

R__LOAD_LIBRARY(libfun4all.so)

#include <triggerbitselection/DstPhotonMBDSelection.h>
R__LOAD_LIBRARY(libtriggerbitselection.so)

void Fun4All_DstPhotonMBDSelection(int nevents = 1, const std::string &dst_triggered_name = "DST_CALOFITTING.root")
{

  // Print the date
  std::time_t t = std::time(0); // get time now
  std::tm* now = std::localtime(&t);
  std::cout << "date = " << (now->tm_year + 1900) << "-" << (now->tm_mon + 1) << "-" << now->tm_mday << std::endl;

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  // se->Verbosity(verbosity);
  recoConsts *rc = recoConsts::instance();
  //===============
  // conditions DB flags
  //===============
  pair<int, int> runseg = Fun4AllUtils::GetRunSegment(dst_triggered_name);
  int runnumber = runseg.first;
  //cout << "run number = " << runnumber << endl;

  // global tag
  // cf gmattson, apply the tower mask from the CDB
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  // // 64 bit timestamp
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  rc->set_IntFlag("RUNNUMBER", runnumber);
  
  Fun4AllInputManager *in = new Fun4AllDstInputManager("DST_TRIGGER");
  in->AddFile(dst_triggered_name);
  se->registerInputManager(in);

  std::string dstname = std::regex_replace(dst_triggered_name, std::regex("DST_CALOFITTING"), "DST_PHOTON_MBD_TRIGGER_LIST");
  std::cout << "dstname = " << dstname << std::endl;
  Fun4AllOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", dstname);
  // Only keep DST#Sync in the output node
  out->StripNode("14001");
  out->StripRunNode("TriggerRunInfo");
  out->StripRunNode("RunHeader");
  out->StripRunNode("CdbUrl");
  out->StripRunNode("Flags");
  out->SplitLevel(99);
  se->registerOutputManager(out);

  ///////////////////
  // analysis modules
  DstPhotonMBDSelection *selection = new DstPhotonMBDSelection("photonmodulename", runnumber);
  se->registerSubsystem(selection);


  // Only keep the relevant nodes: Event number + GL1 info (must be done before fileopen)
  se->BranchSelect("*",0);
  se->BranchSelect("DST#Sync", 1);
  se->BranchSelect("DST#PacketsKeep#14001", 1);
  se->run(nevents);
  se->End();
  se->PrintTimer();
  delete se;

  TFile* f_done_signal = new TFile("DONE.root","recreate");
  std::cout << "All done!" << std::endl;
  gSystem->Exit(0);
}

