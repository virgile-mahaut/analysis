#include <fun4all/Fun4AllUtils.h>
#include <G4_ActsGeom.C>
#include <G4_Global.C>
#include <G4_Magnet.C>
#include <G4_Mbd.C>
#include <GlobalVariables.C>
#include <Trkr_Clustering.C>
#include <Trkr_LaserClustering.C>
#include <Trkr_Reco.C>
#include <Trkr_RecoInit.C>
#include <Trkr_TpcReadoutInit.C>
#include <QA.C>
#include <Calo_Calib.C>

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>

#include <mvtxrawhitqa/MvtxRawHitQA.h>
#include <inttrawhitqa/InttRawHitQA.h>
#include <tpcqa/TpcRawHitQA.h>
#include <trackingqa/InttClusterQA.h>
#include <trackingqa/MicromegasClusterQA.h>
#include <trackingqa/MvtxClusterQA.h>
#include <trackingqa/TpcClusterQA.h>
#include <trackingqa/SiliconSeedsQA.h>
#include <trackingqa/TpcSeedsQA.h>
#include <trackingqa/TrackFittingQA.h>
#include <trackingqa/TpcSiliconQA.h>
#include <trackingqa/VertexQA.h>
#include <trackreco/PHActsTrackProjection.h>

#include <caloreco/CaloGeomMapping.h>
#include <caloreco/RawClusterBuilderTemplate.h>

#include <kfparticle_sphenix/KFParticle_sPHENIX.h>

#include <stdio.h>

#include "HF_selections_SQM26.C"

using namespace HeavyFlavorReco;

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libcdbobjects.so)

bool IsListFile = false;
bool IsDstFile = false;
bool DoUnpacking = false;
bool DoClustering = false;
bool DoSeeding = false;
bool DoFitting = false;

void CheckDstType(const std::string inputDST)
{
    if (inputDST.find("DST_STREAMING_EVENT") != std::string::npos)
    {
        DoUnpacking = true;
        DoClustering = true;
        DoSeeding = true;
        DoFitting = true;
    }
    else if (inputDST.find("DST_TRKR_CLUSTER") != std::string::npos)
    {
        DoUnpacking = false;
        DoClustering = false;
        DoSeeding = true;
        DoFitting = true;
    }
    else if (inputDST.find("DST_TRKR_SEED") != std::string::npos)
    {
        DoUnpacking = false;
        DoClustering = false;
        DoSeeding = false;
        DoFitting = true;
    }
    else if (inputDST.find("DST_TRKR_TRACKS") != std::string::npos)
    {
        DoUnpacking = false;
        DoClustering = false;
        DoSeeding = false;
        DoFitting = false;
    }
    return;
}

std::string GetRunInfo(const std::string &filename)
{
    const std::vector<std::string> runspecies_patterns = {"run2pp", "run3pp", "run3auau", "run3oo"};

    boost::char_separator<char> sep("_");
    boost::tokenizer<boost::char_separator<char>> tok(filename, sep);

    for (const auto &t : tok)
    {
        for (const auto &pattern : runspecies_patterns)
        {
            if (t == pattern)
            {
                return t;
            }
        }
    }

    std::cout << "Cannot extract run info from filename!" << std::endl;
    return "";
}

void Fun4All_HF_SQM26(const int nEvents = 50,                                                                                                                                                                                            //
                const std::string inputDST = "/sphenix/lustre01/sphnxpro/production/run3oo/physics/ana534_2025p010_v001/DST_TRKR_SEED/run_00082400_00082500/DST_TRKR_SEED_run3oo_ana534_2025p010_v001-00082420-00663.root", //
                const std::string outDir = "./",                                                                                                                                                                                    //
                const int nSkip = 0,                                                                                                                                                                                                //
                const bool convertSeeds = false)                                                                                                                                                                                    //
{
    auto se = Fun4AllServer::instance();
    se->Verbosity(1);
    auto rc = recoConsts::instance();

    if (inputDST.find(".list") != std::string::npos)
    {
        IsListFile = true;
        std::cout << "Input is list file" << std::endl;
    }
    if (inputDST.find(".root") != std::string::npos)
    {
        IsDstFile = true;
        std::cout << "Input is dst file" << std::endl;
    }
    if (!IsListFile && !IsDstFile)
    {
        std::cout << "Check your input! Exit" << std::endl;
        gSystem->Exit(0);
    }

    int runnumber = std::numeric_limits<int>::quiet_NaN();
    int segment = std::numeric_limits<int>::quiet_NaN();
    std::string runspecies = "";

    if (IsListFile)
    {
        auto *hitsin = new Fun4AllDstInputManager("InputManager");
        hitsin->AddListFile(inputDST);
        se->registerInputManager(hitsin);

        std::ifstream ifs(inputDST);
        std::string filepath;
        int i = 0;
        while (std::getline(ifs, filepath))
        {
            std::cout << "Adding DST with filepath: " << filepath << std::endl;
            if (i == 0)
            {
                std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(filepath);
                runnumber = runseg.first;
                segment = runseg.second;
                runspecies = GetRunInfo(filepath);
                CheckDstType(filepath);
                break;
            }
            i++;
        }
    }

    if (IsDstFile)
    {
        std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(inputDST);
        runnumber = runseg.first;
        segment = runseg.second;
        runspecies = GetRunInfo(inputDST);

        CheckDstType(inputDST);

        auto *hitsin = new Fun4AllDstInputManager("InputManager");
        hitsin->fileopen(inputDST);
        se->registerInputManager(hitsin);
    }

    rc->set_IntFlag("RUNNUMBER", runnumber);
    rc->set_IntFlag("RUNSEGMENT", segment);
    rc->set_uint64Flag("TIMESTAMP", runnumber);

    Enable::CDB = true;
    rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");
    rc->set_uint64Flag("TIMESTAMP", runnumber);

    std::stringstream nice_runnumber;
    nice_runnumber << std::setw(8) << std::setfill('0') << to_string(runnumber);

    int rounded_up = 100 * (std::ceil((float)(runnumber+0.1) / 100));
    std::stringstream nice_rounded_up;
    nice_rounded_up << std::setw(8) << std::setfill('0') << to_string(rounded_up);

    int rounded_down = 100 * (std::floor((float)(runnumber+0.1) / 100));
    std::stringstream nice_rounded_down;
    nice_rounded_down << std::setw(8) << std::setfill('0') << to_string(rounded_down);

    std::stringstream nice_segment;
    nice_segment << std::setw(5) << std::setfill('0') << to_string(segment);

    std::stringstream nice_skip;
    nice_skip << std::setw(5) << std::setfill('0') << to_string(nSkip);

    std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");

    Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
    ingeo->AddFile(geofile);
    se->registerInputManager(ingeo);

    // set flags
    TRACKING::pp_mode = true;

    Enable::MVTX_APPLYMISALIGNMENT = true;
    ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;

    // distortion calibration mode
    /*
     * set to true to enable residuals in the TPC with
     * TPC clusters not participating to the ACTS track fit
     */
    G4TRACKING::SC_CALIBMODE = false;

    TpcReadoutInit(runnumber);
    // these lines show how to override the drift velocity and time offset values set in TpcReadoutInit
    // G4TPC::tpc_drift_velocity_reco = 0.0073844; // cm/ns
    // TpcClusterZCrossingCorrection::_vdrift = G4TPC::tpc_drift_velocity_reco;
    // G4TPC::tpc_tzero_reco = -5*50;  // ns
    std::cout << " run: " << runnumber << " samples: " << TRACKING::reco_tpc_maxtime_sample << " pre: " << TRACKING::reco_tpc_time_presample << " vdrift: " << G4TPC::tpc_drift_velocity_reco << std::endl;

    G4TPC::REJECT_LASER_EVENTS = true;
    G4TPC::ENABLE_MODULE_EDGE_CORRECTIONS = true;
    // Flag for running the tpc hit unpacker with zero suppression on
    TRACKING::tpc_zero_supp = true;

    // to turn on the default static corrections, enable the two lines below
    G4TPC::ENABLE_STATIC_CORRECTIONS = true;
    G4TPC::USE_PHI_AS_RAD_STATIC_CORRECTIONS = false;

    // to turn on the average corrections, enable the three lines below
    // note: these are designed to be used only if static corrections are also applied
    G4TPC::ENABLE_AVERAGE_CORRECTIONS = true;
    G4TPC::USE_PHI_AS_RAD_AVERAGE_CORRECTIONS = false;
    // to use a custom file instead of the database file:
    G4TPC::average_correction_filename = CDBInterface::instance()->getUrl("TPC_LAMINATION_FIT_CORRECTION");

    G4MAGNET::magfield_rescale = 1;
    TrackingInit();

    output_dir = "./"; // Top dir of where the output nTuples will be written
    trailer = "_" + nice_runnumber.str() + "_" + nice_segment.str() + "_" + nice_skip.str() + ".root";

    if (DoUnpacking)
    {
        for (int felix = 0; felix < 6; felix++)
        {
            Mvtx_HitUnpacking(std::to_string(felix));
        }
        for (int server = 0; server < 8; server++)
        {
            Intt_HitUnpacking(std::to_string(server));
        }
        std::ostringstream ebdcname;
        if (runspecies == "run3pp" || runspecies == "run3auau")
        {
            for (int ebdc = 0; ebdc < 24; ebdc++)
            {
                for (int endpoint = 0; endpoint < 2; endpoint++)
                {
                    ebdcname.str("");
                    if (ebdc < 10)
                    {
                        ebdcname << "0";
                    }
                    ebdcname << ebdc << "_" << endpoint;
                    Tpc_HitUnpacking(ebdcname.str());
                }
            }
        }
        if (runspecies == "run2pp")
        {
            for (int ebdc = 0; ebdc < 24; ebdc++)
            {
                ebdcname.str("");
                if (ebdc < 10)
                {
                    ebdcname << "0";
                }
                ebdcname << ebdc;
                Tpc_HitUnpacking(ebdcname.str());
            }
        }

        Micromegas_HitUnpacking();
    }

    if (DoClustering)
    {
        Mvtx_Clustering();

        Intt_Clustering();

        Tpc_LaserEventIdentifying();

        TPC_Clustering_run2pp();

        Micromegas_Clustering();

        Reject_Laser_Events();
    }

    if (DoSeeding)
    {
        Tracking_Reco_TrackSeed_run2pp();

        auto converter = new TrackSeedTrackMapConverter("SiliconSeedConverter");
        // Default set to full SvtxTrackSeeds. Can be set to
        // SiliconTrackSeedContainer or TpcTrackSeedContainer
        converter->setTrackSeedName("SiliconTrackSeedContainer");
        converter->setTrackMapName("SiliconSvtxTrackMap");
        converter->setFieldMap(G4MAGNET::magfield_tracking);
        converter->Verbosity(0);
        se->registerSubsystem(converter);

        auto finder = new PHSimpleVertexFinder("SiliconVertexFinder");
        finder->Verbosity(0);
        finder->setDcaCut(0.1);
        finder->setTrackPtCut(0.1);
        finder->setBeamLineCut(1);
        finder->setTrackQualityCut(1000000000);
        finder->setNmvtxRequired(3);
        finder->setOutlierPairCut(0.1);
        finder->setTrackMapName("SiliconSvtxTrackMap");
        finder->setVertexMapName("SiliconSvtxVertexMap");
        se->registerSubsystem(finder);

        auto convertertpc = new TrackSeedTrackMapConverter("TpcSeedConverter");
        // Default set to full SvtxTrackSeeds. Can be set to
        // SiliconTrackSeedContainer or TpcTrackSeedContainer
        convertertpc->setTrackSeedName("TpcTrackSeedContainer");
        convertertpc->setTrackMapName("TpcSvtxTrackMap");
        convertertpc->setFieldMap(G4MAGNET::magfield_tracking);
        convertertpc->Verbosity(0);
        se->registerSubsystem(convertertpc);

        auto findertpc = new PHSimpleVertexFinder("TpcSimpleVertexFinder");
        findertpc->Verbosity(0);
        findertpc->setDcaCut(0.5);
        findertpc->setTrackPtCut(0.2);
        findertpc->setBeamLineCut(1);
        findertpc->setTrackQualityCut(1000000000);
        // findertpc->setNmvtxRequired(3);
        findertpc->setRequireMVTX(false);
        findertpc->setOutlierPairCut(0.1);
        findertpc->setTrackMapName("TpcSvtxTrackMap");
        findertpc->setVertexMapName("TpcSvtxVertexMap");
        se->registerSubsystem(findertpc);
    }

    if (DoFitting)
    {
        Tracking_Reco_TrackMatching_run2pp();

        G4TRACKING::convert_seeds_to_svtxtracks = convertSeeds;
        std::cout << "Converting to seeds : " << G4TRACKING::convert_seeds_to_svtxtracks << std::endl;
        /*
         * Either converts seeds to tracks with a straight line/helix fit
         * or run the full Acts track kalman filter fit
         */
        if (G4TRACKING::convert_seeds_to_svtxtracks)
        {
            auto *converter = new TrackSeedTrackMapConverter;
            // Default set to full SvtxTrackSeeds. Can be set to
            // SiliconTrackSeedContainer or TpcTrackSeedContainer
            converter->setTrackSeedName("SvtxTrackSeedContainer");
            converter->setFieldMap(G4MAGNET::magfield_tracking);
            converter->Verbosity(0);
            se->registerSubsystem(converter);
        }
        else
        {
            Tracking_Reco_TrackFit_run2pp();
        }

        // vertexing and propagation to vertex
        Tracking_Reco_Vertex_run2pp();
    }

    output_dir = outDir;

    if (run_pipi_reco)
        create_hf_directories(pipi_reconstruction_name, pipi_output_dir, pipi_output_reco_file);
    if (run_Kpi_reco)
        create_hf_directories(Kpi_reconstruction_name, Kpi_output_dir, Kpi_output_reco_file);
    if (run_KK_reco)
        create_hf_directories(KK_reconstruction_name, KK_output_dir, KK_output_reco_file);
    if (run_ppi_reco)
        create_hf_directories(ppi_reconstruction_name, ppi_output_dir, ppi_output_reco_file);
    if (run_Lambdapi_reco)
        create_hf_directories(Lambdapi_reconstruction_name, Lambdapi_output_dir, Lambdapi_output_reco_file);

    if (run_pipi_reco || run_Kpi_reco || run_KK_reco || run_ppi_reco || run_Lambdapi_reco)// || run_ee_reco)
        init_kfp_dependencies();

    if (run_pipi_reco)
        reconstruct_pipi_mass();
    if (run_Kpi_reco)
        reconstruct_Kpi_mass();
    if (run_KK_reco)
        reconstruct_KK_mass();
    if (run_ppi_reco)
        reconstruct_ppi_mass();
    if (run_Lambdapi_reco)
        reconstruct_Lambdapi_mass();

    se->skip(nSkip);
    se->run(nEvents);
    se->End();
    se->PrintTimer();

    if (run_pipi_reco)
        end_kfparticle(pipi_output_reco_file, pipi_output_dir);
    if (run_Kpi_reco)
        end_kfparticle(Kpi_output_reco_file, Kpi_output_dir);
    if (run_KK_reco)
        end_kfparticle(KK_output_reco_file, KK_output_dir);
    if (run_ppi_reco)
        end_kfparticle(ppi_output_reco_file, ppi_output_dir);
    if (run_Lambdapi_reco)
        end_kfparticle(Lambdapi_output_reco_file, Lambdapi_output_dir);

    delete se;

    std::cout << "Finished" << std::endl;
    gSystem->Exit(0);
}
