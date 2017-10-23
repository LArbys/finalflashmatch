#include <iostream>
#include <string>
#include <vector>

// ROOT
#include "TChain.h"
#include "TBranch.h"

// larlite
#include "DataFormat/opflash.h"

// larcv
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "DataFormat/EventImage2D.h"
#include "Reco3D/AStarTracker.h" // Adrien's 3D tracker

// larlitecv
#include "Base/DataCoordinator.h"
#include "TaggerContourTools/BMTCV.h" // Contour Analysis
#include "FlashMatchInterface/GeneralFlashMatchAlgo.h"

#include "FinalFlashMatch.h"

int main(int nargs, char** argv ) {

  // ====================================================
  // Dev code for final selection.
  // We're given everything we've produced
  //  -- cosmic tags
  //  -- contouring
  //  -- candidate vertex
  //  -- 3D reco for candidate tracks
  //  -- flash matching tools
  // Want another factor of 100 for rejection
  // ====================================================

  // Input files
  // showerqual file: tree with reco shower info (filtered)
  // track reco file:  tree with muon or proton info (filtered)
  // shower reco file: tree with larlite info: flash and mctruth (filtered)
  // larlite opreco file: tree with flashes (unfiltered)

  // Setup config
  // ------------
  std::string cfg_file = "ffmatch.cfg";

  // tagger routine configuration
  larcv::PSet cfg  = larcv::CreatePSetFromFile( cfg_file );
  larcv::PSet pset = cfg.get<larcv::PSet>("FFMatch");


  
  // Setup input
  // ------------
  
  // Opreco data
  larlitecv::DataCoordinator dataco_opreco; // this should load opreco and ssnet file, then we have fucking everything
  dataco_opreco.configure( cfg_file, "StorageManager", "IOManager", "FFMatch" );
  dataco_opreco.add_inputfile( "1e1p_example/larlite_opreco_104036068.root",   "larlite"   );
  dataco_opreco.add_inputfile( "1e1p_example/ssnetout_larcv_104036068.root",   "larcv"   );    
  dataco_opreco.initialize();

  // shower reco data
  larlitecv::DataCoordinator dataco_showerreco;
  dataco_showerreco.add_inputfile( "1e1p_example/shower_reco_out_104036068.root",   "larlite"   );
  dataco_showerreco.initialize();

  // proton 3d reco
  TChain ttracker("_recoTree");
  std::vector<int> *preco_status = 0;
  TBranch *breco_status = 0;
  float proton_endpt[3];
  float proton_startpt[3];
  ttracker.SetBranchAddress( "ProtonStartPoint_X", &proton_startpt[0] );
  ttracker.SetBranchAddress( "ProtonStartPoint_Y", &proton_startpt[1] );
  ttracker.SetBranchAddress( "ProtonStartPoint_Z", &proton_startpt[1] );    
  ttracker.SetBranchAddress( "ProtonEndPoint_X", &proton_endpt[0] );
  ttracker.SetBranchAddress( "ProtonEndPoint_Y", &proton_endpt[1] );
  ttracker.SetBranchAddress( "ProtonEndPoint_Z", &proton_endpt[1] );
  ttracker.SetBranchAddress( "Reco_goodness_v",  &preco_status );
  ttracker.Add( "1e1p_example/tracker_anaout_104036068.root" );  

  // shower info
  TChain tevshower("fEventTree_nueshowers");
  int nshowers;
  tevshower.SetBranchAddress("n_recoshowers", &nshowers);
  tevshower.Add( "1e1p_example/showerqualsingle_104036068.root" );  

  // shower reco
  TChain tshower("fShowerTree_nueshowers");
  float shower_dir[3];
  float shower_pos[3];
  float shower_len;
  tshower.SetBranchAddress("reco_x",     &shower_pos[0]);
  tshower.SetBranchAddress("reco_y",     &shower_pos[1]);
  tshower.SetBranchAddress("reco_z",     &shower_pos[2]);
  tshower.SetBranchAddress("reco_dcosx", &shower_dir[0]);
  tshower.SetBranchAddress("reco_dcosy", &shower_dir[1]);
  tshower.SetBranchAddress("reco_dcosz", &shower_dir[2]);
  tshower.Add("1e1p_example/showerqualsingle_104036068.root");
  tshower.BuildIndex( "event", "subrun" );

  // Setup the FlashMatch Interface
  larcv::PSet genflash_pset = pset.get<larcv::PSet>("GeneralFlashMatchAlgo");
  larlitecv::GeneralFlashMatchAlgoConfig genflash_cfg = larlitecv::GeneralFlashMatchAlgoConfig::MakeConfigFromPSet( genflash_pset );
  larlitecv::GeneralFlashMatchAlgo genflashmatch( genflash_cfg );

  // Setup the tracker

  // Setup the contour
  
  // use shower reco as driver
  int nentries = dataco_showerreco.get_nentries("larlite");

  for (int ientry=0; ientry<nentries; ientry++) {

    // -----------------------
    // Load the entry
    // -----------------------
    dataco_showerreco.goto_entry(ientry,std::string("larlite"));
    int run, subrun, event;
    dataco_showerreco.get_id( run, subrun, event );
    std::cout << "Entry #" << ientry << ": (" << run << "," << subrun << "," << event << ")" << std::endl;    

    // sync up the opreco tree
    dataco_opreco.goto_event( run, subrun, event, std::string("larlite") );

    // load the track reco file
    ttracker.GetEntry(ientry);

    // load shower reco tree
    tevshower.GetEntry(ientry);

    // load opreco event container
    larlite::event_opflash* ev_opflash = (larlite::event_opflash*)dataco_opreco.get_larlite_data( larlite::data::kOpFlash, "simpleFlashBeam" );
    if ( !ev_opflash )
      std::cout << "  error loading opflash info" << std::endl;

    // load raw image
    larcv::EventImage2D* ev_image2d = (larcv::EventImage2D*)dataco_opreco.get_larcv_data( larcv::kProductImage2D, "modimg" );
    const std::vector<larcv::Image2D>& img_v = ev_image2d->Image2DArray();

    // load badch image
    larcv::EventImage2D* ev_badch   = (larcv::EventImage2D*)dataco_opreco.get_larcv_data( larcv::kProductImage2D, "gapchs" );
    const std::vector<larcv::Image2D>& badch_v = ev_badch->Image2DArray();

    // load tagger image
    larcv::EventImage2D* ev_tagger  = (larcv::EventImage2D*)dataco_opreco.get_larcv_data( larcv::kProductImage2D, "combinedtags" );
    const std::vector<larcv::Image2D>& tagger_v = ev_tagger->Image2DArray();

    // ----------------------
    // start analysis
    // ----------------------
    
    std::cout << "  number showers: " << nshowers << std::endl;
    std::cout << "  number of beam-window flashes: " << ev_opflash->size() << std::endl;

    if ( nshowers==0 ) {
      std::cout << "  no shower. continue." << std::endl;
      continue;
    }

    ULong_t showerbytes = tshower.GetEntryWithIndex( event, subrun );
    if ( showerbytes==0 ) {
      std::cout << "  could not load shower info for (event,subrun)=(" << event << "," << subrun << ")" << std::endl;
      continue;
    }

    // perform analysis

    // make the proton
    float maxstepsize = 0.3;
    std::vector<std::vector<float> > proton_track;
    

    
  }

  
  return 0;
}
