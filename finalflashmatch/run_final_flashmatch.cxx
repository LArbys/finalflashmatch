#include <iostream>
#include <string>
#include <vector>
#include <sstream>

// ROOT
#include "TFile.h"
#include "TChain.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TCanvas.h"

// larlite
#include "DataFormat/opflash.h"
#include "DataFormat/vertex.h"
#include "DataFormat/cluster.h"
#include "DataFormat/pfpart.h"
#include "DataFormat/shower.h"

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

  if ( nargs!=2 ) {
    std::cout << "usage: ./run_final_flashmatch [cfg file]" << std::endl;
    return -1;
  }
  
  // Setup config
  // ------------
  std::string cfg_file = argv[1];

  // tagger routine configuration
  larcv::PSet cfg  = larcv::CreatePSetFromFile( cfg_file );
  larcv::PSet pset = cfg.get<larcv::PSet>("FFMatch");


  
  // Setup input
  // ------------
  
  // Opreco data
  larlitecv::DataCoordinator dataco_opreco; // this should load opreco and ssnet file, then we have fucking everything
  std::string inputlist_opreco = pset.get<std::string>("LArLiteOpRecoFilelist",  "larlite" );
  dataco_opreco.add_inputfile( "extbnb_example2/larlite_opreco_55010101.root",   "larlite" );
  dataco_opreco.initialize();

  // shower reco data: the driver
  larlitecv::DataCoordinator dataco_showerreco;
  std::string inputlist_showerreco = pset.get<std::string>("LArLiteShowerFilelist");
  std::string inputlist_ssnet      = pset.get<std::string>("LArCVSSNetFilelist");  
  dataco_showerreco.set_filelist( inputlist_showerreco,  "larlite"   );
  dataco_showerreco.set_filelist( inputlist_ssnet,       "larcv"   );
  dataco_showerreco.initialize();


  // Setup config
  // ------------
  
  // proton 3d reco
  TChain ttracker("_recoTree");
  std::vector<int> *preco_status = 0;
  std::vector<double> *pmuon_energy = 0;
  std::vector<double> *pproton_energy = 0;
  double proton_endpt[3];
  double proton_startpt[3];
  ttracker.SetBranchAddress( "ProtonStartPoint_X", &proton_startpt[0] );
  ttracker.SetBranchAddress( "ProtonStartPoint_Y", &proton_startpt[1] );
  ttracker.SetBranchAddress( "ProtonStartPoint_Z", &proton_startpt[2] );    
  ttracker.SetBranchAddress( "ProtonEndPoint_X", &proton_endpt[0] );
  ttracker.SetBranchAddress( "ProtonEndPoint_Y", &proton_endpt[1] );
  ttracker.SetBranchAddress( "ProtonEndPoint_Z", &proton_endpt[2] );
  ttracker.SetBranchAddress( "Reco_goodness_v",  &preco_status );
  ttracker.SetBranchAddress( "E_muon_v",         &pmuon_energy );
  ttracker.SetBranchAddress( "E_proton_v",       &pproton_energy );

  std::string inputlist_tracker = pset.get<std::string>("TrackerInputlist");
  std::ifstream tracker_flist( inputlist_tracker );
  char trackerbuf[2056];
  std::vector<std::string> tracker_flist_v;  
  do {
    tracker_flist >> trackerbuf;
    std::string tmp = trackerbuf;
    if ( !tmp.empty() )
      tracker_flist_v.push_back( tmp );
    else
      break;
  } while( tracker_flist.good() && !tracker_flist.eof() );
      
  for ( auto const& tracker_file : tracker_flist_v )
    ttracker.Add( tracker_file.c_str() );
  ttracker.BuildIndex( "run", "event" );
  

  // laod input shower list
  std::string inputlist_showerqual = pset.get<std::string>("ShowerqualInputlist");
  std::ifstream showerqual_flist( inputlist_showerqual );
  char showerqualbuf[2056];
  std::vector<std::string> showerqual_flist_v;  
  do {
    showerqual_flist >> showerqualbuf;
    std::string tmp = showerqualbuf;
    if ( !tmp.empty() )
      showerqual_flist_v.push_back( tmp );
    else
      break;
  } while( showerqual_flist.good() && !showerqual_flist.eof() );
        
  
  // shower info
  TChain tevshower("fEventTree_nueshowers");
  int nshowers;
  tevshower.SetBranchAddress("n_recoshowers", &nshowers);
  for ( auto const& showerqual_file : showerqual_flist_v )
    tevshower.Add( showerqual_file.c_str() );
  tevshower.BuildIndex("run","event");
  

  // shower reco
  TChain tshower("fShowerTree_nueshowers");
  double shower_dir[3];
  double shower_pos[3];
  double shower_energy[3];
  double shower_len;
  double shower_mclen = -1.0;
  double shower_mcenergy = -1.0;
  tshower.SetBranchAddress("reco_x",     &shower_pos[0]);
  tshower.SetBranchAddress("reco_y",     &shower_pos[1]);
  tshower.SetBranchAddress("reco_z",     &shower_pos[2]);
  tshower.SetBranchAddress("reco_dcosx", &shower_dir[0]);
  tshower.SetBranchAddress("reco_dcosy", &shower_dir[1]);
  tshower.SetBranchAddress("reco_dcosz", &shower_dir[2]);
  tshower.SetBranchAddress("reco_energy_U", &shower_energy[0]);
  tshower.SetBranchAddress("reco_energy_V", &shower_energy[1]);
  tshower.SetBranchAddress("reco_energy_Y", &shower_energy[2]);  
  tshower.SetBranchAddress("shower_len", &shower_len);
  tshower.SetBranchAddress("mc_length",  &shower_mclen);
  tshower.SetBranchAddress("mc_energy",  &shower_mcenergy);

  for ( auto const& showerqual_file : showerqual_flist_v )
    tshower.Add( showerqual_file.c_str() );
  tshower.BuildIndex("run","event");
  
  // =========================
  // output file
  // -------------------------
  std::string outfile = pset.get<std::string>("OutputFilepath");
  TFile out( outfile.c_str(),"recreate");
  TTree outtree("ffmatch","Final Flash Match result");
  float best_chi2;
  float best_data_totpe;
  float hypo_totpe;
  float hypo_pe[32];
  float hypo_protonpe[32];
  float hypo_showerpe[32];  
  float data_pe[32];
  float vtxpos[3];
  outtree.Branch("vtxpos",vtxpos,"vtxpos[3]/F");
  outtree.Branch("chi2",&best_chi2,"chi2/F");
  outtree.Branch("data_totpe",&best_data_totpe,"data_totpe/F");
  outtree.Branch("hypo_totpe",&hypo_totpe,"hypo_totpe/F");
  outtree.Branch("hypo_pe",hypo_pe,"hypo_pe[32]/F");
  outtree.Branch("data_pe",data_pe,"data_pe[32]/F");
  outtree.Branch("hypo_protonpe",hypo_protonpe,"hypo_protonpe[32]/F");
  outtree.Branch("hypo_showerpe",hypo_showerpe,"hypo_showerpe[32]/F");


  // --------------
  // Algo Setup
  // --------------
  
  // Setup the FlashMatch Interface
  larcv::PSet genflash_pset = pset.get<larcv::PSet>("GeneralFlashMatchAlgo");
  larlitecv::GeneralFlashMatchAlgoConfig genflash_cfg = larlitecv::GeneralFlashMatchAlgoConfig::MakeConfigFromPSet( genflash_pset );
  larlitecv::GeneralFlashMatchAlgo genflashmatch( genflash_cfg );
  float shower_correction_factor = pset.get<float>("ShowerCorrectionFactor");
  bool hasMC  = pset.get<bool>("HasMC");

  // Setup the tracker
  larcv::AStarTracker astracker_algo;
  std::string _spline_file=pset.get<std::string>("TrackerSplineFile"); // example: "../LArCV/app/Reco3D/Proton_Muon_Range_dEdx_LAr_TSplines.root";
  astracker_algo.SetSplineFile(_spline_file);
  astracker_algo.initialize();
  astracker_algo.SetCompressionFactors(1,6);
  astracker_algo.SetVerbose(0);


  
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
    larcv::EventImage2D* ev_image2d = (larcv::EventImage2D*)dataco_showerreco.get_larcv_data( larcv::kProductImage2D, "modimg" );
    const std::vector<larcv::Image2D>& img_v = ev_image2d->Image2DArray();

    // load badch image
    larcv::EventImage2D* ev_badch   = (larcv::EventImage2D*)dataco_showerreco.get_larcv_data( larcv::kProductImage2D, "gapchs" );
    const std::vector<larcv::Image2D>& badch_v = ev_badch->Image2DArray();

    // load tagger image
    larcv::EventImage2D* ev_tagger  = (larcv::EventImage2D*)dataco_showerreco.get_larcv_data( larcv::kProductImage2D, "combinedtags" );
    const std::vector<larcv::Image2D>& tagger_v = ev_tagger->Image2DArray();

    // load vertex/pf-particle
    larlite::event_vertex*  ev_vertex  = (larlite::event_vertex*) dataco_showerreco.get_larlite_data( larlite::data::kVertex,     "dl" );
    larlite::event_pfpart*  ev_pfpart  = (larlite::event_pfpart*) dataco_showerreco.get_larlite_data( larlite::data::kPFParticle, "dl" );    
    larlite::event_cluster* ev_cluster = (larlite::event_cluster*)dataco_showerreco.get_larlite_data( larlite::data::kCluster,    "dl" );
    larlite::event_shower*  ev_shower  = (larlite::event_shower*) dataco_showerreco.get_larlite_data( larlite::data::kShower,     "dl" );
    larlite::event_shower*  ev_shreco  = (larlite::event_shower*) dataco_showerreco.get_larlite_data( larlite::data::kShower,     "showerreco" );        
        
    std::cout << "  number showers: " << nshowers << std::endl;
    std::cout << "  number of beam-window flashes: " << ev_opflash->size() << std::endl;

    if ( nshowers==0 ) {
      std::cout << "  no shower. continue." << std::endl;
      continue;
    }

    if ( !ev_pfpart->empty() )
      std::cout << "  has pfparticle: " << ev_pfpart->size() << std::endl;
    if ( !ev_cluster->empty() )
      std::cout << "  has cluster: " << ev_cluster->size() << std::endl;
    if ( ev_vertex->size() )
      std::cout << "  has vertex: " << ev_vertex->size() << std::endl;
    if ( ev_shreco->size() )
      std::cout << "  has shower (showerreco): " << ev_shreco->size() << std::endl;
    

    ULong_t showerbytes = tshower.GetEntryWithIndex( event, subrun );
    if ( showerbytes==0 ) {
      std::cout << "  could not load shower info for (event,subrun)=(" << event << "," << subrun << ")" << std::endl;
      continue;
    }

    // ----------------------
    // start analysis
    // ----------------------

    const larlite::vertex& vtx = ev_vertex->front();
    std::cout << "  vertex: (" << vtx.X() << "," << vtx.Y() << "," << vtx.Z() << ")" << std::endl;

    // Light Yields
    // proton: 19200
    // electron: 20000
    // muon: 24000
    
    // make the proton qcluster
    flashana::QCluster_t qinteraction;
    flashana::QCluster_t qshower;
    flashana::QCluster_t qproton;

    // get energy
    std::cout << "  proton energy: ";
    for ( auto& protonenergy : *pproton_energy ) {
      std::cout << " " << protonenergy;
    }
    std::cout << std::endl;
    
    // get positions
    float maxstepsize = 0.3;
    std::cout << "  proton "
	      << " start=(" << proton_startpt[0] << "," << proton_startpt[1] << "," << proton_startpt[2] << ") -> "
	      << " end=(" << proton_endpt[0] << "," << proton_endpt[1] << "," << proton_endpt[2] << ")"
	      << std::endl;
    float pdir[3] = {0};
    float plen = 0.;
    for (int i=0; i<3; i++) {
      pdir[i] = proton_endpt[i]-proton_startpt[i];
      plen += pdir[i]*pdir[i];
    }
    plen = sqrt(plen);
    for (int i=0; i<3; i++)
      pdir[i] /= plen;

    // get energy from range
    double proton_energy_fromrange = astracker_algo.GetEnergy("proton",plen); // MeV
    double proton_dedx = proton_energy_fromrange/plen;

    std::cout << "  proton energy from range: " << proton_energy_fromrange << " MeV"
	      << " length=" << plen << " cm"
	      << " de/dx=" << proton_dedx
	      << std::endl;
      
    
    int nsteps = plen/maxstepsize+1;
    float step = plen/float(nsteps);
    for (int istep=0; istep<=nsteps; istep++) {
      float pos[3];
      for (int i=0; i<3; i++)
	pos[i] = proton_startpt[i] + (step*istep)*pdir[i];
      float numphotons = proton_dedx*step*19200;
      flashana::QPoint_t qpt( pos[0], pos[1], pos[2], numphotons );
      qinteraction.push_back( qpt );
      qproton.push_back( qpt );
    }
    
    
    // make shower qcluser: each point corresponds to the number of photons
    const larlite::shower& shreco = ev_shreco->front();
    std::cout << "  shower: "
	      << " dir=(" << shreco.Direction().X() << "," << shreco.Direction().Y() << "," << shreco.Direction().Z() << ")"
	      << " dir2=(" << shower_dir[0] << "," << shower_dir[1] << "," << shower_dir[2] << ") "
	      << " length=" << shreco.Length() << "/" << shower_len
	      << " energy=(" << shower_energy[0] << "," << shower_energy[1] << "," << shower_energy[2] << ")";
    if ( hasMC ) {
      std::cout << " mclength=" << shower_mclen
		<< " mcenergy=" << shower_mcenergy;
    }
    std::cout << std::endl;
    
    if ( shower_energy[2]<1.0 )
      shower_energy[2] = 0.5*(shower_energy[0]+shower_energy[1]);
    shower_energy[2] *= 2.0;
    shower_len = shower_energy[2]/2.4; // MeV / (MeV/cm)

    std::cout << " shower: used energy=" << shower_energy[2] << " used length=" << shower_len << std::endl;
    
    nsteps = shower_len/maxstepsize+1;
    step = plen/float(nsteps);
    for (int istep=0; istep<=nsteps; istep++) {
      double pos[3];
      pos[0] = vtx.X() + (step*istep)*shower_dir[0];
      pos[1] = vtx.Y() + (step*istep)*shower_dir[1];
      pos[2] = vtx.Z() + (step*istep)*shower_dir[2];      
      float numphotons = (shower_correction_factor*step)*20000;
      flashana::QPoint_t qpt( pos[0], pos[1], pos[2], numphotons );
      qinteraction.push_back( qpt );
      qshower.push_back( qpt );
    }


    // convert opflash into flash_t
    std::vector<flashana::Flash_t> dataflash_v = genflashmatch.MakeDataFlashes( *ev_opflash );

    
    // make flash hypothesis
    flashana::Flash_t hypo        = genflashmatch.GenerateUnfittedFlashHypothesis( qinteraction );
    flashana::Flash_t hypo_proton = genflashmatch.GenerateUnfittedFlashHypothesis( qproton );
    flashana::Flash_t hypo_shower = genflashmatch.GenerateUnfittedFlashHypothesis( qshower );    

    // make flash hist
    std::stringstream hname_hypo;
    hname_hypo << "hflash_" << run << "_" << event << "_" << subrun << "_hypo";
    std::stringstream hname_hypo_proton;
    hname_hypo_proton << "hflash_" << run << "_" << event << "_" << subrun << "_hypoproton";
    std::stringstream hname_hypo_shower;
    hname_hypo_shower << "hflash_" << run << "_" << event << "_" << subrun << "_hyposhower";

    TH1D flashhist_hypo( hname_hypo.str().c_str(),"",32,0,32);
    TH1D flashhist_hypo_proton( hname_hypo_proton.str().c_str(),"",32,0,32);
    TH1D flashhist_hypo_shower( hname_hypo_shower.str().c_str(),"",32,0,32);        
    float maxpe_hypo = 0.;
    float petot_hypo = 0.;
    for (int i=0; i<32; i++){
      flashhist_hypo.SetBinContent( i+1, hypo.pe_v[i] );
      flashhist_hypo_proton.SetBinContent(i+1, hypo_proton.pe_v[i] );
      flashhist_hypo_shower.SetBinContent(i+1, hypo_shower.pe_v[i] );      
      if ( maxpe_hypo<hypo.pe_v[i] )
	maxpe_hypo = hypo.pe_v[i];
      petot_hypo += hypo.pe_v[i];
    }// num of pmts
    maxpe_hypo /= petot_hypo;
    flashhist_hypo.Scale(1.0/petot_hypo);
    flashhist_hypo_proton.Scale(1.0/petot_hypo);
    flashhist_hypo_shower.Scale(1.0/petot_hypo);    
    flashhist_hypo.SetLineColor(kRed);
    flashhist_hypo_proton.SetLineColor(kCyan);
    flashhist_hypo_shower.SetLineColor(kMagenta);    
    
    std::vector<TH1D*> flashhist_data_v;
    std::vector<float> petot_data_v;    
    float maxpe_data = 0.;

    for ( size_t i=0; i<ev_opflash->size(); i++) {
      std::stringstream hname_data;
      hname_data << "hflash_" << run << "_" << event << "_" << subrun << "_data" << i;
      TH1D* flashhist_data = new TH1D( hname_data.str().c_str(),"",32,0,32);
      float datatot = 0.;
      for (int ipmt=0; ipmt<32; ipmt++) {
	flashhist_data->SetBinContent(ipmt+1,dataflash_v[i].pe_v[ipmt]);
	datatot += dataflash_v[i].pe_v[ipmt];
	//std::cout << "data [" << ipmt << "] " << dataflash_v[i].pe_v[ipmt] << std::endl;
	if ( dataflash_v[i].pe_v[ipmt]>maxpe_data )
	  maxpe_data = dataflash_v[i].pe_v[ipmt];
	flashhist_data->SetLineColor(kBlack);
      }
      flashhist_data->Scale( 1.0/datatot );
      flashhist_data_v.push_back( flashhist_data );
      maxpe_data /= datatot;
      petot_data_v.push_back( datatot );
    }// end of flashes
    
    TCanvas cflash(hname_hypo.str().c_str(),"",800,600);
    if ( maxpe_hypo<maxpe_data )
      flashhist_hypo.GetYaxis()->SetRangeUser(0,maxpe_data*1.1);
    flashhist_hypo.Draw("hist");
    flashhist_hypo_proton.Draw("histsame");
    flashhist_hypo_shower.Draw("histsame");
    for ( auto const& phist : flashhist_data_v ) {
      phist->Draw("same");
    }
    hname_hypo << ".png";

    cflash.SaveAs( hname_hypo.str().c_str() );


    // calculate simple chi2
    best_chi2 = -1;
    best_data_totpe = -1;
    int best_data = 0;
    int idata = -1;
    for ( auto& phist : flashhist_data_v ) {
      idata++;
      float chi2 = 0.;
      for (int ipmt=0; ipmt<32; ipmt++) {
	float pefrac_data = phist->GetBinContent(ipmt+1);
	float pefrac_hypo = flashhist_hypo.GetBinContent(ipmt+1);

	float pe_data = pefrac_data*petot_data_v[idata];
	float pefrac_data_err = sqrt(pe_data)/petot_data_v[idata];

	float diff = pefrac_hypo - pefrac_data;
	
	std::cout << "[" << ipmt << "] diff=" << diff << " hypo=" << pefrac_hypo << " data=" << pefrac_data << std::endl;
	
	// i know this is all fubar
	if ( pefrac_data_err>0 )
	  chi2 += (diff*diff)/pefrac_data;
	else if (pefrac_data_err==0.0 && pefrac_hypo>0){
	  chi2 += (diff*diff)/pefrac_hypo;
	}
	
      }
      if ( best_chi2<0 || best_chi2>chi2 ) {
	best_chi2 = chi2;
	best_data_totpe = petot_data_v[idata];
	best_data = idata;
      }
    }

    for ( auto& phist : flashhist_data_v )
      delete phist;

    outtree.Fill();
  }

  out.Write();
  
  return 0;
}
