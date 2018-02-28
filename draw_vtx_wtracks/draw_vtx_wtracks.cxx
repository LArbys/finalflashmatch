#include <iostream>
#include <string>
#include <sstream>

// larcv
#include "DataFormat/EventImage2D.h"

// larlite
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

// larlitecv
#include "Base/DataCoordinator.h"

// opencv
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
// larcv opencv utils
#include "CVUtil/CVUtil.h"

int main( int nargs, char** argv ) {

  // draw event, labeling vertex and tracks for 1mu1p events
  // input:
  //  supera/ssnetout file: to get image
  //  opreco larlite file:  to get flash information
  //  vtx intermediate file: to get vertex informtion
  //  run,subrun,event: to ID event
  //  vtx_id: to pick out vertex information
  //  outdir: location to draw image
  //
  // output:
  //  3 pngs per event showing information on the planes

  // We use larlitecv as our interface to the larcv/larlite information
  // We use regular ROOT as interface to vtx file
  // We use opencv to draw
  // The runtime environment assumed with be a dllee_unified container

  // parse arguments
  std::string supera = argv[1];    // ssnetout
  //std::string interf = argv[2];    // intermediate file
  std::string tracker = argv[2];   // stage 2 tracker output
  int run    = std::atoi(argv[3]); 
  int subrun = std::atoi(argv[4]);
  int event  = std::atoi(argv[5]);
  int vtxid  = std::atoi(argv[6]);
  std::string outdir = argv[7];
  
  // for testing purposes
  // std::string supera  = "/home/twongjirad/working/data/larbys/mcc8v6/testdb/bnb5e19/run5925subrun0195/ssnetout-larcv-Run005925-SubRun000195.root";
  // std::string tracker = "/home/twongjirad/working/data/larbys/mcc8v6/testdb/bnb5e19/run5925subrun0195/tracker_reco_59250195.root";
  // std::string interf  = "/home/twongjirad/working/data/larbys/final_files/v1/mcc8v6_bnb5e19_test4_inter.root";
  // int run = 5925;
  // int subrun = 197;
  // int event = 9868;
  // int vtxid = 0;

  // parameters
  float fthreshold = 10.0;
  
  // open the files

  // larcv/larlite
  larlitecv::DataCoordinator dataco;
  dataco.add_inputfile( supera, "larcv" );
  dataco.add_inputfile( tracker, "larlite" );
  dataco.initialize();

  dataco.goto_event( run, subrun, event, "larcv" );
  
  // get the data we need

  // original image
  larcv::EventImage2D* ev_img = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "modimg" );
  larcv::EventImage2D* ev_uburn[3];
  ev_uburn[0] = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "uburn_plane0" );
  ev_uburn[1] = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "uburn_plane1" );
  ev_uburn[2] = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "uburn_plane2" );
  const std::vector<larcv::Image2D>& img_v  = ev_img->Image2DArray();
  const std::vector<larcv::Image2D>* uburn_v[3];
  for (int p=0; p<3; p++) {
    uburn_v[p] = &(ev_uburn[p]->Image2DArray());
  }

  // the reconstructed tracks
  larlite::event_track* ev_tracker = (larlite::event_track*) dataco.get_larlite_data( larlite::data::kTrack, "trackReco" );

  // the reconstucted vertices
  larlite::event_vertex* ev_vertex = (larlite::event_vertex*)dataco.get_larlite_data( larlite::data::kVertex, "trackReco" );

  // that's it! let's make the picture

  const larutil::Geometry* geo = larutil::Geometry::GetME();
  const float cm_per_tick = larutil::LArProperties::GetME()->DriftVelocity()*0.5; // cm/us x 0.5 us/tick
  std::vector<cv::Mat> cvimgs_v;    
  int p=0;
  for ( auto const& img : img_v ) {
    // get cvmat grey-scale version
    //cv::Mat cvimg = larcv::as_mat_greyscale2bgr( img, fthreshold, 200.0 );
    cv::Mat im_gray = larcv::as_gray_mat( img, fthreshold, 256.0, 1.0 );
    cv::Mat cvimg;
    cv::applyColorMap(im_gray, cvimg, cv::COLORMAP_JET);

    // draw tracks on the image
    int itrack = -1;
    for ( auto const& track : *ev_tracker ) {
      itrack++;
      int npts = track.NumberTrajectoryPoints();
      for ( int ipt=0; ipt<npts; ipt++ ) {
	auto const& posv = track.LocationAtPoint( ipt );
	float tick = posv.X()/cm_per_tick + 3200.0;
	Double_t dpos[3] = { posv.X(), posv.Y(), posv.Z() };
	int row    = img.meta().row(tick);
	int wire   = geo->WireCoordinate( dpos, img.meta().plane() );
	int col    = img.meta().col(wire);

	cv::circle( cvimg, cv::Point(col,row), 1, cv::Scalar(255,255,255,255), 1 );
      }
    }

    // draw the vertex
    int ivertex=-1;
    for ( auto const& vertex : *ev_vertex ) {
      ivertex++;
      Double_t dpos[3] = { vertex.X(), vertex.Y(), vertex.Z() };
      float tick = dpos[0]/cm_per_tick + 3200.0;
      int row    = img.meta().row(tick);
      int wire   = geo->WireCoordinate( dpos, img.meta().plane() );
      int col    = img.meta().col(wire);      

      if (ivertex==vtxid )
	cv::circle( cvimg, cv::Point(col,row), 4, cv::Scalar(0,0,255,255), 1 );
      else
	cv::circle( cvimg, cv::Point(col,row), 4, cv::Scalar(255,255,0,255), 1 );
    }
    
    // 
    // 	const larcv::Image2D& badch = badch_v[p];
    // 	for (int i=0; i<(int)badch.meta().cols(); i++) {
    // 	  if ( badch.pixel(10,i)>0 ) {
    // 	    for (int r=0; r<badch.meta().rows(); r++) {
    // 	      cvimg.at<cv::Vec3b>(cv::Point(i,r))[0] = 10;
    // 	      cvimg.at<cv::Vec3b>(cv::Point(i,r))[1] = 10;
    // 	      cvimg.at<cv::Vec3b>(cv::Point(i,r))[2] = 10;	      
    // 	    }
    // 	  }
    // 	}
    std::stringstream ss;
    ss << outdir << "/trackerimg_run" << run << "_subrun" << subrun << "_event" << event << "_plane" << p << ".png";
    cv::imwrite( ss.str(), cvimg );
    cvimgs_v.emplace_back(std::move(cvimg));
    p++;
  }

  // save to png

  
  // done

  return 0;
};
