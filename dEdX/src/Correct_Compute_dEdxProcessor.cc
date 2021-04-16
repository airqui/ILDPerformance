#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>

#include <marlin/Global.h>
#include <marlin/Processor.h>
#include <marlin/AIDAProcessor.h>

#include <lcio.h>
#include <EVENT/LCCollection.h>
#include <EVENT/Track.h>
#include <IMPL/TrackImpl.h>

#include "TCanvas.h"
#include "TImage.h"
#include "TStyle.h"
#include "TMath.h"

#include "Correct_Compute_dEdxProcessor.hh"
#include <signal.h>

Correct_Compute_dEdxProcessor aCorrect_Compute_dEdxProcessor ;

Correct_Compute_dEdxProcessor::Correct_Compute_dEdxProcessor()
  : Processor("Correct_Compute_dEdxProcessor") {
  
  // Processor description
  _description = "Correct_Correct_Compute_dEdxProcessor: Makes a hard angular-based correction of dEdx for all the Tracks in the event. ATTENTION: this processor rewrites the MarlinTrk Collection and it is to be used only for simulations produced with ILCSoft from v02 to v02-02-01" ;
  
  registerInputCollection(LCIO::TRACK,
			  "LDCTrackCollection",
			  "LDC track collection name",
			  _LDCTrackCollection,
			  std::string("MarlinTrkTracks"));

  registerProcessorParameter( "Write_dEdx",
			      "If set, the calculated dEdx value will be rewritten in the its corresponding track (default: false).",
			      _writedEdx,
			      bool(false));

  std::vector<float> _newpar = {  1.00e+00,
				  1.20e-03,
				  1.84e-05};
  registerProcessorParameter( "AngularCorrectionParameters",
			      "parameter for new angular correction dedx= uncorrected_dedx  / f, with f= pol2(theta)",
			      _par,
			      _newpar);

  std::vector<float> _oldpar = {0.635762, -0.0573237};
  registerProcessorParameter( "OutdatedAngularCorrectionParameters",
			      "parameters used in the outdated angular correction outdated_dedx= uncorrected_dedx / f, with f = [0] / ( [0]+[1]*cos(theta)*cos(theta) )",
			      _parv02_02_01,
			      _oldpar);

} 

void Correct_Compute_dEdxProcessor::init() { 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  // it's usually a good idea to
  printParameters();
  
}

void Correct_Compute_dEdxProcessor::processRunHeader( LCRunHeader* ) { 
} 

void Correct_Compute_dEdxProcessor::processEvent( LCEvent * evt ) { 

  //fill values
  if (_writedEdx)
    {
      _LDCCol = evt->getCollection( _LDCTrackCollection ) ;
      int nTrkCand = _LDCCol->getNumberOfElements();
      
      for (int iTRK=0;iTRK<nTrkCand;++iTRK) {
	
	TrackImpl * trkCand = (TrackImpl*) _LDCCol->getElementAt( iTRK );
	
	float dedx=trkCand->getdEdx();
	float dedx_error=trkCand->getdEdxError();
	float trklambda = trkCand->getTanLambda();
	//(trkCand->getTanLambda()/ fabs(trkCand->getTanLambda()) )* sqrt(trkCand->getTanLambda()*trkCand->getTanLambda()/(1.0+trkCand->getTanLambda()*trkCand->getTanLambda()));
	
	std::pair<float,float> new_dedx = getExtraCorrection(dedx, dedx_error, trklambda);
	
	streamlog_out(DEBUG) << "Original dEdx: " <<dedx <<" Error: "<<dedx_error <<std::endl;
	streamlog_out(DEBUG) << "NeW dEdx: " <<new_dedx.first <<" Error: "<<new_dedx.second <<std::endl;
	
	//fill values
	trkCand->setdEdx(new_dedx.first);
	trkCand->setdEdxError(new_dedx.second);
      }
    }   else 
    {
      streamlog_out(ERROR) << " Why do you use this processor and not re-write dEdx ?? Check your steering file." <<std::endl;
      raise(SIGSEGV);
    }
}

void Correct_Compute_dEdxProcessor::check( LCEvent * ) { 
}

void Correct_Compute_dEdxProcessor::end() { 
}


std::pair<float,float> Correct_Compute_dEdxProcessor::getExtraCorrection(float dedx, float dedx_error, float trklambda){

  float trkcos = sqrt(trklambda*trklambda/(1.0+trklambda*trklambda));

  
  // // DEFAULT CORRECTION APPLIED FOR THE MC2020 production obtained using versions v02-02-01 and previous
  // //cal. hit dep.
  // double f1 = 1 + std::exp(-nHit/_ncorrpar);
  // //cal. polar angle dep.
  // // double c=1.0/sqrt(1.0-trkcos*trkcos);
  // // double f2=1.0/(1.0-0.08887*std::log(c)); 

  // //cal. polar angle dep.   20160702
  // //double c = std::acos(trkcos);
  // //if(c>3.141592/2.0) c= 3.141592-c;
  // //double f2 = 1.0/std::pow(c, 0.0703);

  // //change polar angle dep. 20170901
  // double f2 = _acorrpar[0] / (_acorrpar[0] + _acorrpar[1] * trkcos * trkcos);
  // return dedx/(f1*f2);
  

  // FIRST WE UNDO THE CORRECTION APPLIED IN COMPUTE_DEDXPROCESSOR
  // Default values in v02-02-01 for the angular correction
  //std::vector<float> _acorrparv02_02_01 = {0.635762, -0.0573237};
  // set in the steering file
  double f2 = 1.0 / (1.0 + _parv02_02_01[1] * trkcos * trkcos/_parv02_02_01[0]);

  // ****************************************
  // NEW CORRECTION (to be applied to uncorrected dEdx)
  //Parametrization by A. Irles 2021/04/09
  //Using single particle mc2020 samples v02-02-01
  //Analyzed with dEdxAnalyser, fitting histogram NormLambdaFullAll_1
  float theta = fabs(atan(trklambda)*180./M_PI);
  double f3 = _par[0] / (_par[0] + _par[1] * theta  + _par[2] * pow(theta,2) );

  std::pair<float,float> ret;
  if(std::isnan(dedx))
    {
    ret.first =0;
    ret.second=0;
    streamlog_out(DEBUG4) << "NAN!" << std::endl;
  }
  else
  {
    ret.first = dedx*f2*f3;
    ret.second= dedx_error*f2*f3;
  }

  return ret;

}
