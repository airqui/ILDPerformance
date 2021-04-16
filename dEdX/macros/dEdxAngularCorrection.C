#include <TPaveStats.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <TFitResult.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystemFile.h"


void dEdxAngularCorrection() {


  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetOptStat(0);

  TString input_filename = "../tests_adrian/condor_newpar/SingleParticle_MarlinTrks_pol2Correction.root";
  //"../tests_adrian/condor_removecorrection/SingleParticleDST_MarlinTrks_removeCorrection.root";
  // "../steer/rootfiles_2fhad/nominal_removecorrection.root";//condor_PFOs_removecorrection/SingleParticleDST_PFOs_removeCorrection.root";

  TFile *f = new TFile(input_filename);

  TH1F*  dEdx_deviation_lambda_dep = (TH1F*)f->Get("NormLambdaFullAll_1");

  TF1 *f1= new TF1("f1","[0]+[1]*x+[2]*x*x",0,78);
  TF1 *f2= new TF1("f2","[0]+[1]*x*x",0,78);
  TF1 *f3= new TF1("f3","[0]+[1]*cos(TMath::Pi()*x/180.)",0,78);

  TCanvas* c_test = new TCanvas("c_test","c_test",800,800);
  c_test->cd(1);
  dEdx_deviation_lambda_dep->GetXaxis()->SetTitle("lambda/deg ");
  dEdx_deviation_lambda_dep->GetYaxis()->SetTitle("#frac{dE}{dx} deviation");
  dEdx_deviation_lambda_dep->GetYaxis()->SetRangeUser(0.8,1.2);
  dEdx_deviation_lambda_dep->SetMarkerColor(1);
  dEdx_deviation_lambda_dep->SetMarkerStyle(21);
  dEdx_deviation_lambda_dep->SetLineColor(1);
  dEdx_deviation_lambda_dep->SetContour(5);
 

  dEdx_deviation_lambda_dep->Fit(f1,"QREM0");
  dEdx_deviation_lambda_dep->Fit(f2,"QREM0");
  dEdx_deviation_lambda_dep->Fit(f3,"QREM0");
  dEdx_deviation_lambda_dep->Draw("p");
  f1->SetLineColor(2);
  f1->Draw("same");
  f2->SetLineColor(4);
  f2->Draw("same");
  f3->SetLineColor(3);
  f3->Draw("same");

  TLegend * leg = new TLegend(0.2,0.8,0.9,0.9);
  leg->AddEntry(f1,"FIT OPTION 1: [0]+[1]*x+[2]*x*x, with x=theta in degrees","l");
  leg->AddEntry(f2,"FIT OPTION 2: [0]+[1]*x, with x=theta in degrees","l");
  leg->AddEntry(f3,"FIT OPTION 3: [0]+[1]*cos(TMath::Pi()*x/180.), with x=theta in degrees","l");
  leg->Draw();

  // ------------------------
  float a= f1->GetParameter(0);
  float b= f1->GetParameter(1);
  float c= f1->GetParameter(2);
  b/=a;
  c/=a;
  std::cout.setf(ios::scientific);
  // std::cout.setf(ios::showpoint);
  std::cout.precision(2);
  std::cout<<"FIT OPTION 1: [0]+[1]*x+[2]*x*x, with x=theta in degrees"<<std::endl;
  std::cout<<"<parameter name=\"AngularCorrectionParameters\" type=\"FloatVec\">"<<1.0<<" "<< b<<" "<<c<<" </parameter>"<<std::endl;

  
  //----  
  a= f2->GetParameter(0);
  b= f2->GetParameter(1);
  c= f2->GetParameter(2);
  b/=a;
  c/=a;
  std::cout.setf(ios::scientific);
  // std::cout.setf(ios::showpoint);
  std::cout.precision(2);
  std::cout<<"FIT OPTION 2: [0]+[1]*x*x, with x=theta in degrees"<<std::endl;
  std::cout<<"<parameter name=\"AngularCorrectionParameters\" type=\"FloatVec\">"<<1.0<<" "<< b<<" </parameter>"<<std::endl;

  //----  
  a= f3->GetParameter(0);
  b= f3->GetParameter(1);
  b/=a;
  std::cout.setf(ios::scientific);
  // std::cout.setf(ios::showpoint);
  std::cout.precision(2);
  std::cout<<"FIT OPTION 3: [0]+[1]*cos(TMath::Pi()*x/180.), with x=theta in degrees"<<std::endl;
  std::cout<<"<parameter name=\"AngularCorrectionParameters\" type=\"FloatVec\">"<<1.0<<" "<< b<<" </parameter>"<<std::endl;

}
