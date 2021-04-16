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


Double_t myBB(Double_t *x, Double_t *par){
  Float_t xx = x[0];
  Float_t bg = xx/par[0];
  Float_t bb = sqrt(bg*bg/(1.+bg*bg));
  Float_t tmax= par[3]*pow(bg,2.);
  
  //Double_t f = 1.35e-7*(0.5*par[1]*log(par[2]*pow(bg,2.)*tmax) - par[4]*bb*bb - par[5]*bg/2.)/(bb*bb);
  Double_t f = 1.35e-7*(0.5*par[1]*log(par[2]*pow(bg,2.)*tmax) - par[4]*bb*bb - par[5]*bg/2.)/(bb*bb);
  return f;
}

void dEdxPlotandFit() {


  gStyle->SetPadLeftMargin(0.2);

  TString input_filename = "../tests_adrian/condor_newpar/SingleParticle_MarlinTrks_pol2Correction.root";
  // "../steer/rootfiles_2fhad/nominal_removecorrection.root";//condor_PFOs_removecorrection/SingleParticleDST_PFOs_removeCorrection.root";

  TFile *f = new TFile(input_filename);

  TH2F*  dEdx_proton= (TH2F*)f->Get("MydEdxAnalyser/BBProtons");
  TH2F*  dEdx_muon= (TH2F*)f->Get("MydEdxAnalyser/BBMuons");
  TH2F*  dEdx_kaon= (TH2F*)f->Get("MydEdxAnalyser/BBKaons");
  TH2F*  dEdx_pion= (TH2F*)f->Get("MydEdxAnalyser/BBPions");
  TH2F*  dEdx_electron= (TH2F*)f->Get("MydEdxAnalyser/BBElectrons");

  TH1F*  dEdx_mean_proton= (TH1F*)f->Get("BBProtons_1");
  TH1F*  dEdx_mean_muon= (TH1F*)f->Get("BBMuons_1");
  TH1F*  dEdx_mean_kaon= (TH1F*)f->Get("BBKaons_1");
  TH1F*  dEdx_mean_pion= (TH1F*)f->Get("BBPions_1");
  TH1F*  dEdx_mean_electron= (TH1F*)f->Get("BBElectrons_1");

  TF1 *fproton = new TF1("fproton",myBB,0.3,100,6);
  fproton->FixParameter(0,0.938272);//the mass
  fproton->SetParameter(1,7.96117e-02);
  fproton->SetParameter(2,4.13335e+03);
  fproton->SetParameter(3,1.13577e+06);
  fproton->SetParameter(4,1.80555e-01);
  fproton->SetParameter(5,-3.15083e-04);

  TF1 *fmuon = new TF1("fmuon",myBB,0.1,100,6);
  fmuon->FixParameter(0,0.105658);//the mass
  fmuon->SetParameter(1,7.96117e-02);
  fmuon->SetParameter(2,4.13335e+03);
  fmuon->SetParameter(3,1.13577e+06);
  fmuon->SetParameter(4,1.80555e-01);
  fmuon->SetParameter(5,-3.15083e-04);

  TF1 *fpion = new TF1("fpion",myBB,0.1,100,6);
  fpion->FixParameter(0,0.139570);//the mass
  fpion->SetParameter(1,7.96117e-02);
  fpion->SetParameter(2,4.13335e+03);
  fpion->SetParameter(3,1.13577e+06);
  fpion->SetParameter(4,1.80555e-01);
  fpion->SetParameter(5,-3.15083e-04);


  TF1 *felectron = new TF1("felectron",myBB,0.1,100,6);
  felectron->FixParameter(0,0.000510998);//the mass
  felectron->SetParameter(1,-2.40638e-03);
  felectron->SetParameter(2,7.10337e-01);
  felectron->SetParameter(3,2.87718e-01);
  felectron->SetParameter(4,2.87718e-01);
  felectron->FixParameter(5,0);

  TF1 *fkaon = new TF1("fkaon",myBB,0.1,100,6);
  fkaon->FixParameter(0,0.49367);
  fkaon->SetParameter(1,7.96117e-02);
  fkaon->SetParameter(2,4.13335e+03);
  fkaon->SetParameter(3,1.13577e+06);
  fkaon->SetParameter(4,1.80555e-01);
  fkaon->SetParameter(5,-3.15083e-04);



  TCanvas* c_test = new TCanvas("c_test","c_test",1200,800);
  c_test->Divide(3,2);
  c_test->cd(1);
  gPad->SetLogx();
  dEdx_proton->GetXaxis()->SetTitle("momentum ");
  dEdx_proton->GetYaxis()->SetTitle("#frac{dE}{dx}");
  dEdx_proton->Draw("colz");
  dEdx_mean_proton->SetLineColor(4);
  dEdx_mean_proton->SetLineWidth(2);
  dEdx_mean_proton->Draw("lsame");
  dEdx_mean_proton->Fit("fproton","QERM");

  c_test->cd(2);
  gPad->SetLogx();
  dEdx_muon->GetXaxis()->SetTitle("momentum ");
  dEdx_muon->GetYaxis()->SetTitle("#frac{dE}{dx}");
  dEdx_muon->Draw("colz");
  dEdx_mean_muon->SetLineColor(4);
  dEdx_mean_muon->SetLineWidth(2);
  dEdx_mean_muon->Draw("lsame");
  dEdx_mean_muon->Fit("fmuon","QERM");

  c_test->cd(3);
  gPad->SetLogx();
  dEdx_pion->GetXaxis()->SetTitle("momentum ");
  dEdx_pion->GetYaxis()->SetTitle("#frac{dE}{dx}");
  dEdx_pion->Draw("colz");
  dEdx_mean_pion->SetLineColor(4);
  dEdx_mean_pion->SetLineWidth(2);
  dEdx_mean_pion->Draw("lsame");
  dEdx_mean_pion->Fit("fpion","QERM");

  c_test->cd(4);
  gPad->SetLogx();
  dEdx_electron->GetXaxis()->SetTitle("momentum ");
  dEdx_electron->GetYaxis()->SetTitle("#frac{dE}{dx}");
  dEdx_electron->Draw("colz");
  dEdx_mean_electron->SetLineColor(4);
  dEdx_mean_electron->SetLineWidth(2);
  dEdx_mean_electron->Draw("lsame");
  dEdx_mean_electron->Fit("felectron","QERM");

  c_test->cd(5);
  gPad->SetLogx();
  dEdx_kaon->GetXaxis()->SetTitle("momentum ");
  dEdx_kaon->GetYaxis()->SetTitle("#frac{dE}{dx}");
  dEdx_kaon->Draw("colz");
  dEdx_mean_kaon->SetLineColor(4);
  dEdx_mean_kaon->SetLineWidth(2);
  dEdx_mean_kaon->Draw("lsame");
  dEdx_mean_kaon->Fit("fkaon","QERM");

  std::cout<<"#### Parameters for LikelihoodPIDProcessor "<<std::endl;
  std::cout<<" <!--dE/dx parameters for each particle-->"<<std::endl;
  std::cout<<"      <parameter name=\"dEdxParameter_electron\" type=\"std::vector< std::float >\"> "
	   << felectron->GetParameter(1)<<" "<<felectron->GetParameter(2)<<" "<<felectron->GetParameter(3)<<" "<<felectron->GetParameter(4)<<" "<<felectron->GetParameter(5)
	   <<" </parameter>"<<std::endl;
  std::cout<<"      <parameter name=\"dEdxParameter_muon\" type=\"std::vector< std::float >\"> "
	   << fmuon->GetParameter(1)<<" "<<fmuon->GetParameter(2)<<" "<<fmuon->GetParameter(3)<<" "<<fmuon->GetParameter(4)<<" "<<fmuon->GetParameter(5)
	   <<" </parameter>"<<std::endl;
  std::cout<<"      <parameter name=\"dEdxParameter_pion\" type=\"std::vector< std::float >\"> "
	   << fpion->GetParameter(1)<<" "<<fpion->GetParameter(2)<<" "<<fpion->GetParameter(3)<<" "<<fpion->GetParameter(4)<<" "<<fpion->GetParameter(5)
	   <<" </parameter>"<<std::endl;
  std::cout<<"      <parameter name=\"dEdxParameter_kaon\" type=\"std::vector< std::float >\"> "
	   << fkaon->GetParameter(1)<<" "<<fkaon->GetParameter(2)<<" "<<fkaon->GetParameter(3)<<" "<<fkaon->GetParameter(4)<<" "<<fkaon->GetParameter(5)
	   <<" </parameter>"<<std::endl;
  std::cout<<"      <parameter name=\"dEdxParameter_proton\" type=\"std::vector< std::float >\"> "
	   << fproton->GetParameter(1)<<" "<<fproton->GetParameter(2)<<" "<<fproton->GetParameter(3)<<" "<<fproton->GetParameter(4)<<" "<<fproton->GetParameter(5)
	   <<" </parameter>"<<std::endl;

}
