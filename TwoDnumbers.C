#include "CutFlowOK.C"

Double_t GetUpp(Double_t n){

  if (n==0 || n==1) return 1.14;
  else return TMath::Sqrt(n);

}

  

std::vector<int> masses = {5, 10, 20, 30, 40, 50, 60, 70, 80, 85};


std::pair<TH2F*, TH2F*> TwoDnumbers(Int_t dd0cut = 8, TString sample = "signal", Bool_t Draw = true, TString AnalysisResPath = "../MyExternalAnalysis/results/skimmed/", Int_t RunOnN = -1,  Bool_t Scaled = false, Int_t jalg = 2, TString analysis_opt="> d2d dmm anymass1L2M window [2,0.2]"){
  // formulas: atals simple

  // opt: same as CutFlowOK.C
  
  // V[0][n] = n-th x coordinate of the curve
  // V[1][n] = n-th y coordinate of the curve


  
  std::vector<TString> mps = {"m","p"};

  
  TH2F* H = new TH2F("Hfilled","Signal. D_{0,#mu} > 1 mm",19,0-2.5,90+2.5,64,-12,-4);
  TH2F* Hblank = new TH2F("Hblank","Hblank",19,0-2.5,90+2.5,64,-12,-4);

  H->GetXaxis()->SetTitle("M(HN) (GeV)");
  H->GetYaxis()->SetTitle("Log(U^{2})");

  

  std::map<int, std::vector<double>> MapSample;
  if (sample != "signal") MapSample =  CutFlowOK(sample,-1,"n/a",AnalysisResPath,RunOnN,false,jalg,analysis_opt);
 

  Int_t file_counter = 0;
  for (int m : masses){
    for (TString mp : mps){
      for (int unit = 0; unit < 9; unit++){
	for (int decimal = 0; decimal < 6; decimal += 5){


	  TString lt = Form("%s%dp%d", mp.Data(), unit, decimal);
	  

	  TString FileNameToCheck = Form("%s%s", AnalysisResPath.Data(), AnalysisResults("signal",Form("%d",m),lt).Data());
      

	  if (gSystem->AccessPathName(FileNameToCheck)) continue;
	  file_counter++;
     

	  Int_t myid = dcut_id(dd0cut);

	  Double_t bin_content = -1.;

	  if (sample == "signal") MapSample = CutFlowOK("signal",m,lt,AnalysisResPath,RunOnN,false,jalg,analysis_opt);

	  bin_content = 1.*MapSample[m][myid];

	  

	  //MapSample[m][0] is Nnocut * SF;
	  if (Scaled) bin_content *= xsec(sample, Form("%d",m), lt) * LUMI / MapSample[m][0]; 

	  Double_t U = Coupling(Form("%d",m),lt);
	  Double_t LogU2 = 2.*TMath::Log10(U);

	  cout<<"TwoDnumbers.C:: Bin_content "<<sample<<" ("<<m<<","<<LogU2<<") = "<<bin_content<<endl;

	 
	  
	  if (bin_content == 0) Hblank->Fill(1.*m, LogU2, 1e8);
	  else H->Fill(1.*m, LogU2, bin_content);
	  
	 
	}
      }
    }
  }

  if (file_counter==0) cout<<"TwoDnumbers.C:: NO FILES FOUND!"<<endl;

  if (Draw){

    auto c = new TCanvas();

    //H->GetZaxis()->SetRangeUser(1e-6,1); //efficiency
    //H->GetZaxis()->SetRangeUser(1e-4,1e5); //weights
    //H->GetZaxis()->SetRangeUser(1e-10,1e-1); // xsec
    //H->GetZaxis()->SetRangeUser(1e-1,1e5); // number of events
    H->Draw("colz");
    //H->SetMarkerSize(1.6);
    //H->SetTitle("Significance");
    Hblank->SetLineColor(1);
    Hblank->Draw("same BOX");

    

    gStyle->SetOptStat(0);
    //gPad->SetLogz();
    //gStyle->SetPalette(kBlackBody);
    gStyle->SetPalette(kBlueRedYellow);

    c->SaveAs("temp.png");
    c->SaveAs("temp.pdf");
  }

  

  return std::pair<TH2F*, TH2F*>{H, Hblank};

}


void DrawSamples(Int_t dd0cut = 8, TString AnalysisResPath = "../MyExternalAnalysis/results/skimmed/", Int_t RunOnN = -1, Bool_t Scaled = false, Int_t jalg = 2, TString analysis_opt = "< d2d dsigma anymass1L2M"){


  TCanvas *c = new TCanvas("c","c",2500,2200);

  c->Divide(3,3);

  TH2F* Hsig, *Hsig_;
  TH2F* HZbb, *HZbb_;
  TH2F* HZcc, *HZcc_;
  TH2F* HZuds, *HZuds_;
  TH2F* HZmumu, *HZmumu_;
  TH2F* HZtautau, *HZtautau_;
  TH2F* Hmunuqq, *Hmunuqq_;

  // TH1F *hZbb, *hZcc, *hZuds, *hZmumu, *hZtautau, *hmunuqq;

  std::pair<TH2F*, TH2F*> Pair;

  std::vector<double> Values;


  
  c->cd(1);
  TString sample = "signal";
  
  Pair = TwoDnumbers(dd0cut, sample, false, AnalysisResPath, RunOnN, Scaled, jalg, analysis_opt);
  Hsig = (TH2F*)Pair.first->Clone("Hsig");
  Hsig_ = (TH2F*)Pair.second->Clone("Hsig_");
  Hsig->SetTitle(sample);
  Hsig->GetZaxis()->SetRangeUser(1e-1,1e5);
  Hsig->Draw("colz");
  Hsig_->Draw("same BOX");
  gStyle->SetOptStat(0);
  gPad->SetLogz();
  gStyle->SetPalette(kColorPrintableOnGrey);
  


  

  c->cd(2);
  sample = "Zbb";

  Pair = TwoDnumbers(dd0cut, sample, false, AnalysisResPath, RunOnN, Scaled, jalg, analysis_opt);
  HZbb = (TH2F*)Pair.first->Clone("HZbb");
  HZbb_ = (TH2F*)Pair.second->Clone("HZbb_");
  HZbb->SetTitle(sample);
  HZbb->GetZaxis()->SetRangeUser(1e-1,1e5);
  //HZbb->Draw("colz");
  //HZbb_->Draw("same BOX");
  //gStyle->SetOptStat(0);
  //gPad->SetLogz();
  //gStyle->SetPalette(kColorPrintableOnGrey);
  TH1F* hZbb = new TH1F("hZbb","hZbb",masses.size(),0,masses.size());
  for (int i=0; i<masses.size(); i++) {
    hZbb->GetXaxis()->SetBinLabel(i+1,Form("%d",masses.at(i)));
    double content = 0.;
    for (int j = 1; j<HZbb->GetYaxis()->GetNbins(); j++) content = TMath::Max(content, HZbb->GetBinContent( HZbb->GetXaxis()->FindBin(masses.at(i)) , j));
    hZbb->SetBinContent(i+1,content);
  }
  hZbb->Draw("hist text");
  hZbb->SetLineColor(1);
  gStyle->SetOptStat(0);
  gPad->SetLogy();
  

  
  c->cd(3);
  sample = "Zcc";
  
  Pair = TwoDnumbers(dd0cut, sample, false, AnalysisResPath, RunOnN, Scaled, jalg, analysis_opt);
  HZcc = (TH2F*)Pair.first->Clone("HZcc");
  HZcc_ = (TH2F*)Pair.second->Clone("HZcc_");
  HZcc->SetTitle(sample);
  HZcc->GetZaxis()->SetRangeUser(1e-1,1e5);
  //HZcc->Draw("colz");
  //HZcc_->Draw("same BOX");
  //gStyle->SetOptStat(0);
  //gPad->SetLogz();
  //gStyle->SetPalette(kColorPrintableOnGrey);
  TH1F* hZcc = new TH1F("hZcc","hZcc",masses.size(),0,masses.size());
  for (int i=0; i<masses.size(); i++) {
    hZcc->GetXaxis()->SetBinLabel(i+1,Form("%d",masses.at(i)));
    double content = 0.;
    for (int j = 1; j<HZcc->GetYaxis()->GetNbins(); j++) content = TMath::Max(content, HZcc->GetBinContent( HZcc->GetXaxis()->FindBin(masses.at(i)) , j));
    hZcc->SetBinContent(i+1,content);
  }
  hZcc->Draw("hist text");
  hZcc->SetLineColor(1);
  gStyle->SetOptStat(0);
  gPad->SetLogy();


  c->cd(4);
  sample = "Zuds";
  
  Pair = TwoDnumbers(dd0cut, sample, false, AnalysisResPath, RunOnN, Scaled, jalg, analysis_opt);
  HZuds = (TH2F*)Pair.first->Clone("HZuds");
  HZuds_ = (TH2F*)Pair.second->Clone("HZuds_");
  HZuds->SetTitle(sample);
  HZuds->GetZaxis()->SetRangeUser(1e-1,1e5);
  //HZuds->Draw("colz");
  //HZuds_->Draw("same BOX");
  //gStyle->SetOptStat(0);
  //gPad->SetLogz();
  //gStyle->SetPalette(kColorPrintableOnGrey);
  TH1F* hZuds = new TH1F("hZuds","hZuds",masses.size(),0,masses.size());
  for (int i=0; i<masses.size(); i++) {
    hZuds->GetXaxis()->SetBinLabel(i+1,Form("%d",masses.at(i)));
    double content = 0.;
    for (int j = 1; j<HZuds->GetYaxis()->GetNbins(); j++) content = TMath::Max(content, HZuds->GetBinContent( HZuds->GetXaxis()->FindBin(masses.at(i)) , j));
    hZuds->SetBinContent(i+1,content);
  }
  hZuds->Draw("hist text");
  hZuds->SetLineColor(1);
  gStyle->SetOptStat(0);
  gPad->SetLogy();


  c->cd(5);
  sample = "Ztautau";
  
  Pair = TwoDnumbers(dd0cut, sample, false, AnalysisResPath, RunOnN, Scaled, jalg, analysis_opt);
  HZtautau = (TH2F*)Pair.first->Clone("HZtautau");
  HZtautau_ = (TH2F*)Pair.second->Clone("HZtautau_");
  HZtautau->SetTitle(sample);
  HZtautau->GetZaxis()->SetRangeUser(1e-1,1e5);
  //HZtautau->Draw("colz");
  //HZtautau_->Draw("same BOX");
  //gStyle->SetOptStat(0);
  //gPad->SetLogz();
  //gStyle->SetPalette(kColorPrintableOnGrey);
  TH1F* hZtautau = new TH1F("hZtautau","hZtautau",masses.size(),0,masses.size());
  for (int i=0; i<masses.size(); i++) {
    hZtautau->GetXaxis()->SetBinLabel(i+1,Form("%d",masses.at(i)));
    double content = 0.;
    for (int j = 1; j<HZtautau->GetYaxis()->GetNbins(); j++) content = TMath::Max(content, HZtautau->GetBinContent( HZtautau->GetXaxis()->FindBin(masses.at(i)) , j));
    hZtautau->SetBinContent(i+1,content);
  }
  hZtautau->Draw("hist text");
  hZtautau->SetLineColor(1);
  gStyle->SetOptStat(0);
  gPad->SetLogy();



  c->cd(6);
  sample = "Zmumu";
  
  Pair = TwoDnumbers(dd0cut, sample, false, AnalysisResPath, RunOnN, Scaled, jalg, analysis_opt);
  HZmumu = (TH2F*)Pair.first->Clone("HZmumu");
  HZmumu_ = (TH2F*)Pair.second->Clone("HZmumu_");
  HZmumu->SetTitle(sample);
  HZmumu->GetZaxis()->SetRangeUser(1e-1,1e5);
  //HZmumu->Draw("colz");
  //HZmumu_->Draw("same BOX");
  //gStyle->SetOptStat(0);
  //gPad->SetLogz();
  //gStyle->SetPalette(kColorPrintableOnGrey);
  TH1F* hZmumu = new TH1F("hZmumu","hZmumu",masses.size(),0,masses.size());
  for (int i=0; i<masses.size(); i++) {
    hZmumu->GetXaxis()->SetBinLabel(i+1,Form("%d",masses.at(i)));
    double content = 0.;
    for (int j = 1; j<HZmumu->GetYaxis()->GetNbins(); j++) content = TMath::Max(content, HZmumu->GetBinContent( HZmumu->GetXaxis()->FindBin(masses.at(i)) , j));
    hZmumu->SetBinContent(i+1,content);
  }
  hZmumu->Draw("hist text");
  hZmumu->SetLineColor(1);
  gStyle->SetOptStat(0);
  gPad->SetLogy();
  
  

  c->cd(7);
  sample = "munuqq";
  
  Pair = TwoDnumbers(dd0cut, sample, false, AnalysisResPath, RunOnN, Scaled, jalg, analysis_opt);
  Hmunuqq = (TH2F*)Pair.first->Clone("Hmunuqq");
  Hmunuqq_ = (TH2F*)Pair.second->Clone("Hmunuqq_");
  Hmunuqq->SetTitle(sample);
  Hmunuqq->GetZaxis()->SetRangeUser(1,1e5);
  //Hmunuqq->Draw("colz");
  //Hmunuqq_->Draw("same BOX");
  //gStyle->SetOptStat(0);
  //gPad->SetLogz();
  //gStyle->SetPalette(kColorPrintableOnGrey);
  TH1F* hmunuqq = new TH1F("hmunuqq","hmunuqq",masses.size(),0,masses.size());
  for (int i=0; i<masses.size(); i++) {
    hmunuqq->GetXaxis()->SetBinLabel(i+1,Form("%d",masses.at(i)));
    double content = 0.;
    for (int j = 1; j<Hmunuqq->GetYaxis()->GetNbins(); j++) content = TMath::Max(content, Hmunuqq->GetBinContent( Hmunuqq->GetXaxis()->FindBin(masses.at(i)) , j));
    hmunuqq->SetBinContent(i+1,content);
  }
  hmunuqq->Draw("hist text");
  hmunuqq->SetLineColor(1);
  gStyle->SetOptStat(0);
  gPad->SetLogy();
  
  



  c->SaveAs("temp.png");
  c->SaveAs("temp.pdf");


}
