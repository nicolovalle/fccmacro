#include "CutFlowOK.C"

Double_t GetUpp(Double_t n){

  if (n==0 || n==1) return 1.14;
  else return TMath::Sqrt(n);

}

  

std::vector<int> masses = {5, 10, 20, 30, 40, 50, 60, 70, 80, 85};


std::pair<TH2F*, TH2F*> TwoDnumbers(Int_t dd0cut = 0, TString sample = "signal", Bool_t Draw = true, TString AnalysisResPath = "../MyExternalAnalysis/results/skimmed/", Int_t RunOnN = -1,  Bool_t Scaled = false, Int_t jalg = 2, TString analysis_opt="> d2d dmm extra anymass1L2M window [2,0.1]"){
  // formulas: atals simple

  // opt: same as CutFlowOK.C
  
  // V[0][n] = n-th x coordinate of the curve
  // V[1][n] = n-th y coordinate of the curve



  std::vector<TString> mps = {"m","p"};

  
  TH2F* H = new TH2F("Hfilled","Signal. Event selection + r_{vert}>0.5  efficiency.",19,0-2.5,90+2.5,64,-12,-4);
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


	  /* MANUAL IMPLEMENTATION OF D0 CUT ABOVE 70 GEV
          if (m>=70){
	    MapSample = CutFlowOK(sample,m,lt,AnalysisResPath,RunOnN,false,jalg, "< d2d dsigma extra anymass1L2M window [2,0.1]"); myid = dcut_id(8);
	  }
	  */

	  bin_content = 1.*MapSample[m][myid] / 10000.;

	 

	  

	  //MapSample[m][0] is Nnocut * SF;
	  if (Scaled) bin_content *= xsec(sample, Form("%d",m), lt) * LUMI / MapSample[m][0]; 

	  Double_t U = Coupling(Form("%d",m),lt);
	  Double_t LogU2 = 2.*TMath::Log10(U);

	  cout<<"TwoDnumbers.C:: Bin_content "<<sample<<" ("<<m<<","<<LogU2<<") = "<<bin_content<<endl;

	  if (m==40 && LogU2 > -7.5) bin_content=0; 
	  
	  if (bin_content == 0) Hblank->Fill(1.*m, LogU2, 1e8);
	  else H->Fill(1.*m, LogU2, bin_content);
	  
	 
	}
      }
    }
  }

  if (file_counter==0) cout<<"TwoDnumbers.C:: NO FILES FOUND!"<<endl;

  if (Draw){

    auto c = new TCanvas();

    H->GetZaxis()->SetRangeUser(1e-6,0.65); //efficiency
    //H->GetZaxis()->SetRangeUser(1e-4,1e5); //weights
    //H->GetZaxis()->SetRangeUser(1e-10,1e-1); // xsec
    //H->GetZaxis()->SetRangeUser(1e-1,1e5); // number of events
    H->Draw("colz");
    //H->SetMarkerSize(1.6);
    //H->SetTitle("Significance");
    Hblank->SetLineColor(1);
    Hblank->Draw("same BOX");
    //gPad->SetLogz();

    

    gStyle->SetOptStat(0);
    //gPad->SetLogz();
    //gStyle->SetPalette(kBlackBody);
    gStyle->SetPalette(kBlueRedYellow);

    c->SaveAs("temp.png");
    c->SaveAs("temp.pdf");
  }

  

  return std::pair<TH2F*, TH2F*>{H, Hblank};

}

//void SliceSignal


void DrawSamples(Int_t dd0cut = 0, TString AnalysisResPath = "../MyExternalAnalysis/results/skimmed/", Int_t RunOnN = -1, Bool_t Scaled = false, Int_t jalg = 2, TString analysis_opt = "> d2d dmm extra anymass1L2M window [2,0.1]"){

  
  TCanvas *c = new TCanvas("c","c",2500,2200);

  c->Divide(3,3);

  TH2F* Hsig, *Hsig_;
  TH2F* HZbb, *HZbb_;
  TH2F* HZcc, *HZcc_;
  TH2F* HZuds, *HZuds_;
  TH2F* HZud, *HZud_;
  TH2F* HZss, *HZss_;
  TH2F* HZmumu, *HZmumu_;
  TH2F* HZtautau, *HZtautau_;
  TH2F* Hmunuqq, *Hmunuqq_;

  Int_t textsize = 3;

  // TH1F *hZbb, *hZcc, *hZuds, *hZmumu, *hZtautau, *hmunuqq;

  std::pair<TH2F*, TH2F*> Pair;

  std::vector<double> Values;

  Int_t panelcounter = 1;

  
  c->cd(panelcounter); panelcounter++;
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
  


  

  c->cd(panelcounter); panelcounter++;
  sample = "Zbb";

  Pair = TwoDnumbers(dd0cut, sample, false, AnalysisResPath, RunOnN, Scaled, jalg, analysis_opt);
  HZbb = (TH2F*)Pair.first->Clone(Form("H%s (weight: %g)",sample.Data(),Weight(sample)));
  HZbb_ = (TH2F*)Pair.second->Clone("HZbb_");
  HZbb->SetTitle(sample);
  HZbb->GetZaxis()->SetRangeUser(1e-1,1e5);
  //HZbb->Draw("colz");
  //HZbb_->Draw("same BOX");
  //gStyle->SetOptStat(0);
  //gPad->SetLogz();
  //gStyle->SetPalette(kColorPrintableOnGrey);
  TH1F* hZbb = new TH1F("hZbb",Form("H%s (weight: %g)",sample.Data(),Weight(sample)),masses.size(),0,masses.size());
  for (int i=0; i<masses.size(); i++) {
    hZbb->GetXaxis()->SetBinLabel(i+1,Form("%d",masses.at(i)));
    double content = 0.;
    for (int j = 1; j<HZbb->GetYaxis()->GetNbins(); j++) content = TMath::Max(content, HZbb->GetBinContent( HZbb->GetXaxis()->FindBin(masses.at(i)) , j));
    hZbb->SetBinContent(i+1,content);
  }
  hZbb->GetXaxis()->SetTitle("M_{HN} mass bin");
  hZbb->SetMarkerSize(textsize);
  hZbb->Draw("hist text");
  hZbb->SetLineColor(1);
  gStyle->SetOptStat(0);
  gPad->SetLogy();
  

  
  c->cd(panelcounter); panelcounter++;
  sample = "Zcc";
  
  Pair = TwoDnumbers(dd0cut, sample, false, AnalysisResPath, RunOnN, Scaled, jalg, analysis_opt);
  HZcc = (TH2F*)Pair.first->Clone(Form("H%s (weight: %g)",sample.Data(),Weight(sample)));
  HZcc_ = (TH2F*)Pair.second->Clone("HZcc_");
  HZcc->SetTitle(sample);
  HZcc->GetZaxis()->SetRangeUser(1e-1,1e5);
  //HZcc->Draw("colz");
  //HZcc_->Draw("same BOX");
  //gStyle->SetOptStat(0);
  //gPad->SetLogz();
  //gStyle->SetPalette(kColorPrintableOnGrey);
  TH1F* hZcc = new TH1F("hZcc",Form("H%s (weight: %g)",sample.Data(),Weight(sample)),masses.size(),0,masses.size());
  for (int i=0; i<masses.size(); i++) {
    hZcc->GetXaxis()->SetBinLabel(i+1,Form("%d",masses.at(i)));
    double content = 0.;
    for (int j = 1; j<HZcc->GetYaxis()->GetNbins(); j++) content = TMath::Max(content, HZcc->GetBinContent( HZcc->GetXaxis()->FindBin(masses.at(i)) , j));
    hZcc->SetBinContent(i+1,content);
  }
  hZcc->GetXaxis()->SetTitle("M_{HN} mass bin");
  hZcc->SetMarkerSize(textsize);
  hZcc->Draw("hist text");
  hZcc->SetLineColor(1);
  gStyle->SetOptStat(0);
  gPad->SetLogy();


  if (PRODUCTION == "Spring2021"){
  	c->cd(panelcounter); panelcounter++;
  	sample = "Zuds";
  	
  	Pair = TwoDnumbers(dd0cut, sample, false, AnalysisResPath, RunOnN, Scaled, jalg, analysis_opt);
  	HZuds = (TH2F*)Pair.first->Clone(Form("H%s (weight: %g)",sample.Data(),Weight(sample)));
  	HZuds_ = (TH2F*)Pair.second->Clone("HZuds_");
  	HZuds->SetTitle(sample);
  	HZuds->GetZaxis()->SetRangeUser(1e-1,1e5);
  	//HZuds->Draw("colz");
  	//HZuds_->Draw("same BOX");
  	//gStyle->SetOptStat(0);
  	//gPad->SetLogz();
  	//gStyle->SetPalette(kColorPrintableOnGrey);
  	TH1F* hZuds = new TH1F("hZuds",Form("H%s (weight: %g)",sample.Data(),Weight(sample)),masses.size(),0,masses.size());
  	for (int i=0; i<masses.size(); i++) {
  	  hZuds->GetXaxis()->SetBinLabel(i+1,Form("%d",masses.at(i)));
  	  double content = 0.;
  	  for (int j = 1; j<HZuds->GetYaxis()->GetNbins(); j++) content = TMath::Max(content, HZuds->GetBinContent( HZuds->GetXaxis()->FindBin(masses.at(i)) , j));
  	  hZuds->SetBinContent(i+1,content);
  	}
        hZuds->GetXaxis()->SetTitle("M_{HN} mass bin");
        hZuds->SetMarkerSize(textsize);
  	hZuds->Draw("hist text");
  	hZuds->SetLineColor(1);
  	gStyle->SetOptStat(0);
  	gPad->SetLogy();
  }

  if (PRODUCTION == "Winter2023"){
        c->cd(panelcounter); panelcounter++;
  	sample = "Zud";
  	
  	Pair = TwoDnumbers(dd0cut, sample, false, AnalysisResPath, RunOnN, Scaled, jalg, analysis_opt);
  	HZud = (TH2F*)Pair.first->Clone(Form("H%s (weight: %g)",sample.Data(),Weight(sample)));
  	HZud_ = (TH2F*)Pair.second->Clone("HZud_");
  	HZud->SetTitle(sample);
  	HZud->GetZaxis()->SetRangeUser(1e-1,1e5);
  	//HZud->Draw("colz");
  	//HZud_->Draw("same BOX");
  	//gStyle->SetOptStat(0);
  	//gPad->SetLogz();
  	//gStyle->SetPalette(kColorPrintableOnGrey);
  	TH1F* hZud = new TH1F("hZud",Form("H%s (weight: %g)",sample.Data(),Weight(sample)),masses.size(),0,masses.size());
  	for (int i=0; i<masses.size(); i++) {
  	  hZud->GetXaxis()->SetBinLabel(i+1,Form("%d",masses.at(i)));
  	  double content = 0.;
  	  for (int j = 1; j<HZud->GetYaxis()->GetNbins(); j++) content = TMath::Max(content, HZud->GetBinContent( HZud->GetXaxis()->FindBin(masses.at(i)) , j));
  	  hZud->SetBinContent(i+1,content);
  	}
        hZud->GetXaxis()->SetTitle("M_{HN} mass bin");
        hZud->SetMarkerSize(textsize);
  	hZud->Draw("hist text");
  	hZud->SetLineColor(1);
  	gStyle->SetOptStat(0);
  	gPad->SetLogy();

	c->cd(panelcounter); panelcounter++;
  	sample = "Zss";
  	
  	Pair = TwoDnumbers(dd0cut, sample, false, AnalysisResPath, RunOnN, Scaled, jalg, analysis_opt);
  	HZss = (TH2F*)Pair.first->Clone(Form("H%s (weight: %g)",sample.Data(),Weight(sample)));
  	HZss_ = (TH2F*)Pair.second->Clone("HZss_");
  	HZss->SetTitle(sample);
  	HZss->GetZaxis()->SetRangeUser(1e-1,1e5);
  	//HZss->Draw("colz");
  	//HZss_->Draw("same BOX");
  	//gStyle->SetOptStat(0);
  	//gPad->SetLogz();
  	//gStyle->SetPalette(kColorPrintableOnGrey);
  	TH1F* hZss = new TH1F("hZss",Form("H%s (weight: %g)",sample.Data(),Weight(sample)),masses.size(),0,masses.size());
  	for (int i=0; i<masses.size(); i++) {
  	  hZss->GetXaxis()->SetBinLabel(i+1,Form("%d",masses.at(i)));
  	  double content = 0.;
  	  for (int j = 1; j<HZss->GetYaxis()->GetNbins(); j++) content = TMath::Max(content, HZss->GetBinContent( HZss->GetXaxis()->FindBin(masses.at(i)) , j));
  	  hZss->SetBinContent(i+1,content);
  	}
        hZss->GetXaxis()->SetTitle("M_{HN} mass bin");
        hZss->SetMarkerSize(textsize);
  	hZss->Draw("hist text");
  	hZss->SetLineColor(1);
  	gStyle->SetOptStat(0);
  	gPad->SetLogy();
    
  }


  c->cd(panelcounter); panelcounter++;
  sample = "Ztautau";
  
  Pair = TwoDnumbers(dd0cut, sample, false, AnalysisResPath, RunOnN, Scaled, jalg, analysis_opt);
  HZtautau = (TH2F*)Pair.first->Clone(Form("H%s (weight: %g)",sample.Data(),Weight(sample)));
  HZtautau_ = (TH2F*)Pair.second->Clone("HZtautau_");
  HZtautau->SetTitle(sample);
  HZtautau->GetZaxis()->SetRangeUser(1e-1,1e5);
  //HZtautau->Draw("colz");
  //HZtautau_->Draw("same BOX");
  //gStyle->SetOptStat(0);
  //gPad->SetLogz();
  //gStyle->SetPalette(kColorPrintableOnGrey);
  TH1F* hZtautau = new TH1F("hZtautau",Form("H%s (weight: %g)",sample.Data(),Weight(sample)),masses.size(),0,masses.size());
  for (int i=0; i<masses.size(); i++) {
    hZtautau->GetXaxis()->SetBinLabel(i+1,Form("%d",masses.at(i)));
    double content = 0.;
    for (int j = 1; j<HZtautau->GetYaxis()->GetNbins(); j++) content = TMath::Max(content, HZtautau->GetBinContent( HZtautau->GetXaxis()->FindBin(masses.at(i)) , j));
    hZtautau->SetBinContent(i+1,content);
  }
  hZtautau->GetXaxis()->SetTitle("M_{HN} mass bin");
  hZtautau->SetMarkerSize(textsize);
  hZtautau->Draw("hist text");
  hZtautau->SetLineColor(1);
  gStyle->SetOptStat(0);
  gPad->SetLogy();



  c->cd(panelcounter); panelcounter++;
  sample = "Zmumu";
  
  Pair = TwoDnumbers(dd0cut, sample, false, AnalysisResPath, RunOnN, Scaled, jalg, analysis_opt);
  HZmumu = (TH2F*)Pair.first->Clone(Form("H%s (weight: %g)",sample.Data(),Weight(sample)));
  HZmumu_ = (TH2F*)Pair.second->Clone("HZmumu_");
  HZmumu->SetTitle(sample);
  HZmumu->GetZaxis()->SetRangeUser(1e-1,1e5);
  //HZmumu->Draw("colz");
  //HZmumu_->Draw("same BOX");
  //gStyle->SetOptStat(0);
  //gPad->SetLogz();
  //gStyle->SetPalette(kColorPrintableOnGrey);
  TH1F* hZmumu = new TH1F("hZmumu",Form("H%s (weight: %g)",sample.Data(),Weight(sample)),masses.size(),0,masses.size());
  for (int i=0; i<masses.size(); i++) {
    hZmumu->GetXaxis()->SetBinLabel(i+1,Form("%d",masses.at(i)));
    double content = 0.;
    for (int j = 1; j<HZmumu->GetYaxis()->GetNbins(); j++) content = TMath::Max(content, HZmumu->GetBinContent( HZmumu->GetXaxis()->FindBin(masses.at(i)) , j));
    hZmumu->SetBinContent(i+1,content);
  }
  hZmumu->GetXaxis()->SetTitle("M_{HN} mass bin");
  hZmumu->SetMarkerSize(textsize);
  hZmumu->Draw("hist text");
  hZmumu->SetLineColor(1);
  gStyle->SetOptStat(0);
  gPad->SetLogy();
  
  

  c->cd(panelcounter); panelcounter++;
  sample = "munuqq";
  
  Pair = TwoDnumbers(dd0cut, sample, false, AnalysisResPath, RunOnN, Scaled, jalg, analysis_opt);
  Hmunuqq = (TH2F*)Pair.first->Clone(Form("H%s (weight: %g)",sample.Data(),Weight(sample)));
  Hmunuqq_ = (TH2F*)Pair.second->Clone("Hmunuqq_");
  Hmunuqq->SetTitle(sample);
  Hmunuqq->GetZaxis()->SetRangeUser(1,1e5);
  //Hmunuqq->Draw("colz");
  //Hmunuqq_->Draw("same BOX");
  //gStyle->SetOptStat(0);
  //gPad->SetLogz();
  //gStyle->SetPalette(kColorPrintableOnGrey);
  TH1F* hmunuqq = new TH1F("hmunuqq",Form("H%s (weight: %g)",sample.Data(),Weight(sample)),masses.size(),0,masses.size());
  for (int i=0; i<masses.size(); i++) {
    hmunuqq->GetXaxis()->SetBinLabel(i+1,Form("%d",masses.at(i)));
    double content = 0.;
    for (int j = 1; j<Hmunuqq->GetYaxis()->GetNbins(); j++) content = TMath::Max(content, Hmunuqq->GetBinContent( Hmunuqq->GetXaxis()->FindBin(masses.at(i)) , j));
    hmunuqq->SetBinContent(i+1,content);
  }
  hmunuqq->GetYaxis()->SetRangeUser(0.5,1.e5);
  hmunuqq->GetXaxis()->SetTitle("M_{HN} mass bin");
  hmunuqq->SetMarkerSize(textsize);
  hmunuqq->Draw("hist text");
  hmunuqq->SetLineColor(1);
  gStyle->SetOptStat(0);
  gPad->SetLogy();
  
  



  c->SaveAs("temp.png");
  c->SaveAs("temp.pdf");


}
