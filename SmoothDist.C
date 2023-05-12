
#include "getvalues.C"

std::map<int,double> SmoothDist(bool BeforeCuts = true, Double_t d0cut=100, TString analysis_opt="> d2d dmm anymass1L2M"){

  Int_t mass = 40;
  Int_t jalg = 2;
  TString lt = "m3p5";
  Long64_t drawN = (Long64_t) -1;
  Bool_t ScalePlots = true;
  
  TString xlabel = "M_{vis} (GeV/c^{2})";
  TString ylabel = "E_{miss} (GeV)";

  OBS_ID obs1 = omass_selection;
  OBS_ID obs2 = oemiss_selection;

  if (!BeforeCuts){
    obs1 = omass_after_dcut;
    obs2 = oemiss_after_dcut;
  }
  
  TCanvas *c1 = new TCanvas("c1","c1",0,0,1000,700); 

  std::pair<std::vector<Double_t>, Double_t> munuqq1 = getvalues(obs1, "munuqq", mass, "n/a", drawN, d0cut, jalg, analysis_opt);
  std::pair<std::vector<Double_t>, Double_t> Zmumu1 = getvalues(obs1, "Zmumu", mass, "n/a", drawN, d0cut, jalg, analysis_opt);
  std::pair<std::vector<Double_t>, Double_t> Ztautau1 = getvalues(obs1, "Ztautau", mass, "n/a", drawN, d0cut, jalg, analysis_opt);
  std::pair<std::vector<Double_t>, Double_t> Zuds1 = getvalues(obs1, "Zuds",mass, "n/a", drawN, d0cut, jalg, analysis_opt);
  std::pair<std::vector<Double_t>, Double_t> Zbb1 = getvalues(obs1, "Zbb",mass, "n/a", drawN, d0cut, jalg, analysis_opt);
  std::pair<std::vector<Double_t>, Double_t> Zcc1 = getvalues(obs1, "Zcc",mass, "n/a", drawN, d0cut, jalg, analysis_opt);

  std::pair<std::vector<Double_t>, Double_t> munuqq2 = getvalues(obs2, "munuqq", mass, "n/a", drawN, d0cut, jalg, analysis_opt);
  std::pair<std::vector<Double_t>, Double_t> Zmumu2 = getvalues(obs2, "Zmumu", mass, "n/a", drawN, d0cut, jalg, analysis_opt);
  std::pair<std::vector<Double_t>, Double_t> Ztautau2 = getvalues(obs2, "Ztautau", mass, "n/a", drawN, d0cut, jalg, analysis_opt);
  std::pair<std::vector<Double_t>, Double_t> Zuds2 = getvalues(obs2, "Zuds",mass, "n/a", drawN, d0cut, jalg, analysis_opt);
  std::pair<std::vector<Double_t>, Double_t> Zbb2 = getvalues(obs2, "Zbb",mass, "n/a", drawN, d0cut, jalg, analysis_opt);
  std::pair<std::vector<Double_t>, Double_t> Zcc2 = getvalues(obs2, "Zcc",mass, "n/a", drawN, d0cut, jalg, analysis_opt);

  

  TH2F *H2 = new TH2F("Granular2D","Granular2D",200,0,100,120,0,60);
  TH2F *H2S = new TH2F("Smoothed2D","Smoothed2D",200,0,100,120,0,60);

  for (int i=0; i < munuqq1.first.size(); i++){
    H2->Fill(munuqq1.first[i],munuqq2.first[i], Weight("munuqq","n/a","n/a",munuqq1.second));
  }

  for (int i=0; i < Ztautau1.first.size(); i++){
    H2->Fill(Ztautau1.first[i],Ztautau2.first[i], Weight("Ztautau","n/a","n/a",Ztautau1.second));
  }

  for (int i=0; i < Zmumu1.first.size(); i++){
    H2->Fill(Zmumu1.first[i],Zmumu2.first[i], Weight("Zmumu","n/a","n/a",Zmumu1.second));
  }

  for (int i=0; i < Zbb1.first.size(); i++){
    H2->Fill(Zbb1.first[i],Zbb2.first[i], Weight("Zbb","n/a","n/a",Zbb1.second));
  }

  for (int i=0; i < Zcc1.first.size(); i++){
    H2->Fill(Zcc1.first[i],Zcc2.first[i], Weight("Zcc","n/a","n/a",Zcc1.second));
  }

  for (int i=0; i < Zuds1.first.size(); i++){
    H2->Fill(Zuds1.first[i],Zuds2.first[i], Weight("Zuds","n/a","n/a",Zuds1.second));
  }
  

  
  
  for (int i=4; i<H2S->GetNbinsX()-4; i++){
    for (int j=0; j<H2S->GetNbinsY()-4; j++){

      //double content = H2->GetBinContent(i,j) / 49.;

      double content = 0.;
      for (int b1 = -3; b1<=3; b1++)
	for (int b2 = -3; b2<=3; b2++)
	  content += (H2->GetBinContent(i-b1, j-b2)/49. + H2->GetBinContent(i+b1, j+b2)/49.);


      H2S->SetBinContent(i,j,content);
	

    }
  }

  

  H2->GetXaxis()->SetTitle(xlabel);
  H2->GetYaxis()->SetTitle(ylabel);
  H2->Scale(0.012);
  H2->Draw("colz");
  gPad->SetLogz();
  gStyle->SetOptStat(0);

  
  cout<<"*** Integral "<<H2->Integral()<<endl;


  double zmass = 91.1786;

  


  Double_t scalefactor = 0.0102;

  

  std::map<int, double> toret;
  toret.clear();

  for (int imass : std::vector<int>{5,10,20,30,40,50,60,70,80,85}){
    Double_t precoil = (zmass*zmass - 1.* imass*imass ) / (2 * zmass);

    double counter = 0.;
    for (int i=0; i<H2->GetNbinsX(); i++){
      for (int j=0; j< H2->GetNbinsY(); j++){

	if (TMath::Abs( H2->GetXaxis()->GetBinCenter(i) - imass) <= 4. && TMath::Abs( H2->GetYaxis()->GetBinCenter(j) - precoil) <= 3.5 )
	  counter += H2->GetBinContent(i,j) * scalefactor;
      }
    }

    toret[imass] = counter;
  }

  return toret;
}

  
  
  
  /*

  TH1F *h_munuqq = new TH1F("munuqq",Form("#mu#nuqq ",munuqq.first.size()),nbin,bmin,bmax);
  if (drawstack) h_munuqq->SetFillColor(13);
  h_munuqq->SetLineColor(13);
 
  TH1F *h_zmumu = new TH1F("zmumu",Form("Z#rightarrow#mu#mu ",Zmumu.first.size()),nbin,bmin,bmax);
  if (drawstack) h_zmumu->SetFillColor(9);
  h_zmumu->SetLineColor(9);
 
  TH1F *h_ztautau = new TH1F("ztautau",Form("Z#rightarrow#tau#tau ",Ztautau.first.size()),nbin,bmin,bmax);
  if (drawstack) h_ztautau->SetFillColor(8);
  h_ztautau->SetLineColor(8);
  
  TH1F *h_zuds = new TH1F("zuds",Form("Z#rightarrowu/d/s ",Zuds.first.size()),nbin,bmin,bmax);
  if (drawstack) h_zuds->SetFillColor(6);
  h_zuds->SetLineColor(6);
  
  TH1F *h_zbb = new TH1F("zbb",Form("Z#rightarrowbb ",Zbb.first.size()),nbin,bmin,bmax);
  if (drawstack) h_zbb->SetFillColor(4);
  h_zbb->SetLineColor(4);
 
  TH1F *h_zcc = new TH1F("zcc",Form("Z#rightarrowcc ",Zcc.first.size()),nbin,bmin,bmax);
  if (drawstack) h_zcc->SetFillColor(2);
  h_zcc->SetLineColor(2);



  TH1F *h_signalM = new TH1F(Form("signal%d",mass),Form("signal %d ",mass,signalM.first.size()),nbin,bmin,bmax);
  h_signalM->SetLineWidth(3);
  h_signalM->SetLineColor(1);

  TH1F *h_signal30 = new TH1F(Form("signal%d",30),Form("signal %d ",30,signal30.first.size()),nbin,bmin,bmax);
  h_signal30->SetLineWidth(3);
  h_signal30->SetLineColor(30);
  h_signal30->SetLineStyle(9);

  TH1F *h_signal50 = new TH1F(Form("signal%d",50),Form("signal %d ",50,signal50.first.size()),nbin,bmin,bmax);
  h_signal50->SetLineWidth(3);
  h_signal50->SetLineColor(40);
  h_signal50->SetLineStyle(9);

  TH1F *h_signal70 = new TH1F(Form("signal%d",70),Form("signal %d ",70,signal70.first.size()),nbin,bmin,bmax);
  h_signal70->SetLineWidth(3);
  h_signal70->SetLineColor(44);
  h_signal70->SetLineStyle(9);

  TH1F *h_signal20 = new TH1F(Form("signal%d",20),Form("signal %d ",20,signal20.first.size()),nbin,bmin,bmax);
  h_signal20->SetLineWidth(3);
  h_signal20->SetLineColor(30);
  h_signal20->SetLineStyle(9);
  
  
  

  for (double ie : munuqq.first) h_munuqq->Fill(ie);
  cout<<"filled 1"<<endl;
  for (double ie : Zmumu.first) h_zmumu->Fill(ie);
  cout<<"filled 2"<<endl;
  for (double ie : Ztautau.first) h_ztautau->Fill(ie);
  cout<<"filled 3"<<endl;
  for (double ie : Zuds.first) h_zuds->Fill(ie);
  cout<<"filled 4"<<endl;
  for (double ie : Zbb.first) h_zbb->Fill(ie);
  cout<<"filled 5"<<endl;
  for (double ie : Zcc.first) h_zcc->Fill(ie);
  cout<<"filled 6"<<endl;

  for (double ie : signalM.first) h_signalM->Fill(ie);
  cout<<"filled 7: "<<signalM.first.size()<<" events"<<endl;


  for (double ie : signal30.first) h_signal30->Fill(ie);
  cout<<"filled 7: "<<signal30.first.size()<<" events"<<endl;


  for (double ie : signal50.first) h_signal50->Fill(ie);
  cout<<"filled 7: "<<signal50.first.size()<<" events"<<endl;


  for (double ie : signal70.first) h_signal70->Fill(ie);
  cout<<"filled 7: "<<signal70.first.size()<<" events"<<endl;

  for (double ie : signal20.first) h_signal20->Fill(ie);
  cout<<"filled 7: "<<signal20.first.size()<<" events"<<endl;


  
  if (ScalePlots){
    h_munuqq->Scale(Weight("munuqq","n/a","n/a",munuqq.second));
    h_zmumu->Scale(Weight("Zmumu","n/a","n/a",Zmumu.second));
    h_ztautau->Scale(Weight("Ztautau","n/a","n/a",Ztautau.second));
    h_zuds->Scale(Weight("Zuds","n/a","n/a",Zuds.second));
    h_zbb->Scale(Weight("Zbb","n/a","n/a",Zbb.second));
    h_zcc->Scale(Weight("Zcc","n/a","n/a",Zcc.second));
    h_signalM->Scale(Weight("signal",Form("%d",mass),lt,signalM.second));
    h_signal30->Scale(Weight("signal",Form("%d",30),"m3p0",signal30.second));
    h_signal50->Scale(Weight("signal",Form("%d",50),"m3p5",signal50.second));
    h_signal70->Scale(Weight("signal",Form("%d",50),"m4p5",signal70.second));
    h_signal20->Scale(Weight("signal",Form("%d",20),"m4p5",signal20.second));
  }
  

  auto hs = new THStack("hs",Form("%s;%s;%s","",xlabel.Data(),ScalePlots?"Events x Weight":"Events"));;
  hs->Add(h_munuqq);
  hs->Add(h_zmumu);
  hs->Add(h_ztautau);
  hs->Add(h_zuds);
  hs->Add(h_zbb);
  hs->Add(h_zcc);
 
  hs->Draw(drawstack ? "HIST" : "HIST nostack");

  hs->SetMinimum(1);

  
 
  //h_signalM->Draw("HIST same");
  h_signal20->Draw("HIST same");
  //h_signal50->Draw("HIST same");
  h_signal70->Draw("HIST same");

  gPad->SetLogy();

  auto leg = c1->BuildLegend();
  leg->SetMargin(0.25);
  leg->SetBorderSize(0);
  

  if (vline1 != -999){
    auto line1 = new TLine(vline1, 1, vline1, hs->GetMaximum());
    line1->SetLineWidth(2);
    line1->SetLineStyle(2);
    line1->SetLineColor(1);
    line1->Draw("same");
  }

  if (vline2 != -999){
    auto line2 = new TLine(vline2, 1, vline2, hs->GetMaximum());
    line2->SetLineWidth(2);
    line2->SetLineStyle(2);
    line2->SetLineColor(1);
    line2->Draw("same");
  }

  c1->SaveAs("temp.pdf");
  c1->SaveAs("temp.png");
  
 
}
  
  */
