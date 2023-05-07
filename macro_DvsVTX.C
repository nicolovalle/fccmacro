
#include "getvalues.C"

void macro_DvsVTX(){

  TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600); 

  Int_t sigmass = 70;
  TString siglt = "m5p0";

  OBS_ID obs = oVtxXY; Double_t lowbin = 1e-2;  Double_t upbin = 10; 

  //OBS_ID obs = od0; Double_t lowbin = 1e-1;  Double_t upbin = 100;

  
  

  std::pair<std::vector<Double_t>, Double_t> Dsignal = getvalues(obs, "signal", sigmass, siglt, -1);
  std::pair<std::vector<Double_t>, Double_t> DZbb = getvalues(obs, "Zbb", 80, "n/a", 1.e5);
  std::pair<std::vector<Double_t>, Double_t> DZcc = getvalues(obs, "Zcc", 80, "n/a", 7.e4); // proper ratio according to weights.

  TH1F *h_Dsignal = new TH1F("h_Dsignal",Form("Signal %d GeV / Z #rightarrow bb,cc",sigmass),1000,lowbin,upbin);
  TH1F *h_Dbkg = new TH1F("h_Dbkg","Z #rightarrow bb/cc",1000,lowbin,upbin);


  for (double ie : Dsignal.first) h_Dsignal->Fill(ie);
  for (double ie : DZbb.first) h_Dbkg->Fill(ie);
  for (double ie : DZcc.first) h_Dbkg->Fill(ie);


  h_Dsignal->Scale(Weight("signal",Form("%d",sigmass),siglt));
  h_Dbkg->Scale(Weight("Zbb","n/a","n/a",DZbb.second));

  TH1F *C_Dsignal = (TH1F*)h_Dsignal->GetCumulative();
  TH1F *C_Dbkg = (TH1F*)h_Dbkg->GetCumulative();

  C_Dsignal->SetLineColor(2); C_Dsignal->SetLineWidth(2);
  C_Dbkg->SetLineColor(4); C_Dbkg->SetLineWidth(2);

  C_Dsignal->Divide(C_Dbkg);
  
  C_Dsignal->Draw();
  //C_Dbkg->Draw("same");

  C_Dsignal->GetYaxis()->SetRangeUser(1e-9,1e-7);

  gPad->SetLogy();
  gPad->SetLogx();
  gStyle->SetOptStat(0);

  //auto leg = c1->BuildLegend();
  //leg->SetMargin(0.25);
  //leg->SetBorderSize(0);
  
  
  /*

  Int_t mass = 40;
  Int_t jalg = 2;
  TString lt = "m3p5";
  Long64_t drawN = (Long64_t) 1.e5;
  Bool_t ScalePlots = true;
  TString xlabel;
  if (obsID == oVtxXY) xlabel = "sqrt(Vtx_{x}^{2}+Vtx_{y}^{2}) (mm)";
  if (obsID == od0 || obsID == od0sel) xlabel = "D_{0,#mu}/#sigma_{D_{0}}";
  if (obsID == ocosjj) xlabel = "cos(j,j)";
  if (obsID == ocospmissmu) xlabel = "cos(#mu,p_{miss})";
  if (obsID == oMAXcosjmu) xlabel = "cos(j,#mu)";
  if (obsID == oMINcosjmu) xlabel = "min cos(j,#mu)";
  if (obsID == oMINEjet) xlabel = "min E(j) (GeV)";
  if (obsID == omtot) xlabel = "M_{tot = vis+miss} (GeV/c^{2})";
  if (obsID == omass_selection) xlabel = "M_{vis} (GeV/c^{2})";
  if (obsID == oemiss_selection) xlabel = "E_{miss} (GeV)";
  
  
  TCanvas *c1 = new TCanvas("c1","c1",0,0,1000,700); 

  std::pair<std::vector<Double_t>, Double_t> munuqq = getvalues(obsID, "munuqq", mass, "n/a", drawN, 8, jalg);
  std::pair<std::vector<Double_t>, Double_t> Zmumu = getvalues(obsID, "Zmumu", mass, "n/a", drawN, 8, jalg);
  std::pair<std::vector<Double_t>, Double_t> Ztautau = getvalues(obsID, "Ztautau", mass, "n/a", drawN, 8, jalg);
  std::pair<std::vector<Double_t>, Double_t> Zuds = getvalues(obsID, "Zuds",mass, "n/a", drawN, 8, jalg);
  std::pair<std::vector<Double_t>, Double_t> Zbb = getvalues(obsID, "Zbb",mass, "n/a", drawN, 8, jalg);
  std::pair<std::vector<Double_t>, Double_t> Zcc = getvalues(obsID, "Zcc",mass, "n/a", drawN, 8, jalg);

  
  std::pair<std::vector<Double_t>, Double_t> signalM = getvalues(obsID, "signal", mass, lt, -1, 8, jalg);

  
  std::pair<std::vector<Double_t>, Double_t> signal30 = getvalues(obsID, "signal", 30, "m3p0", -1, 8, jalg);
  std::pair<std::vector<Double_t>, Double_t> signal50 = getvalues(obsID, "signal", 50, "m3p5", -1, 8, jalg);
  std::pair<std::vector<Double_t>, Double_t> signal70 = getvalues(obsID, "signal", 70, "m4p0", -1, 8, jalg);

  
  std::pair<std::vector<Double_t>, Double_t> signal20 = getvalues(obsID, "signal", 20, "m2p0", -1, 8, jalg);
  

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
  
  */
}

  
