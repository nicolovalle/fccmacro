
#include "getvalues.C"

void DrawTH1(TString obsID="d0sliding", Int_t nbin=50, Double_t bmin=0, Double_t bmax=100, Bool_t drawstack = true){

  TCanvas *c1 = new TCanvas("c1","c1",0,0,600,400); 

  std::pair<std::vector<Double_t>, Double_t> munuqq = getvalues(obsID, "munuqq", 1.e5);
  std::pair<std::vector<Double_t>, Double_t> Zmumu = getvalues(obsID, "Zmumu", 1.e5);
  std::pair<std::vector<Double_t>, Double_t> Ztautau = getvalues(obsID, "Ztautau", 1.e5);
  std::pair<std::vector<Double_t>, Double_t> Zuds = getvalues(obsID, "Zuds",1.e5);
  std::pair<std::vector<Double_t>, Double_t> Zbb = getvalues(obsID, "Zbb",1.e5);
  std::pair<std::vector<Double_t>, Double_t> Zcc = getvalues(obsID, "Zcc",1.e5);

  //std::pair<std::vector<Double_t>, Double_t> signal30 = getvalues(obsID, "signal", -1, 30);
  std::pair<std::vector<Double_t>, Double_t> signal50 = getvalues(obsID, "signal", -1, 50);
  //std::pair<std::vector<Double_t>, Double_t> signal70 = getvalues(obsID, "signal", -1, 70);

  

  TH1F *h_munuqq = new TH1F("munuqq","munuqq",nbin,bmin,bmax);
  if (drawstack) h_munuqq->SetFillColor(kRed);
  h_munuqq->SetLineColor(kRed);
 
  TH1F *h_zmumu = new TH1F("zmumu","zmumu",nbin,bmin,bmax);
  if (drawstack) h_zmumu->SetFillColor(kBlue);
  h_zmumu->SetLineColor(kBlue);
 
  TH1F *h_ztautau = new TH1F("ztautau","ztautau",nbin,bmin,bmax);
  if (drawstack) h_ztautau->SetFillColor(32);
  h_ztautau->SetLineColor(32);
  
  TH1F *h_zuds = new TH1F("zuds","zuds",nbin,bmin,bmax);
  if (drawstack) h_zuds->SetFillColor(kOrange);
  h_zuds->SetLineColor(kOrange);
  
  TH1F *h_zbb = new TH1F("zbb","zbb",nbin,bmin,bmax);
  if (drawstack) h_zbb->SetFillColor(kGreen);
  h_zbb->SetLineColor(kGreen);
 
  TH1F *h_zcc = new TH1F("zcc","zcc",nbin,bmin,bmax);
  if (drawstack) h_zcc->SetFillColor(kViolet);
  h_zcc->SetLineColor(kViolet);



  TH1F *h_signal50 = new TH1F("signal50","signal50",nbin,bmin,bmax);
  h_signal50->SetLineWidth(3);
  h_signal50->SetLineColor(1);
  
  
  /*
  for (int i=0; i<=Weight("munuqq");i++) for (double ie : munuqq) h_munuqq->Fill(ie);
  cout<<"filled 1"<<endl;
  for (int i=0; i<=Weight("Zmumu");i++) for (double ie : Zmumu) h_zmumu->Fill(ie);
  cout<<"filled 2"<<endl;
  for (int i=0; i<=Weight("Ztautau");i++) for (double ie : Ztautau) h_ztautau->Fill(ie);
  cout<<"filled 3"<<endl;
  for (int i=0; i<=Weight("Zuds");i++) for (double ie : Zuds) h_zuds->Fill(ie);
  cout<<"filled 4"<<endl;
  for (int i=0; i<=Weight("Zbb");i++) for (double ie : Zbb) h_zbb->Fill(ie);
  cout<<"filled 5"<<endl;
  for (int i=0; i<=Weight("Zcc");i++) for (double ie : Zcc) h_zcc->Fill(ie);
  cout<<"filled 6"<<endl;

  for (int i=0; i<=Weight("signal");i++) for (double ie : signal50) h_signal50->Fill(ie);
  cout<<"filled 7"<<endl;
  */

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

  // for (int i=0; i<=Weight("signal","50","n/a",signal50.second);i++)
    for (double ie : signal50.first) h_signal50->Fill(ie);
  cout<<"filled 7"<<endl;


  h_munuqq->Scale(Weight("munuqq","n/a","n/a",munuqq.second));
  h_zmumu->Scale(Weight("Zmumu","n/a","n/a",Zmumu.second));
  h_ztautau->Scale(Weight("Ztautau","n/a","n/a",Ztautau.second));
  h_zuds->Scale(Weight("Zuds","n/a","n/a",Zuds.second));
  h_zbb->Scale(Weight("Zbb","n/a","n/a",Zbb.second));
  h_zcc->Scale(Weight("Zcc","n/a","n/a",Zcc.second));
  h_signal50->Scale(Weight("signal","50","n/a",signal50.second));

  

  auto hs = new THStack("hs","");
  hs->Add(h_munuqq);
  hs->Add(h_zmumu);
  hs->Add(h_ztautau);
  hs->Add(h_zuds);
  hs->Add(h_zbb);
  hs->Add(h_zcc);

  hs->Draw(drawstack ? "HIST" : "HIST nostack");

  h_signal50->Draw("HIST same");

  gPad->SetLogy();

  c1->BuildLegend();
  
  
}

  
