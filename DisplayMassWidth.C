#include "getvalues.C"

double Zmass = 91.1876;



using namespace std;

void DisplayMassWidth(){

  
  
  double de = 0.01;

  double precoil50 = (Zmass*Zmass - 50*50)/(2.*Zmass);
  double precoil20 = (Zmass*Zmass - 20*20)/(2.*Zmass);


  pair<vector<double>, double> sig50e, sig50m, sig20e, sig20m;
  

  sig50m = getvalues(omass_selection,  "signal", 50, "m3p5", -1, 0, 2, "> d2d dsigma anymass1L2M window [99,99]");
  sig20m = getvalues(omass_selection,  "signal", 20, "m3p5", -1, 0, 2, "> d2d dsigma anymass1L2M window [99,99]");
  sig50e = getvalues(oemiss_selection, "signal", 50, "m3p5", -1, 0, 2, "> d2d dsigma anymass1L2M window [99,99]");
  sig20e = getvalues(oemiss_selection, "signal", 20, "m3p5", -1, 0, 2, "> d2d dsigma anymass1L2M window [99,99]");

  int N50 = sig50m.first.size();
  int N20 = sig20m.first.size();

  double dM50 =0, dE50 = 0, dM20 = 0, dE20 = 0;


  int counter = 0;

  for (double x = 0;; x+=de){
    counter = 0;
    for (double s : sig50m.first) if (TMath::Abs( s - 50.) <= x) counter++;
    if (counter >= 0.68*N50) { dM50 = x; break;}
  }

  for (double x = 0;; x+=de){
    counter = 0;
    for (double s : sig20m.first) if (TMath::Abs( s - 20.) <= x) counter++;
    if (counter >= 0.68*N20) { dM20 = x; break;}
  }

  for (double x = 0;; x+=de){
    counter = 0;
    for (double s : sig50e.first) if (TMath::Abs( s - precoil50) <= x) counter++;
    if (counter >= 0.68*N50) { dE50 = x; break;}
  }

  for (double x = 0;; x+=de){
    counter = 0;
    for (double s : sig20e.first) if (TMath::Abs( s - precoil20) <= x) counter++;
    if (counter >= 0.68*N20) { dE20 = x; break;}
  }


  cout<<"dM50= "<<dM50<<endl;
  cout<<"dE50= "<<dE50<<endl;
  cout<<"dM20= "<<dM20<<endl;
  cout<<"dE20= "<<dE20<<endl;
  

  double epsilon = 0; //-0.33;
 
  TH1F *h50m = new TH1F("h50m","",60,46+epsilon,54+epsilon);
  TH1F *h50msig = new TH1F("h50msig","; M_{vis} (GeV/c^{2}); counts",60,46+epsilon,54+epsilon);
  
  TH1F *h20m = new TH1F("h20m","h20m",30,40,60);
  TH1F *h50e = new TH1F("h50e","h50e",30,40,60);
  TH1F *h20e = new TH1F("h20e","h20e",30,40,60);

  
  for (double s: sig50m.first) h50m->Fill(s+epsilon);
  for (double s: sig20m.first) h20m->Fill(s);
  for (double s: sig50e.first) h50e->Fill(s);
  for (double s: sig20e.first) h20e->Fill(s);

  for (int i=1; i< h50msig->GetNbinsX(); i++){
    if (TMath::Abs(h50msig->GetBinCenter(i) - 50 ) <= dM50) h50msig->SetBinContent(i, h50m->GetBinContent(i));
  }

  h50msig->SetTitle("Visible mass distribution for signal sample with M_{HN} = 50 GeV/c^{2}");

  h50m->SetLineColor(1);
  h20m->SetLineColor(1);

  h50e->SetLineColor(2);
  h20e->SetLineColor(2);


 
  h50msig->SetFillStyle(3356);
  h50msig->SetFillColor(1);
  h50msig->SetLineColor(0);
  h50msig->Draw("HIST ");
   h50m->Draw("HIST same");
   //h20e->Draw("same HISO");
  //h20m->Draw("same HISTO");
  //h20m->Draw("same HISTO");

   gStyle->SetOptStat(0);
  
  /*

 




  TGraph *g = new TGraph(nstep,x,y);

  g->SetTitle(Form("Centered on mass M = %d",(int)mass));
  g->GetXaxis()->SetTitle("#DeltaE");
  g->GetYaxis()->SetTitle("Number of events with |m-M|<#DeltaE and |p-P|<#DeltaE");

  g->Draw("APL");

  TLine *l1 = new TLine(2* 0.2 * TMath::Sqrt(mass), 0, 2* 0.2 * TMath::Sqrt(mass), y[nstep-1]);
  TLine *l2 = new TLine(2* 0.3 * TMath::Sqrt(mass), 0, 2* 0.3 * TMath::Sqrt(mass), y[nstep-1]);
  l1->Draw("same");
  l2->Draw("same");
  
  
  */
    
  
  

}
