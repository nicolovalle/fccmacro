#include "getvalues.C"

double Zmass = 91.1876;

double UseWeight(TString opt){
  
  if (opt == "Zbb") return  Weight("Zbb");
  if (opt == "Zcc") return Weight("Zcc");
  if (opt == "Zss") return Weight("Zss");
  if (opt == "Zud") return Weight("Zud");
  if (opt == "Zuds") return Weight("Zuds");
  if (opt == "Zmumu") return Weight("Zmumu");
  if (opt == "Ztautau") return Weight("Ztautau");
  
  if (opt == "munuqq") return Weight("munuqq");
  return 0;
}

using namespace std;

void ScanWidth(){

 
  int mass = 40;
  double de = 0.2;
  const int nstep = 30;

  double precoil = (Zmass*Zmass - mass*mass)/(2.*Zmass);


  pair<vector<double>, double> Zbb, Zcc, Zuds, Zud, Zss, Zmumu, Ztautau, munuqq;
  pair<vector<double>, double> Zbb2, Zcc2, Zuds2, Zud2, Zss2, Zmumu2, Ztautau2, munuqq2;

  Zbb = getvalues(omass_encoded_dcut, "Zbb", -1, "n/a", -1, 8, 2, "< d2d dsigma anymass1L2M window [0,0]");
  Zcc = getvalues(omass_encoded_dcut, "Zcc", -1, "n/a", -1, 8, 2, "< d2d dsigma anymass1L2M window [0,0]");
  if (PRODUCTION == "Spring2021") Zuds = getvalues(omass_encoded_dcut, "Zuds", -1, "n/a", -1, 8, 2, "< d2d dsigma anymass1L2M window [0,0]");
  if (PRODUCTION == "Winter2023") Zud = getvalues(omass_encoded_dcut, "Zud", -1, "n/a", -1, 8, 2, "< d2d dsigma anymass1L2M window [0,0]");
  if (PRODUCTION == "Winter2023") Zss = getvalues(omass_encoded_dcut, "Zss", -1, "n/a", -1, 8, 2, "< d2d dsigma anymass1L2M window [0,0]");
  Zmumu = getvalues(omass_encoded_dcut, "Zmumu", -1, "n/a", -1, 8, 2, "< d2d dsigma anymass1L2M window [0,0]");
  Ztautau = getvalues(omass_encoded_dcut, "Ztautau", -1, "n/a", -1, 8, 2, "< d2d dsigma anymass1L2M window [0,0]");
  munuqq = getvalues(omass_encoded_dcut, "munuqq", -1, "n/a", -1, 8, 2, "< d2d dsigma anymass1L2M window [0,0]");

  Zbb2 = getvalues(oemiss_encoded_dcut, "Zbb", -1, "n/a", -1, 8, 2, "< d2d dsigma anymass1L2M window [0,0]");
  Zcc2 = getvalues(oemiss_encoded_dcut, "Zcc", -1, "n/a", -1, 8, 2, "< d2d dsigma anymass1L2M window [0,0]");
  if (PRODUCTION == "Spring2021") Zuds2 = getvalues(oemiss_encoded_dcut, "Zuds", -1, "n/a", -1, 8, 2, "< d2d dsigma anymass1L2M window [0,0]");
  if (PRODUCTION == "Winter2023") Zud2 = getvalues(oemiss_encoded_dcut, "Zud", -1, "n/a", -1, 8, 2, "< d2d dsigma anymass1L2M window [0,0]");
  if (PRODUCTION == "Winter2023") Zss2 = getvalues(oemiss_encoded_dcut, "Zss", -1, "n/a", -1, 8, 2, "< d2d dsigma anymass1L2M window [0,0]");
  Zmumu2 = getvalues(oemiss_encoded_dcut, "Zmumu", -1, "n/a", -1, 8, 2, "< d2d dsigma anymass1L2M window [0,0]");
  Ztautau2 = getvalues(oemiss_encoded_dcut, "Ztautau", -1, "n/a", -1, 8, 2, "< d2d dsigma anymass1L2M window [0,0]");
  munuqq2 = getvalues(oemiss_encoded_dcut, "munuqq", -1, "n/a", -1, 8, 2, "< d2d dsigma anymass1L2M window [0,0]");
  
  double x[nstep], y[nstep];

  for (int i=0; i<nstep; i++) {
    x[i] = de*(i+1);
    y[i] = 0;
  }


  vector<double> V1, V2;
  double W;

  ////////////
  V1 = Zbb.first;
  V2 = Zbb2.first;
  W = UseWeight("Zbb");
  for (int iv = 0; iv < V1.size(); iv++){
    double im = V1.at(iv);
    double ie = V2.at(iv);
    if (im>0) continue;
    for (int i=0; i<nstep; i++){
      if ( TMath::Abs( TMath::Abs(im) - mass) < x[i] && TMath::Abs(TMath::Abs(ie) - precoil) < x[i]) y[i]+=W;
    }
  }


  //////////////
  V1 = Zcc.first;
  V2 = Zcc2.first;
  W = UseWeight("Zcc");
  for (int iv = 0; iv < V1.size(); iv++){
    double im = V1.at(iv);
    double ie = V2.at(iv);
    if (im>0) continue;
    for (int i=0; i<nstep; i++){
      if ( TMath::Abs( TMath::Abs(im) - mass) < x[i] && TMath::Abs(TMath::Abs(ie) - precoil) < x[i]) y[i]+=W;
    }
  }

  if (PRODUCTION == "Spring2021"){
    //////////////
    V1 = Zuds.first;
    V2 = Zuds2.first;
    W = UseWeight("Zuds");
    for (int iv = 0; iv < V1.size(); iv++){
      double im = V1.at(iv);
      double ie = V2.at(iv);
      if (im>0) continue;
      for (int i=0; i<nstep; i++){
	if ( TMath::Abs( TMath::Abs(im) - mass) < x[i] && TMath::Abs(TMath::Abs(ie) - precoil) < x[i]) y[i]+=W;
      }
    }
  }


  if (PRODUCTION == "Winter2023"){
    //////////////
    V1 = Zud.first;
    V2 = Zud2.first;
    W = UseWeight("Zud");
    for (int iv = 0; iv < V1.size(); iv++){
      double im = V1.at(iv);
      double ie = V2.at(iv);
      if (im>0) continue;
      for (int i=0; i<nstep; i++){
	if ( TMath::Abs( TMath::Abs(im) - mass) < x[i] && TMath::Abs(TMath::Abs(ie) - precoil) < x[i]) y[i]+=W;
      }
    }
    //////////////
    V1 = Zss.first;
    V2 = Zss2.first;
    W = UseWeight("Zss");
    for (int iv = 0; iv < V1.size(); iv++){
      double im = V1.at(iv);
      double ie = V2.at(iv);
      if (im>0) continue;
      for (int i=0; i<nstep; i++){
	if ( TMath::Abs( TMath::Abs(im) - mass) < x[i] && TMath::Abs(TMath::Abs(ie) - precoil) < x[i]) y[i]+=W;
      }
    }
  }




  //////////////
  V1 = Zmumu.first;
  V2 = Zmumu2.first;
  W = UseWeight("Zmumu");
  for (int iv = 0; iv < V1.size(); iv++){
    double im = V1.at(iv);
    double ie = V2.at(iv);
    if (im>0) continue;
    for (int i=0; i<nstep; i++){
      if ( TMath::Abs( TMath::Abs(im) - mass) < x[i] && TMath::Abs(TMath::Abs(ie) - precoil) < x[i]) y[i]+=W;
    }
  }


  //////////////
  V1 = Ztautau.first;
  V2 = Ztautau2.first;
  W = UseWeight("Ztautau");
  for (int iv = 0; iv < V1.size(); iv++){
    double im = V1.at(iv);
    double ie = V2.at(iv);
    if (im>0) continue;
    for (int i=0; i<nstep; i++){
      if ( TMath::Abs( TMath::Abs(im) - mass) < x[i] && TMath::Abs(TMath::Abs(ie) - precoil) < x[i]) y[i]+=W;
    }
  }

  //////////////
  V1 = munuqq.first;
  V2 = munuqq2.first;
  W = UseWeight("munuqq");
  for (int iv = 0; iv < V1.size(); iv++){
    double im = V1.at(iv);
    double ie = V2.at(iv);
    if (im>0) continue;
    for (int i=0; i<nstep; i++){
      if ( TMath::Abs( TMath::Abs(im) - mass) < x[i] && TMath::Abs(TMath::Abs(ie) - precoil) < x[i]) y[i]+=W;
    }
  }





  TGraph *g = new TGraph(nstep,x,y);

  g->SetTitle(Form("Centered on mass M = %d",(int)mass));
  g->GetXaxis()->SetTitle("#DeltaE");
  g->GetYaxis()->SetTitle("Number of events with |m-M|<#DeltaE and |p-P|<#DeltaE");

  g->Draw("APL");

  TLine *l1 = new TLine(2* 0.2 * TMath::Sqrt(mass), 0, 2* 0.2 * TMath::Sqrt(mass), y[nstep-1]);
  TLine *l2 = new TLine(2* 0.3 * TMath::Sqrt(mass), 0, 2* 0.3 * TMath::Sqrt(mass), y[nstep-1]);
  l1->Draw("same");
  l2->Draw("same");
  
  
  
    
  
  

}
