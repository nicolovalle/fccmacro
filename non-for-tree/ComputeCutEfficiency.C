#include "lumisettings.h"

Double_t GetUpper(Double_t x);
float GetContentByLabel(TH1F *h, TString lab);

std::vector<Double_t> ComputeCutEfficiency(TString HNMass="50", TString path="../MyExternalAnalysis/results/"){

  TString slidingcut = "sliding2";
  Double_t bcut = 8;
  
  std::vector<Double_t> Npreselection;
  std::vector<Double_t> Nselection;
  std::vector<Double_t> NslidingNobcut;
  std::vector<Double_t> Nsliding;

  std::vector<Double_t> UppBkg95;
  std::vector<Double_t> Weight;

  

  for (TString ID : id){
    cout<<"Reading "<<ID<<" for HN_mass = "<<HNMass<<endl;
    TString opt = (ID == "signal") ? HNMass : ID;
    TFile *f = new TFile(Form("%s%s",path.Data(),AnalysisResults(opt).Data()));
    TH1F *h = (TH1F*) f->Get("Nevents");
    Npreselection.push_back(GetContentByLabel(h, "preselection"));
    Nselection.push_back(GetContentByLabel(h, "selection2"));
    NslidingNobcut.push_back(GetContentByLabel(h, Form("%s_m%s",slidingcut.Data(),HNMass.Data())));
    if (bcut<0){
      Nsliding.push_back(GetContentByLabel(h, Form("%s_m%s",slidingcut.Data(),HNMass.Data())));
    }
    else{
      TH1F *h2 = (TH1F*) f->Get(Form("%s_m%s/bs_mu",slidingcut.Data(),HNMass.Data()));
      Nsliding.push_back(h2->Integral(0,(int)(bcut*10)));
    }
    h->Delete();
    f->Close();
    f->Delete();
  }

  for (int i=0; i<id.size(); i++){
    Weight.push_back(xsec(id.at(i),HNMass)*LUMI/Nnocut[i]);
    UppBkg95.push_back(GetUpper(Nsliding[i])*Weight[i]);
  }

  for (int i=0;i<id.size();i++){
    cout<<id.at(i)<<" - HN_mass="<<HNMass<<endl;
    cout<<Nnocut.at(i)<<"\t"<<Npreselection.at(i)<<"\t"<<Nselection.at(i)<<"\t"<<Nsliding.at(i)<<"\t"<<endl<<endl;
  }

  cout<<endl<<"============================="<<endl;

  for (int i=0;i<id.size();i++){
    cout<<id.at(i)<<" - HN_mass="<<HNMass<<endl;
    cout<<Nnocut.at(i)<<"\t"<<Npreselection.at(i)*100./Nnocut.at(i)<<"\t"<<Nselection.at(i)*100./Nnocut.at(i)<<"\t"<<Nsliding.at(i)*100./Nnocut.at(i)<<"\t"<<endl<<endl;
  }

  cout<<endl<<"=============================";
  cout<<endl<<"============================="<<endl;

  Double_t S = Nsliding[0]*Weight[0];
  
  for (int i=0;i<id.size();i++){
    
    
    cout<<id.at(i)<<" - HN_mass="<<HNMass<<endl;
    cout<<"Generated: "<<Nnocut[i]<<endl;
    cout<<"Weight: "<<Weight[i]<<endl;
    cout<<"Cut eff: "<<1. - Nsliding[i]/Nnocut[i]<<endl;
    cout<<Nsliding[i]<<" (< "<<GetUpper(Nsliding[i])<<")  --->   < "<<UppBkg95[i]<<endl;
    if (id.at(i) != "signal") cout<<"S/sqrt(S+B) > "<< S / TMath::Sqrt(S + UppBkg95[i])<<endl<<endl;
    
  }

  cout<<endl<<"============================="<<endl;

  std::vector<Double_t> toreturn;
  Double_t Rtotbkg = 0;
  for (int i=1; i<id.size();i++) Rtotbkg += UppBkg95.at(i) * UppBkg95.at(i);
  Rtotbkg = TMath::Sqrt(Rtotbkg);
  Double_t Rsignal = 1.*S;
  Double_t Rbkg_trustLept = 0;
  for (int i=1; i<id.size();i++) {
    if (id.at(i) == "Zmumu" || id.at(i) == "Ztautau") Rbkg_trustLept += Nsliding[i]*Nsliding[i]*Weight[i]*Weight[i];
    else Rbkg_trustLept += UppBkg95.at(i) * UppBkg95.at(i);
  }
  Rbkg_trustLept = TMath::Sqrt(Rbkg_trustLept);
  

  cout<<"+ Signal:       "<<(int)Rsignal/Weight[0]<<" x "<<Weight[0]<<" = "<<Rsignal<<endl;
  cout<<"+ Bkg1:          "<<Rtotbkg<<endl;
  cout<<"    + S/sqrt(S+B) = "<<Rsignal/TMath::Sqrt(Rsignal + Rtotbkg)<<endl;
  

  Double_t DesiredRatio = 2.;
  Double_t ScaleFactor = (DesiredRatio * Rsignal * TMath::Sqrt(4.* Rtotbkg + DesiredRatio * DesiredRatio) + DesiredRatio * DesiredRatio * Rsignal) / (2. * Rsignal * Rsignal);
  cout<<"    + To get significance = "<<DesiredRatio<<", the factor needed on signal is: "<<ScaleFactor<<endl;
  cout<<"    + U^2 = "<<std::scientific<<1.e-4 * ScaleFactor <<endl<<std::defaultfloat;
  toreturn.push_back(ScaleFactor);

  cout<<"+ Bkg2:      "<<Rbkg_trustLept<<endl;
  cout<<"    + S/sqrt(S+B) = "<<Rsignal/TMath::Sqrt(Rsignal + Rbkg_trustLept)<<endl;
  ScaleFactor = (DesiredRatio * Rsignal * TMath::Sqrt(4.* Rbkg_trustLept + DesiredRatio * DesiredRatio) + DesiredRatio * DesiredRatio * Rsignal) / (2. * Rsignal * Rsignal);
  cout<<"    + To get significance = "<<DesiredRatio<<", the factor needed on signal is: "<<ScaleFactor<<endl;
  cout<<"    + U^2 = "<<std::scientific<<1.e-4 * ScaleFactor <<endl<<std::defaultfloat;
  toreturn.push_back(ScaleFactor);
  

  cout<<endl<<"*****************************"<<endl;

  //TString row;
  for (int i=0; i<id.size(); i++){
    TString row = " | ";
    row += id.at(i) + " | " + Form("%7.1f M",Nnocut[i]*1.e-6) + " | " + Form("%7.1f", xsec(id.at(i),HNMass)*LUMI/Nnocut[i]);
    row += " | " + TString((1. - NslidingNobcut[i]/Nnocut[i] < 0.999) ? Form("%2.1f%%",100. - 100.*NslidingNobcut[i]/Nnocut[i]) : Form("%d evt left",(int)NslidingNobcut[i]));
    row += " | " + TString((1. - Nsliding[i]/Nnocut[i] < 0.999) ? Form("%2.1f%%",100. - 100.*Nsliding[i]/Nnocut[i]) : Form("%d evt left",(int)Nsliding[i]));
    row += " | " + TString(Form("%d",(int)(.5+GetUpper(Nsliding[i])*LUMI*xsec(id.at(i),HNMass)/Nnocut[i]))) + " | ";
    cout<<row<<endl;
  }

  // cout<<row<<endl;
    

  return toreturn;

  
}


Double_t CPsum(Int_t N, Int_t trials, Double_t p){

  Double_t toret = 0;
  for (int i=0; i<=N; i++){

    toret += ROOT::Math::binomial_pdf(i,p,trials);
  }

  return toret;
}



Double_t GetUpperCP(Int_t success, Int_t trials, Double_t CL){

  Double_t precision = 1./trials;

  Double_t seed = 1.*success/trials;
  seed = seed + 2.* TMath::Sqrt( seed * (1.-seed) / trials) + precision;

  Double_t direction;

  

  Double_t toret = seed;

  for (int i=0; i<5; i++){

    seed = toret;
    precision = precision / 10.;

    if ( CPsum(success, trials, seed) > 1.-CL) direction = 1.;
    else direction = -1.;

    for (Double_t pp = seed ;; pp = pp + direction*precision){

      toret = pp;

      Double_t newsum = CPsum(success, trials, pp);

      //cout<<"f("<<pp<<") = "<<newsum<<endl;

      if (pp < 0. || pp > 1.) break;

      if ( direction < 0 &&  newsum > 1.-CL) break;

      if ( direction > 0 &&  newsum < 1.-CL) break; 
    }
    
    toret = toret - direction * precision * 0.5;
    cout<<"Approx with precision "<<precision<<"   --> "<<toret<<endl;
  }

  return toret - direction * precision * 0.5;
}

  
  



Double_t GetUpper95(Double_t x){

  if (x==0) return 3;
  if (x==1) return 4.74;
  if (x==2) return 6.3;
  if (x==3) return 7.75;
  if (x==4) return 9.15;
  if (x==5) return 10.51;
  if (x==6) return 11.84;
  if (x==7) return 13.15;
  if (x==8) return 14.43;
  if (x==9) return 15.7;
  return x+1.65*TMath::Sqrt(x);
}
  

Double_t GetUpper(Double_t x){
  return GetUpper95(x);
}

float GetContentByLabel(TH1F *h, TString lab){

  float toret = -1;
  for (int i=1; i <= h->GetNbinsX(); i++){

    if ((TString)h->GetXaxis()->GetBinLabel(i) == lab)
      return h->GetBinContent(i);
  }

  return toret;
}

    

