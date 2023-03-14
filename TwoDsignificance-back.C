#include "results.h"
#include "lumisettings.h"

Double_t AtlasZ(Double_t s, Double_t b, Double_t sig){

  Double_t n = s+b;

  if (sig<1){

    return TMath::Sqrt(2 *( n *TMath::Log(n/b) - (n-b)  ));

    
  }

  

  Double_t FirstLog = n*(b+sig*sig) / (b*b + n*sig*sig);
  FirstLog = TMath::Log(FirstLog);

  Double_t SecondLog = sig*sig*(n-b) / (b* (b+sig*sig));
  SecondLog = TMath::Log( 1 + SecondLog);

  Double_t Rad = n*FirstLog - (b*b / (sig*sig)) * SecondLog;
  Rad = TMath::Sqrt(2*Rad);

  return Rad;
  
}

Double_t GetUpp(Double_t n){

  if (n==0 || n==1) return 1.14;
  else return TMath::Sqrt(n);

}

  

std::vector<std::vector<double>> getarray(TString lifetime){

  if (lifetime == "m1p0") return CutFlowSignal_m1p0;
  if (lifetime == "m1p5") return CutFlowSignal_m1p5;
  if (lifetime == "m2p0") return CutFlowSignal_m2p0;
  if (lifetime == "m2p5") return CutFlowSignal_m2p5;
  if (lifetime == "m3p0") return CutFlowSignal_m3p0;
  if (lifetime == "m3p5") return CutFlowSignal_m3p5;
  if (lifetime == "m4p0") return CutFlowSignal_m4p0;
  if (lifetime == "m4p5") return CutFlowSignal_m4p5;
  if (lifetime == "m5p0") return CutFlowSignal_m5p0;
  if (lifetime == "m5p5") return CutFlowSignal_m5p5;
  if (lifetime == "m6p0") return CutFlowSignal_m6p0;
  if (lifetime == "m6p5") return CutFlowSignal_m6p5;
  if (lifetime == "m7p0") return CutFlowSignal_m7p0;
  
  
}

void TwoDsignificance(){

  std::vector<TString> lifetime = {"m1p0", "m1p5", "m2p0", "m2p5", "m3p0", "m3p5", "m4p0", "m4p5", "m5p0", "m5p5", "m6p0", "m6p5", "m7p0"};

  TH2F* H = new TH2F("H","H",7,15,85,35,-11,-4);

  H->GetXaxis()->SetTitle("M(HN) (GeV)");
  H->GetYaxis()->SetTitle("Log(U^{2})");

  for (int m = 30; m< 71;  m+=20){

    for (TString lt : lifetime){

    

      std::vector<std::vector<double>> vSig = getarray(lt);

      Int_t myid = 12;

      Double_t signal = vSig[(int)(m/10-1)][myid];
      Double_t Zmumu = CutFlowZmumu[(int)(m/10-1)][myid];
      Double_t Ztautau = CutFlowZtautau[(int)(m/10-1)][myid];
      Double_t Zbb = CutFlowZbb[(int)(m/10-1)][myid];
      Double_t Zcc = CutFlowZcc[(int)(m/10-1)][myid];
      Double_t Zuds = CutFlowZuds[(int)(m/10-1)][myid];
      Double_t munuqq = CutFlowmunuqq[(int)(m/10-1)][myid];

      Double_t SigmaBkg = TMath::Sqrt( TMath::Power(GetUpp(Zmumu)*Weight("Zmumu"),2) + TMath::Power(GetUpp(Ztautau)*Weight("Ztautau"),2) + TMath::Power(GetUpp(Zbb)*Weight("Zbb"),2) + TMath::Power(GetUpp(Zcc)*Weight("Zcc"),2) + TMath::Power(GetUpp(Zuds)*Weight("Zuds"),2) + TMath::Power(GetUpp(munuqq)*Weight("munuqq"),2));

      Double_t U = Coupling(Form("%d",m),lt);
      
      Double_t Y = 2.*TMath::Log10(U);
      Double_t X = 1.*m;

      Double_t totsig = signal*Weight("signal", Form("%d",m), lt);
      Double_t totbkg = Zmumu*Weight("Zmumu") + Ztautau * Weight("Ztautau") + Zbb * Weight("Zbb") + Zcc * Weight("Zcc") + Zuds * Weight("Zuds") +  munuqq * Weight("munuqq");

      //Double_t Z = totsig / TMath::Sqrt(totsig + totsig + SigmaBkg);
      Double_t Z = AtlasZ(totsig,totbkg,0);

      cout<<"M "<<m<<"     LT "<<lt<<"    logU2 "<<Y<<"    sig "<<totsig<<"    bkg "<<totbkg<<"   Errbkg "<<SigmaBkg<<"    Z "<<Z<<endl;
      cout<<"-------------------------------------------------"<<endl;
      cout<<signal<<" "<<Zbb<<" "<<Zcc<<" "<<Zuds<<" "<<Zmumu<<" "<<Ztautau<<" "<<munuqq<<endl;
      cout<<"-------------------------------------------------"<<endl;
      
      

      if (Z<100) H->Fill(X,Y,Z);
       

      
      
      
    }
  }

  H->GetZaxis()->SetRangeUser(1e-5,1000);

  
  H->Draw("colz text");
  H->SetMarkerSize(1.6);
  H->SetTitle("Significance (bkg+3 #sigma)");

  gStyle->SetOptStat(0);
  gPad->SetLogz();
  gStyle->SetPalette(kBlackBody);
}
