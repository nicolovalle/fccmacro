#include "results.h"

#include "CutFlowOK.C"

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

  


std::vector<std::vector<double>> TwoDsignificance(Int_t dd0cut = 8, Bool_t Draw = true){

  // V[0][n] = n-th x coordinate of the curve
  // V[1][n] = n-th y coordinate of the curve


  std::vector<int> masses = {2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 85};
  std::vector<TString> mps = {"m","p"};

  //std::vector<TString> lifetime = {"m1p0", "m1p5", "m2p0", "m2p5", "m3p0", "m3p5", "m4p0", "m4p5", "m5p0", "m5p5", "m6p0", "m6p5", "m7p0"};

  TH2F* H = new TH2F("H","H",18,0,90,35,-11,-4);

  H->GetXaxis()->SetTitle("M(HN) (GeV)");
  H->GetYaxis()->SetTitle("Log(U^{2})");

  //X,Y,SIGNIFICANCE:
  std::vector<double> TDPlot_M;
  std::vector<double> TDPlot_U2;
  std::vector<double> TDPlot_Z;

  for (int m : masses){
    for (TString mp : mps){
      for (int unit = 0; unit < 9; unit++){
	for (int decimal = 0; decimal < 6; decimal += 5){


	  TString lt = Form("%s%dp%d", mp.Data(), unit, decimal);

    


	  TString AnalysisResPath = "/home/nvalle/FCC/localcopy/workdir/MyExternalAnalysis/results/skimmed/";

	  TString FileNameToCheck = Form("%s%s", AnalysisResPath.Data(), AnalysisResults("signal",Form("%d",m),lt).Data());
      
	  if (gSystem->AccessPathName(FileNameToCheck)) continue;
     
	  Int_t myid = dd0cut + 4; // 12 means d0Cut = 8

      

	  Double_t signal = CutFlowOK("signal",m,lt,AnalysisResPath)[myid];

	  if (signal == 0){
	    cout<<"From "<<FileNameToCheck<<" the signal has 0 events"<<endl;
	    continue;
	  }

	  /*
	  Double_t Zmumu = CutFlowZmumu[(int)(m/10-1)][myid];
	  Double_t Ztautau = CutFlowZtautau[(int)(m/10-1)][myid];
	  Double_t Zbb = CutFlowZbb[(int)(m/10-1)][myid];
	  Double_t Zcc = CutFlowZcc[(int)(m/10-1)][myid];
	  Double_t Zuds = CutFlowZuds[(int)(m/10-1)][myid];
	  Double_t munuqq = CutFlowmunuqq[(int)(m/10-1)][myid];
	  */
     
	  
	    Double_t Zmumu = CutFlowOK("Zmumu",m,lt,AnalysisResPath)[myid];
	    Double_t Ztautau = CutFlowOK("Ztautau",m,lt,AnalysisResPath)[myid];
	    Double_t Zbb = CutFlowOK("Zbb",m,lt,AnalysisResPath)[myid];
	    Double_t Zcc = CutFlowOK("Zcc",m,lt,AnalysisResPath)[myid];
	    Double_t Zuds = CutFlowOK("Zuds",m,lt,AnalysisResPath)[myid];
	    Double_t munuqq = CutFlowOK("munuqq",m,lt,AnalysisResPath)[myid];
	  

	  Double_t SigmaBkg = TMath::Sqrt( TMath::Power(GetUpp(Zmumu)*Weight("Zmumu"),2) + TMath::Power(GetUpp(Ztautau)*Weight("Ztautau"),2) + TMath::Power(GetUpp(Zbb)*Weight("Zbb"),2) + TMath::Power(GetUpp(Zcc)*Weight("Zcc"),2) + TMath::Power(GetUpp(Zuds)*Weight("Zuds"),2) + TMath::Power(GetUpp(munuqq)*Weight("munuqq"),2));

	  Double_t U = Coupling(Form("%d",m),lt);
      
	  Double_t Y = 2.*TMath::Log10(U);
	  Double_t X = 1.*m;

	  Double_t totsig = signal*Weight("signal", Form("%d",m), lt);
	  Double_t totbkg = Zmumu*Weight("Zmumu") + Ztautau * Weight("Ztautau") + Zbb * Weight("Zbb") + Zcc * Weight("Zcc") + Zuds * Weight("Zuds") +  munuqq * Weight("munuqq");

	  Double_t Z = totsig / TMath::Sqrt(totsig + totbkg);// + SigmaBkg);
	  //Double_t Z = AtlasZ(totsig,totbkg,0);

	  cout<<"M "<<m<<"     LT "<<lt<<"    logU2 "<<Y<<"    sig "<<totsig<<"    bkg "<<totbkg<<"   Errbkg "<<SigmaBkg<<"    Z "<<Z<<endl;
	  cout<<"-------------------------------------------------"<<endl;
	  cout<<signal<<" "<<Zbb<<" "<<Zcc<<" "<<Zuds<<" "<<Zmumu<<" "<<Ztautau<<" "<<munuqq<<endl;
	  cout<<"-------------------------------------------------"<<endl;
      
      

	  if (Z<100 && Z>0.1) {
	    TDPlot_M.push_back(X);
	    TDPlot_U2.push_back(Y);
	    TDPlot_Z.push_back(Z);
	    H->Fill(X,Y,Z);
	  }

	}
      }
    }
  }

      
      

  // DRAWING THE LINE BY LINEAR INTERPOLATION

  std::vector<double> cx;
  std::vector<double> cy;

  
  for (int i : masses){

    std::vector<double> ycoordinates;

    for (int j=0; j < TDPlot_M.size(); j++){
      if (int(TDPlot_M[j]) == i) ycoordinates.push_back(TDPlot_U2[j]);
    }

    std::sort(ycoordinates.begin(), ycoordinates.end());
    

    double spre = 0;
    double spost = 0;
    double upre = -99;
    double upost = -99;
  
    for (double iY : ycoordinates){

      spre = spost;
      upre = upost;

      for (int j=0; j < TDPlot_M.size(); j++){
	if (int(TDPlot_M[j]) == i && TDPlot_U2[j] == iY) spost = TDPlot_Z[j];
      }
      upost = iY;

      if (spre < 2 && spost >= 2){

	cx.push_back(i);

	double interpol = upre + (2.-spre)*((upost-upre)/(spost-spre));

	cy.push_back(interpol);

      }     
    }
  }


  const Int_t npoint = cx.size();
  double vcx[npoint], vcy[npoint];
  for (int i = 0; i<npoint; i++) {vcx[i] = cx[i]; vcy[i] = cy[i];}

  auto gg = new TGraph(npoint,vcx,vcy);
    
   

  H->GetZaxis()->SetRangeUser(1e-5,1000);

  if (Draw){
  
    H->Draw("colz text");
    H->SetMarkerSize(1.6);
    H->SetTitle("Significance");

    //gg->Draw("same");

    gStyle->SetOptStat(0);
    gPad->SetLogz();
    gStyle->SetPalette(kBlackBody);
  }

  std::vector<std::vector<double>> toret;
  toret.push_back(cx);
  toret.push_back(cy);

  return toret;
}


void ScanD0Cut(){

  auto c = new TCanvas("c","c",1000,600);

  Int_t dcut;

  std::vector<std::vector<double>> XY;
  
  dcut = 2;
  double_t x2[50], y2[50];
  XY = TwoDsignificance(dcut,false);
  for (int i=0; i<XY[0].size(); i++){
    x2[i] = XY[0][i];
    y2[i] = XY[1][i];
  }
  auto g2 = new TGraph(XY[0].size(), x2, y2);
  g2->SetLineColor(dcut);
  g2->SetTitle(Form("D0 cut = %d",dcut));


  dcut = 4;
  double_t x4[50], y4[50];
  XY = TwoDsignificance(dcut,false);
  for (int i=0; i<XY[0].size(); i++){
    x4[i] = XY[0][i];
    y4[i] = XY[1][i];
  }
  auto g4 = new TGraph(XY[0].size(), x4, y4);
  g4->SetLineColor(dcut);
  g4->SetTitle(Form("D0 cut = %d",dcut));


  dcut = 6;
  double_t x6[50], y6[50];
  XY = TwoDsignificance(dcut,false);
  for (int i=0; i<XY[0].size(); i++){
    x6[i] = XY[0][i];
    y6[i] = XY[1][i];
  }
  auto g6 = new TGraph(XY[0].size(), x6, y6);
  g6->SetLineColor(dcut);
  g6->SetTitle(Form("D0 cut = %d",dcut));


  dcut = 8;
  double_t x8[50], y8[50];
  XY = TwoDsignificance(dcut,false);
  for (int i=0; i<XY[0].size(); i++){
    x8[i] = XY[0][i];
    y8[i] = XY[1][i];
  }
  auto g8 = new TGraph(XY[0].size(), x8, y8);
  g8->SetLineColor(dcut);
  g8->SetTitle(Form("D0 cut = %d",dcut));


  dcut = 10;
  double_t x10[50], y10[50];
  XY = TwoDsignificance(dcut,false);
  for (int i=0; i<XY[0].size(); i++){
    x10[i] = XY[0][i];
    y10[i] = XY[1][i];
  }
  auto g10 = new TGraph(XY[0].size(), x10, y10);
  g10->SetLineColor(dcut);
  g10->SetTitle(Form("D0 cut = %d",dcut));
  


  dcut = 12;
  double_t x12[50], y12[50];
  XY = TwoDsignificance(dcut,false);
  for (int i=0; i<XY[0].size(); i++){
    x12[i] = XY[0][i];
    y12[i] = XY[1][i];
  }
  auto g12 = new TGraph(XY[0].size(), x12, y12);
  g12->SetLineColor(dcut);
  g12->SetTitle(Form("D0 cut = %d",dcut));


  dcut = 16;
  double_t x16[50], y16[50];
  XY = TwoDsignificance(dcut,false);
  for (int i=0; i<XY[0].size(); i++){
    x16[i] = XY[0][i];
    y16[i] = XY[1][i];
  }
  auto g16 = new TGraph(XY[0].size(), x16, y16);
  g16->SetLineColor(dcut);
  g16->SetTitle(Form("D0 cut = %d",dcut));

  
  dcut = 20;
  double_t x20[50], y20[50];
  XY = TwoDsignificance(dcut,false);
  for (int i=0; i<XY[0].size(); i++){
    x20[i] = XY[0][i];
    y20[i] = XY[1][i];
  }
  auto g20 = new TGraph(XY[0].size(), x20, y20);
  g20->SetLineColor(dcut);
  g20->SetTitle(Form("D0 cut = %d",dcut));


  // g2->Draw("");
  g4->Draw();
  g6->Draw("same");
  
  g8->Draw("same");
  g10->Draw("same");
  g12->Draw("same");
  g16->Draw("same");
  g20->Draw("same");
  

  c->BuildLegend();
}
