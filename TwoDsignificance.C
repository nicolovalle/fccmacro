#include "CutFlowOK.C"


std::map<int,double> SmoothedLL ={{5 , 831.16352}, {10 , 2574.8944}, { 20 , 2986.9197}, { 30 , 3721.1144}, { 40 , 1468.0616}, { 50 , 2259.2084}, { 60 , 7043.7216}, { 70 , 28927.867}, { 80 , 192722.72}, { 85 , 167167.23 }};

std::map<int,double> SmoothedLL_noRoll = {{ 5 , 447.50697},{ 10 , 1350.2468},{ 20 , 1506.8090},{ 30 , 1753.5736},{ 40 , 670.10691},{ 50 , 1148.5794},{ 60 , 3377.6229},{ 70 , 12411.982},{ 80 , 81230.329},{ 85 , 81761.084 }};

std::map<int,double> WeightedLL =  {{ 5 , 0.0000000}, { 10 , 0.0000000}, { 20 , 15.134251}, { 30 , 369.36671}, { 40 , 261.55994}, { 50 , 1051.4037}, { 60 , 3161.1757}, { 70 , 12383.366}, { 80 , 83187.383}, { 85 , 81564.470 }};


Double_t AtlasZ(Double_t s, Double_t b, Double_t sig){
  // http://cds.cern.ch/record/2736148/files/ATL-PHYS-PUB-2020-025.pdf

  Double_t n = s+b;

  if (b<=1.) return AtlasZ(s,1.01,sig);

  if (sig<1e-6)

    return TMath::Sqrt(2 *( n *TMath::Log(n/b) - (n-b) ));

  
  Double_t FirstLog = n*(b+sig*sig) / (b*b + n*sig*sig);
  FirstLog = TMath::Log(FirstLog);

  Double_t SecondLog = sig*sig*(n-b) / (b* (b+sig*sig));
  SecondLog = TMath::Log( 1. + SecondLog);

  Double_t Rad = n*FirstLog - (b*b / (sig*sig)) * SecondLog;
  Rad = TMath::Sqrt(2.*Rad);

  return Rad;
  
}

Double_t GetUpp(Double_t n){

  if (n==0 || n==1) return 1.14;
  else return TMath::Sqrt(n);

}

  


std::vector<std::vector<double>> TwoDsignificance(Int_t dd0cut = 100, TString formula = "myZ", double addsigmabkg=0., Bool_t Draw = true, TString AnalysisResPath = "../MyExternalAnalysis/results/skimmed/", Int_t jalg = 2, TString analysis_opt="> d2d dmm anymass1L2M window [4,3]"){
  // formulas: atals simple signal myZ

  // opt: same as CutFlowOK.C
  
  // V[0][n] = n-th x coordinate of the curve
  // V[1][n] = n-th y coordinate of the curve


  std::vector<int> masses = {5, 10, 20, 30, 40, 50, 60, 70, 80, 85};
  std::vector<TString> mps = {"m","p"};

  //std::vector<TString> lifetime = {"m1p0", "m1p5", "m2p0", "m2p5", "m3p0", "m3p5", "m4p0", "m4p5", "m5p0", "m5p5", "m6p0", "m6p5", "m7p0"};

  TH2F* H = new TH2F("H","H",91,0-2.5,90+2.5,64,-12,-4);
  //TH2F* HBLACK = new TH2F("HBLACK","HBLACK",91,0-2.5,90+2.5,128,-12,-4);

  

  H->GetXaxis()->SetTitle("M(HN) (GeV)");
  H->GetYaxis()->SetTitle("Log(U^{2})");

  //X,Y,SIGNIFICANCE:
  std::vector<double> TDPlot_M;
  std::vector<double> TDPlot_U2;
  std::vector<double> TDPlot_Z;
  double TargetZ; // interpolating value, chosen accordingly to "formula"

  TDPlot_M.clear();
  TDPlot_U2.clear();
  TDPlot_Z.clear();


 
  std::map<int, std::vector<double>> bkgMapZmumu = CutFlowOK("Zmumu",-1,"n/a",AnalysisResPath,-1,false,jalg,analysis_opt);
  std::map<int, std::vector<double>> bkgMapZtautau = CutFlowOK("Ztautau",-1,"n/a",AnalysisResPath,-1,false,jalg,analysis_opt);
  std::map<int, std::vector<double>> bkgMapZbb = CutFlowOK("Zbb",-1,"n/a",AnalysisResPath,-1,false,jalg,analysis_opt);
  std::map<int, std::vector<double>> bkgMapZcc = CutFlowOK("Zcc",-1,"n/a",AnalysisResPath,-1,false,jalg,analysis_opt);
  std::map<int, std::vector<double>> bkgMapZuds = CutFlowOK("Zuds",-1,"n/a",AnalysisResPath,-1,false,jalg,analysis_opt);
  std::map<int, std::vector<double>> bkgMapmunuqq = CutFlowOK("munuqq",-1,"n/a",AnalysisResPath,-1,false,jalg,analysis_opt);

  ofstream LOG;
  LOG.open("LOG_TwoDsignificance.txt", std::ios_base::app);

  GetAvailableDatapoints(AnalysisResPath);


  for (int ip = 0; ip<AvailableDatapoints.size(); ip++){

    int m = AvailableDatapoints.at(ip).first;

    

    TString lt = AvailableDatapoints.at(ip).second;


    Int_t myid = dcut_id(dd0cut);
    Double_t signal = CutFlowOK("signal",m,lt,AnalysisResPath,-1,false,jalg,analysis_opt)[m][myid];

    Double_t U = Coupling(Form("%d",m),lt);

    Double_t Zmumu = bkgMapZmumu[m][myid];
    Double_t Ztautau = bkgMapZtautau[m][myid];
    Double_t Zbb = bkgMapZbb[m][myid];
    Double_t Zcc = bkgMapZcc[m][myid];
    Double_t Zuds = bkgMapZuds[m][myid];
    Double_t munuqq = bkgMapmunuqq[m][myid];

    Double_t SigmaBkg = TMath::Sqrt( TMath::Power(GetUpp(Zmumu)*Weight("Zmumu"),2) + TMath::Power(GetUpp(Ztautau)*Weight("Ztautau"),2) + TMath::Power(GetUpp(Zbb)*Weight("Zbb"),2) + TMath::Power(GetUpp(Zcc)*Weight("Zcc"),2) + TMath::Power(GetUpp(Zuds)*Weight("Zuds"),2) + TMath::Power(GetUpp(munuqq)*Weight("munuqq"),2));

    Double_t Y = 2.*TMath::Log10(U);
    // if (TMath::IsNaN(Y)) continue; // not needed --> should be protected by AvailableDatapoints
    Double_t X = 1.*m;
    
    Double_t totsig = signal*Weight("signal", Form("%d",m), lt);
    Double_t totbkg = Zmumu*Weight("Zmumu") + Ztautau * Weight("Ztautau") + Zbb * Weight("Zbb") + Zcc * Weight("Zcc") + Zuds * Weight("Zuds") +  munuqq * Weight("munuqq");
    //totbkg = WeightedLL[m];

    Double_t Z;
	  
    if (formula == "simple"){
      Z = totsig / TMath::Sqrt(totsig + totbkg + addsigmabkg*SigmaBkg);
      TargetZ = 2.;
    }
    
    else if (formula == "atlas"){
      Z = AtlasZ(totsig,totbkg + addsigmabkg*SigmaBkg,0);
      TargetZ = 2.;
    }
    
    else if (formula == "signal"){
      Z = totsig;
      TargetZ = 3.;
    }
    
    else if (formula == "myZ"){
      Z = 1.-ROOT::Math::poisson_cdf(int(totbkg + addsigmabkg*SigmaBkg), totbkg+totsig+addsigmabkg*SigmaBkg);
      TargetZ = 0.95;
    }
    
    else{
      cout<<"TwoDsignificance.C:: ERROR - FORMULA NOT RECOGNIZED"<<endl;
    }

    TDPlot_M.push_back(X);
    TDPlot_U2.push_back(Y);
    TDPlot_Z.push_back(Z);

  }


  const int npoints = TDPlot_Z.size();

  Double_t aX[npoints], aY[npoints], aZ[npoints];

  for (int i=0; i<npoints; i++){
    aX[i] = TDPlot_M.at(i);
    aY[i] = TDPlot_U2.at(i);
    aZ[i] = TDPlot_Z.at(i);
    //cout<<aX[i]<<"\t"<<aY[i]<<"\t"<<aZ[i]<<endl;
  }

  
  
  
  /*
  cout<<"TwoDsignificance.C:: interpolating..."<<endl;
  for (int i=1; i<H->GetNbinsX(); i++){
    cout<<i<<"/"<<H->GetNbinsX()<<endl;
    for (int j=1; j<H->GetNbinsY(); j++){
      cout<<".."<<j<<"/"<<H->GetNbinsY()<<endl;
      double ix = H->GetXaxis()->GetBinCenter(i);
      double iy = H->GetYaxis()->GetBinCenter(j);
      double iz = Gr->Interpolate(ix,iy);
      H->Fill(i,j,iz);
    }
  }
  */


  

	  
   

    
  

  

      
      

  // DRAWING THE LINE BY LINEAR INTERPOLATION

  for (int j=0; j < TDPlot_M.size(); j++) TDPlot_U2[j] = TMath::Power(10, TDPlot_U2[j]);

  std::vector<double> cx;
  std::vector<double> cy;

  
  for (int i : masses){

    std::vector<double> ycoordinates;

    for (int j=0; j < TDPlot_M.size(); j++){
      if (int(TDPlot_M[j]) == i) ycoordinates.push_back(TDPlot_U2[j]);
    }

    std::sort(ycoordinates.begin(), ycoordinates.end());
    

    double outvalue = -99;
    double spre = 0;
    double spost = 0;
    double upre =  outvalue;
    double upost = outvalue;
  
    for (double iY : ycoordinates){

      spre = spost;
      upre = upost;

      for (int j=0; j < TDPlot_M.size(); j++){
	if (int(TDPlot_M[j]) == i && TDPlot_U2[j] == iY) spost = TDPlot_Z[j];
      }
      upost = iY;

      if (spre < TargetZ && spost >= TargetZ && upre > outvalue && upost > outvalue){

	cx.push_back(i);

	double interpol = upre + (TargetZ-spre)*((upost-upre)/(spost-spre));

	cy.push_back(TMath::Log10(interpol));

	LOG<<"dcut="<<dd0cut<<") "<<i<<" "<<upre<<" "<<upost<<" "<<interpol<<endl;
	LOG<<"========================================="<<endl;

      }     
    }
  }

  std::vector<std::vector<double>> toret;
  toret.push_back(cx);
  toret.push_back(cy);

  LOG.close();

  for (int j=0; j<cx.size(); j++) cout<<"Interpolation result: ---- M "<<cx[j]<<" -> "<<cy[j]<<endl;

  // Drawing part
  

  if (Draw){

    TGraph2D *Gr = new TGraph2D(npoints, aX, aY, aZ);

  auto c = new TCanvas();
  //Gr->SetNpx(20);
  //Gr->SetNpy(60);
  Gr->Draw("colz");
  //double conts[] = {0.5, 0.95, 1.};
  //Gr->SetContour(3,conts);
  //c->Update();

  

    //auto c = new TCanvas();

    const Int_t npoint = cx.size();
    double vcx[npoint], vcy[npoint];
    for (int i = 0; i<npoint; i++) {vcx[i] = cx[i]; vcy[i] = cy[i];}
     
    auto gg = new TGraph(npoint,vcx,vcy);

    if (TargetZ == 2.)
      H->GetZaxis()->SetRangeUser(1e-5,1000);
    else if (TargetZ  == 0.05)
      H->GetZaxis()->SetRangeUser(0,1);
    //H->Draw("colz");
    H->SetMarkerSize(1.6);
    H->SetTitle("Significance");

    //HBLACK->Draw("same BOX");

    gg->Draw("same L");

    gStyle->SetOptStat(0);
    gPad->SetLogz();
    //gStyle->SetPalette(kBlackBody);

    c->SaveAs("temp.png");
  }
  

  
  return toret;
}


void CompareAnalyses(){

  auto c = new TCanvas("c","c",1000,600);

  std::vector<std::vector<double>> XY;

  Int_t dcut = 8;
  TString AnalysisOptions = "> d3d dsigma";

  int colcounter = 1;

  TString formula = "atlas";
  
  Double_t x0[50], y0[50];
  XY = TwoDsignificance(dcut,formula,0.,false,"../MyExternalAnalysis/results/skimmed/",2,"< d2d dsigma anymass1L2M window [1.5,0.25]");
  //XY = TwoDsignificance(100,"atlas",0.,false,"../MyExternalAnalysis/results-V230322/skimmed_loose/",2,"> d2d dmm anymass1L2M");
  for (int i=0; i<XY[0].size(); i++){
    x0[i] = XY[0][i];
    y0[i] = XY[1][i];
  }
  auto g0 = new TGraph(XY[0].size(), x0, y0);
  g0->SetLineColor(colcounter);
  g0->SetTitle("1.5* 25% sqrt(M)");
  g0->SetLineWidth(2);
  colcounter++;
  

  Double_t x1[50], y1[50];
  XY = TwoDsignificance(dcut,formula,0.,false,"../MyExternalAnalysis/results/skimmed/",2,"< d2d dsigma anymass1L2M window [1.5,0.25,1.5,0.35]");
  //XY = TwoDsignificance(100,"atlas",0.,false,"../MyExternalAnalysis/results/",2,"> d2d dmm anymass1L2M");
  for (int i=0; i<XY[0].size(); i++){
    x1[i] = XY[0][i];
    y1[i] = XY[1][i];
  }
  auto g1 = new TGraph(XY[0].size(), x1, y1);
  g1->SetLineColor(colcounter);
  g1->SetLineStyle(7);
  g1->SetTitle("bkg: 35% sqrt(M)");
  g1->SetLineWidth(2);
  colcounter++;


  Double_t x2[50], y2[50];
  XY = TwoDsignificance(dcut,formula,0.,false,"../MyExternalAnalysis/results/skimmed/",2,"< d2d dsigma anymass1L2M window [1.5, 0.25, 1.5, 0.20]");
  //XY = TwoDsignificance(8,"atlas",0.,false,"../MyExternalAnalysis/results-V230322/skimmed_loose/",2,"< d2d dsigma anymass1L2M");
  for (int i=0; i<XY[0].size(); i++){
    x2[i] = XY[0][i];
    y2[i] = XY[1][i];
  }
  auto g2 = new TGraph(XY[0].size(), x2, y2);
  g2->SetLineColor(colcounter);
  g2->SetTitle("bkg: 20% sqrt(M)");
  g2->SetLineWidth(2);
  colcounter++;

  /*
  Double_t x3[50], y3[50];
  XY = TwoDsignificance(8,"atlas",0.,false,"../MyExternalAnalysis/results/",2,"< d2d dsigma anymass1L2M");
  for (int i=0; i<XY[0].size(); i++){
    x3[i] = XY[0][i];
    y3[i] = XY[1][i];
  }
  auto g3 = new TGraph(XY[0].size(), x3, y3);
  g3->SetLineColor(colcounter);
  g3->SetLineStyle(7);
  g3->SetTitle("D0#mu < 8#sigma, sig+irreduc beam spread");
  g3->SetLineWidth(2);
  colcounter++;
  


  Double_t x4[50], y4[50];
  //XY = TwoDsignificance(dcut,"simple",0.,false,"../MyExternalAnalysis/results/skimmed/",2,AnalysisOptions);
  XY = TwoDsignificance(16,"atlas",0.,false,"../MyExternalAnalysis/results-V230322/skimmed_loose/",2,"< d2d dsigma anymass1L2M");
  for (int i=0; i<XY[0].size(); i++){
    x4[i] = XY[0][i];
    y4[i] = XY[1][i];
  }
  auto g4 = new TGraph(XY[0].size(), x4, y4);
  g4->SetLineColor(colcounter);
  g4->SetTitle("D0#mu <16#sigma, no beam spread");
  g4->SetLineWidth(2);
  //colcounter++;


  Double_t x5[50], y5[50];
  XY = TwoDsignificance(16,"atlas",0.,false,"../MyExternalAnalysis/results/",2,"< d2d dsigma anymass1L2M");
  for (int i=0; i<XY[0].size(); i++){
    x5[i] = XY[0][i];
    y5[i] = XY[1][i];
  }
  auto g5 = new TGraph(XY[0].size(), x3, y3);
  g5->SetLineColor(colcounter);
  g5->SetLineStyle(7);
  g5->SetTitle("D0#mu < 16#sigma, sig+irreduc beam spread");
  g5->SetLineWidth(2);
  colcounter++;


  */
  


  auto gax = (TGraph*)g0->Clone("ciao");
  gax->SetMarkerColor(0);
  gax->SetLineColor(0);
  gax->SetTitle("Curve at significance #approx 2");
  gax->GetYaxis()->SetRangeUser(-12,-5);
  gax->GetYaxis()->SetTitle("Log (U^{2})");
  gax->GetXaxis()->SetLimits(0,90);
  gax->GetXaxis()->SetTitle("M_{HN} (GeV)");
  
  gax->Draw("A C");


  
  
  g0->Draw("same");
  g1->Draw("same");
  g2->Draw("same");
  //g3->Draw("same");
  //g4->Draw("same");
  //g5->Draw("same");
 

  

  c->BuildLegend();

  c->SetGridy();
  c->SaveAs("temp.png");


  
}
  


void ScanD0Cut(){

  auto c = new TCanvas("c","c",1000,600);

  Int_t dcut;

  Int_t jalg = 2;

  std::vector<std::vector<double>> XY;
  int colcounter=1;

  TString AnalysisOptions = "< d3d dsigma anymass1L2M";

  dcut = 4;
  double_t x4[50], y4[50];
  XY = TwoDsignificance(dcut,"atlas",0.,false,"../MyExternalAnalysis/results/skimmed/",jalg,AnalysisOptions);
  for (int i=0; i<XY[0].size(); i++){
    x4[i] = XY[0][i];
    y4[i] = XY[1][i];
  }
  auto g4 = new TGraph(XY[0].size(), x4, y4);
  g4->SetLineColor(colcounter);
  g4->SetTitle(Form("Impact par < %d #sigma",dcut));
  g4->SetLineWidth(2);
  colcounter++;

 
  dcut = 8;
  double_t x8[50], y8[50];
  XY = TwoDsignificance(dcut,"atlas",0.,false,"../MyExternalAnalysis/results/skimmed/",jalg,AnalysisOptions);
  for (int i=0; i<XY[0].size(); i++){
    x8[i] = XY[0][i];
    y8[i] = XY[1][i];
  }
  auto g8 = new TGraph(XY[0].size(), x8, y8);
  g8->SetLineColor(colcounter);
  g8->SetTitle(Form("Impact par < %d #sigma",dcut));
  g8->SetLineWidth(2);
  colcounter++;

 
  dcut = 30;
  double_t x12[50], y12[50];
  XY = TwoDsignificance(dcut,"atlas",0.,false,"../MyExternalAnalysis/results/skimmed/",jalg,AnalysisOptions);
  for (int i=0; i<XY[0].size(); i++){
    x12[i] = XY[0][i];
    y12[i] = XY[1][i];
  }
  auto g12 = new TGraph(XY[0].size(), x12, y12);
  g12->SetLineColor(colcounter);
  g12->SetTitle(Form("Impact par < %d #sigma",dcut));
  g12->SetLineWidth(2);
  colcounter++;


  dcut = 100;
  double_t x16[50], y16[50];
  XY = TwoDsignificance(dcut,"atlas",0.,false,"../MyExternalAnalysis/results/skimmed/",jalg,AnalysisOptions);
  for (int i=0; i<XY[0].size(); i++){
    x16[i] = XY[0][i];
    y16[i] = XY[1][i];
  }
  auto g16 = new TGraph(XY[0].size(), x16, y16);
  g16->SetLineColor(colcounter);
  g16->SetTitle(Form("Impact par < %d #sigma",dcut));
  g16->SetLineWidth(2);
  colcounter++;

  
  dcut = 200;
  double_t x20[50], y20[50];
  XY = TwoDsignificance(dcut,"atlas",0.,false,"../MyExternalAnalysis/results/skimmed/",jalg,AnalysisOptions);
  for (int i=0; i<XY[0].size(); i++){
    x20[i] = XY[0][i];
    y20[i] = XY[1][i];
  }
  auto g20 = new TGraph(XY[0].size(), x20, y20);
  g20->SetLineColor(colcounter+1);
  g20->SetTitle(Form("Impact par < %d #sigma",dcut));
  g20->SetLineWidth(2);
  colcounter++;


  

  auto gax = (TGraph*)g4->Clone("ciao");
  gax->SetMarkerColor(0);
  gax->SetLineColor(0);
  gax->SetTitle("Curve at significance #approx 2");
  gax->GetYaxis()->SetRangeUser(-14,-4);
  gax->GetYaxis()->SetTitle("Log (U^{2})");
  gax->GetXaxis()->SetLimits(0,90);
  gax->GetXaxis()->SetTitle("M_{HN} (GeV)");
  
  gax->Draw("A C");


  
  
  g4->Draw("same");
  g8->Draw("same");
  g12->Draw("same");
  g16->Draw("same");
  g20->Draw("same");


  
  
  c->BuildLegend();

  c->SetGridy();
  c->SaveAs("temp.png");
}

