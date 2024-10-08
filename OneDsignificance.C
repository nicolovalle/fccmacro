#include "CutFlowOK.C"





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

  


std::vector<std::vector<double>> OneDsignificance(Int_t m = 5, Int_t dd0cut = 8, TString formula = "myZ", double addsigmabkg=0., Bool_t Draw = true, TString AnalysisResPath = "../MyExternalAnalysis/results/skimmed/", Int_t jalg = 2, TString analysis_opt="< d2d dsigma anymass1L2M window [2,0.2]"){
  // formulas: atals simple

  // opt: same as CutFlowOK.C
  
  // V[0][n] = n-th x coordinate of the curve
  // V[1][n] = n-th y coordinate of the curve


  
  std::vector<TString> mps = {"m","p"};

  //std::vector<TString> lifetime = {"m1p0", "m1p5", "m2p0", "m2p5", "m3p0", "m3p5", "m4p0", "m4p5", "m5p0", "m5p5", "m6p0", "m6p5", "m7p0"};

  

  //LOG(COUPLING^2), SIGNIFICANCE
  std::vector<double> TR_U2;
  std::vector<double> TR_Z;

  TR_U2.clear();
  TR_Z.clear();
  


 
  std::map<int, std::vector<double>> bkgMapZmumu = CutFlowOK("Zmumu",m,"n/a",AnalysisResPath,-1,false,jalg,analysis_opt,true);
  std::map<int, std::vector<double>> bkgMapZtautau = CutFlowOK("Ztautau",m,"n/a",AnalysisResPath,-1,false,jalg,analysis_opt,true);
  std::map<int, std::vector<double>> bkgMapZbb = CutFlowOK("Zbb",m,"n/a",AnalysisResPath,-1,false,jalg,analysis_opt,true);
  std::map<int, std::vector<double>> bkgMapZcc = CutFlowOK("Zcc",m,"n/a",AnalysisResPath,-1,false,jalg,analysis_opt,true);
  std::map<int, std::vector<double>> bkgMapZuds = CutFlowOK("Zuds",m,"n/a",AnalysisResPath,-1,false,jalg,analysis_opt,true);
  std::map<int, std::vector<double>> bkgMapmunuqq = CutFlowOK("munuqq",m,"n/a",AnalysisResPath,-1,false,jalg,analysis_opt,true);

  GetAvailableDatapoints(AnalysisResPath);
 
    
    
	for (int ip = 0; ip<AvailableDatapoints.size(); ip++){

	  if (AvailableDatapoints.at(ip).first != m) continue;
	  TString lt = AvailableDatapoints.at(ip).second;


     
	  Int_t myid = dcut_id(dd0cut);


	  Double_t signal = CutFlowOK("signal",m,lt,AnalysisResPath,-1,false,jalg,analysis_opt,true)[m][myid];
	  

	  Double_t U = Coupling(Form("%d",m),lt);
	  
	  if (signal == 0){
	    cout<<"Signal m="<<m<<" lt="<<lt<<" has 0 events"<<endl;
	    //continue;
	  }

	
     
	  
	  Double_t Zmumu = bkgMapZmumu[m][myid];
	  Double_t Ztautau = bkgMapZtautau[m][myid];
	  Double_t Zbb = bkgMapZbb[m][myid];
	  Double_t Zcc = bkgMapZcc[m][myid];
	  Double_t Zuds = bkgMapZuds[m][myid];
	  Double_t munuqq = bkgMapmunuqq[m][myid];
	  

	  Double_t SigmaBkg = TMath::Sqrt( TMath::Power(GetUpp(Zmumu)*Weight("Zmumu"),2) + TMath::Power(GetUpp(Ztautau)*Weight("Ztautau"),2) + TMath::Power(GetUpp(Zbb)*Weight("Zbb"),2) + TMath::Power(GetUpp(Zcc)*Weight("Zcc"),2) + TMath::Power(GetUpp(Zuds)*Weight("Zuds"),2) + TMath::Power(GetUpp(munuqq)*Weight("munuqq"),2));

	 
      
	  Double_t Y = 2.*TMath::Log10(U);
	  Double_t X = 1.*m;

	  Double_t totsig = signal*Weight("signal", Form("%d",m), lt);
	  Double_t totbkg = Zmumu*Weight("Zmumu") + Ztautau * Weight("Ztautau") + Zbb * Weight("Zbb") + Zcc * Weight("Zcc") + Zuds * Weight("Zuds") +  munuqq * Weight("munuqq");

	  Double_t Z, TargetZ;

	  if (formula == "simple"){
      Z = totsig / TMath::Sqrt(totsig + totbkg + addsigmabkg*SigmaBkg);
      TargetZ = 2.;
    }
    
    else if (formula == "atlas"){
      if (totsig == 0) Z = 0;
      else Z = AtlasZ(totsig,totbkg + addsigmabkg*SigmaBkg,0);
      TargetZ = 2.;
    }
    
    else if (formula == "signal"){
      Z = totsig;
      TargetZ = 3.;
    }
    
    else if (formula == "myCL"){
      Z = 1.-ROOT::Math::poisson_cdf(int(totbkg + addsigmabkg*SigmaBkg), totbkg+totsig+addsigmabkg*SigmaBkg);
      TargetZ = 1.-0.05;
    }

    else if (formula == "myZ"){
      
      double alpha = ROOT::Math::poisson_cdf(int(totbkg + addsigmabkg*SigmaBkg), totbkg+totsig+addsigmabkg*SigmaBkg);
      if (alpha < 1.e-20) Z = 0.-ROOT::Math::gaussian_quantile(1.e-20/2.,1.);
      else Z = 0.-ROOT::Math::gaussian_quantile(alpha/2.,1.);
      if (Z<0) Z = 0;
      if (totsig == 0) Z = 0;
      TargetZ = 2.;
    }
    
    else{
      cout<<"TwoDsignificance.C:: ERROR - FORMULA NOT RECOGNIZED"<<endl;
    }

	  cout<<"M "<<m<<"     LT "<<lt<<"    logU2 "<<Y<<"    sig "<<totsig<<"    bkg "<<totbkg<<"   Errbkg "<<SigmaBkg<<"    Z "<<Z<<endl;
	  //LOG<<"-------------------------------------------------"<<endl;
	  cout<<"Events: "<<signal<<" "<<Zbb<<" "<<Zcc<<" "<<Zuds<<" "<<Zmumu<<" "<<Ztautau<<" "<<munuqq<<endl;
	  //LOG<<"-------------------------------------------------"<<endl;
      

	  TR_U2.push_back(Y);
	  TR_Z.push_back(Z);

	 

	}
  
 

  if (Draw){

    auto c = new TCanvas();

    const Int_t npoint = TR_U2.size();
    double vcx[npoint], vcy[npoint];
    for (int i = 0; i<npoint; i++) {vcx[i] = TR_U2[i]; vcy[i] = TR_Z[i];}
     
    auto gg = new TGraph(npoint,vcx,vcy);

    gg->SetLineColor(2);

    gg->Draw("A L");
    gg->GetXaxis()->SetTitle("Log(U^{2})");
    gg->GetYaxis()->SetTitle("Significance");

    gStyle->SetOptStat(0);
    //gPad->SetLogy();
    gPad->SetGridy();
    gPad->SetGridx();
    

    c->SaveAs("temp.png");
  }

  
  return std::vector<std::vector<double>>{TR_U2, TR_Z};
}

void CutScan(){

  auto c = new TCanvas("c","c",1000,600);

  auto hs = new THStack("hs","");

  std::vector<std::vector<double>> VV;

  std::vector<TGraph*> grv;

  double x[1000][50], y[1000][50];

  int icut = 0;
  int npoint;
  TString title;
  std::vector<TString> titv;
  TGraph* gr;

  std::vector<double> cjj_scan = {-0.8, -0.88, -0.96, -0.75, -0.70};
  std::vector<double> mjm_scan = {-0.98, -0.85, -1.00};
  std::vector<double> Mjm_scan = {0.8, 0.5, 0.96};

  for (double cjj : cjj_scan){
    for (double mjm : mjm_scan){
      for (double Mjm : Mjm_scan){

	TString opp = Form("[%1.2f,%1.2f,%1.2f]",cjj,mjm,Mjm);
	title = opp;

	cout<<"RUNNING CUT VARIATION: "<<opp<<endl;
	VV = OneDsignificance(80, 8,"atlas", 0., false, "../MyExternalAnalysis/results/skimmed/", 2, " < d2d dsigma cutvariation "+opp);
  
	npoint = VV[0].size();
	for (int i=0;i<npoint;i++) {x[icut][i] = VV[0][i]; y[icut][i]=VV[1][i];}
	gr = new TGraph(npoint,x[icut],y[icut]);
	grv.push_back((TGraph*)gr->Clone(title));
	titv.push_back(title);
	icut++;
	
      }
    }
  }


  Double_t xlimit[2]={-12,-3};
  Double_t ylimit[2]={1e-3,1e2};
  
  grv.at(0)->GetXaxis()->SetLimits(xlimit[0],xlimit[1]);
  grv.at(0)->GetYaxis()->SetRangeUser(ylimit[0],ylimit[1]);
  grv.at(0)->SetLineWidth(2);
  grv.at(0)->SetLineColor(1);
  grv.at(0)->SetTitle(titv.at(0));
  grv.at(0)->Draw("A C");
  

  
  for (int i=1; i<grv.size(); i++){
    grv.at(i)->SetLineColor(i+1);
    grv.at(i)->SetTitle(titv.at(i));
    grv.at(i)->Draw("same C");
  }


  Double_t line[2][50];
  for (int i=0; i<50; i++){line[0][i]=xlimit[0]+1.*i*(xlimit[1]-xlimit[0])/50.; line[1][i]=2.; }
  auto grl = new TGraph(50,line[0],line[1]);

  grl->SetLineColor(1);
  grl->SetLineStyle(9);
  grl->Draw("same C");

  
  gPad->SetLogy();
  c->BuildLegend();
  c->SaveAs("temp.png");
}
  
  

  
/*
void CutScan(){

  auto c = new TCanvas("c","c",1000,600);

  std::vector<std::vector<double>> XY;

  Int_t dcut = 8;
  TString AnalysisOptions = "> d3d dsigma jalg";

  int colcounter = 1;
  
  Double_t x0[50], y0[50];
  //XY = TwoDsignificance(dcut,"simple",0.,false,"../MyExternalAnalysis/results/skimmed/",0,AnalysisOptions);
  XY = TwoDsignificance(100,"simple",0.,false,"../MyExternalAnalysis/results/skimmed_extraloose/",2,"< d2d dmm anymass1L2M");
  for (int i=0; i<XY[0].size(); i++){
    x0[i] = XY[0][i];
    y0[i] = XY[1][i];
  }
  auto g0 = new TGraph(XY[0].size(), x0, y0);
  g0->SetLineColor(colcounter);
  g0->SetTitle("D0#mu < 1mm");
  g0->SetLineWidth(2);
  colcounter++;
  

  Double_t x1[50], y1[50];
  //XY = TwoDsignificance(dcut,"simple",0.,false,"../MyExternalAnalysis/results/skimmed/",2,AnalysisOptions);
  XY = TwoDsignificance(100,"simple",0.,false,"../MyExternalAnalysis/results/skimmed_extraloose/",2,"> d2d dmm anymass1L2M");
  for (int i=0; i<XY[0].size(); i++){
    x1[i] = XY[0][i];
    y1[i] = XY[1][i];
  }
  auto g1 = new TGraph(XY[0].size(), x1, y1);
  g1->SetLineColor(colcounter);
  g1->SetTitle("D0#mu > 1mm");
  g1->SetLineWidth(2);
  colcounter++;


  Double_t x2[50], y2[50];
  //XY = TwoDsignificance(dcut,"simple",0.,false,"../MyExternalAnalysis/results/skimmed/",2,AnalysisOptions);
  XY = TwoDsignificance(0,"simple",0.,false,"../MyExternalAnalysis/results/skimmed_extraloose/",2,"> d2d dmm anymass1L2M");
  for (int i=0; i<XY[0].size(); i++){
    x2[i] = XY[0][i];
    y2[i] = XY[1][i];
  }
  auto g2 = new TGraph(XY[0].size(), x2, y2);
  g2->SetLineColor(colcounter);
  g2->SetTitle("no D0#mu cut");
  g2->SetLineWidth(2);
  colcounter++;


  auto gax = (TGraph*)g0->Clone("ciao");
  gax->SetMarkerColor(0);
  gax->SetLineColor(0);
  gax->SetTitle("Curve at significance #approx 2");
  gax->GetYaxis()->SetRangeUser(-12,-6);
  gax->GetYaxis()->SetTitle("Log (U^{2})");
  gax->GetXaxis()->SetLimits(0,90);
  gax->GetXaxis()->SetTitle("M_{HN} (GeV)");
  
  gax->Draw("A C");


  
  
  g0->Draw("same");
  g1->Draw("same");
  g2->Draw("same");
 

  

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
  gax->GetYaxis()->SetRangeUser(-15,-6);
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

*/
