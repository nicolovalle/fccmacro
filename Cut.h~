// BRANCHES of the eventsTree

TChain *TREE = new TChain("eventsTree");

Long64_t oNEntries;
Int_t oNReco;
Int_t oNMuon;
std::vector<Int_t>* oNJet = 0;
std::vector<Double_t>* oPxJet1 = 0;
std::vector<Double_t>* oPyJet1 = 0;
std::vector<Double_t>* oPzJet1 = 0;
std::vector<Double_t>* oEJet1 = 0;
std::vector<Double_t>* oPxJet2 = 0;
std::vector<Double_t>* oPyJet2 = 0;
std::vector<Double_t>* oPzJet2 = 0;
std::vector<Double_t>* oEJet2 = 0;
Double_t oPxMu;
Double_t oPyMu;
Double_t oPzMu;
Double_t oEMu;
Double_t oPxMiss;
Double_t oPyMiss;
Double_t oPzMiss;
Double_t oEMiss;
Double_t oPxVis;
Double_t oPyVis;
Double_t oPzVis;
Double_t oEVis;
Double_t oMuD0;
Double_t oMuZ0;
Double_t oMuD0sig;
Double_t oMuZ0sig;
Double_t oCosPmiss;
Double_t oCosPmissMu;
Double_t oCosjj;
Double_t oCosj1Mu;
Double_t oCosj2Mu;
Double_t oMjjMu;
Double_t oMVis;
Double_t oMjjMuPmiss;
Double_t oNtracks;
Double_t oVtx_x;
Double_t oVtx_y;
Double_t oVtx_z;
Double_t oVtx_chi2;
Double_t oVtx_ntrk;


void TREESetBranch(){
  TREE->SetBranchAddress("NEntries_tchain",&oNEntries);
  TREE->SetBranchAddress("NReco",&oNReco);
  TREE->SetBranchAddress("NMuon",&oNMuon);
  TREE->SetBranchAddress("NJet",&oNJet);
  TREE->SetBranchAddress("PxJet1",&oPxJet1);
  TREE->SetBranchAddress("PyJet1",&oPyJet1);
  TREE->SetBranchAddress("PzJet1",&oPzJet1);
  TREE->SetBranchAddress("EJet1",&oEJet1);
  TREE->SetBranchAddress("PxJet2",&oPxJet2);
  TREE->SetBranchAddress("PyJet2",&oPyJet2);
  TREE->SetBranchAddress("PzJet2",&oPzJet2);
  TREE->SetBranchAddress("EJet2",&oEJet2);
  TREE->SetBranchAddress("PxMu",&oPxMu);
  TREE->SetBranchAddress("PyMu",&oPyMu);
  TREE->SetBranchAddress("PzMu",&oPzMu);
  TREE->SetBranchAddress("EMu",&oEMu);
  TREE->SetBranchAddress("PxMiss",&oPxMiss);
  TREE->SetBranchAddress("PyMiss",&oPyMiss);
  TREE->SetBranchAddress("PzMiss",&oPzMiss);
  TREE->SetBranchAddress("EMiss",&oEMiss);
  TREE->SetBranchAddress("PxVis",&oPxVis);
  TREE->SetBranchAddress("PyVis",&oPyVis);
  TREE->SetBranchAddress("PzVis",&oPzVis);
  TREE->SetBranchAddress("EVis",&oEVis);
  TREE->SetBranchAddress("MuD0",&oMuD0);
  TREE->SetBranchAddress("MuZ0",&oMuZ0);
  TREE->SetBranchAddress("MuD0sig",&oMuD0sig);
  TREE->SetBranchAddress("MuZ0sig",&oMuZ0sig);
  TREE->SetBranchAddress("CosPmiss",&oCosPmiss);
  TREE->SetBranchAddress("CosPmissMu",&oCosPmissMu);
  TREE->SetBranchAddress("Cosjj",&oCosjj);
  TREE->SetBranchAddress("Cosj1Mu",&oCosj1Mu);
  TREE->SetBranchAddress("Cosj2Mu",&oCosj2Mu);
  TREE->SetBranchAddress("MjjMu",&oMjjMu);
  TREE->SetBranchAddress("MVis",&oMVis);
  TREE->SetBranchAddress("MjjMuPmiss",&oMjjMuPmiss);
  TREE->SetBranchAddress("Ntracks",&oNtracks);
  TREE->SetBranchAddress("Vtx_x",&oVtx_x);
  TREE->SetBranchAddress("Vtx_y",&oVtx_y);
  TREE->SetBranchAddress("Vtx_z",&oVtx_z);
  TREE->SetBranchAddress("Vtx_chi2",&oVtx_chi2);
  TREE->SetBranchAddress("Vtx_ntrk",&oVtx_ntrk);
}

// Variable filled at each Fill by BUILD_DERIVATE

TLorentzVector lvj1;
TLorentzVector lvj2;
TLorentzVector lvmu;
TLorentzVector lvmiss;
TLorentzVector lvvis;

void BUILD_DERIVATE(int jalg){


   lvj1 = TLorentzVector(oPxJet1->at(jalg), oPyJet1->at(jalg), oPzJet1->at(jalg), oEJet1->at(jalg));
   lvj2 = TLorentzVector(oPxJet2->at(jalg), oPyJet2->at(jalg), oPzJet2->at(jalg), oEJet2->at(jalg));
   lvmu = TLorentzVector(oPxMu, oPyMu, oPzMu, oEMu);
   lvmiss = TLorentzVector(oPxMiss, oPyMiss, oPzMiss, oEMiss);
   lvvis = TLorentzVector(oPxVis, oPyVis, oPzVis, oEVis);
  
}


// Medium Mass: the one used at the beginning for any mass
Bool_t SELECTION_MM_2JET(int jalg){


  Double_t cosjj = TMath::Cos(lvj1.Angle(lvj2.Vect()));
  Double_t cosj1mu = TMath::Cos(lvj1.Angle(lvmu.Vect()));
  Double_t cosj2mu = TMath::Cos(lvj2.Angle(lvmu.Vect()));
    
  Bool_t sel1 = (TMath::Abs(lvmiss.CosTheta()) < 0.94);
  Bool_t sel2 = (TMath::Cos(lvmiss.Angle(lvmu.Vect())) < 0.80);
  Bool_t sel3 = (oEJet1->at(jalg) >= 3. && oEJet2->at(jalg) >= 3.);
  Bool_t sel41 = cosjj > -0.8;
  Bool_t sel42 = cosjj < 0.98;
  //Bool_t sel5_old = (TMath::Min(lvj1.M(), lvj1.M()) > 1.8 && TMath::Min(lvj1.M2(), lvj1.M2()) > 0);
  Bool_t sel6 = (TMath::Max( cosj1mu, cosj2mu ) < 0.8);

  Bool_t s21 = (TMath::Min(lvj1.M(), lvj2.M()) > 0.2 && TMath::Min(lvj1.M2(), lvj2.M2()) > 0);
  Bool_t s22 = (TMath::Min( cosj1mu, cosj2mu ) > -0.98);
  Bool_t s23 = ((lvj1 + lvj2 + lvmiss + lvmu).M() > 80.);

  Bool_t r1 = (oNReco > 0);

  return sel1 && sel2 && sel3 && sel41 && sel42 && sel6 && s21 && s22 && s23 && r1;
 
}

// Low mass 1 jet
Bool_t SELECTION_LM_1JET(int jalg){

  // they are exactly the same in case of 1 jet
  Double_t cosjmu1 = TMath::Cos(lvj1.Angle(lvmu.Vect()));
  Double_t cosjmu2 = TMath::Cos(lvj2.Angle(lvmu.Vect()));
  Double_t cosjmu = TMath::Max(cosjmu1, cosjmu2);
      
  Bool_t sel1 = (TMath::Abs(lvmiss.CosTheta()) < 0.94);
  Bool_t sel2 = (TMath::Cos(lvmiss.Angle(lvmu.Vect())) < 0.50);
  Bool_t sel3 = (oEJet1->at(jalg) >= 3.);
   
  Bool_t sel4 = (-0.5 < cosjmu && cosjmu < 0.96);
  Bool_t sel5 = ((lvvis+lvmiss).M() > 80);

  Bool_t r1 = (oNReco > 0);

  return sel1 && sel2 && sel3 && sel4 && sel5 && r1;

}


Bool_t SELECTION_LM_2JET(int jalg){

  Bool_t sel1 = SELECTION_LM_1JET(jalg);

  Double_t cosjj = TMath::Cos(lvj1.Angle(lvj2.Vect()));
  Bool_t sel2 = (cosjj > -0.8);

  return sel1 && sel2;
}



// SOME TEMPORARY CUT VARIATION

Bool_t SEL_MM_2J(int jalg, double_t CUTcosjj, double CUTmincosjmu, double CUTmaxcosjmu){


  Double_t cosjj = TMath::Cos(lvj1.Angle(lvj2.Vect()));
  Double_t cosj1mu = TMath::Cos(lvj1.Angle(lvmu.Vect()));
  Double_t cosj2mu = TMath::Cos(lvj2.Angle(lvmu.Vect()));
    
  Bool_t sel1 = (TMath::Abs(lvmiss.CosTheta()) < 0.94);
  Bool_t sel2 = (TMath::Cos(lvmiss.Angle(lvmu.Vect())) < 0.80);
  Bool_t sel3 = (oEJet1->at(jalg) >= 3. && oEJet2->at(jalg) >= 3.);
  Bool_t sel41 = cosjj > CUTcosjj;
  Bool_t sel42 = cosjj < 0.98;
  //Bool_t sel5_old = (TMath::Min(lvj1.M(), lvj1.M()) > 1.8 && TMath::Min(lvj1.M2(), lvj1.M2()) > 0);
  Bool_t sel6 = (TMath::Max( cosj1mu, cosj2mu ) < CUTmaxcosjmu);

  Bool_t s21 = (TMath::Min(lvj1.M(), lvj2.M()) > 0.2 && TMath::Min(lvj1.M2(), lvj2.M2()) > 0);
  Bool_t s22 = (TMath::Min( cosj1mu, cosj2mu ) > CUTmincosjmu);
  Bool_t s23 = ((lvj1 + lvj2 + lvmiss + lvmu).M() > 80.);

  return sel1 && sel2 && sel3 && sel41 && sel42 && sel6 && s21 && s22 && s23;
 
 
}


std::vector<TString> ExplodeString(TString s, TString token){

  std::vector<TString> toret;
  int index=0;
  TString s2 = s;

  while (s2.Index(token)>0){
    index = s2.Index(token);
    TString is = s2(0,index);
    toret.push_back(is);
    s2 = (TString)(s2(index+1,s2.Length()));
  }
  TString is = s2;
  toret.push_back(s2);

  return toret;
  
  
}

std::vector<float> GetFloatArray(TString opt){
  int i1 = opt.Index("[");
  int i2 = opt.Index("]");
  if (i1<0 || i2<0 ) {cout<<"Cut.h:: ERROR - BAD [] OPTION DECODING"<<endl; return std::vector<float>{};}

  TString s = opt(i1+1, i2-i1-1);

  std::vector<TString> V = ExplodeString(s,",");
  
  std::vector<float> toret;
  for (TString iopt : V){
    toret.push_back(std::stof(iopt.Data()));
  }

  return toret;
  
		  
}
  

  






