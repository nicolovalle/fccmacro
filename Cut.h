// BRANCHES of the eventsTree

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

  return sel1 && sel2 && sel3 && sel41 && sel42 && sel6 && s21 && s22 && s23;
 
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

  return sel1 && sel2 && sel3 && sel4 && sel5;

}


Bool_t SELECTION_LM_2JET(int jalg){

  Bool_t sel1 = SELECTION_LM_1JET(jalg);

  Double_t cosjj = TMath::Cos(lvj1.Angle(lvj2.Vect()));
  Bool_t sel2 = (cosjj > -0.8);

  return sel1 && sel2;
}
  

  






