
#include "lumisettings.h"

Double_t Muon_mass = 0.105658;
Double_t Z_mass = 91.1876;



std::pair<std::vector<Double_t>, Double_t> getvalues(TString obsID, TString opt="signal", Long64_t RunOnN = -1, Int_t mass=50, Double_t d0cut=8, TString dir="../MyExternalAnalysis/results/"){
  // !! With the last option = true, the getvalues[0] is the scale factor!!
  // available obsID: mtot, d0[selection], d0sliding[selection]
  // If sliding is not needed, you can use any mass.


  TString Dir = dir;

  if (obsID == "d0" || obsID == "d0sliding")
    Dir = Dir+"/skimmed/";

  TString fname=Form("%s%s",Dir.Data(), AnalysisResults(opt,Form("%d",mass)).Data());

  cout<<"getvalues.C:: Opening file "<<fname<<endl;

  std::vector<Double_t> toret;
  toret.clear();
  
  TFile *File = new TFile(fname);
  TTree *T = (TTree*)File->Get("eventsTree");


  Long64_t oNEntries;
  Int_t oNReco;
  Int_t oNMuon;
  Int_t oNJet;
  Double_t oPxJet1;
  Double_t oPyJet1;
  Double_t oPzJet1;
  Double_t oEJet1;
  Double_t oPxJet2;
  Double_t oPyJet2;
  Double_t oPzJet2;
  Double_t oEJet2;
  Double_t oPxMu;
  Double_t oPyMu;
  Double_t oPzMu;
  Double_t oEMu;
  Double_t oPxMiss;
  Double_t oPyMiss;
  Double_t oPzMiss;
  Double_t oEMiss;
  Double_t oMuD0;
  Double_t oMuZ0;
  Double_t oMuD0sig;
  Double_t oMuZ0sig;


  T->SetBranchAddress("NEntries_tchain",&oNEntries);
  T->SetBranchAddress("NReco",&oNReco);
  T->SetBranchAddress("NMuon",&oNMuon);
  T->SetBranchAddress("NJet",&oNJet);
  T->SetBranchAddress("PxJet1",&oPxJet1);
  T->SetBranchAddress("PyJet1",&oPyJet1);
  T->SetBranchAddress("PzJet1",&oPzJet1);
  T->SetBranchAddress("EJet1",&oEJet1);
  T->SetBranchAddress("PxJet2",&oPxJet2);
  T->SetBranchAddress("PyJet2",&oPyJet2);
  T->SetBranchAddress("PzJet2",&oPzJet2);
  T->SetBranchAddress("EJet2",&oEJet2);
  T->SetBranchAddress("PxMu",&oPxMu);
  T->SetBranchAddress("PyMu",&oPyMu);
  T->SetBranchAddress("PzMu",&oPzMu);
  T->SetBranchAddress("EMu",&oEMu);
  T->SetBranchAddress("PxMiss",&oPxMiss);
  T->SetBranchAddress("PyMiss",&oPyMiss);
  T->SetBranchAddress("PzMiss",&oPzMiss);
  T->SetBranchAddress("EMiss",&oEMiss);
  T->SetBranchAddress("MuD0",&oMuD0);
  T->SetBranchAddress("MuZ0",&oMuZ0);
  T->SetBranchAddress("MuD0sig",&oMuD0sig);
  T->SetBranchAddress("MuZ0sig",&oMuZ0sig);

  Long64_t EntriesTree = T->GetEntries();

  Double_t nPreselection = 0, nOneMuon = 0, nSelection = 0;
  Int_t nSliding[81], nBcut[25][81];

  for (int i=0;i<81;i++) for (int j=0; j<25; j++) nSliding[i] = nBcut[j][i] = 0;


  Long64_t NN = EntriesTree;
  Double_t SF = 1.;
  if (RunOnN >= 0 && RunOnN < EntriesTree) {NN = RunOnN; SF = 1.*NN/EntriesTree;}
  cout<<"getvalues.C:: Running on "<<NN<<" entries out of "<<EntriesTree<<endl;
  cout<<"getvalues.C:: Scale factor is "<<SF<<endl;

  nOneMuon = NN;


  for (Long64_t i = 0; i < NN; i++){

    if (i%1000000 == 0) cout<<i<<" / "<<NN<<endl;

    T->GetEntry(i);

    nPreselection = oNEntries * SF;


    if (oNMuon != 1  || oNJet != 2 || oEMiss < 0) {nPreselection--; nOneMuon--; continue;}

    
    TLorentzVector lvj1(oPxJet1, oPyJet1, oPzJet1, oEJet1);
    TLorentzVector lvj2(oPxJet2, oPyJet2, oPzJet2, oEJet2);
    TLorentzVector lvmu(oPxMu, oPyMu, oPzMu, oEMu);
    TLorentzVector lvmiss(oPxMiss, oPyMiss, oPzMiss, oEMiss);

    if (obsID == "mtot") {
      toret.push_back((lvj1 + lvj2 + lvmiss + lvmu).M());
      continue;
    }

    Double_t cosjj = TMath::Cos(lvj1.Angle(lvj2.Vect()));
    Double_t cosj1mu = TMath::Cos(lvj1.Angle(lvmu.Vect()));
    Double_t cosj2mu = TMath::Cos(lvj2.Angle(lvmu.Vect()));
    
    Bool_t sel1 = (TMath::Abs(lvmiss.CosTheta()) < 0.94);
    if (!sel1) continue;
    Bool_t sel2 = (TMath::Cos(lvmiss.Angle(lvmu.Vect())) < 0.80);
    if (!sel2) continue;
    Bool_t sel3 = (oEJet1 >= 3. && oEJet2 >= 3.);
    if (!sel3) continue;
    Bool_t sel4 = (cosjj > -0.8 && cosjj < 0.98);
    if (!sel4) continue;
    //Bool_t sel5_old = (TMath::Min(lvj1.M(), lvj1.M()) > 1.8 && TMath::Min(lvj1.M2(), lvj1.M2()) > 0);
    Bool_t sel6 = (TMath::Max( cosj1mu, cosj2mu ) < 0.8);
    if (!sel6) continue;

    Bool_t s21 = (TMath::Min(lvj1.M(), lvj2.M()) > 0.2 && TMath::Min(lvj1.M2(), lvj2.M2()) > 0);
    if (!s21) continue;
    Bool_t s22 = (TMath::Min( cosj1mu, cosj2mu ) > -0.98);
    if (!s22) continue;
    Bool_t s23 = ((lvj1 + lvj2 + lvmiss + lvmu).M() > 80.);
    if (!s23) continue;

    Bool_t selection = sel1 && sel2 && sel3 && sel4 && sel6 && s21 && s22 && s23;

   
    if (selection){

      
      if (obsID == "d0") toret.push_back(TMath::Abs(oMuD0sig));

      Double_t im = mass;
      Double_t pRecoil = (Z_mass*Z_mass - 1.* im * im ) / (2 * Z_mass);
      Bool_t slid1 = (TMath::Abs(lvmiss.P() - pRecoil) <= 3.5);
      Bool_t slid2 = (TMath::Abs( (lvj1+lvj2+lvmu).M() - 1.*im) <= 4.);

      if (slid1 && slid2){

	if (obsID == "d0sliding") toret.push_back(TMath::Abs(oMuD0sig));
      }
      
    }

    

    

  
  }

  
  cout<<fname<<" READ "<<endl;
  return std::pair<std::vector<Double_t>, Double_t>{toret, SF};
}
