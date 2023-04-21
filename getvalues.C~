
#include "lumisettings.h"
#include "Cut.h"

Double_t Muon_mass = 0.105658;
Double_t Z_mass = 91.1876;



std::pair<std::vector<Double_t>, Double_t> getvalues(TString obsID, TString opt="signal", Int_t mass = 50, TString lifetime = "n/a", Long64_t RunOnN = -1, Double_t d0cut=8, Int_t jalg = 2, TString analysis_opt = "< d2d dsigma anymass1L2M", TString dir="../MyExternalAnalysis/results/"){
  
  // getvalue.second is the scale factor
  // If sliding is not needed, you can use any mass.

  TREE->Reset();

  TString Dir = dir;

  if (obsID == "d0" || obsID == "d0sliding" || obsID == "mass1mm")
    Dir = Dir+"/skimmed/";

  TString fname=Form("%s%s",Dir.Data(), AnalysisResults(opt,Form("%d",mass),lifetime).Data());

  cout<<"getvalues.C:: Opening file "<<fname<<endl;

  std::vector<Double_t> toret;
  toret.clear();
  

  TREE->Add(fname);
  TREESetBranch();
  
  Long64_t EntriesTree = TREE->GetEntries();

  

  Double_t nPreselection = 0, nOneMuon = 0, nSelection = 0;
 

  Long64_t NN = EntriesTree;
  Double_t SF = 1.;
  if (RunOnN >= 0 && RunOnN < EntriesTree) {NN = RunOnN; SF = 1.*NN/EntriesTree;}
  cout<<"getvalues.C:: Running on "<<NN<<" entries out of "<<EntriesTree<<endl;
  cout<<"getvalues.C:: Scale factor is "<<SF<<endl;

  nOneMuon = NN;


  for (Long64_t i = 0; i < NN; i++){

    if (i%1000000 == 0) cout<<i<<" / "<<NN<<endl;

    TREE->GetEntry(i);

    nPreselection = oNEntries * SF;

    
    BUILD_DERIVATE(jalg);

    if (obsID == "mtot") {
      toret.push_back((lvvis + lvmiss).M());
      continue;
    }

    if (obsID == "cosjj"){
      toret.push_back(TMath::Cos(lvj1.Angle(lvj2.Vect())));
      continue;
    }

    if (obsID == "MAXcosjmu"){
      toret.push_back(TMath::Max(oCosj1Mu, oCosj2Mu));
      continue;
    }

    if (obsID == "MINcosjmu"){
      toret.push_back(TMath::Min(oCosj1Mu, oCosj2Mu));
    }

    Bool_t selection = false;

    if (analysis_opt.Contains("anymass1L2M")){

      if (oNJet->at(jalg) == 2) selection = SELECTION_MM_2JET(jalg);
      if (oNJet->at(jalg) == 1) selection = SELECTION_LM_1JET(jalg);
      
    }

    if (analysis_opt.Contains("anymass1L2L")){

      if (oNJet->at(jalg) == 2) selection = SELECTION_LM_2JET(jalg);
      if (oNJet->at(jalg) == 1) selection = SELECTION_LM_1JET(jalg);
      
    }
    
   

    

    
   
    if (selection){

      if (obsID == "mass1mm"){

	if (oMuD0 > 1.) toret.push_back(lvvis.M());
      }

      
      if (obsID == "d0") toret.push_back(TMath::Abs(oMuD0sig));

      Double_t im = mass;
      Double_t pRecoil = (Z_mass*Z_mass - 1.* im * im ) / (2 * Z_mass);
      Bool_t slid1 = (TMath::Abs(lvmiss.P() - pRecoil) <= 3.5);
      Bool_t slid2 = (TMath::Abs( lvvis.M() - 1.*im) <= 4.);

      if (obsID == "d0sliding" && slid1 && slid2) toret.push_back(TMath::Abs(oMuD0sig));
      

      if (obsID == "mass_after_dcut"){

	Bool_t cut_condition = false;

	if (analysis_opt.Contains("d3d") and analysis_opt.Contains("dsigma"))
	  cut_condition =  (TMath::Sqrt(oMuD0sig*oMuD0sig + oMuZ0sig*oMuZ0sig) < 1.*d0cut);

	else if (analysis_opt.Contains("d2d") && analysis_opt.Contains("dsigma"))
	  cut_condition =  (TMath::Abs(oMuD0sig) < 1.*d0cut);

	else if (analysis_opt.Contains("d2d") && analysis_opt.Contains("dmm"))
	  cut_condition = (TMath::Abs(oMuD0) < 1.*d0cut/100.);

	else cout<<"CutFlowOK.C:: ERROR - OPTIONS NOT SUPPORTED"<<endl;

	if (analysis_opt.Contains(">")) cut_condition = !cut_condition;

	if (cut_condition) toret.push_back(lvvis.M());
	

	
      }



      
      
    }

    

    

  
  }

  
  cout<<fname<<" READ "<<endl;
  return std::pair<std::vector<Double_t>, Double_t>{toret, SF};
}
