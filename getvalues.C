
#include "lumisettings.h"
#include "Cut.h"

Double_t Muon_mass = 0.105658;
Double_t Z_mass = 91.1876;

enum OBS_ID {
  od0,      // oMuD0sig after 1 muon selection
  od0sel,   // oMuD0sig after event selection driven by analysis_opt and jalg
  od0sliding,
  omass1mm,
  omtot,
  ocosjj,
  ocospmiss,
  ocospmissmu,
  oMAXcosjmu,
  oMINcosjmu,
  oMINEjet,
  omass_selection, // Visible mass after event selection (no sliding cuts)
  oemiss_selection, // Emiss after event selection (no sliding cuts)
  omass_after_dcut,
  oVtxXY,   // Vtx distance to 0 on XY after 1 muon selection
  oVtxXYZ,   // Vtx distance to 0 in 3D after 1 muon selection
  oVtxXYsliding,   // Vtx distance to 0 on XY after selection driven by analysis_opt and jalg + sliding cuts
  oVtxXYZsliding   // Vtx distance to 0 in 3D after selection driven by analysis_opt and jalg + sliding cuts
};



std::pair<std::vector<Double_t>, Double_t> getvalues(OBS_ID obsID, TString opt="signal", Int_t mass = 50, TString lifetime = "n/a", Long64_t RunOnN = -1, Double_t d0cut=8, Int_t jalg = 2, TString analysis_opt = "< d2d dsigma anymass1L2M", TString dir="../MyExternalAnalysis/results/"){
  
  // getvalue.second is the scale factor
  // If sliding is not needed, you can use any mass.

  TREE->Reset();

  TString Dir = dir;

  if (obsID == od0sel || obsID == od0sliding || obsID == omass1mm ||
      obsID == oVtxXYsliding || obsID == oVtxXYZsliding || obsID == omass_selection)
    Dir = Dir+"/skimmed/"; // remember to check the target of the symlink!

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

    //cout<<"Njet: "<<oNJet->at(0)<<" "<<oNJet->at(1)<<" "<<oNJet->at(2)<<endl;

    //if (oNJet->at(jalg) != 1) continue;
    if (obsID == od0){
      toret.push_back(TMath::Abs(oMuD0sig));
      continue;
    }

    if (obsID == omtot) {
      toret.push_back((lvvis + lvmiss).M());
      continue;
    }

    if (obsID == ocospmiss){
      toret.push_back(TMath::Cos(lvmiss.Theta()));
      continue;
    }

    if (obsID == ocospmissmu){
      toret.push_back(TMath::Cos(lvmiss.Angle(lvmu.Vect())));
      continue;
    }

    if (obsID == ocosjj){
      toret.push_back(TMath::Cos(lvj1.Angle(lvj2.Vect())));
      continue;
    }

    if (obsID == oMAXcosjmu){
      toret.push_back(TMath::Max(oCosj1Mu, oCosj2Mu));
      continue;
    }

    if (obsID == oMINcosjmu){
      toret.push_back(TMath::Min(oCosj1Mu, oCosj2Mu));
      continue;
    }

    if (obsID == oMINEjet){
      toret.push_back(TMath::Min(oEJet1->at(jalg),oEJet2->at(jalg)));
      continue;
    }

    if (obsID == oVtxXY){
      toret.push_back(TMath::Sqrt(oVtx_x*oVtx_x + oVtx_y*oVtx_y));
      continue;
    }

    if (obsID == oVtxXYZ){
      toret.push_back(TMath::Sqrt(oVtx_x*oVtx_x + oVtx_y*oVtx_y + oVtx_z*oVtx_z));
      continue;
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
      
      
      
      if (obsID == omass1mm){
	if (oMuD0 > 1.) toret.push_back(lvvis.M());
	continue;
      }

      
      if (obsID == od0sel) {
	toret.push_back(TMath::Abs(oMuD0sig));
	continue;
      }

      Double_t im = mass;
      Double_t pRecoil = (Z_mass*Z_mass - 1.* im * im ) / (2 * Z_mass);
      Bool_t slid1 = (TMath::Abs(lvmiss.P() - pRecoil) <= 3.5);
      Bool_t slid2 = (TMath::Abs( lvvis.M() - 1.*im) <= 4.);


      if (obsID == omass_selection){
	toret.push_back(lvvis.M());
	continue;
      }

      if (obsID == oemiss_selection){
	toret.push_back(oEMiss);
	continue;
      }

      if (slid1 && slid2){
	if (obsID == od0sliding) toret.push_back(TMath::Abs(oMuD0sig));
	if (obsID == oVtxXYsliding) toret.push_back(TMath::Sqrt(oVtx_x*oVtx_x + oVtx_y*oVtx_y));
	if (obsID == oVtxXYZsliding) toret.push_back(TMath::Sqrt(oVtx_x*oVtx_x + oVtx_y*oVtx_y + oVtx_z*oVtx_z));
	continue;
      }
								      
      

      if (obsID == omass_after_dcut){

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

  
  cout<<"getvalues.C:: "<<fname<<" READ "<<endl;
  
  cout<<"getvalues.C:: returning array with "<<toret.size()<<" elements"<<endl<<endl;
  return std::pair<std::vector<Double_t>, Double_t>{toret, SF};
}
