#include "lumisettings.h"
#include "Cut.h"

Double_t Muon_mass = 0.105658;
Double_t Z_mass = 91.1876;


std::vector<int> possible_masses = {1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100};

std::vector<int> possible_dcut = {0,1,2,3,4,5,6,7,8,10,12,16,20,25,30,35,40,45,50,75,100,125,150,175,200,225,250,275,300,350,400,500}; // sigma units

int dcut_id(int dcut){
  // used to return the correct index in CutFlowOK[M] vector
  std::vector<int>::iterator itr = std::find(possible_dcut.begin(), possible_dcut.end(), dcut);
  if (itr != possible_dcut.cend()) return std::distance(possible_dcut.begin(), itr) + 4;
  else return -1;
}

std::map<int, std::vector<double>> CutFlowOK(TString opt="signal", Int_t mass=80, TString lifetime = "n/a", TString dir="../MyExternalAnalysis/results/", Long64_t RunOnN = -1, Bool_t CutByCutFlow = false, Int_t jalg = 2, TString analysis_opt="< d2d dsigma anymass1L2M window [2,0.25]", Bool_t FixedMass = false){


  // analysis options available:
  // ">" or "<" for the type of cut in impact parameter
  // "dsigma" or "dmm" to cut in number of sigma or in unit 10-5 m
  // "d2d" or "d3d" to cut on D0 or D0+Z0
  // "fixedwindow" to set \pm 4 and \pm 3.5 in the Mvis, Emiss window, or "window [A,B]" where A is in unit of B*sqrt(mass), or "window [Asig,Bsig,Abkg,Bbkg]" using two different parameterizations for signal and background.
  // type of analysis:
  ///// "anymass1L2M" : mass independent, it uses old analysis for all cases with 2 Jets and LM-1j analysis for all cases with 1 jet
  ///// "anymass1L2L": mass independent, it uses LM-1j analysis or LM-2j analyses according to number of jets 
  
  // mass is used only to open the corresponding file. The output will contain the cut for all the possible_masses unless "FixedMass = true"

  // To be read like this: CutFlowOK[M] where M(int) is the Mass in GeV. CutFlowOK[M] is a vector:
  
  //// Standard ouptut:
  //// 0: Nnocut(opt,mass,lifetime)  <--- written in lumisettings.h, not read from the file!
  //// 1: (Tree entries && nMuon==1)  2: selection  3: sliding(mass)
  //// 4...24:  sliding(mass) & d0cut = index-4
    
  //// When CutByCutFlow is set: <-- NOT SUPPORTED
  //// 0: Nnocut(opt,mass,lifetime) <--- written in lumisettings.h, not read from the file!
  //// 1: entries of the tree,   2: one muon selection (should be equal to 1)
  //// 3: 1& cos(pmiss)<0.94
  //// 4: 1& cos(pmiss,mu)<0.8
  //// 5: 1& ejet>3
  //// 6: 1& cosjj
  //// 7:
  //// 8: 1& cos(j,mu)<0.8
  //// 9: 1& Mj > 0.2 & M2j > 0
  //// 10: 1& cos(jmu) > -0.98
  //// 11: 1& Mtot > 80
  //// 12: selection


  std::vector<int> here_possible_masses = possible_masses;
  if (FixedMass) here_possible_masses = std::vector<int>{mass};
 

  TREE->Reset();
  
  TString fname=Form("%s%s",dir.Data(), AnalysisResults(opt,Form("%d",mass),lifetime).Data());

  if (fname.Contains(".root")){
    TREE->Add(fname);
    cout<<"CutFlowOK:: Running on a single file: "<<fname<<endl;
  }

  else{
    TString rootfile;
    std::ifstream infile(fname);
    while(!infile.eof()){
      infile >> rootfile;
      if (rootfile=="") continue;
      TREE->Add(rootfile);
      cout<<"CutFlowOK:: "<<rootfile<<" added to the chain."<<endl;
    }     
  }

  
  

  std::map<int,std::vector<double>> toret;
  toret.clear();

  for (int im : here_possible_masses)
    toret[im] = std::vector<double>{};
  
  TREESetBranch();
  

  Long64_t EntriesTree = TREE->GetEntries();

  Double_t nPreselection = 0, nOneMuon = 0, nSelection = 0;
  Int_t nSliding[200], nBcut[500][200]; //first index: Dcut (in sigma). Second index analysis mass (in gev)
  Int_t nCutByCut[50][200]; //first index: see above. Second index: analysis mass

  for (int i=0;i<200;i++) for (int j=0; j<25; j++) nSliding[i] = nBcut[j][i] = 0;
  for (int i=0;i<50;i++) for (int j=0; j<200;j ++) nCutByCut[i][j] = 0;


  

  Long64_t NN = EntriesTree;
  Double_t SF = 1.;
  if (RunOnN >= 0 && RunOnN < EntriesTree) {NN = RunOnN; SF = 1.*NN/EntriesTree;}
  cout<<"CutFlowOK.C:: Running on "<<NN<<" entries out of "<<EntriesTree<<endl;
  cout<<"CutFlowOK.C:: Scale factor is "<<SF<<endl;

  nOneMuon = NN;

  for (Long64_t i = 0; i < NN; i++){

    if (i%1000000 == 0) cout<<i<<" / "<<NN<<endl;

    
    TREE->GetEntry(i);

    nPreselection = oNEntries * SF;

    for (int j=0; j<200; j++) nCutByCut[0][j] ++;

    if (oNMuon != 1  || oNJet->at(jalg) < 1 || oEMiss < 0) {nPreselection--; nOneMuon--; continue;}

    for (int j=0; j<200; j++) nCutByCut[1][j] ++;


    // this builds the lorentz vectors
    BUILD_DERIVATE(jalg);
    

    Bool_t selection = false;

    if (analysis_opt.Contains("cutvariation")){


      std::vector<float> vcuts = GetFloatArray(analysis_opt);
      double CUTcosjj = vcuts.at(0);
      double CUTmincosjmu = vcuts.at(1);
      double CUTmaxcosjmu = vcuts.at(2);

      if (oNJet->at(jalg) == 2) selection = SEL_MM_2J(jalg,CUTcosjj,CUTmincosjmu,CUTmaxcosjmu);
      if (oNJet->at(jalg) == 1) selection = SELECTION_LM_1JET(jalg);
      
    }

    else selection = SELECTION_STRING_KINE(jalg,analysis_opt);
    
    if (!selection) continue;
    
    
      
   
    if (selection){

      for (int j=0; j<200; j++) nCutByCut[11][j]++;

      nSelection ++;

      
      for (int im : here_possible_masses){

	Bool_t slidingsel = SELECTION_STRING_SLIDING(im, analysis_opt, opt);
	

	if (slidingsel){
	  nSliding[im]++;
	  for (int id0 = 0; id0 < possible_dcut.size(); id0++){

	    int dcut = possible_dcut[id0];

	    Bool_t cut_condition = SELECTION_STRING_DCUT(dcut, analysis_opt);
	    
	    if (cut_condition) nBcut[id0][im]++;
	    
	  } // loop on id0
	}

	
      }
    }

  
  }

  if (!CutByCutFlow){
    for (int im : here_possible_masses){
      toret[im].push_back(1.*Nnocut(opt,im,lifetime)*SF);
      toret[im].push_back(1.*nOneMuon);
      toret[im].push_back(1.*nSelection);
      toret[im].push_back(1.*nSliding[im]);
      for (int id0 = 0; id0 < possible_dcut.size(); id0++)
	toret[im].push_back(1.*nBcut[id0][im]);
    }
  }
  if (CutByCutFlow){
    for (int im : here_possible_masses){
      toret[im].push_back(1.*Nnocut(opt,im,lifetime)*SF);
      for (int cbc = 0; cbc<12; cbc++)
	toret[im].push_back(nCutByCut[cbc][im]);
    }
  }
  
  

  cout << "Preselection: "<<nPreselection<<endl;
  cout << "One muon    : "<<nOneMuon<<endl;
  cout << "Selection   : "<<nSelection<<endl;

  

  return toret;
}

      
   
