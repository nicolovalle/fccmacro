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

std::map<int, std::vector<double>> CutFlowOK(TString opt="signal", Int_t mass=80, TString lifetime = "n/a", TString dir="../MyExternalAnalysis/results/", Long64_t RunOnN = -1, Bool_t CutByCutFlow = false, Int_t jalg = 2, TString analysis_opt="< d2d dsigma anymass1L2M"){

  // CUT BY CUT FLOW NOT IMPLEMENTED YET FOR LESS THAN 2 JETS

  // analysis options available:
  // ">" or "<" for the type of cut in impact parameter
  // "dsigma" or "dmm" to cut in number of sigma or in unit 10-5 m
  // "d2d" or "d3d" to cut on D0 or D0+Z0
  // type of analysis:
  ///// "anymass1L2M" : mass independent, it uses old analysis for all cases with 2 Jets and LM-1j analysis for all cases with 1 jet
  ///// "anymass1L2L": mass independent, it uses LM-1j analysis or LM-2j analyses according to number of jets 
  
  // mass is used only to open the corresponding file. The output will contain the cut for all the possible_masses

  // To be read like this: CutFlowOK[M] where M(int) is the Mass in GeV. CutFlowOK[M] is a vector:
  
  //// Standard ouptut:
  //// 0: Nnocut(opt,lifetime)  <--- written in lumisettings.h, not read from the file!
  //// 1: nOneMuon  2: selection  3: sliding(mass)
  //// 4...24:  sliding(mass) & d0cut = index-4
    
  //// When CutByCutFlow is set:
  //// 0: Nnocut(opt,lifetime) <--- written in lumisettings.h, not read from the file!
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
 
    

  TChain *T = new TChain("eventsTree");
  //TString fname = "../MyExternalAnalysis/results/AnalysisResultsfilelist_Zbb_highstat_pt7.txt.root";
  TString fname=Form("%s%s",dir.Data(), AnalysisResults(opt,Form("%d",mass),lifetime).Data());

  if (fname.Contains(".root")){
    T->Add(fname);
    cout<<"CutFlowOK:: Running on a single file: "<<fname<<endl;
  }

  else{
    TString rootfile;
    std::ifstream infile(fname);
    while(!infile.eof()){
      infile >> rootfile;
      if (rootfile=="") continue;
      T->Add(rootfile);
      cout<<"CutFlowOK:: "<<rootfile<<" added to the chain."<<endl;
    }     
  }

  
  

  std::map<int,std::vector<double>> toret;
  toret.clear();

  for (int im : possible_masses)
    toret[im] = std::vector<double>{};
  
  
  //TFile *File = new TFile(fname);
  //TTree *T = (TTree*)File->Get("eventsTree");

  
 

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
  T->SetBranchAddress("PxVis",&oPxVis);
  T->SetBranchAddress("PyVis",&oPyVis);
  T->SetBranchAddress("PzVis",&oPzVis);
  T->SetBranchAddress("EVis",&oEVis);
  T->SetBranchAddress("MuD0",&oMuD0);
  T->SetBranchAddress("MuZ0",&oMuZ0);
  T->SetBranchAddress("MuD0sig",&oMuD0sig);
  T->SetBranchAddress("MuZ0sig",&oMuZ0sig);

  Long64_t EntriesTree = T->GetEntries();

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

    
    T->GetEntry(i);

    nPreselection = oNEntries * SF;

    for (int j=0; j<200; j++) nCutByCut[0][j] ++;

    if (oNMuon != 1  || oNJet->at(jalg) < 1 || oEMiss < 0) {nPreselection--; nOneMuon--; continue;}

    for (int j=0; j<200; j++) nCutByCut[1][j] ++;


    // this builds the lorentz vectors
    BUILD_DERIVATE(jalg);
    

    Bool_t selection = false;

    if (analysis_opt.Contains("anymass1L2M")){

      if (oNJet->at(jalg) == 2) selection = SELECTION_MM_2JET(jalg);
      if (oNJet->at(jalg) == 1) selection = SELECTION_LM_1JET(jalg);
      
    }


    if (analysis_opt.Contains("anymass1L2L")){

      if (oNJet->at(jalg) == 2) selection = SELECTION_LM_2JET(jalg);
      if (oNJet->at(jalg) == 1) selection = SELECTION_LM_1JET(jalg);
      
    }
    
    if (!selection) continue;
    
    
      
   
    if (selection){

      for (int j=0; j<200; j++) nCutByCut[11][j]++;

      nSelection ++;

      
      for (int im : possible_masses){

	Double_t pRecoil = (Z_mass*Z_mass - 1.* im * im ) / (2 * Z_mass);
	Bool_t slid1 = (TMath::Abs(lvmiss.P() - pRecoil) <= 3.5);
	Bool_t slid2 = (TMath::Abs( (lvvis).M() - 1.*im) <= 4.);

	if (slid1 && slid2){
	  nSliding[im]++;
	  for (int id0 = 0; id0 < possible_dcut.size(); id0++){
	    Bool_t cut_condition = false;
	    
	    if (analysis_opt.Contains("d3d") and analysis_opt.Contains("dsigma"))
	      cut_condition =  (TMath::Sqrt(oMuD0sig*oMuD0sig + oMuZ0sig*oMuZ0sig) < 1.*possible_dcut[id0]);
	    else if (analysis_opt.Contains("d2d") && analysis_opt.Contains("dsigma"))
	      cut_condition =  (TMath::Abs(oMuD0sig) < 1.*possible_dcut[id0]);
	    else if (analysis_opt.Contains("d2d") && analysis_opt.Contains("dmm"))
	      cut_condition = (TMath::Abs(oMuD0) < 1.*possible_dcut[id0]/100.);
	    else
	      cout<<"CutFlowOK.C:: ERROR - OPTIONS NOT SUPPORTED"<<endl;
	    
	    if (analysis_opt.Contains(">")) cut_condition = !cut_condition;
	    
	    if (cut_condition) nBcut[id0][im]++;
	    
	  } // loop on id0
	}

	
      }
    }

  
  }

  if (!CutByCutFlow){
    for (int im : possible_masses){
      toret[im].push_back(1.*Nnocut(opt,lifetime)*SF);
      toret[im].push_back(1.*nOneMuon);
      toret[im].push_back(1.*nSelection);
      toret[im].push_back(1.*nSliding[im]);
      for (int id0 = 0; id0 < possible_dcut.size(); id0++)
	toret[im].push_back(1.*nBcut[id0][im]);
    }
  }
  if (CutByCutFlow){
    for (int im : possible_masses){
      toret[im].push_back(1.*Nnocut(opt,lifetime)*SF);
      for (int cbc = 0; cbc<12; cbc++)
	toret[im].push_back(nCutByCut[cbc][im]);
    }
  }
  
  

  cout << "Preselection: "<<nPreselection<<endl;
  cout << "One muon    : "<<nOneMuon<<endl;
  cout << "Selection   : "<<nSelection<<endl;

  

  return toret;
}





//TString maketablestring

TString maketableline1(TString A1, int A2, int B1, TString C1, double D1, double E1, double F1, double G1, double H1,  double I1,  double J1, double K1,  double L1,  double M1,  double N1, double O1,  TString bkgcolor="#ffffff" , TString color="#990000" ){

/*
     A             B           C             D           E        F          G          H            I               J               K             L            M                N              O
| type        |   Mana  |  lifetime   |  log(U2)    |  xsec   |  Ngen   |  weight  |   Nsel    |   Nsliding   |   Ndcut = 4    |  Ndcut = 8 |  Ndcut = 20  | Ndcut  = 50 | Ndcut  = 100 | Ndcut  = 200 |    
| M if signal |         |             |             |         |  wghted |          |  wghted   | wghted       | wghted         | wghted     | wghted       | wghted      | wghted       | wghted       |   
   
*/
  

  TString A = Form("<td> %s <br> M=%d GeV</td>", A1.Data(), A2);
  if (A2<0) A = Form("<td> %s </td>",A1.Data());
  TString B = Form("<td> %d </td>", B1);
  TString C = Form("<td> %s </td>", C1.Data());
  TString D;
  if (D1>0) D = Form("<td> %3.3f </td>", 2.*TMath::Log10(D1));
  if (D1<=0) D = Form("<td> </td>");	      
  TString E = Form("<td> %f </td>", E1);
  TString F = Form("<td> %10.2f k <br> <ccc style=\"color:%s;\"> %10.2f k </ccc></td>", F1/1000., color.Data(), (F1/1000.)*G1);
  TString G = Form("<td> %10.2f </td>", G1);
  TString H = Form("<td> %10.2f <br> <ccc style=\"color:%s;\"> %10.2f </ccc></td>", H1, color.Data(), H1*G1);
  TString I = Form("<td> %10.2f <br> <ccc style=\"color:%s;\"> %10.2f </ccc></td>", I1, color.Data(), I1*G1);
  TString J = Form("<td> %10.2f <br> <ccc style=\"color:%s;\"> %10.2f </ccc></td>", J1, color.Data(), J1*G1);
  TString K = Form("<td> %10.2f <br> <ccc style=\"color:%s;\"> %10.2f </ccc></td>", K1, color.Data(), K1*G1);
  TString L = Form("<td> %10.2f <br> <ccc style=\"color:%s;\"> %10.2f </ccc></td>", L1, color.Data(), L1*G1);
  TString M = Form("<td> %10.2f <br> <ccc style=\"color:%s;\"> %10.2f </ccc></td>", M1, color.Data(), M1*G1);
  TString N = Form("<td> %10.2f <br> <ccc style=\"color:%s;\"> %10.2f </ccc></td>", N1, color.Data(), N1*G1);
  TString O = Form("<td> %10.2f <br> <ccc style=\"color:%s;\"> %10.2f </ccc></td>", O1, color.Data(), O1*G1);

  return "<tr align=\"center\" style=\"background-color:"+bkgcolor+";\">"+A+B+C+D+E+F+G+H+I+J+K+L+M+N+O+"</tr>";
  
  
}

TString maketableline2(TString a1, Int_t a2, TString b, double c, std::vector<double> v, TString bkgcolor="#ffffff"){

  /*

    A          B          C       D       E       F             G                  H              I        J1     J2    K               L                    M                   N
  | type  | lifetife | log(U2) | Ngen | ExtAna | OneMuon | cos(pmiss)<0.94 | cos(pmiss,mu)<0.8 | ejet>3 |cosjj | cosjj | cos(j,mu)<0.8 | Mj > 0.2 & M2j > 0 | cos(jmu) > -0.98 | Mtot > 80 | 
   */

  TString A = Form("<td> %s <br> M=%d GeV</td>", a1.Data(),a2);
  if (a2<0) A = Form("<td> %s </td>", a1.Data());
  TString B = Form("<td> %s </td>", b.Data());
  TString C;
  if (c>0)  C = Form("<td> %3.3f </td>", 2.*TMath::Log10(c));
  if (c<=0) C = Form("<td> </td>");

  TString toret = A+B+C;

  for (int i: std::vector<int>{0,1}){
    toret = toret+ Form("<td> %10.0f </td>",v.at(i));
  }
  toret = toret+ Form("<td> %10.0f <br> <ccc style=\"color:#000099\"> %3.2f%% </ccc></td>",v.at(2), 100.*v.at(2)/v.at(0));
  double comb_eff = 1;
  for (int i: std::vector<int>{3,4,5,6,7,8,9,10,11}){
    double efff = v.at(i)/v.at(2);
    comb_eff *= efff;
    TString pcol = efff < 0.5 ? "#a300a3" : efff < 0.8 ? "#ff1d1d" : "#009900";
    toret = toret+ Form("<td> %10.0f <br> <ccc style=\"color:%s\"> %3.2f%% </ccc></td>",v.at(i), pcol.Data(), 100.*efff);
  }

  toret = toret+ Form("<td> %10.0f <br> <ccc style=\"color:#009900; font-size:small\"> %3.2f (<em>%3.2f</em>) %% </ccc></td>",v.at(12), 100.*v.at(12)/v.at(2), 100.*comb_eff );

  return "<tr align=\"center\" style=\"background-color:"+bkgcolor+";\">"+toret+"</tr>";

}


void makeHTMLtable(TString outfile="./summary.html", int iopt = 1, TString AnalysisResPath = "../MyExternalAnalysis/results/skimmed/", Int_t jalg = 2,  Bool_t ProcessOnlySignal = false, Int_t RunOnN = -1, Bool_t upload = false, TString comments_on_top="jalg==2 - asking D0 > dcut in sigma unit"){

  // iopt. 1: summary and b cut     2: selection cut flow

  TString Headers;
  if (iopt == 1){
    Headers =  "<tr align=\"center\" style=\"background-color: #ff00ff;\"><td>Sample</td><td>M_ana</td><td>Log(ctau)</td><td>log(U2)</td><td>xsec/pb</td><td>Ngen</td><td>Weight</td><td>EvtSel</td><td>KinSel</td><td>b=4sig</td><td>b=8sig</td><td>b=20sig</td><td>b=50sig</td><td>b=100sig</td><td>b=200sig</td></tr>";
  }
  if (iopt == 2){
    Headers = "<tr align=\"center\" style=\"background-color: #ff00ff;\"><td>Sample</td><td>Log(ctau)</td><td>log(U2)</td><td>Ngen</td><td>MyAnalysis</td><td>one muon</td><td style=\"background-color: #ffb299;\">cos(pmiss)<.94</td><td style=\"background-color: #ffb299;\">cos(pmiss,mu)<.8</td><td style=\"background-color: #ffb299;\">Ejet>3</td><td style=\"background-color: #ffb299;\">cosjj>-.8</td><td style=\"background-color: #ffb299;\">cosjj<.98</td><td style=\"background-color: #ffb299;\">cos(j,mu)<.8</td><td style=\"background-color: #ffb299;\">Mj>.2</td><td style=\"background-color: #ffb299;\">cos(j,mu)>-.98</td><td style=\"background-color: #ffb299;\">M>80</td><td>evt.sel.</td></tr>";
  }

  ofstream Outfile;
  Outfile.open(outfile);
  Outfile<<"<br>"<<endl<<comments_on_top<<"<br>"<<endl;
  if (iopt == 1) Outfile<<"<br>Red numbers are multiplied by the weight<br><br>"<<endl;
  Outfile<<"<table border=\"1\">"<<endl;
  
  

  std::vector<double> V;

  std::set<int> enabled_masses;

  
    
  for ( int m : possible_masses ){
    bool PrintHeader = true;
    for ( TString mp : std::vector<TString>{"p","m"}){
      for (int unit = 0; unit < 9; unit++){
	for (int decimal = 0; decimal < 6; decimal += 5){

	  TString mstr = Form("%d",m);
	  TString lt = Form("%s%dp%d", mp.Data(), unit, decimal);

	  TString lt_for_tab = Form("%s%d.%d", mp=="m"?"-":"+", unit, decimal);


	  TString FileNameToCheck = Form("%s%s", AnalysisResPath.Data(), AnalysisResults("signal",Form("%d",m),lt).Data());

	  if (gSystem->AccessPathName(FileNameToCheck)) continue;

	  if (PrintHeader) {Outfile<<Headers<<endl; PrintHeader = false;}

	  enabled_masses.insert(m);
	  //4 8 20 50 100 200
	  if (iopt == 1) V = CutFlowOK("signal",m,lt,AnalysisResPath,-1,false,jalg,"> d2d dsigma jalg")[m];
	  if (iopt == 2) V = CutFlowOK("signal",m,lt,AnalysisResPath,-1,true,jalg,"> d2d dsigma jalg")[m];

	  Double_t w = Weight("signal",mstr,lt);
	  
 
	  
	  TString line;

	  //                                    A1      A2 B  C1           D                    E                       F1  G1, H    I      J                K              L                M               N                   O
	  if (iopt == 1)  line = maketableline1("Signal",m, m, lt_for_tab, Coupling(mstr,lt), xsec("signal",mstr,lt), V[0],  w, V[2], V[3], V[dcut_id(4)], V[dcut_id(8)], V[dcut_id(20)], V[dcut_id(50)], V[dcut_id(100)], V[dcut_id(200)]);
	  

	  if (iopt == 2) line = maketableline2("Signal",m,lt_for_tab,Coupling(mstr,lt), V);
	  

	  
	  Outfile<<line;
	  Outfile<<endl;
	  
	}
      }
    }
  }

  std::vector<TString> SomeColors = {"#ffbf00", "#9fe2bf", "#6495ed", "ffee00", "#ffbf00", "#9fe2bf", "#6495ed", "ffee00", "#ffbf00", "#9fe2bf", "#6495ed", "ffee00", "#ffbf00", "#9fe2bf", "#6495ed", "ffee00"};

  // mass is "enabled" is there is at least one signal file with that mass
  int cindex = -1;

  if (!ProcessOnlySignal){
    for (int m : enabled_masses){
         
    
    cindex++;

    if (iopt == 2 && cindex > 0) continue;

    Outfile<<Headers<<endl;
    
    for (TString bkg : std::vector<TString>{"Zbb","Zcc","Zuds","Zmumu","Ztautau","munuqq"} ){

      TString mstr = Form("%d",m);
      Double_t w = Weight(bkg);

      if (iopt == 1) V = CutFlowOK(bkg, m, "n/a", AnalysisResPath, RunOnN, false, jalg,"> d2d dsigma jalg")[m];
      if (iopt == 2) V = CutFlowOK(bkg, m, "n/a", AnalysisResPath, RunOnN, true, jalg,"> d2d dsigma jalg")[m];


      TString line;
      

				   
      
      if (iopt == 1)  line = maketableline1(bkg, -1, m,  "",        -1,                xsec(bkg,mstr),         V[0],  w, V[2], V[3], V[dcut_id(4)], V[dcut_id(8)], V[dcut_id(20)], V[dcut_id(50)], V[dcut_id(100)], V[dcut_id(200)], SomeColors[cindex]);

                       

      if (iopt == 2) line = maketableline2(bkg,-1,"",-1,V,"#ffffcc");
      Outfile<<line;
      Outfile<<endl;
      
      
    }
    }
  }

  Outfile<<"</table>"<<endl;
  Outfile.close();

  if (upload) gSystem->Exec("scp *.html nvalle@lxplus.cern.ch:/eos/user/n/nvalle/www/fcc/.");
    

}



