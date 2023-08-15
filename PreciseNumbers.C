
#include "CutFlowOK.C"

/*
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
*/

TString RoundWithUnits(double n){
  TString toret;
  if (n < 1000 && (int)n == n) toret = Form("%d",(int)n);
  else if (n < 1) toret = Form("%g",n);
  else if (n < 10) toret = Form("%.2f",n);
  else if (n < 1000) toret = Form("%.1f",n);
  else if (n < 1e6) toret = Form("%.1fk",n/1000.);
  else if (n < 1e9) toret = Form("%.1fM",n/1.e6);
  else if (n < 1e12) toret = Form("%.1fB",n/1.e9);
  else toret = Form("%g",n);
  return toret;
}

TString BkgLatex(TString opt){
  TString toret = opt;
  if (opt == "Zbb") toret = "$Z\\to bb$";
  if (opt == "Zcc") toret = "$Z\\to cc$";
  if (opt == "Zss") toret = "$Z\\to ss$";
  if (opt == "Zud") toret = "$Z\\to u/d$";
  if (opt == "Zuds") toret = "$Z\\to u/d/s$";
  if (opt == "Zmumu") toret = "$Z\\to \\mu\\mu$";
  if (opt == "Ztautau") toret = "$Z\\to \\tau\\tau$";
  if (opt == "munuqq") toret = "$\\mu\\nu qq$";

  return toret;
}

TString SampLatex(TString opt, int mass){

  if (opt != "signal") return BkgLatex(opt);
  else return (TString)Form("$%d \\egev$",mass);

}

void PreciseNumbers(TString opt="signal", Int_t mass=80, TString lifetime = "n/a", TString dir="../MyExternalAnalysis/results/"){

  // analysis options:
  //  only1L ,  only2M  ,  combined_selection anymass1L2M,  withdcut anymass1L2M,


  
  Int_t jalg = 2;

  Int_t Counters[100];
  for (int i=0;i<100;i++) Counters[i] = 0;

  Int_t PreciseNum1[100], PreciseNum2[100]; // 1 and 2 jets
  for (int i=0;i<100;i++) PreciseNum1[i] = PreciseNum2[i] = 0;

 

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

  
  
  
  TREESetBranch();
  

  Long64_t EntriesTree = TREE->GetEntries();

  

  Long64_t NN = EntriesTree;
  Double_t SF = 1.;
  cout<<"CutFlowTable.C:: Running on "<<NN<<" entries out of "<<EntriesTree<<endl;
  //cout<<"CutFlowTable.C:: oNEntries = "<<oNEntries<<endl;
  cout<<"CutFlowTable.C:: Scale factor is "<<SF<<endl;

  

  Double_t Ngen = Nnocut(opt,Form("%d",mass),lifetime) * SF;
  Int_t NOneMuon = NN;
  Int_t NTwoJets = 0;
  Int_t NOneJet = 0;

  PreciseNum1[0] = PreciseNum2[0] = Ngen;
  PreciseNum1[1] = PreciseNum2[1] = NN;
  
 

  for (Long64_t i = 0; i < NN; i++){

    

    if (i%1000000 == 0) cout<<i<<" / "<<NN<<endl;

    TREE->GetEntry(i);

    if (i==0)
      cout<<"CutFlowTable.C:: oNEntries = "<<oNEntries<<endl;

    if (oNMuon != 1  || oNJet->at(jalg) < 1 || oEMiss < 0) {Ngen -= SF; NOneMuon--;  continue;} // this conditions should be always satifies in the root tree
    if (oNElectron > 0) continue;

    PreciseNum1[2]++;
    PreciseNum2[2]++;

    

    

    // this builds the lorentz vectors
    BUILD_DERIVATE(jalg);



   
    
   

    if (oNJet->at(jalg) == 2){


      PreciseNum2[3]++;

      NTwoJets++;
      
      //                       0             1           2        3         4            5            6          7
      // Header: two_jets, cos(pmiss), cos(pmiss,mu),  Ej+Mj,  cos(jj),  cos(j,mu)<,  cos(j,mu)>,   Mtot,    COMBINED
      // All divided by the number of events with OneMuon && TwoJets (but the first, which is NTwoJets/NOneMuon)

      Double_t cosjj = TMath::Cos(lvj1.Angle(lvj2.Vect()));
      Double_t cosj1mu = TMath::Cos(lvj1.Angle(lvmu.Vect()));
      Double_t cosj2mu = TMath::Cos(lvj2.Angle(lvmu.Vect()));
    
      Bool_t sel1 = (TMath::Abs(lvmiss.CosTheta()) < 0.94);
      if (!sel1) continue;
      PreciseNum2[4]++;
      
      Bool_t sel2 = (TMath::Cos(lvmiss.Angle(lvmu.Vect())) < 0.80);
      if (!sel2) continue;
      PreciseNum2[5]++;
      
      Bool_t sel3 = (oEJet1->at(jalg) >= 3. && oEJet2->at(jalg) >= 3.);
      if (!sel3) continue;
      PreciseNum2[6]++;
      
      Bool_t sel41 = cosjj > -0.8;
      if (!sel41) continue;
      PreciseNum2[7]++;
      
      Bool_t sel42 = cosjj < 0.98;
      if (!sel41) continue;
      PreciseNum2[8]++;
      
      Bool_t sel6 = (TMath::Max( cosj1mu, cosj2mu ) < 0.8);
      if (!sel6) continue;
      PreciseNum2[9]++;

    
      Bool_t s22 = (TMath::Min( cosj1mu, cosj2mu ) > -0.98);
      if (!s22) continue;
      PreciseNum2[10]++;

      Bool_t s21 = (TMath::Min(lvj1.M(), lvj2.M()) > 0.2 && TMath::Min(lvj1.M2(), lvj2.M2()) > 0);
      if (!s21) continue;
      PreciseNum2[11]++;
      

      
      Bool_t s23 = ((lvj1 + lvj2 + lvmiss + lvmu).M() > 80.);
      if (!s23) continue;
      PreciseNum2[12]++;

      
    }


    if (oNJet->at(jalg) == 1){

     

      PreciseNum1[3]++;
      
      NOneJet++;
      
      //                       0             1           2        3         4            5            6          
      // Header: one_jets, cos(pmiss), cos(pmiss,mu),  Ej+Mj, cos(j,mu)<,  cos(j,mu)>,   Mtot,    COMBINED
      // All divided by the number of events with OneMuon && OneJets (but the first, which is NOneJets/NOneMuon)

      Double_t cosjmu1 = TMath::Cos(lvj1.Angle(lvmu.Vect()));
      Double_t cosjmu2 = TMath::Cos(lvj2.Angle(lvmu.Vect()));
      Double_t cosjmu = TMath::Max(cosjmu1, cosjmu2);
      
      Bool_t sel1 = (TMath::Abs(lvmiss.CosTheta()) < 0.94);
      if (!sel1) continue;
      PreciseNum1[4]++;
      
      Bool_t sel2 = (TMath::Cos(lvmiss.Angle(lvmu.Vect())) < 0.50);
      if (!sel2) continue;
      PreciseNum1[5]++;
      
      Bool_t sel3 = (oEJet1->at(jalg) >= 3.);
      if (!sel3) continue;
      PreciseNum1[6]++;
   
      Bool_t sel4low = (-0.5 < cosjmu);
      if (!sel4low) continue;
      PreciseNum1[7]++;
      
      Bool_t sel4up = (cosjmu < 0.96);
      if (!sel4up) continue;
      PreciseNum1[8]++;
      
      Bool_t sel5 = ((lvvis+lvmiss).M() > 80);
      if (!sel5) continue;
      PreciseNum1[9]++;

    

 
    }


   




    
  } // loop over events



  cout<<"==================="<<endl;
  cout<<opt<<endl;
  cout<<"==================="<<endl;
  cout<<"Ngen                   \t "<<PreciseNum1[0]<<endl;
  cout<<"Filter & 1 mu          \t "<<PreciseNum1[1]<<endl;
  cout<<"& nlep==1              \t "<<PreciseNum1[2]<<endl;
  cout<<"------------------"<<endl;
  cout<<"Filter, nlep=1, 2 jets \t "<<PreciseNum2[3]<<endl;
  cout<<"& |cos th_miss|<0.94   \t "<<PreciseNum2[4]<<endl;
  cout<<"& cos_mu_emiss < 0.8   \t "<<PreciseNum2[5]<<endl;
  cout<<"& Ejet >= 3            \t "<<PreciseNum2[6]<<endl;
  cout<<"& cosjj > -0.8         \t "<<PreciseNum2[7]<<endl;
  cout<<"& cosjj < 0.98         \t "<<PreciseNum2[8]<<endl;
  cout<<"& max_cosjmu < 0.8     \t "<<PreciseNum2[9]<<endl;
  cout<<"& min_cosjmu > -0.98   \t "<<PreciseNum2[10]<<endl;
  cout<<"& Mj>0.2 and M2j>0     \t "<<PreciseNum2[11]<<endl;
  cout<<"& Mvis > 80            \t "<<PreciseNum2[12]<<endl;
  cout<<"------------------"<<endl;
  cout<<"Filter, nlep=1, 1 jet  \t "<<PreciseNum1[3]<<endl;
  cout<<"& |cos th_miss|<0.94   \t "<<PreciseNum1[4]<<endl;
  cout<<"& cos_mu_emiss < 0.5   \t "<<PreciseNum1[5]<<endl;
  cout<<"& Ejet >= 3            \t "<<PreciseNum1[6]<<endl;
  cout<<"& cosjmu > -0.5        \t "<<PreciseNum1[7]<<endl;
  cout<<"& cosjmu < 0.96        \t "<<PreciseNum1[8]<<endl;
  cout<<"& Mvis > 80            \t "<<PreciseNum1[9]<<endl;
 


  exit(0);

  
			 
}

    
