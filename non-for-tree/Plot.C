#include "lumisettings.h"

void Plot(TString cut_histo="preselection Mjjmu", TString HNMass="50", TString path="../MyExternalAnalysis/results/", Bool_t Scale=true){

  

  TCanvas *c1 = new TCanvas("c1","c1",0,0,600,400);

  THStack *HS = new THStack("HS","");

  TH1F *g;

  Int_t counter = 0;

  TString ttitle, thname;

  TString objname_in_file = cut_histo;
  
  if (objname_in_file.Contains("sliding"))
    objname_in_file.ReplaceAll(" ",Form("_m%s/",HNMass.Data()));
  else
    objname_in_file.ReplaceAll(" ","/");
  
  auto legend = new TLegend();
  
  for (TString ID : id){

    
    TString opt = (ID == "signal") ? HNMass : ID;
    TString filename = (TString)Form("%s%s",path.Data(),AnalysisResults(opt).Data());
    cout<<"Opening "<<filename<<endl;

    TFile *f = new TFile(filename);

    g = (TH1F*)f->Get(objname_in_file)->Clone(ID);

    ttitle = Form("%s - %s",g->GetTitle(),objname_in_file.Data());
    thname = ttitle;
    
    g->SetDirectory(nullptr);
    g->SetLineColor(counter+1);
    g->SetTitle(ID);
    if (Scale){
      Double_t factor = xsec(ID,HNMass)*LUMI / Nnocut(ID);
      g->Scale(factor);
      g->SetTitle(Form("%s x %7.2f",ID.Data(),factor));
      //g->SetDrawOption("C");
    }
    

    HS->Add(g);
    legend->AddEntry(g,g->GetTitle(),"l");

    f->Close();
    f->Delete();

    counter ++;
  }


  cout<<"------------"<<thname<<endl;
  HS->SetTitle(ttitle);
  
  HS->Draw("c nostack");
  legend->Draw("same");
  //c1->BuildLegend();

  c1->SaveAs(Form("%s.root",objname_in_file.ReplaceAll("/","-").Data()));

 


}
