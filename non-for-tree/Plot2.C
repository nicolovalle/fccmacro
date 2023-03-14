#include "lumisettings.h"

void Plot2(TString cut_histo="selection2 Mjjmu", TString HNMass="50", TString path="../MyExternalAnalysis/results/", Bool_t Scale=true){

  cout<<"14"<<endl;

  TCanvas *c1 = new TCanvas("c1","c1",0,0,600,400);

  THStack *HS = new THStack();

  TH1F *hsig;

  TH1F *g;

  Int_t counter = 0;

  TString ttitle, thname;

  TString objname_in_file = cut_histo;
  
  if (objname_in_file.Contains("sliding"))
    objname_in_file.ReplaceAll(" ",Form("_m%s/",HNMass.Data()));
  else
    objname_in_file.ReplaceAll(" ","/");
  
  auto legend = new TLegend();
  
  std::vector<TString> id2  = {"signal","munuqq","Ztautau","Zmumu","Zuds","Zbb","Zcc"};

  std::vector<TH1F*> VH;
  VH.clear();

  for (TString ID : id2){

    
    
    
    TString opt = (ID == "signal") ? HNMass : ID;
    TString filename = (TString)Form("%s%s",path.Data(),AnalysisResults(opt).Data());
    cout<<"Opening "<<filename<<endl;

    TFile *f = new TFile(filename);

    VH.push_back( (TH1F*)f->Get(objname_in_file)->Clone(ID));

    ttitle = Form("%s - %s",VH.at(VH.size()-1)->GetTitle(),objname_in_file.Data());
    thname = ttitle;
    
    VH.at(VH.size()-1)->SetDirectory(nullptr);
    VH.at(VH.size()-1)->SetFillColor(counter+1);
    VH.at(VH.size()-1)->SetTitle(ID);
    VH.at(VH.size()-1)->Rebin(5);
    if (Scale){
      Double_t factor = xsec(ID,HNMass)*LUMI / Nnocut(ID);
      VH.at(VH.size()-1)->Scale(factor);
      VH.at(VH.size()-1)->SetTitle(Form("%s x %7.2f",ID.Data(),factor));
      //VH.at(VH.size()-1)->SetDrawOption("C");
    }
    

    if (ID == "signal") {
       hsig = (TH1F*) VH.at(VH.size()-1)->Clone("Signal");
       legend->AddEntry(VH.at(VH.size()-1),VH.at(VH.size()-1)->GetTitle(),"l");
       //hsig->Scale(1.e6);
    }
    else {
      HS->Add(VH.at(VH.size()-1));
      legend->AddEntry(VH.at(VH.size()-1),VH.at(VH.size()-1)->GetTitle(),"l");
    }
    

    f->Close();
    f->Delete();

    counter ++;
  }


  cout<<"------------"<<thname<<endl;
  // HS->SetTitle(ttitle);
  
  //HS->Draw("HIST");
  hsig->Draw("HIST same");
  gPad->SetLogy();
  
  //legend->Draw("same");
  //c1->BuildLegend();

  //c1->SaveAs(Form("%s.root",objname_in_file.ReplaceAll("/","-").Data()));

 


}
