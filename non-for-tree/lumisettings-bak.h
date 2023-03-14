Double_t LUMI = 1.5e8;

Double_t
Nsignal =  1.e5, //0 
  NZbb =    9.7e8, // 1.e7, //1
  NZcc =     9.7e8, // 1.e7, //2
  NZmumu =   1.e7, //3
  NZtautau = 1.e7, //4
  NZuds =    9.3e8, // 1.e8, //5
  Nmunuqq =  5.e5; //6 

std::map<TString, Double_t> xs_vs_m = {
  {"10",0.3665},
  {"20",0.3397},
  {"30",0.3045},
  {"40",0.2598},
  {"50",0.2061},
  {"60",0.1455},
  {"70",0.08358},
  {"80",0.03137}
};

  
std::vector<TString> id = {"signal","Zbb","Zcc","Zmumu","Ztautau","Zuds","munuqq"};


std::vector<Double_t> Nnocut = {Nsignal, NZbb, NZcc, NZmumu, NZtautau, NZuds, Nmunuqq};

TString AnalysisResults(TString opt){

  
  TString toret;
  if      (opt == "Zbb")     toret = "AnalysisResults_Zbb_highstat.root";
  else if (opt == "Zcc")     toret = "AnalysisResults_Zcc_highstat.root";
  else if (opt == "Zuds")    toret = "AnalysisResults_Zuds_highstat.root";
  else if (opt == "Zmumu")   toret = "AnalysisResults-Zmumu.root";
  else if (opt == "Ztautau") toret = "AnalysisResults-Ztautau.root";
  else if (opt == "munuqq")  toret = "AnalysisResults-munuqq.root";
  else                       toret = Form("AnalysisResults-signal-M-%s.root",opt.Data());

  return toret;
}


Double_t xsec(TString id, TString HNMass = "50"){

  if (id == "Zbb")     return 6645.46;
  if (id == "Zcc")     return 5215.46;
  if (id == "Zmumu")   return 1462.09;
  if (id == "Ztautau") return 1476.58;
  if (id == "Zuds")    return 18616.5;
  if (id == "munuqq")  return 0.003192;
  if (id == "signal")  return xs_vs_m[HNMass];
  return -1.;
}
