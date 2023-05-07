Double_t LUMI = 1.5e8;
//Double_t LUMI = 3.3e3;

Double_t
  Nsignal =  1.e5, //0
  NsignalLT = 1.e4,
  NZbb =    9.8e8, // 1.e7, //1
  NZcc =     9.9e8, // 1.e7, //2
  NZuds =    1.e9, // 1.e8, //5
  NZmumu =   1.e7, //3
  NZtautau = 1.e7, //4
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

  
//std::vector<TString> id = {"signal","Zbb","Zcc","Zuds","Zmumu","Ztautau","munuqq"};



Double_t Nnocut(TString opt, TString HNMass="any", TString lifetime="n/a"){

  if (opt == "signal" && lifetime == "m0p5" && HNMass == "40") return 5000.;
  
  if (opt == "Zbb")     return NZbb;
  else if (opt == "Zcc")     return NZcc;
  else if (opt == "Zmumu")   return NZmumu;
  else if (opt == "Ztautau") return NZtautau;
  else if (opt == "Zuds")    return NZuds;
  else if (opt == "munuqq")  return Nmunuqq;
  
  if (opt == "signal" && lifetime != "n/a") return NsignalLT;
  else if (opt=="signal") return Nsignal;

  return -1;
  
 
}
Double_t Nnocut(TString opt, int HNMass, TString lifetime="n/a"){
  return Nnocut(opt,Form("%d",HNMass),lifetime);
}



TString AnalysisResults(TString opt, TString HNMass = "50", TString lifetime="n/a" ){

  
  TString toret="FILE-NOT-FOUND";
  if      (opt == "Zbb")     toret = "AnalysisResults-Zbb_highstat.root";
  else if (opt == "Zcc")     toret = "AnalysisResults-Zcc_highstat.root";
  else if (opt == "Zuds")    toret = "AnalysisResults-Zuds_highstat.root";
  else if (opt == "Zmumu")   toret = "AnalysisResults-Zmumu.root";
  else if (opt == "Ztautau") toret = "AnalysisResults-Ztautau.root";
  else if (opt == "munuqq")  toret = "AnalysisResults-munuqq.root";
  else if (opt == "signal" && lifetime == "n/a")  toret = Form("AnalysisResults-signal-M-%s.root",HNMass.Data());
  else if (lifetime == "n/a") toret = Form("AnalysisResults-signal-M-%s.root",opt.Data());
  else if (lifetime != "n/a") toret = Form("AnalysisResults-signal_10k_%s_%s.root",HNMass.Data(),lifetime.Data());

  return toret;
}
TString AnalysisResults(TString opt, int HNMass = 50, TString lifetime="n/a" ){
  return AnalysisResults(opt,Form("%d",HNMass),lifetime);
}



Double_t xsec(TString opt, TString HNMass = "50", TString lifetime="n/a"){

  if (opt == "Zbb")           return 6645.46;
  else if (opt == "Zcc")      return 5215.46;
  else if (opt == "Zmumu")    return 1462.09;
  else if (opt == "Ztautau")  return 1476.58;
  else if (opt == "Zuds")     return 18616.5;
  else if (opt == "munuqq")   return 0.003192;
  else if (opt == "signal" && lifetime == "n/a") return xs_vs_m[HNMass];

  else{

    TString spec = Form("HNL_%s_%s",HNMass.Data(),lifetime.Data());

    /* FIRST VALUES USED
    if (spec == "HNL_30_m1p0") return 4.198388e-07;
    if (spec == "HNL_30_m1p5") return 1.3276468e-06;
    if (spec == "HNL_30_m2p0") return 4.198399e-06;
    if (spec == "HNL_30_m2p5") return 1.3276468e-05;
    if (spec == "HNL_30_m3p0") return 4.198388e-05;
    if (spec == "HNL_30_m3p5") return 0.00013276468;
    if (spec == "HNL_30_m4p0") return 0.0004198399;
    if (spec == "HNL_30_m4p5") return 0.0013276468;
    if (spec == "HNL_30_m5p0") return 0.004198388;
    if (spec == "HNL_30_m5p5") return 0.013276468;
    if (spec == "HNL_30_m6p0") return 0.04198399;
    if (spec == "HNL_30_m6p5") return 0.13276468;
    if (spec == "HNL_30_m7p0") return 0.4198388;
    if (spec == "HNL_50_m1p0") return 1.82648e-08;
    if (spec == "HNL_50_m1p5") return 5.776428e-08;
    if (spec == "HNL_50_m2p0") return 1.82648e-07;
    if (spec == "HNL_50_m2p5") return 5.77644e-07;
    if (spec == "HNL_50_m3p0") return 1.82648e-06;
    if (spec == "HNL_50_m3p5") return 5.776428e-06;
    if (spec == "HNL_50_m4p0") return 1.82648e-05;
    if (spec == "HNL_50_m4p5") return 5.77644e-05;
    if (spec == "HNL_50_m5p0") return 0.000182648;
    if (spec == "HNL_50_m5p5") return 0.0005776428;
    if (spec == "HNL_50_m6p0") return 0.00182648;
    if (spec == "HNL_50_m6p5") return 0.00577644;
    if (spec == "HNL_50_m7p0") return 0.0182648;
    if (spec == "HNL_70_m1p0") return 9.2784699423e-10;
    if (spec == "HNL_70_m1p5") return 2.93409692649e-09;
    if (spec == "HNL_70_m2p0") return 9.2784669402e-09;
    if (spec == "HNL_70_m2p5") return 2.93409592649e-08;
    if (spec == "HNL_70_m3p0") return 9.2784699423e-08;
    if (spec == "HNL_70_m3p5") return 2.93409692649e-07;
    if (spec == "HNL_70_m4p0") return 9.2784669402e-07;
    if (spec == "HNL_70_m4p5") return 2.93409592649e-06;
    if (spec == "HNL_70_m5p0") return 9.2784699423e-06;
    if (spec == "HNL_70_m5p5") return 2.93409692649e-05;
    if (spec == "HNL_70_m6p0") return 9.2784669402e-05;
    if (spec == "HNL_70_m6p5") return 0.000293409592649;
    if (spec == "HNL_70_m7p0") return 0.00092784699423;
    */

    if (spec == "HNL_10_m0p5") return 5.026031e-05;
    if (spec == "HNL_10_m1p0") return 0.00015893683;
    if (spec == "HNL_10_m1p5") return 0.0005026031;
    if (spec == "HNL_10_m2p0") return 0.0015893683;
    if (spec == "HNL_10_m2p5") return 0.005026031;
    if (spec == "HNL_10_m3p0") return 0.015893683;
    if (spec == "HNL_10_m3p5") return 0.05026031;
    if (spec == "HNL_10_m4p0") return 0.15893683;
    if (spec == "HNL_10_m4p5") return 0.5026031;
    if (spec == "HNL_10_m5p0") return 1.5893683;
    if (spec == "HNL_10_p0p0") return 1.5893683e-05;
    if (spec == "HNL_10_p0p5") return 5.026031e-06;
    if (spec == "HNL_10_p1p0") return 1.5893683e-06;
    if (spec == "HNL_10_p1p5") return 5.026031e-07;
    if (spec == "HNL_20_m0p5") return 1.3288316e-06;
    if (spec == "HNL_20_m1p0") return 4.202092e-06;
    if (spec == "HNL_20_m1p5") return 1.3288317e-05;
    if (spec == "HNL_20_m2p0") return 4.202092e-05;
    if (spec == "HNL_20_m2p5") return 0.00013288316;
    if (spec == "HNL_20_m3p0") return 0.0004202092;
    if (spec == "HNL_20_m3p5") return 0.0013288317;
    if (spec == "HNL_20_m4p0") return 0.004202092;
    if (spec == "HNL_20_m4p5") return 0.013288316;
    if (spec == "HNL_20_m5p0") return 0.04202092;
    if (spec == "HNL_20_m5p5") return 0.13288317;
    if (spec == "HNL_20_m6p0") return 0.4202092;
    if (spec == "HNL_20_m6p5") return 1.3288316;
    if (spec == "HNL_20_p0p0") return 4.202092e-07;
    if (spec == "HNL_20_p0p5") return 1.3288317e-07;
    if (spec == "HNL_20_p1p0") return 4.202092e-08;
    if (spec == "HNL_2_m1p0") return 0.67516377;
    if (spec == "HNL_30_m0p5") return 1.4580685e-07;
    if (spec == "HNL_30_m1p0") return 4.610882e-07;
    if (spec == "HNL_30_m1p5") return 1.4580685e-06;
    if (spec == "HNL_30_m2p0") return 4.610892e-06;
    if (spec == "HNL_30_m2p5") return 1.4580685e-05;
    if (spec == "HNL_30_m3p0") return 4.610882e-05;
    if (spec == "HNL_30_m3p5") return 0.00014580685;
    if (spec == "HNL_30_m4p0") return 0.0004610892;
    if (spec == "HNL_30_m4p5") return 0.0014580685;
    if (spec == "HNL_30_m5p0") return 0.004610882;
    if (spec == "HNL_30_m5p5") return 0.014580685;
    if (spec == "HNL_30_m6p0") return 0.04610892;
    if (spec == "HNL_30_m6p5") return 0.14580685;
    if (spec == "HNL_30_m7p0") return 0.4610882;
    if (spec == "HNL_30_p0p0") return 4.610892e-08;
    if (spec == "HNL_30_p0p5") return 1.4580685e-08;
    if (spec == "HNL_40_m0p5") return 2.720248e-08;
    if (spec == "HNL_40_m1p0") return 8.602117e-08;
    if (spec == "HNL_40_m1p5") return 2.720248e-07;
    if (spec == "HNL_40_m2p0") return 8.602118e-07;
    if (spec == "HNL_40_m2p5") return 2.720248e-06;
    if (spec == "HNL_40_m3p0") return 8.602117e-06;
    if (spec == "HNL_40_m3p5") return 2.720248e-05;
    if (spec == "HNL_40_m4p0") return 8.602118e-05;
    if (spec == "HNL_40_m4p5") return 0.0002720248;
    if (spec == "HNL_40_m5p0") return 0.0008602117;
    if (spec == "HNL_40_m5p5") return 0.002720248;
    if (spec == "HNL_40_m6p0") return 0.008602118;
    if (spec == "HNL_40_m6p5") return 0.02720248;
    if (spec == "HNL_40_m7p0") return 0.08602117;
    if (spec == "HNL_40_p0p0") return 8.602118e-09;
    if (spec == "HNL_50_m1p0") return 2.002058e-08;
    if (spec == "HNL_50_m1p5") return 6.330149e-08;
    if (spec == "HNL_50_m2p0") return 2.001758e-07;
    if (spec == "HNL_50_m2p5") return 6.33015e-07;
    if (spec == "HNL_50_m3p0") return 2.001748e-06;
    if (spec == "HNL_50_m3p5") return 6.330149e-06;
    if (spec == "HNL_50_m4p0") return 2.001758e-05;
    if (spec == "HNL_50_m4p5") return 6.33015e-05;
    if (spec == "HNL_50_m5p0") return 0.0002001748;
    if (spec == "HNL_50_m5p5") return 0.0006330149;
    if (spec == "HNL_50_m6p0") return 0.002001758;
    if (spec == "HNL_50_m6p5") return 0.00633015;
    if (spec == "HNL_50_m7p0") return 0.02002058;
    if (spec == "HNL_5_m0p5") return 0.001885223;
    if (spec == "HNL_5_m1p0") return 0.005961641;
    if (spec == "HNL_5_m1p5") return 0.01885223;
    if (spec == "HNL_5_m2p0") return 0.05961641;
    if (spec == "HNL_5_m2p5") return 0.1885223;
    if (spec == "HNL_5_m3p0") return 0.5961641;
    if (spec == "HNL_5_m3p5") return 1.885223;
    if (spec == "HNL_5_p0p0") return 0.0005961641;
    if (spec == "HNL_5_p0p5") return 0.0001885223;
    if (spec == "HNL_5_p1p0") return 5.961641e-05;
    if (spec == "HNL_5_p1p5") return 1.885223e-05;
    if (spec == "HNL_5_p2p0") return 5.961641e-06;
    if (spec == "HNL_5_p2p5") return 1.885223e-06;
    if (spec == "HNL_60_m1p5") return 1.5333045e-08;
    if (spec == "HNL_60_m2p0") return 4.848728e-08;
    if (spec == "HNL_60_m2p5") return 1.5333046e-07;
    if (spec == "HNL_60_m3p0") return 4.848738e-07;
    if (spec == "HNL_60_m3p5") return 1.5333045e-06;
    if (spec == "HNL_60_m4p0") return 4.848738e-06;
    if (spec == "HNL_60_m4p5") return 1.5333046e-05;
    if (spec == "HNL_60_m5p0") return 4.848738e-05;
    if (spec == "HNL_60_m5p5") return 0.00015333045;
    if (spec == "HNL_60_m6p0") return 0.0004848738;
    if (spec == "HNL_60_m6p5") return 0.0015333046;
    if (spec == "HNL_60_m7p0") return 0.004848748;
    if (spec == "HNL_70_m2p0") return 9.9566775336e-09;
    if (spec == "HNL_70_m2p5") return 3.1484036525e-08;
    if (spec == "HNL_70_m3p0") return 9.9566805346e-08;
    if (spec == "HNL_70_m3p5") return 3.1484236525e-07;
    if (spec == "HNL_70_m4p0") return 9.9566775336e-07;
    if (spec == "HNL_70_m4p5") return 3.1484036525e-06;
    if (spec == "HNL_70_m5p0") return 9.9566805346e-06;
    if (spec == "HNL_70_m5p5") return 3.1484236525e-05;
    if (spec == "HNL_70_m6p0") return 9.9566775336e-05;
    if (spec == "HNL_70_m6p5") return 0.00031484036525;
    if (spec == "HNL_70_m7p0") return 0.00099566805346;
    if (spec == "HNL_80_m3p0") return 1.0628312e-08;
    if (spec == "HNL_80_m3p5") return 3.3609825e-08;
    if (spec == "HNL_80_m4p0") return 1.0628313e-07;
    if (spec == "HNL_80_m4p5") return 3.3609825e-07;
    if (spec == "HNL_80_m5p0") return 1.0628312e-06;
    if (spec == "HNL_80_m5p5") return 3.3609825e-06;
    if (spec == "HNL_80_m6p0") return 1.0628313e-05;
    if (spec == "HNL_80_m6p5") return 3.3609825e-05;
    if (spec == "HNL_80_m7p0") return 0.00010628312;
    if (spec == "HNL_85_m4p0") return 2.1419744e-08;
    if (spec == "HNL_85_m4p5") return 6.773503e-08;
    if (spec == "HNL_85_m5p0") return 2.1419743e-07;
    if (spec == "HNL_85_m5p5") return 6.773503e-07;
    if (spec == "HNL_85_m6p0") return 2.1419744e-06;
    if (spec == "HNL_85_m6p5") return 6.773503e-06;
    if (spec == "HNL_85_m7p0") return 2.1419743e-05;


    return -1;
    
  }
}
Double_t xsec(TString opt, int HNMass = 50, TString lifetime="n/a"){
  return xsec(opt,Form("%d",HNMass),lifetime);
}

Double_t Weight(TString opt, TString HNMass="specify only for signal", TString lifetime="n/a", Double_t custom_weight = 1.){

  return xsec(opt,HNMass,lifetime)*LUMI/(Nnocut(opt,HNMass,lifetime) * custom_weight);
 
}


Double_t Coupling(TString HNMass, TString lifetime){

  TString spec = Form("HNL_%s_%s",HNMass.Data(),lifetime.Data());

  /* FIRST VALUES
  if (spec == "HNL_30_m1p0") return 1.175106e-05;
  if (spec == "HNL_30_m1p5") return 2.089667e-05;
  if (spec == "HNL_30_m2p0") return 3.716012e-05;
  if (spec == "HNL_30_m2p5") return 6.608108e-05;
  if (spec == "HNL_30_m3p0") return 1.175106e-04;
  if (spec == "HNL_30_m3p5") return 2.089667e-04;
  if (spec == "HNL_30_m4p0") return 3.716012e-04;
  if (spec == "HNL_30_m4p5") return 6.608108e-04;
  if (spec == "HNL_30_m5p0") return 1.175106e-03;
  if (spec == "HNL_30_m5p5") return 2.089667e-03;
  if (spec == "HNL_30_m6p0") return 3.716012e-03;
  if (spec == "HNL_30_m6p5") return 6.608108e-03;
  if (spec == "HNL_30_m7p0") return 1.175106e-02;
  if (spec == "HNL_50_m1p0") return 2.974684e-06;
  if (spec == "HNL_50_m1p5") return 5.289819e-06;
  if (spec == "HNL_50_m2p0") return 9.406776e-06;
  if (spec == "HNL_50_m2p5") return 1.672788e-05;
  if (spec == "HNL_50_m3p0") return 2.974684e-05;
  if (spec == "HNL_50_m3p5") return 5.289819e-05;
  if (spec == "HNL_50_m4p0") return 9.406776e-05;
  if (spec == "HNL_50_m4p5") return 1.672788e-04;
  if (spec == "HNL_50_m5p0") return 2.974684e-04;
  if (spec == "HNL_50_m5p5") return 5.289819e-04;
  if (spec == "HNL_50_m6p0") return 9.406776e-04;
  if (spec == "HNL_50_m6p5") return 1.672788e-03;
  if (spec == "HNL_50_m7p0") return 2.974684e-03;
  if (spec == "HNL_70_m1p0") return 1.053186e-06;
  if (spec == "HNL_70_m1p5") return 1.872860e-06;
  if (spec == "HNL_70_m2p0") return 3.330468e-06;
  if (spec == "HNL_70_m2p5") return 5.922503e-06;
  if (spec == "HNL_70_m3p0") return 1.053186e-05;
  if (spec == "HNL_70_m3p5") return 1.872860e-05;
  if (spec == "HNL_70_m4p0") return 3.330468e-05;
  if (spec == "HNL_70_m4p5") return 5.922503e-05;
  if (spec == "HNL_70_m5p0") return 1.053186e-04;
  if (spec == "HNL_70_m5p5") return 1.872860e-04;
  if (spec == "HNL_70_m6p0") return 3.330468e-04;
  if (spec == "HNL_70_m6p5") return 5.922503e-04;
  if (spec == "HNL_70_m7p0") return 1.053186e-03;
  */

  if (spec == "HNL_10_m0p5") return 1.116472e-04;
  if (spec == "HNL_10_m1p0") return 1.985399e-04;
  if (spec == "HNL_10_m1p5") return 3.530594e-04;
  if (spec == "HNL_10_m2p0") return 6.278383e-04;
  if (spec == "HNL_10_m2p5") return 1.116472e-03;
  if (spec == "HNL_10_m3p0") return 1.985399e-03;
  if (spec == "HNL_10_m3p5") return 3.530594e-03;
  if (spec == "HNL_10_m4p0") return 6.278383e-03;
  if (spec == "HNL_10_m4p5") return 1.116472e-02;
  if (spec == "HNL_10_m5p0") return 1.985399e-02;
  if (spec == "HNL_10_p0p0") return 6.278383e-05;
  if (spec == "HNL_10_p0p5") return 3.530594e-05;
  if (spec == "HNL_10_p1p0") return 1.985399e-05;
  if (spec == "HNL_10_p1p5") return 1.116472e-05;
  if (spec == "HNL_20_m0p5") return 1.884013e-05;
  if (spec == "HNL_20_m1p0") return 3.350301e-05;
  if (spec == "HNL_20_m1p5") return 5.957771e-05;
  if (spec == "HNL_20_m2p0") return 1.059458e-04;
  if (spec == "HNL_20_m2p5") return 1.884013e-04;
  if (spec == "HNL_20_m3p0") return 3.350301e-04;
  if (spec == "HNL_20_m3p5") return 5.957771e-04;
  if (spec == "HNL_20_m4p0") return 1.059458e-03;
  if (spec == "HNL_20_m4p5") return 1.884013e-03;
  if (spec == "HNL_20_m5p0") return 3.350301e-03;
  if (spec == "HNL_20_m5p5") return 5.957771e-03;
  if (spec == "HNL_20_m6p0") return 1.059458e-02;
  if (spec == "HNL_20_m6p5") return 1.884013e-02;
  if (spec == "HNL_20_p0p0") return 1.059458e-05;
  if (spec == "HNL_20_p0p5") return 5.957771e-06;
  if (spec == "HNL_20_p1p0") return 3.350301e-06;
  if (spec == "HNL_2_m1p0") return 1.396039e-02;
  if (spec == "HNL_30_m0p5") return 6.608108e-06;
  if (spec == "HNL_30_m1p0") return 1.175106e-05;
  if (spec == "HNL_30_m1p5") return 2.089667e-05;
  if (spec == "HNL_30_m2p0") return 3.716012e-05;
  if (spec == "HNL_30_m2p5") return 6.608108e-05;
  if (spec == "HNL_30_m3p0") return 1.175106e-04;
  if (spec == "HNL_30_m3p5") return 2.089667e-04;
  if (spec == "HNL_30_m4p0") return 3.716012e-04;
  if (spec == "HNL_30_m4p5") return 6.608108e-04;
  if (spec == "HNL_30_m5p0") return 1.175106e-03;
  if (spec == "HNL_30_m5p5") return 2.089667e-03;
  if (spec == "HNL_30_m6p0") return 3.716012e-03;
  if (spec == "HNL_30_m6p5") return 6.608108e-03;
  if (spec == "HNL_30_m7p0") return 1.175106e-02;
  if (spec == "HNL_30_p0p0") return 3.716012e-06;
  if (spec == "HNL_30_p0p5") return 2.089667e-06;
  if (spec == "HNL_40_m0p5") return 3.089376e-06;
  if (spec == "HNL_40_m1p0") return 5.493774e-06;
  if (spec == "HNL_40_m1p5") return 9.769465e-06;
  if (spec == "HNL_40_m2p0") return 1.737284e-05;
  if (spec == "HNL_40_m2p5") return 3.089376e-05;
  if (spec == "HNL_40_m3p0") return 5.493774e-05;
  if (spec == "HNL_40_m3p5") return 9.769465e-05;
  if (spec == "HNL_40_m4p0") return 1.737284e-04;
  if (spec == "HNL_40_m4p5") return 3.089376e-04;
  if (spec == "HNL_40_m5p0") return 5.493774e-04;
  if (spec == "HNL_40_m5p5") return 9.769465e-04;
  if (spec == "HNL_40_m6p0") return 1.737284e-03;
  if (spec == "HNL_40_m6p5") return 3.089376e-03;
  if (spec == "HNL_40_m7p0") return 5.493774e-03;
  if (spec == "HNL_40_p0p0") return 1.737284e-06;
  if (spec == "HNL_50_m1p0") return 2.979364e-06;
  if (spec == "HNL_50_m1p5") return 5.298141e-06;
  if (spec == "HNL_50_m2p0") return 9.421575e-06;
  if (spec == "HNL_50_m2p5") return 1.675419e-05;
  if (spec == "HNL_50_m3p0") return 2.979364e-05;
  if (spec == "HNL_50_m3p5") return 5.298141e-05;
  if (spec == "HNL_50_m4p0") return 9.421575e-05;
  if (spec == "HNL_50_m4p5") return 1.675419e-04;
  if (spec == "HNL_50_m5p0") return 2.979364e-04;
  if (spec == "HNL_50_m5p5") return 5.298141e-04;
  if (spec == "HNL_50_m6p0") return 9.421575e-04;
  if (spec == "HNL_50_m6p5") return 1.675419e-03;
  if (spec == "HNL_50_m7p0") return 2.979364e-03;
  if (spec == "HNL_5_m0p5") return 6.821484e-04;
  if (spec == "HNL_5_m1p0") return 1.213050e-03;
  if (spec == "HNL_5_m1p5") return 2.157143e-03;
  if (spec == "HNL_5_m2p0") return 3.836002e-03;
  if (spec == "HNL_5_m2p5") return 6.821484e-03;
  if (spec == "HNL_5_m3p0") return 1.213050e-02;
  if (spec == "HNL_5_m3p5") return 2.157143e-02;
  if (spec == "HNL_5_p0p0") return 3.836002e-04;
  if (spec == "HNL_5_p0p5") return 2.157143e-04;
  if (spec == "HNL_5_p1p0") return 1.213050e-04;
  if (spec == "HNL_5_p1p5") return 6.821484e-05;
  if (spec == "HNL_5_p2p0") return 3.836002e-05;
  if (spec == "HNL_5_p2p5") return 2.157143e-05;
  if (spec == "HNL_60_m1p5") return 3.109415e-06;
  if (spec == "HNL_60_m2p0") return 5.529409e-06;
  if (spec == "HNL_60_m2p5") return 9.832835e-06;
  if (spec == "HNL_60_m3p0") return 1.748553e-05;
  if (spec == "HNL_60_m3p5") return 3.109415e-05;
  if (spec == "HNL_60_m4p0") return 5.529409e-05;
  if (spec == "HNL_60_m4p5") return 9.832835e-05;
  if (spec == "HNL_60_m5p0") return 1.748553e-04;
  if (spec == "HNL_60_m5p5") return 3.109415e-04;
  if (spec == "HNL_60_m6p0") return 5.529409e-04;
  if (spec == "HNL_60_m6p5") return 9.832835e-04;
  if (spec == "HNL_60_m7p0") return 1.748553e-03;
  if (spec == "HNL_70_m2p0") return 3.330468e-06;
  if (spec == "HNL_70_m2p5") return 5.922503e-06;
  if (spec == "HNL_70_m3p0") return 1.053186e-05;
  if (spec == "HNL_70_m3p5") return 1.872860e-05;
  if (spec == "HNL_70_m4p0") return 3.330468e-05;
  if (spec == "HNL_70_m4p5") return 5.922503e-05;
  if (spec == "HNL_70_m5p0") return 1.053186e-04;
  if (spec == "HNL_70_m5p5") return 1.872860e-04;
  if (spec == "HNL_70_m6p0") return 3.330468e-04;
  if (spec == "HNL_70_m6p5") return 5.922503e-04;
  if (spec == "HNL_70_m7p0") return 1.053186e-03;
  if (spec == "HNL_80_m3p0") return 5.707272e-06;
  if (spec == "HNL_80_m3p5") return 1.014912e-05;
  if (spec == "HNL_80_m4p0") return 1.804798e-05;
  if (spec == "HNL_80_m4p5") return 3.209435e-05;
  if (spec == "HNL_80_m5p0") return 5.707272e-05;
  if (spec == "HNL_80_m5p5") return 1.014912e-04;
  if (spec == "HNL_80_m6p0") return 1.804798e-04;
  if (spec == "HNL_80_m6p5") return 3.209435e-04;
  if (spec == "HNL_80_m7p0") return 5.707272e-04;
  if (spec == "HNL_85_m4p0") return 1.087665e-05;
  if (spec == "HNL_85_m4p5") return 1.934173e-05;
  if (spec == "HNL_85_m5p0") return 3.439499e-05;
  if (spec == "HNL_85_m5p5") return 6.116391e-05;
  if (spec == "HNL_85_m6p0") return 1.087665e-04;
  if (spec == "HNL_85_m6p5") return 1.934173e-04;
  if (spec == "HNL_85_m7p0") return 3.439499e-04;


  return -1;
  
}
Double_t Coupling(int HNMass, TString lifetime){
  return Coupling(Form("%d",HNMass), lifetime);
}


void DrawCol()
{
   Int_t i,n;
   Double_t x,y;
   TLatex *l;

   TGraph *g = (TGraph*)gPad->GetListOfPrimitives()->FindObject("Graph");
   n = g->GetN();
   TMarker *m;
   for (i=1; i<n; i++) {
      g->GetPoint(i,x,y);
      m = new TMarker(x,y,(i%17)+51);
      m->SetMarkerSize(1.5);
      // m->SetMarkerStyle((i%17)+51);;
      m->Paint();
   }
}



void MapPoint(){

  Double_t x[1000], y[10000];
  Int_t npoint = 0;

  for (int m = 0; m < 95; m++){
    for (TString mp : std::vector<TString>{"m","p"}){
      for (int unit = 0; unit < 9; unit++){
	for (int decimal = 0; decimal < 6; decimal += 5){

	  TString lt = Form("%s%dp%d", mp.Data(), unit, decimal);

	  if (Coupling(Form("%d",m), lt) > 0){
	    x[npoint] = 1.*m;
	    //y[npoint] = ( mp == "m" ? -1. : 1.)*(unit+0.1*decimal);
	    y[npoint] = 2.*TMath::Log10(Coupling(Form("%d",m), lt));
	    npoint++;
	  }
	  
	}
      }
    }
  }

  cout<<"N points = "<<npoint<<endl;
  
  auto g = new TGraph(npoint,x,y);
  g->SetMarkerStyle(8);
  g->SetMarkerSize(1);

  TExec *ex = new TExec("ex","DrawCol();");

   g->GetListOfFunctions()->Add(ex);
  g->GetYaxis()->SetRangeUser(-12,-2);
   //g->GetYaxis()->SetRangeUser(-8,3);
  g->GetXaxis()->SetLimits(0,90);
  g->GetXaxis()->SetTitle("M_{HN} (GeV)");
  g->GetYaxis()->SetTitle("Log (U^{2})");
  //g->GetYaxis()->SetTitle("Log (decay length / m)");
  g->SetMarkerColor(0);
  g->SetMarkerSize(0);
  g->SetTitle("");
  g->Draw("A P");


}


