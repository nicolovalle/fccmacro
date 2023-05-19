double eeta(double a, double b, double c, double d){

  TLorentzVector l;
  l.SetPxPyPzE(a,b,c,d);
  return l.Eta();
}
