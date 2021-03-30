#include "Calculations.h"


//Pick one
//#define Carbon
#define Lead

string targ;

void CompareCross(double th0){

#ifdef Carbon
E_min = 900;
n_E = 28;
n_Th = 71;
ReadTable("carbon12.txt");
LoadHor("c12_FSU.dat",E_min,n_E,n_Th,th0);
const int Z = 6;
const int A = 12;
targ = "C12";
#endif

#ifdef Lead
E_min = 550;
n_E = 35;
n_Th = 66;
ReadTable("pb208.txt");
LoadHor("lead.dat",E_min,n_E,n_Th,th0);
const int Z = 82;
const int A = 208;
targ = "Pb208";
#endif


  int m = energy.size();

  //[0] is for Horowitz, [1] - Mott*FF2
  TGraph *gr[2];

  vector<double> Mott, qsq, FF, sigma;
  double thisE, thisMott, thisQsq; 

  for(int i = 0; i < n_E; i++){ 
   
   thisE = E_min + double(i)*E_step;
   //Calculating the Mott cross section and Q^2
   thisMott = CrossSection(thisE/1000, Z, th0);
   thisQsq = Qsq(thisE/1000,th0,A);   

   Mott.push_back(thisMott);
   qsq.push_back(thisQsq);
  


  }     


  //Now we have Qsq. We now need to calculate FF^2. I will use a linear interpolation
  //Recall q is vector from the FF^2 table 

 //for(int j = 0; j < qsq.size(); j++){ 
  
  for(int i = 0; i < qsq.size(); i++){
   
   //Searches for the first value of qsq which is larger or equal to actual value that we compute    
   vector<double>::iterator iterQsq = upper_bound(q.begin(),q.end(),qsq.at(i));  
   
   int indxQsq = iterQsq - q.begin();    

   if(indxQsq == 0)indxQsq = 1; //use min value
   if(indxQsq <= 0 || indxQsq >= q.size()) return 0;


   double q1, q2, f1, f2;

   q1 = q.at(indxQsq);
   q2 = q.at(indxQsq-1);
   
   f1 = FFSq.at(indxQsq);
   f2 = FFSq.at(indxQsq-1);

//   cout << q1 << "  " << f1 << "  " << q2 << "   " << f2 << endl;


   FF.push_back(Interpolate(qsq.at(i),q2,f2,q1,f1));

  sigma.push_back(FF.at(i)*Mott.at(i));
   


  }


  gr[0] = new TGraph(m, &energy[0], &xs1[0]);
  gr[1] = new TGraph(m, &energy[0],&sigma[0]);
  int color[2] = {2,4};


  for(int i = 0; i < 2; i++){ 
   gr[i]->SetMarkerStyle(20);
   gr[i]->SetMarkerColor(color[i]);
   gr[i]->SetTitle(Form("#frac{d#sigma}{d#Omega} vs. E for Fixed Angle = %.1f^{o},%s",targ.c_str(),th0)); 
   gr[i]->GetXaxis()->SetTitle("E (GeV/c)");
   gr[i]->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega} (mbarns/str)");

 } 
 
  TCanvas *c = new TCanvas();
  gr[0]->Draw("AP");
  gr[1]->Draw("Psame");
  c->SaveAs(Form("CrossFor%s_%.1f.jpg",targ.c_str(),th0)); 



}
