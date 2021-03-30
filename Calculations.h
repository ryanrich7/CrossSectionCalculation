#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

//Relevant constants
const double hbarc = 0.19738; //Units of GeV-fm, 1 fm = 10^-15 m
const double d2r = TMath::Pi()/180;
const double alpha = pow(137,-1);
const double amu = 0.93149432;

double CrossSection(double, int, double);
void ReadTable(string );
double Qsq(double, double, int);
void LoadHor(string, int, int, int);
double Interpolate(double, double, double, double, double);

//For form factor tables from HAPLOG 3009
vector<double> q, FFSq;

//For Horowitz table
vector<double> energy,COMth, Xs, Asym;
//This is to get energy dependence of cross section for a given angle
vector<double> xs1;


int E_min, n_E, n_Th;
int E_step = 50;

int m; 

//E in GeV, th0 in degrees
double CrossSection(double E, int Z, double th0){

  //Mott cross section -- from HAMC formula document
  double Mott;

  double th1 = 0.5*th0*d2r;

  double num = pow(alpha*Z*hbarc*TMath::Cos(th1),2);
  double den = 400*pow(E*TMath::Sin(th1)*TMath::Sin(th1),2);


  Mott = num/den;

  return Mott*1000;//Units of mbars/str to compare with Horowitz 



}

void ReadTable(string filename){


   ifstream data;
   data.open(filename.c_str());

   double Q, FF;

   while( data >> Q >> FF){ q.push_back(pow(Q*hbarc,2)); FFSq.push_back(FF); }


}


double Qsq(double E, double th0, int A){

 double M = A*amu; //target mass in GeV
 double recoil = 1/(1 + (E/M)*(1-TMath::Cos(th0*d2r))); 

 double Eprime = E*recoil;

 double qsq = 2*E*Eprime*(1-TMath::Cos(th0*d2r));


 return qsq;



}


void LoadHor(string Horfile, int E_min, int n_E, int n_Th, double thisAngVal){

 ifstream hor;
 hor.open(Horfile.c_str());
 
 double thisEnergy, COM, xs, asym;
 string dummy;

  for(int i = 0; i < n_E; i++){
     hor >> dummy;
     thisEnergy = E_min + double(i)*E_step;
     energy.push_back(thisEnergy/1000);//In GeV
   for(int j = 0; j < n_Th; j++){ 
     hor >> COM >> xs >> asym;
     COMth.push_back(COM); 
     Xs.push_back(xs); 
     Asym.push_back(asym);
  
    if(COM == thisAngVal) { xs1.push_back(xs); } 


   }
 }


}

double Interpolate(double x, double x0, double y0, double x1, double y1){

   double slope = (y1-y0)/(x1-x0);
   double y = slope*(x-x0)+y0;

   return y;

}

