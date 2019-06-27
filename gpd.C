//////
// ROOT Macro for calculating and plotting gpd, w/ parameters
// Kyle Pressler
// SIWIF
//////
// How to use me.
// From root prompt i.e. root [0], use command,
// .L gpd.C+   (this will rebuild the library in the root interpreter)
// start()     (this calls our starting function, it can be named whatever, here it is start)
// printgpd1D()  (this will print the gpd from 0 to 1 with xi,t defined in arguments
// printgpd3D() (prints to .dat file gpd over all three variables)
//////
// Right now this function only evaluates H at a given point. Once I am sure that it can
// do that correctly I wish to be able to save to an array and file lots of data points to be plotted.
//////

#include "TROOT.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <fstream>

Double_t Integrate(Double_t (*func)(Double_t *,Double_t *), Double_t from, Double_t to, Double_t n, Double_t *x, Double_t *par, Int_t r);// These are function declarations, the functions themselves are defined below.
Double_t hint(Double_t *x, Double_t *par);											 // Function declarations tell the compiler, in this case the ROOT Cling interpreter
Double_t exactint(Double_t *x, Double_t *par);	
Double_t eint(Double_t *x, Double_t *par);
Double_t exact(Double_t *x, Double_t *par);									                 // what kinds of functions will be called, it tells them their names and what arguments
Double_t h(Double_t *x, Double_t *par);											 // they accept. For example, hint() is passed two pointers, *x, and *par. When the function
Double_t e(Double_t *x, Double_t *par);  											 // is called later, it will be given what should be passed along, namely two arrays, x, and para.
 																 // which contain the variables and parameters.
void start() // This is the main function.                                                                                
{
	//Double_t x[4] = {.5,0,0,0}; // Array containing variables, (X,Xi,t,k)
	//Double_t para[3] = {0.367,0.583,0.963}; // Array containing parameters, (m,Mx,Mlambda)

	return;
}

Double_t Integrate(Double_t (*func)(Double_t *,Double_t *), Double_t from, Double_t to, Double_t n, Double_t *x, Double_t *par, Int_t r) // This function will integrate any math function of four variables and n-many parameters
{																// using Simpson's method.
	Double_t result;													
	Double_t h = (to - from) / n;
 	Double_t sum1 = 0.0;
	Double_t sum2 = 0.0;
	Double_t y[4] = {x[0],x[1],x[2],x[3]};
	
	for(Int_t i = 0;i < n ;i++)
	{
		y[r] = (from + h * i + h / 2.0);
		sum1 += func(y,par);
	}
	for(Int_t i = 1;i < n;i++)
	{
		y[r] = (from + h * i + h / 2.0);
		sum2 += func(y,par);
	}
	
	y[r] = from;
	Double_t a = func(y,par);
	y[r] = to;
	Double_t b = func(y,par);

	result = h / 6.0 * (a + b + 4.0 * sum1 + 2.0 * sum2); 
	
	return result;
}

Double_t hint(Double_t *x, Double_t *par) // This function defines the "integrand" to be integrated in the GPD H.
{
	Float_t xx = x[0]; // Take the values of the variables out of the array x[-]
	Float_t ss = x[1]; // and put them into new float variables named xx,ss,etc.
	Float_t tt = x[2]; //
	Float_t kk = x[3]; //

	Float_t m = par[0];						       // These simply define all of the relevant constants used in
	Float_t mx = par[1];						       // the definition of H.
	Float_t mv = par[2];						       // In particular the three parameters are pulled out of the 
	Float_t M = .938;						       // array passed to the function.
	Float_t MU2 = xx*pow(M,2) - (xx*pow((1-xx),-1)*pow(mx,2)) - pow(mv,2); //
	Float_t xp = (xx-ss)/(1-ss);					       //
	Float_t Dperp = pow(-1*pow(M,2)*pow(ss,2)*pow((1-ss),-1) - tt,0.5);    //
	Float_t a = (m + M*xx);						       //
	Float_t b = (m + M*xp);						       //

	Double_t D0 = ((1-xx)*(MU2)) - pow(kk,2); // Same as above
	Double_t D1 = (1-xp)*MU2 - pow(kk,2) - pow(Dperp,2) * pow((1-xp),2);
	Double_t D2 = (1-xp)*kk*Dperp;

	Double_t num, denom; // I split the integrand into a numerator and denominator to be able to check that both are outputting the right values.
	
	num = ((((a*b) + pow(kk,2))*D1) + (-2*pow(D2,2)));
	denom = pow(( pow(D1,2) + (-4*pow(D2,2))),1.5)*pow(D0,2);

	//cout << D0 << endl;
	//cout << denom << endl;

	Double_t result;

	result = kk * (num/denom);

	return result;
}

Double_t exactint(Double_t *x, Double_t *par) // Function that defines the integrand of the "Analytically Solved" form of H(X,0,0) given in line 42 of Liuti paper.
{
	Float_t xx = x[0]; 
	Float_t kk = x[3];

	Float_t m = par[0];  // Parameters
	Float_t mx = par[1]; //
	Float_t mv = par[2]; //
	
	Float_t M = .938; 
	Float_t MU2 = xx*pow(M,2) - (xx*pow((1-xx),-1)*pow(mx,2)) - pow(mv,2);
	Float_t a = (m + M*xx);

	Double_t result;
	
	result = kk * (pow(a,2)+pow(kk,2))/pow((1-xx)*MU2-pow(kk,2),4);

	return result;
}

Double_t eint(Double_t *x, Double_t *par)
{
	Float_t xx = x[0];
	Float_t ss = x[1];
	Float_t tt = x[2];
	Float_t kk = x[4];

	Float_t m = par[0];						       // These simply define all of the relevant constants used in
	Float_t mx = par[1];						       // the definition of H.
	Float_t mv = par[2];						       // In particular the three parameters are pulled out of the 
	Float_t M = .938;						       // array passed to the function.
	Float_t MU2 = xx*pow(M,2) - (xx*pow((1-xx),-1)*pow(mx,2)) - pow(mv,2); //
	Float_t xp = (xx-ss)/(1-ss);					       //
	Float_t Dperp = pow(-1*pow(M,2)*pow(ss,2)*pow((1-ss),-1) - tt,0.5);    //
	Float_t a = (m + M*xx);						       //
	
	Double_t D0 = ((1-xx)*(MU2)) - pow(kk,2); // Same as above
	Double_t D1 = (1-xp)*MU2 - pow(kk,2) - pow(Dperp,2) * pow((1-xp),2);
	Double_t D2 = (1-xp)*kk*Dperp;

	Double_t num, denom;
	num = 2 * M * (1-ss) * (-2*M*(xx-xp)*pow(kk,2) - (a)*(1-xp)*D1);
	denom = pow(( pow(D1,2) + (-4*pow(D2,2))),1.5)*pow(D0,2);
	
	Double_t result;

	result = kk * (num/denom);

	return result;
}

Double_t exact(Double_t *x, Double_t *par)
{
	Float_t xx = x[0];	
	
	Double_t N = 1.468;
		
	Double_t h;
	Double_t Reg = pow(xx,-.222);

	h = 4 * N * pow(1-xx,3) * Integrate(exactint,0,10,1000,x,par,3)*Reg;

	return h;
}

Double_t h(Double_t *x, Double_t *par) // Function which defines the total GPD H, compare to eqn. 38
{

	Float_t xx = x[0];
	Float_t ss = x[1];
	Float_t tt = x[2];	
	Float_t alpha = par[3];
	Float_t alphaprime = par[4]*pow((1-xx),par[5]);
	Float_t beta = 10.0;	

	Double_t N = 1.468;
	Double_t E = e(x,par);
	Double_t Regarg = -1*(alpha + (alphaprime*tt) + (beta*tt));

	Double_t Reg = pow(xx,Regarg);

	Double_t h,a,b,c;
	
	a = N * (1-ss/2) * pow(1-xx,3) * pow(1-ss,-2);
	b = Integrate(hint,0,10,1000,x,par,3);
	c = pow(ss,2) * pow((4 * (1-ss)),-1) * E;

	h = -4*((a*b) + c)*Reg;

	//cout << "a = " << a << endl;
	//cout << "b = " << b << endl;	//These are just for testing, output the values of a,b,c to check.
	//cout << "c = " << c << endl;

	return h;
}

Double_t e(Double_t *x, Double_t *par)
{
	Float_t xx = x[0];
	Float_t ss = x[1];
	Float_t tt = x[2];
	Float_t alpha = par[3];
	Float_t alphaprime = par[4]*pow((1-xx),par[5]);
	Float_t beta = 10.0;

	Double_t Regarg = -1*(alpha + (alphaprime*tt) + (beta*tt));
	Double_t Reg = pow(xx,Regarg);
	Double_t N = 1.468;

	Double_t a,b,e;

	a = N * (1-ss/2) * pow(1-xx,3) * pow(1-ss,-1);
	b = Integrate(eint,0,10,1000,x,par,3);

	e = a*b*Reg;

	return e;
}

void printgpd1D(Double_t xi, Double_t t, Int_t n)
{
	ofstream file;
	file.open("output1D.dat");
	Double_t x[4] = {0,0,0,0};
	Double_t para[6] = {0.367,0.583,0.963,0.222,2.443,0.6649};

	Double_t xx,xl,xr;

	xx = 0.0;
	xl = 0.0;
	xr = 1.0;

	Double_t stepx = (xr-xl)/n;

	for(int i=0;i<=n;i++)
	{
		x[0] = xx;
		x[1] = xi;
		x[2] = t;
		x[3] = 0;

		file << xx << "\t" << xi << "\t" << t << "\t" << h(x,para) << endl;

		xx += stepx;
	}
	file.close();
}


void printgpd3D(Int_t n)
{
	ofstream file;
	file.open("output3D.dat");
	Double_t x[4] = {.5,0,0,0};
	Double_t para[6] = {0.367,0.583,0.963,0.222,2.443,0.6649};
	
	Double_t xl,xr,xil,xir,tl,tr,xx,xi,t;
	Double_t stepx,stepxi,stept;
	Double_t M = .938;
	Double_t result[n];
	Long_t count = 1;

	xx = 0.001;
	t = -2.0;
	xi = 0.0;
	xl = 0.0;
	xr = 1.0;
	tl = -2.0;
	tr = 0.0;
	xil = 0.0;
	xir = 1.0;

	stepx = (xr-xl)/n;
	stept = (tr-tl)/n;

	for(int i=0;i<=n;i++)
	{
		for(int i=0;i<=n;i++)
		{
			xir = (pow((pow(t,2))-(4*pow(M,2)*t),.5))/(2*pow(M,2)) + ((t) / (2*pow(M,2)));
			stepxi = (xir-xi)/n;
			while(xi <= xir && xi <= xx)
			{	
				x[0] = xx;
				x[1] = xi;
				x[2] = t;
				x[3] = 0;
				result[i] = h(x,para);
				file << xx << "\t" << xi << "\t" << t << "\t" << result[i] << endl;				
			xi += stepxi;
			cout << count << endl;
			count += 1;
			}
			xi=0;
			t += stept;
		}
		t=-2.0;
		xx += stepx;
	}
	xx=0;
	
	file.close();
}
	
