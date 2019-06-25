//////
// ROOT Macro for calculating and plotting gpd, w/ parameters
// Kyle Pressler
// SIWIF
//////
// How to use me.
// From root prompt i.e. root [0], use command,
// .L gpd.C+   (this will rebuild the library in the root interpreter)
// start()     (this calls our starting function, it can be named whatever, here it is start)
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

Double_t Integrate(Double_t (*func)(Double_t *,Double_t *), Double_t from, Double_t to, Double_t n, Double_t *x, Double_t *par); // These are function declarations, the functions themselves are defined below.
Double_t gpdint(Double_t *x, Double_t *par);											 // Function declarations tell the compiler, in this case the ROOT Cling interpreter
Double_t solved(Double_t *x, Double_t *par);											 // what kinds of functions will be called, it tells them their names and what arguments
Double_t gpd(Double_t *x, Double_t *par);											 // they accept. For example, gpdint() is passed two pointers, *x, and *par. When the function
																 // is called later, it will be given what should be passed along, namely two arrays, x, and para.
void start() // This is the main function.                                                                                 // which contain the variables and parameters.
{
	Double_t x[4] = {.5,0,0,.5}; // Array containing variables, (X,Xi,t,k)
	Double_t para[3] = {0.367,0.583,0.963}; // Array containing parameters, (m,Mx,Mlambda)

	cout << "H(" << x[0] << "," << x[1] << "," << x[2] << ") = " << gpd(x,para) << endl;    // Output gpd H(X,Xi,t)
	cout << "Integrand =  " << gpdint(x,para) << endl; // Output gpd integrand evaluated at H(X,0,0,k)
	cout << "Exact H(" << x[0] << ",0,0) = " << solved(x,para) << endl; // Output solved version H(X,0,0)

	return;
}

Double_t Integrate(Double_t (*func)(Double_t *,Double_t *), Double_t from, Double_t to, Double_t n, Double_t *x, Double_t *par) // This function will integrate any math function of four variables and n-many parameters
{																// using Simpson's method.
	Double_t result;													
	Double_t h = (to - from) / n;
 	Double_t sum1 = 0.0;
	Double_t sum2 = 0.0;
	Double_t y[4] = {x[0],x[1],x[2],x[3]};
	
	for(Int_t i = 0;i < n ;i++)
	{
		y[3] = (from + h * i + h / 2.0);
		sum1 += func(y,par);
	}
	for(Int_t i = 1; i<n;i++)
	{
		y[3] = (from + h * i + h / 2.0);
		sum2 += func(y,par);
	}
	
	y[3] = from;
	Double_t a = func(y,par);
	y[3] = to;
	Double_t b = func(y,par);

	result = h / 6.0 * (a + b + 4.0 * sum1 + 2.0 * sum2); 
	
	return result;
}

Double_t gpdint(Double_t *x, Double_t *par) // This function defines the "integrand" to be integrated in the GPD H.
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

	num = ((a*b+pow(kk,2))*D1 -2*pow(D2,2));
	denom = pow((pow(D1,2)-4*pow(D2,2)),1.5)*pow(D0,2);

	Double_t result;

	result = kk * num/denom;

	return result;
}

Double_t solved(Double_t *x, Double_t *par) // Function that defines the "Analytically Solved" form of H(X,0,0) given in line 43 of Liuti paper.
{
	Float_t xx = x[0]; 

	Float_t m = par[0];  // Parameters
	Float_t mx = par[1]; //
	Float_t mv = par[2]; //

	Float_t N = 1.468; // Constants
	Float_t M = .938; 
	Float_t MU2 = xx*pow(M,2) - (xx*pow((1-xx),-1)*pow(mx,2)) - pow(mv,2);
	Float_t a = (m + M*xx);

	Double_t result;
	
	result = N * (3.14159/12) * (2*pow(a,2)+MU2 / pow(MU2,3))*pow(1-xx,4);

	return result;
}

Double_t gpd(Double_t *x, Double_t *par) // Function which defines the total GPD H, compare to eqn. 38
{

	Float_t xx = x[0];
	Float_t ss = x[1];
		
	Double_t N = 1.468;
	Double_t E = 1;

	Double_t h,a,b,c;
	
	a = N * (1-ss/2) * pow(1-xx,3) * pow(1-ss,-2);
	b = Integrate(gpdint,0,100,1000,x,par);
	c = pow(ss,2) * pow((4 * (1-ss)),-1) * E;
	h = a * (b + c);

	//cout << "a = " << a << endl;
	//cout << "b = " << b << endl;	//These are just for testing, output the values of a,b,c to check.
	//cout << "c = " << c << endl;

	return h;
}
	
