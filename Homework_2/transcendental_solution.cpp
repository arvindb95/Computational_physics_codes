#include <iostream>
#include <math.h>
using namespace std;

double lhs(double delta, double eta){
	return sqrt(pow((delta/eta),2.0)-1.0);
}
double f1(double delta, double eta){
	return lhs(delta, eta) - tan(eta);
}	
double f2(double delta, double eta){
	return lhs(delta, eta) + 1.0/tan(eta);
}
double df1(double delta, double x0, double stepsize){
	return (f1(delta, x0 - 2*stepsize) - 8*f1(delta, x0 - stepsize) + 8*f1(delta, x0 + stepsize) - f1(delta, x0 + 2*stepsize))/(12*stepsize) ;
}
double df2(double delta, double x0, double stepsize){
	        return (f2(delta, x0 - 2*stepsize) - 8*f2(delta, x0 - stepsize) + 8*f2(delta, x0 + stepsize) - f2(delta, x0 + 2*stepsize))/(12*stepsize) ;
}
double NewtonRaphsonf1(double delta, double x0, double stepsize){
	double der = df1(delta, x0, stepsize);
	double x1 = x0 - (f1(delta, x0)/der);
	while ((x0 - x1) >= pow(stepsize, 8.0)){
		x0 = x1;
		der = df1(delta, x0, stepsize);
		x1 = x0 - (f1(delta, x0)/der);
	}
	return x1;	
}
double NewtonRaphsonf2(double delta, double x0, double stepsize){
	double der = df2(delta, x0, stepsize);
	double x1 = x0 - (f2(delta, x0)/der);
	while ((x0 - x1) >= pow(stepsize, 8.0)){
		x0 = x1;
		der = df2(delta, x0, stepsize);
		x1 = x0 - (f2(delta, x0)/der);
		}       
	return x1;
}
int main(){
	const double PI = 3.141592653589793238463;

	double a = 3.0;               // the half width of the well in angstroms 
	double m = 1.0;               // mass of electron in 1 me units
	double V0 = 10.0;             // height/depth of well in eV
	double hbar = 1.0;            // hbar
	double delta = sqrt(2.0*m*V0*pow(a,2.0)/pow(hbar,2.0));
	
	double even_guess = PI/4.0 ;   // Guess for even solution
	double odd_guess = 0.999*PI ;  // Guess for odd solution
	double even_sol = NewtonRaphsonf1(delta, even_guess,pow(10,-3.0)) ; // Obtaining the solution using the Newton Raphson method
	double odd_sol = NewtonRaphsonf2(delta, odd_guess,pow(10,-3.0)) ;  // Obtaining the solution using the Newton Raphson method
	
	cout << "eta for even sol : " << even_sol << endl;
	cout << "eta for odd sol : " << odd_sol << endl;

	return 0;

}
