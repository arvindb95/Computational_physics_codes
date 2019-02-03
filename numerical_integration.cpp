#include <iostream>
#include <math.h>
using namespace std;

double f_l(double u){
	double e = -1.0/3.0;
	return 3.0*pow((1 - pow(u,3.0)),e);
}

double f_r(double s){
	double g = 3.0/2.0;
	double h = -2.0/3.0;
	return g*pow((1 - pow(s,g)),h);

}

int main(){
	double a_l = 0.0;
	double b_l = pow(0.5,1.0/3.0);
	double number_of_steps = pow(10,4);
	double stepsize_l = (b_l - a_l)/number_of_steps;

	double trap_l = 0;
	double simp_l = 0;
	for (double i = a_l + stepsize_l; i < b_l; i += 2.0*stepsize_l){
		trap_l += (stepsize_l/2.0)*(f_l(i-stepsize_l) + 2.0*f_l(i) + f_l(i+stepsize_l));
		simp_l += (stepsize_l/3.0)*(f_l(i-stepsize_l) + 4.0*f_l(i) + f_l(i+stepsize_l));
	}


        double a_r = 0.0;
        double b_r = pow(0.5,2.0/3.0);
        double stepsize_r = (b_r - a_r)/number_of_steps;

        double trap_r = 0;
        double simp_r = 0;
        for (double i = a_r + stepsize_r; i < b_r; i += 2.0*stepsize_r){
		trap_r += (stepsize_r/2.0)*(f_r(i-stepsize_r) + 2.0*f_r(i) + f_r(i+stepsize_r));
		simp_r += (stepsize_r/3.0)*(f_r(i-stepsize_r) + 4.0*f_r(i) + f_r(i+stepsize_r));
	}

	double final_trap = trap_l + trap_r;
	double final_simp = simp_l + simp_r;

	cout << "Step size : (left) : " << stepsize_l << " (right) : "<< stepsize_r << endl;	
	cout << "Trapezoidal method : " << final_trap << endl;
	cout << "Simpsons method : " << final_simp << endl;
	return 0;
}
