#include <iostream>
#include <math.h> // math library
using namespace std;

double square(double x){
	return pow(x,2.0);

}

double firstDirivative_O2(double x0, double stepsize){
	return (square(x0 + stepsize) - square(x0 - stepsize))/(2*stepsize) ;
}


double firstDirivative_O4(double x0, double stepsize){
	        return (square(x0 - 2*stepsize) - 8*square(x0 - stepsize) + 8*square(x0 + stepsize) - square(x0 + 2*stepsize))/(12*stepsize) ;
}

int main(){
	cout << firstDirivative_O2(1,0.1) << endl;
	return 0;
}


