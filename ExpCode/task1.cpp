#include <cmath> 
#include <iostream>

double myexp(double,double);
int main(){
	double x = -10;
	double tol = 0.00000001;
	double calc = myexp(x,tol);
	std::cout.precision(15);
	std::cout << "Approximated value " << calc <<"\n";
	std::cout << "Reference evalue " << exp(x) <<"\n";
	double error = std::abs((exp(x)-calc));
	std::cout <<"error : " << error <<"\n";
	std::cout <<"toler : "<< tol   <<"\n";
	
	if(error < tol){
			std::cout <<"error less than tol"<<"\n";
	}
	if(error >= tol){
			std::cout <<"error larger tol"<<"\n";
	}
}
//Solves taylor series  
//1 + x + x*x/2 + x*x/2*x/3, no use of power of faculty
double myexp(double x, double tol){
	double next_Term = x;
	double output = 1;
	double last = 0;
	bool growing = true;
	//prevF = x; %Error is given for come c between 0 and x. potential error largest when c is x
	//Stop condition when update is smaller what double can store.
	for(int i = 2; ( std::fabs((output-last)) > tol) || growing  ;i++){
		last = output;
		output += next_Term;

        double addFactor = x/i;
		next_Term = next_Term * addFactor;
		//std::cout << next_Term<<"\n";

        if (std::fabs(addFactor) < 1){
            growing = false;
        }
	}
	output += next_Term;
	return output;
}
