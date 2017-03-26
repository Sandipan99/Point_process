#include<iostream>
#include<vector>
#include<cmath>
#include<dlib/optimization.h>
#include "optimization_function.h"
#include<fstream>


void PoissonProcess(std::vector<double> event){
	column_vector lambda(1);
	lambda = 0.5;	
	Poisson_process f(event);
	find_min_using_approximate_derivatives(bfgs_search_strategy(),objective_delta_stop_strategy(),f,lambda,-1);
	std::cout<< "solution: " << lambda << std::endl;
}


void UnivarHawkesProcess(std::vector<double> event){
	column_vector starting_point(3);
	starting_point = 0.5,0.7,1.0;
	Univar_Hawkes f(event);
	find_min_using_approximate_derivatives(bfgs_search_strategy(),objective_delta_stop_strategy(),f,starting_point,-1);
	std::cout<< "solution: " << starting_point << std::endl;	
}

int main(int argc, char* argv[]){
	double a;
	std::vector<double> event;
	std::ifstream in;
	in.open("event_sequence_poisson");
	while(in >> a){
		event.push_back(a);
	}
	in.close();
	//PoissonProcess(event);
	UnivarHawkesProcess(event);
}
