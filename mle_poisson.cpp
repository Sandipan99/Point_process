#include<iostream>
#include<vector>
#include<cmath>
#include<dlib/optimization.h>
#include<fstream>

using namespace dlib;

typedef matrix<double,0,1> column_vector;

class function{	
	private:
		std::vector<double> i_a_t;
		double t_n;
		double sum = 0.0;
		double n;
	public:
		function(const std::vector<double> input){
			i_a_t.push_back(input[0]);
			sum+=input[0];
			for(int i=1;i<(int)input.size();i++){
				i_a_t.push_back(input[i]-input[i-1]);
				sum+=input[i]-input[i-1];
			}
			t_n = input[(int)input.size()-1];
			n = (int)input.size(); 
		}
		double operator()(const column_vector x) const{
			return(x(0)*t_n - n*log(x(0)) - sum);
		}
};


int main(int argc, char* argv[]){
	double a;
	std::vector<double> event;
	column_vector lambda(1);
	lambda = 0.5;
	std::ifstream in;
	in.open("event_sequence_poisson");
	while(in >> a){
		event.push_back(a);
	}
	in.close();
	function f(event);
	find_min_using_approximate_derivatives(bfgs_search_strategy(),objective_delta_stop_strategy(),f,lambda,-1);
	std::cout<< "solution: " << lambda << std::endl;

}
