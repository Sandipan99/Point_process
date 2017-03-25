#include "optimization_function.h"

Poisson_process::Poisson_process(const std::vector<double> input){
	i_a_t.push_back(input[0]);
	sum+=input[0];
	for(int i=1;i<(int)input.size();i++){
		i_a_t.push_back(input[i]-input[i-1]);
		sum+=input[i]-input[i-1];
	}
	t_n = input[(int)input.size()-1];
	n = (int)input.size(); 
}

double Poisson_process::operator()(const column_vector x)const{
	return(x(0)*t_n - n*log(x(0)) - sum);	
}
