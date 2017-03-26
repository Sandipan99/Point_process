#include "optimization_function.h"
#include<cmath>

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

Univar_Hawkes::Univar_Hawkes(const std::vector<double> input){
	for(int i=0;i<n;i++)
		arrival.push_back(input[i]);
	n = (int)input.size();
	t_n = input[n-1];
	
}

double Univar_Hawkes::operator()(const column_vector x)const{
	double s_1 = 0.0, s_2 = 0.0, s_t, u;
	double R[n];
	for(int i=0;i<(int)arrival.size();i++){
		s_1 += exp((-1)*x(2)*(t_n - arrival[i])) - 1;
	}
	s_1 = s_1*x(1)/x(2);
	R[0] = 0;
	for(int i=1;i<n;i++){
		R[i] = exp((-1)*x(2)*(arrival[i]-arrival[i-1]))*(1+R[i-1]);
	}
		
	for(int i=0;i<n;i++)
		s_2 += log(x(0) + x(1)*R[i]);
	return(x(0)*t_n - s_1 - s_2);
}
