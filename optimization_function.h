
#ifndef __OPTIMIZATION_FUNCTION_H_INCLUDED__
#define __OPTIMIZATION_FUNCTION_H_INCLUDED__

#include<vector>
#include<tuple>
#include<map>
#include<dlib/optimization.h>

using namespace dlib;
typedef matrix<double,0,1> column_vector;

class Poisson_process{
	private:
		std::vector<double> i_a_t;
		double t_n;
		double sum = 0.0;
		double n;
	public:
		Poisson_process(std::vector<double>);
		double operator()(const column_vector) const;
};


class Univar_Hawkes{
	private:
		std::vector<double> arrival;
		double t_n;
		int n;

	public:
		Univar_Hawkes(std::vector<double>);
		double operator()(const column_vector) const;
		 
};

/*
class Multivar_Hawkes{
	private:
		std::map<int,double> arrival;
		double t_n;
		int n;
		int dim;
		column_vector beta_m(1);
	public:
		Multivar_Hawkes(std::vector<std::tuple<int,double> >, int);
		double operator()(const column_vector) const; //parametrs to be placed in a single vector with proper ordering maintained 
	
};
*/
#endif
