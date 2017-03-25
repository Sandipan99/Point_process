
#ifndef __OPTIMIZATION_FUNCTION_H_INCLUDED__
#define __OPTIMIZATION_FUNCTION_H_INCLUDED__

#include<vector>
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

};

#endif
