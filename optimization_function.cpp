#include "optimization_function.h"
#include<cmath>
#include<iostream>

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
	n = (int)input.size();
	for(int i=0;i<n;i++)
		arrival.push_back(input[i]);	
	t_n = input[n-1];
	
}

double Univar_Hawkes::operator()(const column_vector x)const{
	double s_1 = 0.0, s_2 = 0.0, s_t, u;
	double R[n];
	for(int i=0;i<(int)arrival.size();i++){
		s_1 += exp((-1)*0.7*(t_n - arrival[i])) - 1;
	}
	s_1 = s_1*x(1)/0.7;
	
	R[0] = 0;
	for(int i=1;i<n;i++){
		R[i] = exp((-1)*0.7*(arrival[i]-arrival[i-1]))*(1+R[i-1]);
	}
		
	for(int i=0;i<n;i++)
		s_2 += log(x(0) + x(1)*R[i]);
	//std::cout << x(0)*t_n + s_1 - s_2 << std::endl;
	return(x(0)*t_n - s_1 - s_2);	
}

/*
Multivar_Hawkes::Multivar_Hawkes(std::vector<std::tuple<int,double> >input, int d){
	int id;
	double val;	
	n = (int)input.size();
	dim = d;
	for(int i=0;i<n;i++){
		id = std::get<0>(input[i]);
		val = std::get<1>(input[i]);
		arrival.push_back(std::make_tuple(id,val));
	}
		
}

double Multivar_Hawkes::operator()(column_vector x)const{
	double fin_res = 0.0;
	int k;	
	for(int i=0;i<dim;i++){    // calculating likelihood for each dimension. final likelihood will be sum of all
		double t_1 = 0.0,s_1 = 0.0,s_2 = 0.0,s_t,s_3=0.0,alpha,beta,t,t_a,t_b;
		int m;
		t_1 = x(i)*t_n;    // first term in the expression....
		for(int j=0;j<dim;j++){
			alpha = x(dim+dim*i+j);
			beta = x(3*dim+dim*i+j);
			s_3 = alpha/beta;
			s_t = 0.0;
			for(k=0;k<n;k++){
				m = std::get<0>(arrival[k]);
				if(m==j){
					t = std::get<1>(arrival[k]);
					s_t+=1-exp((-1)*beta*(t_n-t));
				}
			}
			s_3 = s_3*s_t;
			s_1+=s_3; // second term in the expression....
		}
		std::vector<double> arrival_p;
		for(k=0;k<n;k++){
			m = std::get<0>(arrival[k]);
			if(m==i)
				arrival_p.push_back(std::get<1>(arrival[k]));
		}
		k = (int)arrival_p.size(); 
		// calculating R 
		double R[dim][k];
		double s_12,s_21;
		for(int a=0;a<dim;a++){
			R[a][0] = 0.0;
			beta = x(3*dim+a);
			for(int b=1;b<k;b++){
				t_a = arrival_p[b-1];
				t_b = arrival_p[b];
				s_12 = 0.0;
				for(int l=0;l<n;l++){	
					t = std::get<1>(arrival[l]);
					if((t>t_a)and(t<t_b))
				 		s_12+= exp((-1)*beta*(arrival_p[b]-t));
				}
				R[a][b] = s_12 + exp((-1)*beta*(arrival_p[b]-arrival_p[b-1]));
			}
		}
		for(int a=0;a<k;a++){
			//s_2 = log(x(i);
			s_21 = 0.0;
			for(b=0;b<dim;b++){
				alpha = x(dim+b);
				s_21 += alpha*R[b][a];
			}
			s_2 += log(log(x(i)+s_21)); // third term in the expression	
		}
		fin_res+=t_1+s_1-s_2;			
	}
	return fin_res;	
}
*/
