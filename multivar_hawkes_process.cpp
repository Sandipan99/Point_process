//implements multivariate hawkes process using ogata's modified thinning algorithm
//input dimension,T,mu,alpha,beta

#include<iostream>
#include<stdio.h>
#include<vector>
#include<armadillo>
#include<math.h>
#include<fstream>
#include<time.h>
#include<map>


int present(std::map<int,std::vector<float> > m, int e){
	std::map<int,std::vector<float> >::iterator it = m.find(e);
	if(it!=m.end())
		return 1;
	else	return 0;
}

float calculate_lambda(arma::mat mu,arma::mat alpha,arma::mat beta,std::map<int,std::vector<float> > tau, float s, int dim){
	int i,j;
	std::vector<float>::iterator it;
	float lambda = 0.0;
	for(i=0;i<dim;i++){
		lambda += mu(i);
		for(j=0;j<dim;j++){
			if(present(tau,j)){
				for(it=tau[j].begin();it!=tau[j].end();it++)
					lambda+=alpha(i,j)*exp((-1)*beta(i,j)*(s-*it));
			}
		}
	}	
	return lambda;
}


void MultiVariateHawkes(arma::mat mu,arma::mat alpha,arma::mat beta, int T, int dim){
	float lambda_b,lambda_s,lambda_m;
	arma::vec n(dim,1);
	int i=0,k;
	float s=0.0,w,D,u;
	std::map<int,std::vector<float> > tau;
	std::vector<float>::iterator it;

	std::srand(time(0));

	while(s<T){
		lambda_b = calculate_lambda(mu,alpha,beta,tau,s,dim);

		u = (float)(rand())/RAND_MAX;
		w = (-1)*log(u)/lambda_b;
		s = s+w;
		D =  (float)(rand())/RAND_MAX;

		lambda_s = calculate_lambda(mu,alpha,beta,tau,s,dim);	

		if(D*lambda_b<=lambda_s){
			k=1;
			lambda_m = calculate_lambda(mu,alpha,beta,tau,s,k);
			
			while(D*lambda_b>lambda_m){
				k++;
				lambda_m = calculate_lambda(mu,alpha,beta,tau,s,k);
			}
			
			n(k-1)++;
			tau[k-1].push_back(s);
		}
	}
	std::ofstream out;
	out.open("event_sequence_hawkes_mulv");
	for(i=0;i<dim;i++){
		for(it=tau[i].begin();it!=tau[i].end();it++)
			out << i << " " << *it <<std::endl;
	}
	out.close();
}



int main(int argc, char* argv[]){

	int dim = 2;
	int T = 20,iterate = 1;
	arma::vec mu = {0.1,0.5};
	arma::mat alpha = {{0.1,0.7},{0.5,0.2}};
	arma::mat beta = {{1.2,1.0},{0.8,0.6}};
	
	while(iterate>0){
		MultiVariateHawkes(mu, alpha, beta, T, dim);
		iterate--;
	}
}
