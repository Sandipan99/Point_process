//implements univariate hawkes process follwing ogata's algorithm

//input: intensity distribution parameter (u_p), end time (t)
//output: sequence of events 

#include<stdlib.h>
#include<time.h>
#include<iostream>
#include<vector>
#include<math.h>
#include<fstream>
#include<string>


void UniVariateHawkes(float mu, float alpha, float beta, int T, int type){
	std::vector<float> tau;
	std::vector<float> lambda_t;
	std::vector<float> t_s;
	std::vector<float>::iterator it;
	float lambda_s,lambda_b,u,w,D,s=0.0;
	int i,n=0;
	std::srand(time(NULL));
	lambda_b = mu;
	while(true){
		lambda_b = mu;
		if(tau.size()>0){
			for(it=tau.begin();it!=tau.end();it++)	
				lambda_b+= alpha*exp((-1)*beta*(s-*it));
		}
		u = (float)(rand())/RAND_MAX;
		w = (-1)*log(u)/lambda_b;
		s = s+w;
		D =  (float)(rand())/RAND_MAX;
		lambda_s = mu;
		for(it=tau.begin();it!=tau.end();it++)	
			lambda_s+= alpha*exp((-1)*beta*(s-*it));

		if(D*lambda_b<=lambda_s){
			n+=1;
			if(type==1){			
				if(s<T){
					lambda_t.push_back(lambda_s);
					tau.push_back(s);
				}
				else	break;
			}
			else{
				if(n<=T){
					lambda_t.push_back(lambda_s);
					tau.push_back(s);
				}
				else	break;
				
			}
		}
	}

	std::ofstream out;
	out.open("event_sequence_hawkes_univ");
	for(it=tau.begin();it!=tau.end();it++){
		out<<*it<<std::endl;
	}
	out.close();
	
}

int main(int argc, char* argv[]){
	
	//int T = atoi(argv[1]);
	int events = atoi(argv[1]); //number of events to generate or end time
	float mu = atof(argv[2]);
	float alpha = atof(argv[3]);
	float beta = atof(argv[4]);
	int type = atoi(argv[5]); // determines the type of the simulation 1. simulation with specified end time 2. simulation with specified number of events

	UniVariateHawkes(mu, alpha, beta, events,type);
		
}
