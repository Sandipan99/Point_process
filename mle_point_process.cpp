
#include<iostream>
#include<vector>
#include<cmath>
#include<dlib/optimization.h>
#include "optimization_function.h"
#include<fstream>
#include<string>

/**

@author Sandipan Sikdar

**/


bool present(std::vector<int> m,int i){
	std::vector<int>::iterator it;
	it = std::find(m.begin(),m.end(),i);
	if(it!=m.end())
		return true;
	else	return false;
}

void PoissonProcess(std::vector<double> event){
	column_vector lambda(1);
	lambda = 0.5;	
	Poisson_process f(event);
	find_min_using_approximate_derivatives(bfgs_search_strategy(),objective_delta_stop_strategy(),f,lambda,-1);
	std::cout<< "solution: " << lambda << std::endl;
}


void UnivarHawkesProcess(std::vector<double> event){
	column_vector starting_point(2);
	double b = 0.7;
	while(b<5){
		starting_point = 0.8,0.3;
		Univar_Hawkes f(event,b);
		find_min_using_approximate_derivatives(bfgs_search_strategy(),objective_delta_stop_strategy(),f,starting_point,-1);
		std::cout<< "solution: " << starting_point << std::endl;
		b+=0.25;
		break;
	}	
}

int main(int argc, char* argv[]){
	int a,b,t,n,curr,cnt;
	std::map<int,std::vector<int> > node_arr;
	std::vector<double> event;
	std::vector<int> node;
	std::ifstream in;
	std::string str;
	in.open("event_sequence_hawkes_univ");
	while(in >> t >> a >> b){
		if(!present(node,a))
			node.push_back(a);
		if(!present(node,b))
			node.push_back(b);
		node_arr[a].push_back(t);
		node_arr[b].push_back(t);
	}
	in.close();
	
	std::vector<int>::iterator it;
	for(it=node.begin();it!=node.end();it++){
		n=*it;
		curr = node_arr[n][0];
		cnt = 1;
		for(int i=1;i<(int)node_arr[n].size();i++){
			if(curr!=node_arr[n][i]){
				curr = node_arr[n][i];
				event.push_back(node_arr[n][i]);
			}
		}
		
		UnivarHawkesProcess(event);
		event.clear();
		break;
	}
	//std::cout<<event.size()<<std::endl;
	//PoissonProcess(event);
	//UnivarHawkesProcess(event);
}
