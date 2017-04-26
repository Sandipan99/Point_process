//simulates poisson process.....
//input: end time (t), intensity (lambda)
//output: sequence of times of the occurrence of events

#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<vector>
#include<fstream>
#include<time.h>

int main(int argc, char* argv[]){
	std::vector<float> e;
	int i=0,j=0;
	int t = atoi(argv[1]);
	float u,v,i_a_t;
	float lambda = atof(argv[2]);
	int type = atoi(argv[3]);
	int n=0;
	std::srand(time(0));
	i_a_t = 0.0;
	while(true){
		u = (float)(rand())/RAND_MAX;
		v = (-1)*(1/lambda)*log(1-u);
		i_a_t+=v;	
		if(type==1){   // generate specific number of events 
			if(n<t){
				e.push_back(i_a_t);
				n++;
			}
			else	break;
		}
		else{	      // generate upto a specific time
			if(i_a_t<t)
				e.push_back(i_a_t);
			else	break;
		}
	}
	i++;
	
	std::ofstream out;
	out.open("event_sequence_poisson");
	std::vector<float>::iterator it;
	for(it=e.begin();it!=e.end();it++)
		out<<*it<<std::endl;
	out.close();
}
