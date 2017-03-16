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
	int i = 0,j=0;
	int t = atoi(argv[1]);
	float u,v,i_a_t;
	float lambda = atof(argv[2]);
	int sim = 1,n=0;
	//int events[sim];
	//for(j=0;j<sim;j++)
	//	events[j] = 0;
	//std::cout << t << "," << lambda << std::endl;
	std::srand(time(0));
	while(i<sim){
		i_a_t = 0.0;
		while(i_a_t<t){
			u = (float)(rand())/RAND_MAX;
			v = (-1)*(1/lambda)*log(1-u);
			i_a_t+=v;
			n++;
			e.push_back(i_a_t);
			//events[i]+=1;
			if(n>50)
				break;
		}
		i++;
	}
	
	//for(i=0;i<sim;i++)
	//	std::cout<<events[i]<<std::endl;
	std::ofstream out;
	out.open("event_sequence_poisson");
	std::vector<float>::iterator it;
	for(it=e.begin();it!=e.end();it++)
		out<<*it<<std::endl;
	out.close();
}
