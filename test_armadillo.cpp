#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<armadillo>

int main(int argc, char* argv[]){
	arma::mat A(5,5);
	A.ones();
	arma::mat B = A + A; 	
	std::cout << B.size() << std::endl;	
}
