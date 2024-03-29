
#include <iostream>
#include "GetPot"
#include "Gradiente.hpp"
#include "Parametri.h"
#include <cmath>



// f
double fun(const std::vector<double> &x) 
{
    return x[0]*x[1]+4*x[0]*x[0]*x[0]*x[0]+x[1]*x[1]+3*x[0];     // f(X1,X2)=f(x[0],x[1])
};

// gradiente di f
std::vector<double> dfun(const std::vector<double> &x) 
{
    std::vector<double> V(2);
    V[0]= x[1]+16*x[0]*x[0]*x[0]+3;   //derivata di fun rispetto X1
    V[1]= x[0]+2*x[1];      //derivata di fun rispetto X2
    return V;
};
int main(int argc, char **argv) {
    
    GetPot command_line(argc, argv);	
    int T=command_line("scelta",3);
    int d=command_line("gradiente",1);

    GetPot datafile("data");
    double x_zero= datafile("data/x_iniziale/x_zero",0.0);
    double x_uno= datafile("data/x_iniziale/x_uno",0.0);
    unsigned int k= datafile("data/Parametri/k",1000);
    const double tol_res= datafile("data/Parametri/tol_res",10e-6);
    const double tol_step= datafile("data/Parametri/tol_step",10e-6);
    const double alpha_zero= datafile("data/Parametri/alpha_zero",0.05);
    const double mu= datafile("data/Parametri/mu",0.3);
    const double teta= datafile("data/Parametri/teta",0.04);

    std::vector<double> xzero= {x_zero,x_uno};
    std::vector<double> x_sol(2);
    Parametri P{k,tol_res,tol_step,alpha_zero,xzero,mu,teta};
    

    if (T==1){
	    std::cout<<"scelta Exponential decay"<<std::endl;
	    if(d==1){
		    std::cout<<"usato gradiente esatto"<<std::endl;
		    Gradiente<expDec>(fun,dfun,P,x_sol);}
	    else{
		    std::cout<<"usata differenze finite centrate"<<std::endl;
	            Gradiente<expDec>(fun,derivfun,P,x_sol);}
    }

    if (T==2){
	    std::cout<<"scelta Inverse decay"<<std::endl;
	    if(d==1){
		    std::cout<<"usato gradiente esatto"<<std::endl;
                    Gradiente<invDec>(fun,dfun,P,x_sol);}
            else{
                    std::cout<<"usata differenze finite centrate"<<std::endl;
                    Gradiente<invDec>(fun,derivfun,P,x_sol);}
    }

    if (T==3){
	    std::cout<<"scelta Armijo rule"<<std::endl;
            if(d==1){
		    std::cout<<"usato gradiente esatto"<<std::endl;
		    Gradiente<A>(fun,dfun,P,x_sol);}
	    else{
		    std::cout<<"usata differenze finite centrate"<<std::endl;
		    Gradiente<A>(fun,derivfun,P,x_sol);}
    }

    std::cout<<"la soluzione: x_1="<<x_sol[0]<<" x_2="<<x_sol[1]<<std::endl;
    
    return 0;
}

