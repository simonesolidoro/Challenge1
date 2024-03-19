
#include <iostream>
#include "Gradiente.cpp"
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
int main() {
    std::vector<double> xzero={0,0};
    std::vector<double> x_sol(2);
    unsigned int k=1000;
    const double tol_res=10e-6;
    const double tol_step= 10e-6;
    const double alpha_zero= 0.05;
    const double mu=0.2;
    const double teta=0.04;
    Parametri P{k,tol_res,tol_step,alpha_zero,xzero,mu,teta};
    x_sol= Gradiente(fun,dfun,P);
    std::cout<<"la soluzione: "<<x_sol[0]<<" "<<x_sol[1]<<std::endl;
    return 0;
}
// prova recupero vecchio 
