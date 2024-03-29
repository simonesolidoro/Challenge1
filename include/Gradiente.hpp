#include "functional"
#include "Parametri.h"
#include <vector>
#include <math.h>

//dichiarazioni
enum scelta {expDec, invDec, A};
template <scelta T>
void Gradiente(const std::function<double(const std::vector<double> &)> &, const std::function<std::vector<double>(const std::vector<double> &)> &, const Parametri &, std::vector<double>&);

double norma(const std::vector<double> &x);   // norma di vettore

std::vector<double> sotVet(const std::vector<double> &x, const std::vector<double> &y); // sottrazione membro a membro tra vettori

std::vector<double> Prod_sca_vett(const double c, const std::vector<double> &x); // prodotto tra scalare e vettore

double norma_di_diff(const std::vector<double> &,const std::vector<double> &); // unione di sotVet e norma per fare unico ciclo for

double fun(const std::vector<double> &x); // fun definita in main ridichiarata per usarla in differenze finite
std::vector<double> derivfun(const std::vector<double> &x); // calcola gradiente con metodo differenze finite centrate ed h=0.01

// definizioni
template <scelta T>
void Gradiente(const std::function<double(const std::vector<double> &)> &fun,
               const std::function<std::vector<double>(const std::vector<double> &)> &dfun,
               const Parametri & P,
               std::vector<double> &x_sol){
    std::vector<double> x_old(P.x_zero);  
    std::vector<double> x_new(P.x_zero);  
    double alpha=P.alpha_zero;
    double alpha_zero=P.alpha_zero; //copia di alphazero per non dover ogni volta richiamare P.alpha_zero
    unsigned int k=0;   // contatore iterazioni
    bool no_conv= 1;    // inizializzato vero per garantire prima iterazione
    while (no_conv) {
        //scelta di alpha, (if constexpr valutato at compile time quindi anche se in while rimane efficiente)
        //expoetial decay 
        if constexpr (T == expDec) {
            alpha = alpha_zero * std::exp(-P.mu * k);
        }
        //inverse decay
        if constexpr (T == invDec ) {
            alpha = alpha_zero / (1 + P.mu * k);
        }
        // Armijo

        if constexpr (T == A) {
            alpha = alpha_zero; // ad ogni step ricerca di nuova alpha riparte da alpha_zero
            // Armijo rule, check if the sufficient decrease condition is satisfied
            while ((fun(x_old) - fun(sotVet(x_old, Prod_sca_vett(alpha, dfun(x_old))))) < P.teta * alpha * pow(norma(dfun(x_old)), 2)) {
                alpha *= 0.5;
            }
        }


        //calcolo xk+1
        for (size_t i = 0; i < P.x_zero.size(); i++) {
            x_new[i] = x_old[i] - alpha * dfun(x_old)[i];
        }

        //incrementa numero iterazione per la successiva
        k += 1;

        // check convergenza
        no_conv = k < P.max_it && std::abs(fun(x_new) - fun(x_old)) > P.tol_res && norma_di_diff(x_new, x_old) > P.tol_step;
                   // se no_conv=vero-> nessun criterio di convergenza soddisfatto e quindi rientra nel ciclo while
                   // se almeno un criterio di convergenza è soddisfatto -> no_conv=falso ed esce dal ciclo while

        // aggiorna xk= xk+1
        for (size_t i=0; i<P.x_zero.size();i++){
            x_old[i] = x_new[i];
        }

    }
    // rturn:copia ultima x_new(soluzione) in reference x_sol
    for(size_t i=0; i<x_new.size();i++)
        x_sol[i]=x_new[i];
};

//calcolo gradiente con diff finite
std::vector<double> derivfun(const std::vector<double> &x){
	std::vector<double> V(2);
	const double h=0.01;
	V[0]= ( fun({x[0]+h,x[1]}) - fun({x[0]-h,x[1]}) )/(2*h);   //derivata di fun rispetto X1
	V[1]= ( fun({x[0],x[1]+h}) - fun({x[0],x[1]-h}) )/(2*h);      //derivata di fun rispetto X2
	return V;
};

// funzioni ausiliarie

double norma(const std::vector<double> &x){
    double sum=0.0;
    for (auto i : x){
        sum+=(i*i);
    }
    return sqrt(sum);
}

double norma_di_diff(const std::vector<double> &x, const std::vector<double> &y){
    double sum=0.0;
    double diff=0.0; // variabile ausiliaria per differenza tra x[i] e y[i]
    for (size_t i=0; i<x.size(); i++){
        diff=x[i]-y[i];
        sum+=(diff*diff);
    }
    return sqrt(sum);
}

std::vector<double> sotVet(const std::vector<double> &x, const std::vector<double> &y){
    std::vector<double> V(x.size(),0.0);
    for (size_t i=0; i<x.size(); i++){
        V[i]=x[i]-y[i];
    }
    return V;
}

std::vector<double> Prod_sca_vett(const double c, const std::vector<double> &x){
    std::vector<double> V(x.size(),0.0);
    for (size_t i=0; i<x.size(); i++){
        V[i]=x[i]*c;
    }
    return V;
}

