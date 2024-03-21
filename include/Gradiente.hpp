#include "functional"
#include "Parametri.h"
#include <vector>
#include <math.h>

//dichiarazioni
enum scelta {expDec, invDec, A};
template <scelta T>
void Gradiente(const std::function<double(const std::vector<double> &)> &, const std::function<std::vector<double>(const std::vector<double> &)> &, const Parametri, std::vector<double>&);

double norma(const std::vector<double> &x);   // norma di vettore

std::vector<double> sotVet(const std::vector<double> &x, const std::vector<double> &y); // sottrazione membro a membro tra vettori

std::vector<double> Prod_sca_vett(const double c, const std::vector<double> &x); // prodotto tra scalare e vettore

double norma_di_diff(const std::vector<double> &,const std::vector<double> &); // unione di sotVet e norma per fare unico ciclo for


// definizioni
template <scelta T>
void Gradiente(const std::function<double(const std::vector<double> &)> &fun,
               const std::function<std::vector<double>(const std::vector<double> &)> &dfun,
               const Parametri P,
               std::vector<double> &x_sol){
    std::vector<double> x_old(P.x_zero);  
    std::vector<double> x_new(P.x_zero);  
    double alpha=P.alpha_zero;
    unsigned int k=0;   // contatore iterazioni
    bool no_conv= 1;    // inizializzato vero per garantire prima iterazione
    while (no_conv) {
        //scelta di alpha
        //expoetial decay 
        if constexpr (T == expDec) {
            alpha = P.alpha_zero * std::exp(-P.mu * k);
        }
        //inverse decay
        if constexpr (T == invDec ) {
            alpha = P.alpha_zero / (1 + P.mu * k);
        }
        // Armijo

        if constexpr (T == A) {
            alpha = P.alpha_zero; // ad ogni step ricerca di nuova alpha riparte da alpha_zero
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
}

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

