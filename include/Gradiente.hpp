#include "functional"
#include "Parametri.h"
#include <vector>
#include <math.h>

//dichiarazioni
//@note Try to use namespaces to avoid name clashes
enum scelta {expDec, invDec, A};
//@note you use a vector of vectors for a matricies. It is fine, but rememeber that while an array of arrays 
// stores elements contigous in memory, a vector of vector does not since tthe vectors for each row are stored in different memory locations.
// This can lead to cache misses and slower performance. If you are working with large matricies, consider using a single vector and indexing it manually.
// Also, consider using a library like Eigen for matrix operations.

//@note comment the functions declarations for readeability.
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
            //@note always sset a maximum number of iterations to avoid infinite loops
            // Armijo should converge if f differentiable with Lipschitz continuous gradient
            // but we cannot guarantee that roundoff errors do not cause the loop to run forever
            while ((fun(x_old) - fun(sotVet(x_old, Prod_sca_vett(alpha, dfun(x_old))))) < P.teta * alpha * pow(norma(dfun(x_old)), 2)) {
                alpha *= 0.5;
            }
        }


        //calcolo xk+1
        for (size_t i = 0; i < P.x_zero.size(); i++) {
            x_new[i] = x_old[i] - alpha * dfun(x_old)[i];
        }

        //incrementa numero iterazione per la successiva
        k += 1; //@note use ++k instead of k+=1 for better performance

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
    //@note why not simply x_sol=x_new?
    //In fact here you can also do a move assignment since x_new is not used anymore
    // x_sol=std::move(x_new);
    // and save a useless copy.
    for(size_t i=0; i<x_new.size();i++)
        x_sol[i]=x_new[i];
};

//calcolo gradiente con diff finite
std::vector<double> derivfun(const std::vector<double> &x){
	std::vector<double> V(2);
	const double h=0.01;
    //@note Why assume a bivariate function? You can easily generalize to a multivariate function since you are using a vector.
	V[0]= ( fun({x[0]+h,x[1]}) - fun({x[0]-h,x[1]}) )/(2*h);   //derivata di fun rispetto X1
	V[1]= ( fun({x[0],x[1]+h}) - fun({x[0],x[1]-h}) )/(2*h);      //derivata di fun rispetto X2
	return V;
};

// funzioni ausiliarie

//@note I would have put these utilities in a separate file since they are of general use.
//Note: if you define a (non template) function in a header file, you should declare it inline to avoid multiple definition errors.
//This is because the header file may be included in multiple translation units, and the function is defined in each of them.
//inline double norma(const std::vector<double> &x){

double norma(const std::vector<double> &x){
    double sum=0.0;
    for (auto i : x){
        sum+=(i*i);//@note good use of range based for loop. The paranthesis is not needed here, but it is fine.
    }
    return sqrt(sum);
}

double norma_di_diff(const std::vector<double> &x, const std::vector<double> &y){
    double sum=0.0;
    double diff=0.0; // variabile ausiliaria per differenza tra x[i] e y[i]
    for (size_t i=0; i<x.size(); i++){
        diff=x[i]-y[i];
        sum+=(diff*diff);//@note I see you like putting paranthesis around the expression. It is not needed, but it is fine.
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

//@note Here and before you coud have overloaded the operators insttead of creating free functions
std::vector<double> Prod_sca_vett(const double c, const std::vector<double> &x){
    std::vector<double> V(x.size(),0.0);
    for (size_t i=0; i<x.size(); i++){
        V[i]=x[i]*c;
    }
    return V;
}

