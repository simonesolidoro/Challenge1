
#include "functional"
#include "Parametri.h"
#include <vector>
std::vector<double> Gradiente(const std::function<double(const std::vector<double> &)> &fun,
		              const std::function<std::vector<double>(const std::vector<double> &)> &dfun,
			      const Parametri P);

double norma(const std::vector<double> &x);   // norma di vettore

std::vector<double> sotVet(const std::vector<double> &x, const std::vector<double> &y); // sottrazione membro a membro tra vettori

std::vector<double> Prod_sca_vett(const double c, const std::vector<double> &x); // prodotto tra scalare e vettore

double norma_di_diff(const std::vector<double> &,const std::vector<double> &); // unione di sotVet e norma per fare unico ciclo for

// per permettere scelta metodo per alpha da finire
//template <typename T>
//auto new_alpha(unsigned int k, const double mu, const double alpha_zero, T strategia_alpha, const double d);


