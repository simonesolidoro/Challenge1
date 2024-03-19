

#ifndef PACS_CHALLENGE1_PARAMETRI_H
#define PACS_CHALLENGE1_PARAMETRI_H
#include <vector>

struct Parametri{
    unsigned int max_it;
    const double tol_res;
    const double tol_step;
    const double alpha_zero;
    const std::vector<double> x_zero;
    const double mu;
    const double teta; // teta tra (0, 0.05) di condizione sufficiente di Armijo rule
};

#endif //PACS_CHALLENGE1_PARAMETRI_H
