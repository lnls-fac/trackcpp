#include "auxiliary.h"


double get_magnetic_rigidity(const double energy) {
    double gamma = (energy/1e6) / (electron_rest_energy_MeV);
    double beta  = sqrt(1 - 1/(gamma*gamma));
    double b_rho = beta * energy / light_speed; // [T.m]
    return b_rho;
}
