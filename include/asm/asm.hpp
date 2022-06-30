#ifndef ASM_HPP
#define ASM_HPP

#include <vector>
#include "psr.hpp"

std::vector<double> make_g_vector(std::vector<Correlation>& corrs, double stdev, double correntropy_factor, bool print);

std::vector<double> make_G_matrix(std::vector<Correlation>& corrs, double stdev, double correntropy_factor, bool print);

std::vector<double> solve_system(std::vector<double>& g, std::vector<double>& G, bool print);

#endif