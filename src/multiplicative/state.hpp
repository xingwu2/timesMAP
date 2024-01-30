#ifndef TIMESMAP_MULTIPLICATIVE_STATE_H_
#define TIMESMAP_MULTIPLICATIVE_STATE_H_

/*
    timesMAP - multiplicative polygenic model
    Copyright (C) 2023-2024 Lucas Czech and Xing Wu

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lczech@carnegiescience.edu>
    Xing Wu <wxu@carnegiescience.edu>
    Department of Plant Biology, Carnegie Institution For Science
    260 Panama Street, Stanford, CA 94305, USA
*/

#include <string>
#include <vector>

// =================================================================================================
//     Multiplicative State
// =================================================================================================

/**
 * @brief Collection of all posteriors of the multiplicative model.
 */
struct MultiplicativePosteriors
{
    std::vector<double> alpha;
    std::vector<double> beta;
    std::vector<int> gamma;
    double pi      = 0.0;
    double sigma_1 = 0.0;
    double sigma_e = 0.0;
};

/**
 * @brief Internally cached values of the multiplicative model.
 */
struct MultiplicativeCache
{
    // Intermediate values kept between the sampling update functions to avoid recomputation
    std::vector<double> c_col_sum_squared;
    std::vector<double> circle_product_x_beta;
    std::vector<double> product_c_alpha;
};

/**
 * @brief Convergence statistics computed by the multiplicative model.
 */
struct MultiplicativeStatistics
{
    size_t polygenicity            = 0;
    double genetic_var             = 0.0;
    double pheno_var               = 0.0;
    double large_beta_ratio        = 0.0;
    double large_beta_heritability = 0.0;
    double total_heritability      = 0.0;
};

#endif // include guard
