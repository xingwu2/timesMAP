#ifndef TIMESMAP_MULTIPLICATIVE_HYPER_PARAMS_H_
#define TIMESMAP_MULTIPLICATIVE_HYPER_PARAMS_H_

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

#include "CLI/CLI.hpp"

#include "common/input_data.hpp"

#include <string>
#include <vector>

// =================================================================================================
//     Hyperparameters
// =================================================================================================

/**
 * @brief Collection of all hyper params for the multiplicative model.
 */
struct MultiplicativeHyperparams
{
    // Distribution hyperparameters.
    // pi_b is initialized to ratio num_snps / pi_b_ratio
    double pi_a = 1.0;
    double pi_b = 1.0;
    double pi_b_ratio = 10.0;
    double sigma_1_a = 1.0;
    double sigma_1_b = 1.0;
    double sigma_e_a = 1.0;
    double sigma_e_b = 1.0;

    // Sample validity
    double min_beta = 0.05;
    double large_beta = 0.3;

    // Convergence Criteria
    size_t num_top_betas = 5;
    double geweke_first = 0.1;
    double geweke_last = 0.5;
    size_t geweke_intervals = 20;
    double convergence_max_zscore = 1.5;
};

// -------------------------------------------------------------------------
//     add_multiplicative_hyper_params_cli
// -------------------------------------------------------------------------

/**
 * @brief Add all command line options for the multiplicative hyper params to a given command.
 */
inline void add_multiplicative_hyper_params_cli( CLI::App& app, MultiplicativeHyperparams& hyper )
{
    // ---------------------------------------
    //     Distribution hyperparameters
    // ---------------------------------------

    app.add_option(
        "--pi-a",
        hyper.pi_a,
        "Initial value for the `a` parameter of the beta distribution for pi."
    )->group( "Hyperparameters: Distribution" );

    // pi_b is initialized to ratio num_snps / pi_b_ratio
    // app.add_option(
    //     "--pi-b",
    //     hyper.pi_b,
    //     "Initial value for the `b` parameter of the beta distribution for pi."
    // )->group( "Hyperparameters: Distribution" );

    app.add_option(
        "--pi-b-ratio",
        hyper.pi_b_ratio,
        "Ratio to compute the initial value for the `b` parameter of the beta distribution for pi. "
        "The `b` parameter is then computed as `num_snps / pi_b_ratio`."
    )->group( "Hyperparameters: Distribution" );

    app.add_option(
        "--sigma-1-a",
        hyper.sigma_1_a,
        "Initial value for the `a` parameter of the gamma distribution for sigma_1."
    )->group( "Hyperparameters: Distribution" );

    app.add_option(
        "--sigma-1-b",
        hyper.sigma_1_b,
        "Initial value for the `b` parameter of the gamma distribution for sigma_1."
    )->group( "Hyperparameters: Distribution" );

    app.add_option(
        "--sigma-e-a",
        hyper.sigma_e_a,
        "Initial value for the `a` parameter of the gamma distribution for sigma_e."
    )->group( "Hyperparameters: Distribution" );

    app.add_option(
        "--sigma-e-b",
        hyper.sigma_e_b,
        "Initial value for the `b` parameter of the gamma distribution for sigma_e."
    )->group( "Hyperparameters: Distribution" );

    // ---------------------------------------
    //     Sample validity
    // ---------------------------------------

    app.add_option(
        "--min-beta",
        hyper.min_beta,
        "Minimum value for a beta draw to be considered. Below that, beta is set to 0."
    )->group( "Hyperparameters: Sample Validity" );

    app.add_option(
        "--large-beta",
        hyper.large_beta,
        "Threshold for a beta to be considered large, for the sample statistics computation."
    )->group( "Hyperparameters: Sample Validity" );

    // ---------------------------------------
    //     Convergence Criteria
    // ---------------------------------------

    app.add_option(
        "--num-top-betas",
        hyper.num_top_betas,
        "Test convergence of the top n beta values."
    )->group(
        "Hyperparameters: Convergence Criteria"
    );

    app.add_option(
        "--geweke-first",
        hyper.geweke_first,
        "First fraction of values for the Geweke (1992) z-score to test convergence."
    )->group(
        "Hyperparameters: Convergence Criteria"
    )->check(
        CLI::Range( static_cast<double>(0.0), static_cast<double>(1.0) )
    );

    app.add_option(
        "--geweke-last",
        hyper.geweke_last,
        "Last fraction of values for the Geweke (1992) z-score to test convergence."
    )->group(
        "Hyperparameters: Convergence Criteria"
    )->check(
        CLI::Range( static_cast<double>(0.0), static_cast<double>(1.0) )
    );

    app.add_option(
        "--geweke-intervals",
        hyper.geweke_intervals,
        "Intervals of values for the Geweke (1992) z-score to test convergence."
    )->group(
        "Hyperparameters: Convergence Criteria"
    );

    app.add_option(
        "--max-zscore",
        hyper.convergence_max_zscore,
        "Maximum z-score that all tests have to be below to reach convergence."
    )->group( "Hyperparameters: Convergence Criteria" );
}

// -------------------------------------------------------------------------
//     initialize_hyperparams
// -------------------------------------------------------------------------

/**
 * @brief Init the hyper params using the command line options and data.
 */
inline void initialize_hyperparams( Data const& data, MultiplicativeHyperparams& hyper )
{
    hyper.pi_b = data.num_snps / hyper.pi_b_ratio;
    if( hyper.pi_b <= 1.0 ) {
        throw std::runtime_error( "Too few SNPs, chain would not converge." );
    }
    if( hyper.pi_b <= 10000 ) {
        hyper.pi_b = 10000;
    }
}

#endif // include guard
