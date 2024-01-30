#ifndef TIMESMAP_COMMON_RUN_PARAMS_H_
#define TIMESMAP_COMMON_RUN_PARAMS_H_

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

#include <string>
#include <vector>

// =================================================================================================
//     Run Parameters
// =================================================================================================

/**
 * @brief Collection of all global run params used for general Markov chains.
 */
struct Runparams
{
    // Global params
    size_t num_chains = 1;
    size_t num_threads = 1;
    size_t random_seed = 0;
    std::string output_prefix;
    size_t interation_report_interval = 2500;

    // Chain params
    size_t burnin_iterations = 2000;
    size_t regular_iterations = 10000;
    size_t sanity_iterations = 100;
    size_t test_convergence_start = 10000;
    size_t test_convergence_interval = 1000;
    size_t test_convergence_stop = 100000;
};

/**
 * @brief Set up the command line interface for the global run params.
 */
inline void add_runparam_cli( CLI::App& app, Runparams& run )
{
    // Global params

    app.add_option(
        "--num-chains",
        run.num_chains,
        "Number of chains to run in total."
    )->group( "Run Parameters: Global" );

    app.add_option(
        "--num-threads",
        run.num_threads,
        "Number of threads to use for running the chains."
    )->group( "Run Parameters: Global" );

    app.add_option(
        "--random-seed",
        run.random_seed,
        "Initial random seed to use. Incremented by one for each chain."
    )->group( "Run Parameters: Global" );

    app.add_option(
        "--output-prefix",
        run.output_prefix,
        "Prefix to use for output file names, to distinguish runs."
    )->group( "Run Parameters: Global" );

    // Chain params

    app.add_option(
        "--burnin-iterations",
        run.burnin_iterations,
        "Number of initial iterations that are discarded."
    )->group( "Run Parameters: Chain" );

    app.add_option(
        "--regular-iterations",
        run.regular_iterations,
        "Number of iterations to run the chain."
    )->group( "Run Parameters: Chain" );

    app.add_option(
        "--sanity-iterations",
        run.sanity_iterations,
        "Number of sanity interations"
    )->group( "Run Parameters: Chain" );

    app.add_option(
        "--test-convergence-start",
        run.test_convergence_start,
        "TODO"
    )->group( "Run Parameters: Chain" );

    app.add_option(
        "--test-convergence-interval",
        run.test_convergence_interval,
        "TODO"
    )->group( "Run Parameters: Chain" );

    app.add_option(
        "--test-convergence-stop",
        run.test_convergence_stop,
        "TODO"
    )->group( "Run Parameters: Chain" );
}

#endif // include guard
