#ifndef TIMESMAP_COMMON_SAMPLER_H_
#define TIMESMAP_COMMON_SAMPLER_H_

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

#include "common/run_params.hpp"

#include <random>
#include <string>
#include <vector>

#include "genesis/utils/core/logging.hpp"

// =================================================================================================
//     Sampler
// =================================================================================================

/**
 * @brief Base class for running one Markov chain.
 *
 * The class is meant to be derived from, to implement the actual sampling functions.
 * Here, we only take care of the overall structure, running a chain until convergence.
 *
 * To run multiple chains, several instances of this class can be run in parallel.
 */
class Sampler
{
public:

    // -------------------------------------------------------------------------
    //     Constructors and Rule of Five
    // -------------------------------------------------------------------------

    Sampler(
        Runparams const& run,
        size_t chain_num,
        size_t random_seed
    )
        : run(run)
        , chain_num( chain_num )
        , rand_gen( random_seed )
    {
        // Optional: Automatic random initialization if no seed given
        if( random_seed == 0 ) {
            std::random_device rand_dev;
            rand_gen.seed( rand_dev() );
        }
    }

    virtual ~Sampler() = default;

    Sampler( Sampler const& ) = default;
    Sampler( Sampler&& )      = default;

    Sampler& operator= ( Sampler const& ) = default;
    Sampler& operator= ( Sampler&& )      = default;

    // -------------------------------------------------------------------------
    //     Public Members
    // -------------------------------------------------------------------------

    /**
     * @brief Run the Markov chain.
     */
    void run_chain()
    {
        size_t it = 1;
        size_t max_it = run.burnin_iterations + run.regular_iterations;

        while( it <= max_it && it <= run.test_convergence_stop ) {
            if( it % run.interation_report_interval == 0 ) {
                LOG_MSG1 << "Chain " << chain_num << ", iteration " << it;
            } else {
                LOG_MSG2 << "Chain " << chain_num << ", iteration " << it;
            }

            // Do all sampling updates.
            sampling_step();

            // Compute sanity metrics to see if this is good. If not, discard and start again,
            // using the values that we just computed as the new distributions.
            if( it > run.sanity_iterations && ! valid_sample()) {
                continue;
            }

            // After the burn-in, we keep a trace of all samples
            if( it > run.burnin_iterations ) {
                push_to_trace( it );
            }

            // After a while, start testing for convergence
            if(
                it > run.burnin_iterations + run.test_convergence_start &&
                it % run.test_convergence_interval == 0
            ) {
                LOG_MSG2 << "Testing convergence";
                if( has_converged() ) {
                    break;
                }

                // If we have not converged yet, we run more iterations.
                // We remove a batch from the trace to reduce memory.
                pop_from_trace( run.test_convergence_interval );
                max_it += run.test_convergence_interval;
            }

            ++it;
        }
        if( it >= run.test_convergence_stop ) {
            LOG_INFO << "Chain has not converged after " << it << " iterations";
        }
    }

    // -------------------------------------------------------------------------
    //     Virtual Members
    // -------------------------------------------------------------------------

    virtual void sampling_step() = 0;
    virtual bool valid_sample() = 0;
    virtual bool has_converged() = 0;
    virtual void push_to_trace( size_t iteration ) = 0;
    virtual void pop_from_trace( size_t num_elements ) = 0;

    // -------------------------------------------------------------------------
    //     Protected Member Variables
    // -------------------------------------------------------------------------

protected:

    Runparams run;

    // Keep track of the number of the chain, for user output only
    size_t chain_num;

    // Random number generator for the chain
    // TODO: other engines are faster - if they are good enough for MCMC, we might want to switch
    std::default_random_engine rand_gen;

};

#endif // include guard
