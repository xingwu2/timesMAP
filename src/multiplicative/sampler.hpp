#ifndef TIMESMAP_MULTIPLICATIVE_SAMPLER_H_
#define TIMESMAP_MULTIPLICATIVE_SAMPLER_H_

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

#include "common/sampler.hpp"
#include "multiplicative/hyper_params.hpp"
#include "multiplicative/state.hpp"
#include "multiplicative/trace.hpp"

#include <cassert>
#include <string>
#include <vector>

// =================================================================================================
//     Multiplicative Sampler
// =================================================================================================

class MultiplicativeSampler : public Sampler
{
public:

    // -------------------------------------------------------------------------
    //     Constructors and Rule of Five
    // -------------------------------------------------------------------------

    MultiplicativeSampler(
        Runparams const& run,
        size_t chain_num,
        size_t random_seed,
        Data const& data,
        MultiplicativeHyperparams const& hyper
    );

    virtual ~MultiplicativeSampler() = default;

    MultiplicativeSampler( MultiplicativeSampler const& ) = default;
    MultiplicativeSampler( MultiplicativeSampler&& )      = default;

    MultiplicativeSampler& operator= ( MultiplicativeSampler const& ) = default;
    MultiplicativeSampler& operator= ( MultiplicativeSampler&& )      = default;

    // -------------------------------------------------------------------------
    //     Init
    // -------------------------------------------------------------------------

    void init_posteriors();
    void init_cache();

    // -------------------------------------------------------------------------
    //     Helpers
    // -------------------------------------------------------------------------

    template< class Generator >
    double draw_pi( double a, double b, Generator& gen )
    {
        // Need to get a beta distrib by using two gamma distribs,
        // see https://stackoverflow.com/a/10359049
        std::gamma_distribution<> dist_gamma_a( a, 1.0 );
        std::gamma_distribution<> dist_gamma_b( b, 1.0 );
        auto const draw_a = dist_gamma_a( gen );
        auto const draw_b = dist_gamma_b( gen );
        return draw_a / ( draw_a + draw_b );
    }

    template< class Generator >
    double draw_sigma_std_dev( double a, double b, Generator& gen )
    {
        // see https://stackoverflow.com/a/10359049
        std::gamma_distribution<> dist_gamma( a, 1.0 / b );
        auto const sigma_1_neg2 = dist_gamma( gen );
    	return std::sqrt( 1.0 / sigma_1_neg2 );
    }

    std::vector<double> circle_product_matrix(
        genesis::utils::Matrix<unsigned char> const& x, std::vector<double> const& beta
    );

    // -------------------------------------------------------------------------
    //     Sampling
    // -------------------------------------------------------------------------

    double sample_pi();
    double sample_sigma_1();
    double sample_sigma_e();
    void sample_and_update_alpha();
    void sample_and_update_gamma();
    void sample_and_update_beta();

    virtual void sampling_step() override;

    // -------------------------------------------------------------------------
    //     Convergence
    // -------------------------------------------------------------------------

    virtual bool has_converged() override;

    // -------------------------------------------------------------------------
    //     Stats
    // -------------------------------------------------------------------------

    MultiplicativeStatistics compute_sample_statistics();
    virtual bool valid_sample() override;

    // -------------------------------------------------------------------------
    //     Trace
    // -------------------------------------------------------------------------

    virtual void push_to_trace( size_t iteration ) override
    {
        assert( hyper );
        trace.push_to_trace( iteration, post, stats, hyper );
    }

    virtual void pop_from_trace( size_t num_elements ) override
    {
        trace.pop_from_trace( num_elements );
    }

    MultiplicativeTrace const& get_trace() const
    {
        return trace;
    }

    // -------------------------------------------------------------------------
    //     Private Member Variables
    // -------------------------------------------------------------------------

private:

    // External data pointers, to avoid copies.
    Data const& data;
    MultiplicativeHyperparams const& hyper;

    // Internal data for the current state.
    MultiplicativeCache cache;
    MultiplicativePosteriors post;
    MultiplicativeStatistics stats;

    // Keep track of the chain.
    MultiplicativeTrace trace;

};

#endif // include guard
