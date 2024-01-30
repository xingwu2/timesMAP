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

#include "multiplicative/sampler.hpp"

#include "common/geweke.hpp"

#include "genesis/utils/core/logging.hpp"
#include "genesis/utils/math/moments.hpp"

#include <functional>
#include <numeric>

// =================================================================================================
//     Constructors and Rule of Five
// =================================================================================================

MultiplicativeSampler::MultiplicativeSampler(
    Runparams const& run,
    size_t chain_num,
    size_t random_seed,
    Data const& data,
    MultiplicativeHyperparams const& hyper
)
    : Sampler( run, chain_num, random_seed )
    , data( data )
    , hyper( hyper )
{
    // Boundary checks. This model only works with a single phenotype.
    if( data.y.cols() != 1 ) {
        throw std::invalid_argument( "Input y has more than one column" );
    }

    // Init internal data
    init_posteriors();
    init_cache();
    trace.entries.reserve( run.test_convergence_start );
}

// =================================================================================================
//     Init
// =================================================================================================

// -------------------------------------------------------------------------
//     init_posteriors
// -------------------------------------------------------------------------

void MultiplicativeSampler::init_posteriors()
{
    auto const num_snps  = data.x.cols();
    auto const num_covar = data.c.cols();

    // Init scalars randomly (needed for multi-threaded chains)
    assert( hyper.pi_b > 1.0 );
    post.pi = draw_pi( hyper.pi_a, hyper.pi_b, rand_gen );
    post.sigma_1 = draw_sigma_std_dev( hyper.sigma_1_a, hyper.sigma_1_b, rand_gen );
    post.sigma_e = draw_sigma_std_dev( hyper.sigma_e_a, hyper.sigma_e_b, rand_gen );
    LOG_DBG << "Initial pi: " << post.pi;
    LOG_DBG << "pi_b: " << hyper.pi_b;
    LOG_DBG << "Initial sigma_1: " << post.sigma_1;
    LOG_DBG << "Initial sigma_e: " << post.sigma_e;

    // Init alpha
    post.alpha.resize( num_covar );
    std::uniform_real_distribution<> dist_uniform_0_1( 0.0, 1.0 );
    for( auto& v : post.alpha ) {
        v = dist_uniform_0_1( rand_gen );
    }

    // Init gamma
    post.gamma.resize( num_snps );
    // std::binomial_distribution<> dist_binomial( 1, post.pi );
    std::bernoulli_distribution dist_bernoulli( post.pi );
    for( auto& v : post.gamma ) {
        // v = dist_binomial( rand_gen );
        v = dist_bernoulli( rand_gen );
    }
    //LOG_DBG << "Initial gamma: " << post.gamma;
    // Init beta
    post.beta.resize( num_snps );
    std::normal_distribution<> dist_norm{ 0.0, post.sigma_1 };
    for( size_t i = 0; i < num_snps; ++i ) {
        post.beta[i] = post.gamma[i] == 0 ? 0.0 : dist_norm( rand_gen );
    }
    //LOG_DBG << "Initial beta: " << post.beta;
    // return post;
}

// -------------------------------------------------------------------------
//     init_cache
// -------------------------------------------------------------------------

void MultiplicativeSampler::init_cache()
{
    // Init the caches.
    cache.c_col_sum_squared.resize( data.num_covar, 0.0 );
    cache.product_c_alpha.resize( data.num_indiv, 0.0 );
    //cache.circle_product_x_beta.resize( data.num_indiv, 0.0 );

    // Pre-compute the per-column sum of squares of c.
    // This is constant throughout the run.
    for( size_t c = 0; c < data.num_covar; ++c ) {
        for( size_t i = 0; i < data.num_indiv; ++i ) {
            cache.c_col_sum_squared[c] += data.c(i, c) * data.c(i, c);
        }
    }
    // compute the circle product x_beta
    cache.circle_product_x_beta = circle_product_matrix( data.x, post.beta );

    //compute the c_alpha (matrix multiplication)
     for( size_t i = 0; i < data.num_indiv; ++i ) {
        for( size_t c = 0; c < data.num_covar; ++c ) {
            cache.product_c_alpha[i] += data.c(i, c) * post.alpha[ c ];
        }
    }

    //refresh_cache_cache( data, post, cache );
}

// =================================================================================================
//     Helpers
// =================================================================================================

// -------------------------------------------------------------------------
//     circle_product_matrix
// -------------------------------------------------------------------------

std::vector<double> MultiplicativeSampler::circle_product_matrix(
    genesis::utils::Matrix<unsigned char> const& x, std::vector<double> const& beta
) {
    // Make sure that everything has the right dimensions.
    assert( x.cols() == beta.size() );
    auto result = std::vector<double>( x.rows() );

    // Element-wise multiplication of each row of X with beta,
    // then multiply along the columns to get the product for each row.
    for( size_t r = 0; r < x.rows(); ++r ) {
        result[r] = 1.0;
        for( size_t c = 0; c < x.cols(); ++c ) {
            result[r] *= 1.0 + x( r, c ) * beta[c];
        }
    }
    return result;
}

// =================================================================================================
//     Sampling
// =================================================================================================

// -------------------------------------------------------------------------
//     sample_pi
// -------------------------------------------------------------------------

double MultiplicativeSampler::sample_pi()
{
    // Sum the gammas, and compute new a and b
    size_t sum_gamma_p = 0;
    size_t sum_gamma_n = 0;
    for( auto g : post.gamma ) {
        assert( g == 0 || g == 1 );
        sum_gamma_p += g;
        sum_gamma_n += 1 - g;
    }
    auto const a_new = static_cast<double>( sum_gamma_p ) + hyper.pi_a;
    auto const b_new = static_cast<double>( sum_gamma_n ) + hyper.pi_b;

    // Draw from distribution
    return draw_pi( a_new, b_new, rand_gen );
}

// -------------------------------------------------------------------------
//     sample_sigma_1
// -------------------------------------------------------------------------

double MultiplicativeSampler::sample_sigma_1()
{
    // Compute new a and b
    auto const num_snps = post.gamma.size();
    size_t sum_gamma = 0;
    double sum_beta_gamma_sqr = 0.0;
    assert( num_snps == post.beta.size() );
    for( size_t s = 0; s < num_snps; ++s ) {
        sum_gamma += post.gamma[s];
        sum_beta_gamma_sqr += post.beta[s] * post.beta[s] * post.gamma[s];
    }
    auto const a_new = 0.5 * static_cast<double>( sum_gamma ) + hyper.sigma_1_a;
	auto const b_new = 0.5 * sum_beta_gamma_sqr + hyper.sigma_1_b;

    // Draw from distribution
    return draw_sigma_std_dev( a_new, b_new, rand_gen );
}

// -------------------------------------------------------------------------
//     sample_sigma_e
// -------------------------------------------------------------------------

double MultiplicativeSampler::sample_sigma_e()
{
    assert( data.y.cols() == 1 );
    assert( data.num_indiv == data.y.rows() );
    assert( data.num_indiv == cache.circle_product_x_beta.size() );
    assert( data.num_indiv == cache.product_c_alpha.size() );

    // Compute sum of squared residuals
    double sum_square_resid = 0.0;
    for( size_t i = 0; i < data.num_indiv; ++i ) {
        auto const resid = data.y( i, 0 ) - cache.circle_product_x_beta[i] - cache.product_c_alpha[i];
        sum_square_resid += resid * resid;
    }

    // Compute new a and b
    auto const a_new = static_cast<double>( data.num_indiv ) / 2.0 + hyper.sigma_e_a;
    auto const b_new = sum_square_resid / 2.0 + hyper.sigma_e_b;

    // Draw from distribution
    return draw_sigma_std_dev( a_new, b_new, rand_gen );
}

// -------------------------------------------------------------------------
//     sample_and_update_alpha
// -------------------------------------------------------------------------

void MultiplicativeSampler::sample_and_update_alpha()
{
    assert( data.num_indiv == data.y.rows() );
    assert( data.num_indiv == cache.circle_product_x_beta.size() );
    assert( data.num_indiv == cache.product_c_alpha.size() );
    assert( data.num_covar == post.alpha.size() );
    assert( data.num_covar == cache.c_col_sum_squared.size() );

    auto const sigma_e_neg2 = std::pow( post.sigma_e, -2 );
    double dot_prod = 0.0;

    if( data.num_covar == 1 ) {
        auto const c = 0;
        for( size_t i = 0; i < data.num_indiv; ++i ) {
            auto const y_x_beta = data.y( i, 0 ) - cache.circle_product_x_beta[i];
            dot_prod += y_x_beta * data.c(i, c);
        }

        auto const new_variance = 1.0 / ( cache.c_col_sum_squared[c] * sigma_e_neg2 );
        auto const new_mean = new_variance * dot_prod * sigma_e_neg2;

        // Draw a new alpha
        std::normal_distribution<> dist_normal( new_mean, std::sqrt( new_variance ));
        post.alpha[c] = dist_normal( rand_gen );

        // recompute the c_alpha
        for( size_t i = 0; i < data.num_indiv; ++i ) {
            cache.product_c_alpha[i] = data.c(i, c) * post.alpha[c];
        }

    } else{
        for( size_t c = 0; c < data.num_covar; ++c ) {
            for( size_t i = 0; i < data.num_indiv; ++i ) {
                // the ith individual c_alpha without the cth covariate
                cache.product_c_alpha[i] -= data.c(i, c) * post.alpha[c];
                auto const c_alpha_negc_i = cache.product_c_alpha[i];
                auto const y_negi = data.y( i, 0 ) - c_alpha_negc_i - cache.circle_product_x_beta[i];
                dot_prod += y_negi * data.c(i, c);
            }
            auto const new_variance = 1.0 / ( cache.c_col_sum_squared[c] * sigma_e_neg2 );
            auto const new_mean = new_variance * dot_prod * sigma_e_neg2;

            // Draw a new alpha
            std::normal_distribution<> dist_normal( new_mean, std::sqrt( new_variance ));
            post.alpha[c] = dist_normal( rand_gen );

            // recompute the c_alpha by adding the c*alpha back in
            for( size_t i = 0; i < data.num_indiv; ++i ) {
                cache.product_c_alpha[i] += data.c(i, c) * post.alpha[c];
            }
        }
    }
}
// -------------------------------------------------------------------------
//     sample_and_update_gamma
// -------------------------------------------------------------------------

void MultiplicativeSampler::sample_and_update_gamma()
{
    assert( data.num_snps == post.gamma.size() );

    // Some constants
    auto const sigma_e_neg2 = std::pow( post.sigma_e, -2 );
    auto const sigma_1_neg2 = std::pow( post.sigma_1, -2 );

    // Update gamma for all snps
    for( size_t s = 0; s < data.num_snps; ++s ) {

        // compute the norm of x_beta_neg_s * x[s]
        double norm_x_beta_x = 0.0;
        double dot_prod_residuals = 0.0;

        for( size_t i = 0; i < data.num_indiv; ++i ) {
            auto const x_beta_negs_i = cache.circle_product_x_beta[i] / (1+ data.x(i,s) * post.beta[s]);
            auto const x_beta_neg_x = x_beta_negs_i * data.x(i, s);
            // Update our sums
            norm_x_beta_x += x_beta_neg_x * x_beta_neg_x;

            // update the residual dot product
             auto const residual = data.y( i, 0 ) - cache.product_c_alpha[i] - x_beta_negs_i;
             dot_prod_residuals += residual * x_beta_neg_x;
        }

        auto const variance = 1.0 / ( norm_x_beta_x * sigma_e_neg2 + sigma_1_neg2 );
        auto const mean = variance * dot_prod_residuals * sigma_e_neg2;

        // Compute the intermediate values
        auto const d = norm_x_beta_x * sigma_e_neg2 / sigma_1_neg2 + 1.0;
        auto const f = std::sqrt( 1.0 / d );
        auto const a = f * std::exp( 0.5 * mean * mean / variance );
        auto const b = ( 1.0 - post.pi ) / ( 1.0 - post.pi + post.pi * a );

        // Draw a new gamma
        std::bernoulli_distribution dist_bernoulli( 1.0 - b );
        post.gamma[s] = dist_bernoulli( rand_gen );
    }

}

// -------------------------------------------------------------------------
//     sample_and_update_beta
// -------------------------------------------------------------------------

void MultiplicativeSampler::sample_and_update_beta()
{
    assert( data.num_indiv == data.x.rows() );
    assert( data.num_indiv == cache.circle_product_x_beta.size() );
    assert( data.num_snps  == data.x.cols() );
    assert( data.num_snps  == post.beta.size() );
    assert( data.num_snps  == post.gamma.size() );

    // Some constants
    auto const sigma_e_neg2 = std::pow( post.sigma_e, -2 );
    auto const sigma_1_neg2 = std::pow( post.sigma_1, -2 );

    // Update beta for all snps
    for( size_t s = 0; s < data.num_snps; ++s ) {

        // first remove the sth beta from the circle product
        for( size_t i = 0; i < data.num_indiv; ++i ) {
            auto const circ_prod = data.x(i, s) * post.beta[s] + 1.0;
            cache.circle_product_x_beta[i] /= circ_prod;
        }

        if( post.gamma[s] == 0 ) {
            // if gamma[s] ==0. then cache.circle_product_x_beta doesnt need adjustment
            post.beta[s] = 0.0;
        } else{
            // compute the norm of x_beta_neg_s * x[s]
            double norm_x_beta_x = 0.0;
            double dot_prod_residuals = 0.0;

            for( size_t i = 0; i < data.num_indiv; ++i ) {
                auto const x_beta_neg_x = cache.circle_product_x_beta[i] * data.x(i, s);
                // Update our sums
                norm_x_beta_x += x_beta_neg_x * x_beta_neg_x;

                // update the residual dot product
                auto const residual = data.y( i, 0 ) - cache.product_c_alpha[i] - cache.circle_product_x_beta[i];
                dot_prod_residuals += residual * x_beta_neg_x;
            }
            auto const variance = 1.0 / ( norm_x_beta_x * sigma_e_neg2 + sigma_1_neg2 );
            auto const mean = variance * dot_prod_residuals * sigma_e_neg2;

             // Draw a new beta
            std::normal_distribution<> dist_normal(mean, std::sqrt( variance ));
            post.beta[s] = dist_normal( rand_gen );

            // Ensure minimum beta threshold
            assert( hyper.min_beta > 0.0 );
            if( std::abs( post.beta[s] ) < hyper.min_beta ) {
                // if the updated beta is 0, then the circle_product_x_beta doesnt need adjustment
                post.beta[s] = 0.0;
            } else{
                for( size_t i = 0; i < data.num_indiv; ++i ) {
                    auto const circ_prod = data.x(i, s) * post.beta[s] + 1.0;
                    cache.circle_product_x_beta[i] *= circ_prod;
                }
            }
        }
    }
}

// -------------------------------------------------------------------------
//     sampling_step
// -------------------------------------------------------------------------

void MultiplicativeSampler::sampling_step()
{
    //assert( std::isfinite( hyper.min_sigma_1 ) && hyper.min_sigma_1 > 0.0 );
    post.sigma_1 = sample_sigma_1();

    // no need to do this condition anymore
    // if( post.sigma_1 < hyper.min_sigma_1 ) {
    //     post.sigma_1 = hyper.min_sigma_1;
    //     post.pi = 0.0;
    // } else {
    //     post.pi = sample_pi();
    // }

    post.pi = sample_pi();
    post.sigma_e = sample_sigma_e();
    sample_and_update_gamma();
    sample_and_update_alpha();
    sample_and_update_beta();

    // Compute the stats, which will be needed later for the convergence test
    stats = compute_sample_statistics();
}

// =================================================================================================
//     Convergence
// =================================================================================================

bool MultiplicativeSampler::has_converged()
{
    // We compute the z scores of all our posterior traces of interest.
    // Our trace stores the posterios etc in a way that is not easy to iterate over,
    // so we need to build vectors of the variables first. This is kind of expensive,
    // but we only run this function every 1000 or so iterations, so that's okay for now.
    // Might optimize later.
    double max_zscore = 0.0;

    // Local helper function to select entries from the trace and update the max score.
    auto update_zscore_with_trace_variable = [&](
        MultiplicativeHyperparams const& hyper, MultiplicativeTrace const& trace,
        std::function<double(MultiplicativeTrace::Entry const&)> select_value,
        double& max_zscore
    ) {
        // Get a vector of the variable along all traces.
        auto const trace_size = trace.entries.size();
        std::vector<double> values;
        values.reserve( trace_size );
        for( size_t t = 0; t < trace_size; ++t ) {
            values.push_back( select_value( trace.entries[t] ));
        }

        // Compute the zscore of the values, and if any is larger than the current one, store it.
        auto const zscores = geweke_z_scores(
            values, hyper.geweke_first, hyper.geweke_last, hyper.geweke_intervals
        );
        for( auto z : zscores ) {
            max_zscore = std::max( max_zscore, std::abs( z ));
        }
    };

    // Test alpha convergence for all covariates
    for( size_t c = 0; c < data.num_covar; ++c ) {
        update_zscore_with_trace_variable(
            hyper, trace,
            [&]( MultiplicativeTrace::Entry const& entry ){
                assert( data.num_covar == entry.post.alpha.size() );
                return entry.post.alpha[c];
            },
            max_zscore
        );
    }

    // Test convergence of the top n beta values
    for( size_t b = 0; b < hyper.num_top_betas; ++b ) {
        update_zscore_with_trace_variable(
            hyper, trace,
            [&]( MultiplicativeTrace::Entry const& entry ){
                assert( hyper.num_top_betas == entry.top_betas.size() );
                return entry.top_betas[b];
            },
            max_zscore
        );
    }

    // Test convergence of large beta_ratio
    update_zscore_with_trace_variable(
        hyper, trace,
        [&]( MultiplicativeTrace::Entry const& entry ){
            return entry.stats.large_beta_ratio;
        },
        max_zscore
    );

    // Test convergence of sigma_1
    // update_zscore_with_trace_variable(
    //     hyper, trace,
    //     [&]( MultiplicativeTrace::Entry const& entry ){
    //        return entry.post.sigma_1;
    //     },
    //     max_zscore
    // );

    // Test convergence of sigma_e
    update_zscore_with_trace_variable(
        hyper, trace,
        [&]( MultiplicativeTrace::Entry const& entry ){
            return entry.post.sigma_e;
        },
        max_zscore
    );

    // Test convergence of total_heritability
    update_zscore_with_trace_variable(
        hyper, trace,
        [&]( MultiplicativeTrace::Entry const& entry ){
            return entry.stats.total_heritability;
        },
        max_zscore
    );

    // Finally, test if the largest zscore we found is good for convergence
    if( max_zscore < hyper.convergence_max_zscore ) {
        LOG_INFO << "Convergence has been reached";
        return true;
    }
    return false;
}

// =================================================================================================
//     Stats
// =================================================================================================

MultiplicativeStatistics MultiplicativeSampler::compute_sample_statistics()
{
    assert( data.num_indiv == data.x.rows() );
    assert( data.num_indiv == data.y.rows() );
    assert( data.num_indiv == cache.circle_product_x_beta.size() );
    assert( data.num_indiv == cache.product_c_alpha.size() );
    assert( data.num_snps  == data.x.cols() );
    assert( data.num_snps  == post.beta.size() );
    assert( data.num_snps  == post.gamma.size() );
    assert( std::isfinite( hyper.large_beta ) && hyper.large_beta > 0.0 );

    // Return value
    MultiplicativeStatistics stats;

    // We need to compute some variances, using a one pass algorithm for that
    genesis::utils::Moments geno_moments;
    genesis::utils::Moments pheno_moments;
    genesis::utils::Moments large_beta_moments;

    // Get the large beta count, and sum up the polygenicity
    size_t large_beta_count = 0;
    for( size_t s = 0; s < data.num_snps; ++s ) {
        if( std::abs( post.beta[s] ) > hyper.large_beta ) {
            ++large_beta_count;
        }
        stats.polygenicity += post.gamma[s];
    }

    // Go through the individuals and compute the relevant variances and sums
    for( size_t i = 0; i < data.num_indiv; ++i ) {
        geno_moments.push( cache.circle_product_x_beta[i] );
        pheno_moments.push( data.y( i, 0 ) - cache.product_c_alpha[i] );

        // Compute circle product for the current individual for its moment
        double circ_prod = 1.0;
        for( size_t s = 0; s < data.num_snps; ++s ) {
            if( std::abs( post.beta[s] ) > hyper.large_beta ) {
                circ_prod *= 1.0 + data.x(i, s) * post.beta[s];
            }
        }
        large_beta_moments.push( circ_prod );
    }

    // Compute our final statistics
    stats.genetic_var = geno_moments.variance();
    stats.pheno_var = pheno_moments.variance();
    stats.large_beta_ratio = static_cast<double>( large_beta_count ) / post.beta.size();
    stats.large_beta_heritability = large_beta_moments.variance() / stats.pheno_var;
    stats.total_heritability = stats.genetic_var / stats.pheno_var;
    return stats;
}

bool MultiplicativeSampler::valid_sample()
{
    // Previous alternative checks, not currently used.
    // bool const herit_too_large = stats.large_beta_heritability > 1.0;
    // bool const herit_larger_total = stats.large_beta_heritability > stats.total_heritability;

    // Check that we got a usable sample
    if( stats.total_heritability > 1 or stats.total_heritability < 0.01 ) {
        LOG_DBG << "unrealistic beta sample:";
        LOG_DBG << " - polygenicity            " << stats.polygenicity;
        LOG_DBG << " - genetic_var             " << stats.genetic_var;
        LOG_DBG << " - pheno_var               " << stats.pheno_var;
        LOG_DBG << " - large_beta_ratio        " << stats.large_beta_ratio;
        LOG_DBG << " - large_beta_heritability " << stats.large_beta_heritability;
        LOG_DBG << " - total_heritability      " << stats.total_heritability;
        return false;
    }
    return true;
}
