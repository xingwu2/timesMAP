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

#include "additive/trace.hpp"

#include "genesis/utils/core/logging.hpp"
#include "genesis/utils/io/output_target.hpp"
#include "genesis/utils/math/moments.hpp"
#include "genesis/utils/math/ranking.hpp"
#include "genesis/utils/text/string.hpp"

// =================================================================================================
//     Public Members
// =================================================================================================

void AdditiveTrace::push_to_trace(
    size_t iteration,
    AdditivePosteriors const& post,
    AdditiveStatistics const& stats,
    AdditiveHyperparams const& hyper
) {
    // Prepare a trace entry with copies of the current samplign cache
    AdditiveTrace::Entry entry;
    entry.iteration = iteration;
    entry.post = post;
    entry.stats = stats;

    // For the zscore test, we also want to use the top n absolute beta values of each iteration.
    // We use an efficient selection algorithm that does not need to sort the whole vector.
    // We get the result with the signs still there, but sorted by absolute value, largest first.
    auto abs_greater = []( double l, double r ){
        return std::abs( l ) > std::abs( r );
    };
    entry.top_betas = genesis::utils::n_first_elements(
        post.beta.begin(), post.beta.end(), hyper.num_top_betas, abs_greater
    );

    // Move the entry to the end of the trace
    entries.push_back( std::move( entry ));
}

void AdditiveTrace::pop_from_trace( size_t num_elements )
{
    num_elements = std::min( num_elements, entries.size() );
    entries.erase( entries.begin(), entries.begin() + num_elements );
}

void AdditiveTrace::summarize_trace( Data const& data, std::string const& file_prefix ) const
{
    using namespace genesis::utils;

    // Prepare summary statistics of all posterior variables
    auto mom_alpha = std::vector<Moments>( data.num_covar );
    auto mom_beta  = std::vector<Moments>( data.num_snps );
    auto mom_gamma = std::vector<Moments>( data.num_snps );
    Moments mom_pi;
    Moments mom_sigma_1;
    Moments mom_sigma_e;
    Moments mom_poly;

    // Use the whole trace for the statistics
    for( auto const& entry : entries ) {
        // Push alpha
        assert( data.num_covar == entry.post.alpha.size() );
        for( size_t c = 0; c < data.num_covar; ++c ) {
            mom_alpha[c].push( entry.post.alpha[c] );
        }

        // Push beta
        assert( data.num_snps == entry.post.beta.size() );
        for( size_t s = 0; s < data.num_snps; ++s ) {
            mom_beta[s].push( entry.post.beta[s] );
        }

        // Push gamma
        assert( data.num_snps == entry.post.gamma.size() );
        for( size_t s = 0; s < data.num_snps; ++s ) {
            mom_gamma[s].push( static_cast<double>( entry.post.gamma[s] ));
        }

        // Push others
        mom_pi.push( entry.post.pi );
        mom_sigma_1.push( entry.post.sigma_1 );
        mom_sigma_e.push( entry.post.sigma_e );
        mom_poly.push( entry.stats.polygenicity );
    }

    // Debug print of results
    LOG_DBG << "Trace summary:";
    LOG_DBG << "alpha";
    for( auto const& mom_a : mom_alpha ) {
        LOG_DBG << " - " << mom_a.mean() << " ± " << mom_a.stddev();
    }
    LOG_DBG << "beta";
    for( auto const& mom_b : mom_beta ) {
        LOG_DBG << " - " << mom_b.mean() << " ± " << mom_b.stddev();
    }
    LOG_DBG << "gamma";
    for( auto const& mom_g : mom_gamma ) {
        LOG_DBG << " - " << mom_g.mean() << " ± " << mom_g.stddev();
    }
    LOG_DBG << "pi           " << mom_pi.mean() << " ± " << mom_pi.stddev();
    LOG_DBG << "sigma_1      " << mom_sigma_1.mean() << " ± " << mom_sigma_1.stddev();
    LOG_DBG << "sigma_e      " << mom_sigma_e.mean() << " ± " << mom_sigma_e.stddev();
    LOG_DBG << "polygenicity " << mom_poly.mean() << " ± " << mom_poly.stddev();

    // Also write a trace summary file.
    auto summary_target = to_file( file_prefix + "_trace_summary.csv" );
    (*summary_target) << "value\tmean\tstddev\n";
    (*summary_target) << "pi\t" << mom_pi.mean() << "\t" << mom_pi.stddev() << "\n";
    (*summary_target) << "sigma_1\t" << mom_sigma_1.mean() << "\t" << mom_sigma_1.stddev() << "\n";
    (*summary_target) << "sigma_e\t" << mom_sigma_e.mean() << "\t" << mom_sigma_e.stddev() << "\n";
    (*summary_target) << "polygenicity\t" << mom_poly.mean() << "\t" << mom_poly.stddev() << "\n";
}

void AdditiveTrace::write_trace( std::string const& file_prefix ) const
{
    using namespace genesis::utils;

    // Open file targets. Can make them compressed if we want.
    auto var_target   = to_file( file_prefix + "_trace_var.csv" );
    auto alpha_target = to_file( file_prefix + "_trace_alpha.csv" );
    auto beta_target  = to_file( file_prefix + "_trace_beta.csv" );
    auto gamma_target = to_file( file_prefix + "_trace_gamma.csv" );

    // Write headers for the var target, as that one has some mumbp jumbo columns.
    (*var_target) << "iteration\tpi\tsigma_1\tsigma_e\tpolygenicity\tgenetic_var\tpheno_var\t";
    (*var_target) << "large_beta_ratio\tlarge_beta_heritability\ttotal_heritability\n";

    // Now write the whole trace.
    for( auto const& entry : entries ) {
        (*var_target) << entry.iteration << "\t";
        (*var_target) << entry.post.pi << "\t";
        (*var_target) << entry.post.sigma_1 << "\t";
        (*var_target) << entry.post.sigma_e << "\t";
        (*var_target) << entry.stats.polygenicity << "\t";
        (*var_target) << entry.stats.genetic_var << "\t";
        (*var_target) << entry.stats.pheno_var << "\t";
        (*var_target) << entry.stats.large_beta_ratio << "\t";
        (*var_target) << entry.stats.large_beta_heritability << "\t";
        (*var_target) << entry.stats.total_heritability << "\n";

        join( alpha_target->ostream(), entry.post.alpha, "\t" );
        join( beta_target->ostream(),  entry.post.beta, "\t" );
        join( gamma_target->ostream(), entry.post.gamma, "\t" );
        (*alpha_target) << "\n";
        (*beta_target) << "\n";
        (*gamma_target) << "\n";
    }
}
