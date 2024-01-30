#ifndef TIMESMAP_COMMON_GEWEKE_H_
#define TIMESMAP_COMMON_GEWEKE_H_

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

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <vector>

#include "genesis/utils/math/moments.hpp"

// =================================================================================================
//     Convergence
// =================================================================================================

// -------------------------------------------------------------------------
//     geweke_z_scores
// -------------------------------------------------------------------------

/**
 * @brief Compute the Geweke z-scores for convergence testing.
 *
 * See Geweke (1992) for details. This implementation is inspired by the implementation of the
 * function in pymc2, see https://github.com/pymc-devs/pymc2/blob/6b1b51ddea1a74c50d9a027741252b30810b29e0/pymc/diagnostics.py#L236
 */
inline std::vector<double> geweke_z_scores(
    std::vector<double> const& x, double first = 0.1, double last = 0.5, size_t intervals = 20
) {
    // Input checks
    if(
        ! std::isfinite(first) || first <= 0.0 || first >= 1.0 ||
        ! std::isfinite(last)  || last  <= 0.0 || last  >= 1.0 ||
        first >= last || first + last >= 1.0 || intervals == 0
    ) {
        throw std::invalid_argument(
            "Invalid intervals for Geweke convergence analysis: first=" + std::to_string( first ) +
            ", last=" + std::to_string( last ) + ", intervals=" + std::to_string( intervals )
        );
    }

    // Prepare and check boundary cases
    std::vector<double> zscores;
    if( x.size() < 2 ) {
        return zscores;
    }

    // Last index
    auto const end_idx = x.size() - 1;

    // Start intervals going up to the <last>% of the chain
    auto const last_start_idx = ( 1.0 - last ) * static_cast<double>( end_idx );

    // Loop over intervals
    size_t start_idx = 0;
    while( start_idx < static_cast<size_t>( last_start_idx )) {

        // Compute moments of first slice
        genesis::utils::Moments first_slice_moments;
        auto const first_slice_end_idx = start_idx + static_cast<size_t>(first * (end_idx - start_idx));
        assert( start_idx < x.size() && first_slice_end_idx < x.size() );
        assert( start_idx <= first_slice_end_idx );
        for( size_t i = start_idx; i < first_slice_end_idx; ++i ) {
            first_slice_moments.push( x[i] );
        }

        // Compute moments of last slice
        genesis::utils::Moments last_slice_moments;
        auto const last_slice_start_idx = end_idx - static_cast<size_t>(last * (end_idx - start_idx));
        assert( last_slice_start_idx < x.size() );
        for( size_t i = last_slice_start_idx; i < x.size(); ++i ) {
            last_slice_moments.push( x[i] );
        }

        // Compute z statistic
        auto const sqrt_var_sum = std::sqrt(
            first_slice_moments.variance() + last_slice_moments.variance()
        );
        auto const z = ( first_slice_moments.mean() - last_slice_moments.mean() ) / sqrt_var_sum;

        // Store result. The python implementation that we based this function on
        // also stores the start_idx here, but we do not need that for our case here, see
        // https://github.com/owlas/pymc3/blob/master/pymc3/diagnostics.py
        zscores.push_back( z );

        // Move to the next start index
        start_idx += static_cast<size_t>( last_start_idx / static_cast<double>( intervals - 1 ));
    }
    return zscores;
}


#endif // include guard
