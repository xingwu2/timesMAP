#ifndef TIMESMAP_MULTIPLICATIVE_TRACE_H_
#define TIMESMAP_MULTIPLICATIVE_TRACE_H_

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

#include "common/input_data.hpp"
#include "multiplicative/hyper_params.hpp"
#include "multiplicative/state.hpp"

#include <string>
#include <vector>

// =================================================================================================
//     Multiplicative Trace
// =================================================================================================

/**
 * @brief Helper class to handle the trace of a multiplicative model chain, and compute statistics.
 */
class MultiplicativeTrace
{
public:

    // -------------------------------------------------------------------------
    //     Inner Classes
    // -------------------------------------------------------------------------

    struct Entry
    {
        size_t iteration = 0;
        MultiplicativePosteriors post;
        MultiplicativeStatistics stats;
        std::vector<double> top_betas;
    };

    // -------------------------------------------------------------------------
    //     Constructors and Rule of Five
    // -------------------------------------------------------------------------

    MultiplicativeTrace() = default;
    ~MultiplicativeTrace() = default;

    MultiplicativeTrace( MultiplicativeTrace const& ) = default;
    MultiplicativeTrace( MultiplicativeTrace&& )      = default;

    MultiplicativeTrace& operator= ( MultiplicativeTrace const& ) = default;
    MultiplicativeTrace& operator= ( MultiplicativeTrace&& )      = default;

    // -------------------------------------------------------------------------
    //     Public Members
    // -------------------------------------------------------------------------

    void push_to_trace(
        size_t iteration,
        MultiplicativePosteriors const& post,
        MultiplicativeStatistics const& stats,
        MultiplicativeHyperparams const& hyper
    );

    void pop_from_trace( size_t num_elements );

    void summarize_trace( Data const& data, std::string const& file_prefix ) const;
    void write_trace( std::string const& file_prefix ) const;

    // -------------------------------------------------------------------------
    //     Member Variables
    // -------------------------------------------------------------------------

    std::vector<Entry> entries;

};

#endif // include guard
