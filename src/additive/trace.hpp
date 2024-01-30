#ifndef TIMESMAP_ADDITIVE_TRACE_H_
#define TIMESMAP_ADDITIVE_TRACE_H_

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
#include "additive/hyper_params.hpp"
#include "additive/state.hpp"

#include <string>
#include <vector>

// =================================================================================================
//     Additive Trace
// =================================================================================================

/**
 * @brief Helper class to handle the trace of a additive model chain, and compute statistics.
 */
class AdditiveTrace
{
public:

    // -------------------------------------------------------------------------
    //     Inner Classes
    // -------------------------------------------------------------------------

    struct Entry
    {
        size_t iteration = 0;
        AdditivePosteriors post;
        AdditiveStatistics stats;
        std::vector<double> top_betas;
    };

    // -------------------------------------------------------------------------
    //     Constructors and Rule of Five
    // -------------------------------------------------------------------------

    AdditiveTrace() = default;
    ~AdditiveTrace() = default;

    AdditiveTrace( AdditiveTrace const& ) = default;
    AdditiveTrace( AdditiveTrace&& )      = default;

    AdditiveTrace& operator= ( AdditiveTrace const& ) = default;
    AdditiveTrace& operator= ( AdditiveTrace&& )      = default;

    // -------------------------------------------------------------------------
    //     Public Members
    // -------------------------------------------------------------------------

    void push_to_trace(
        size_t iteration,
        AdditivePosteriors const& post,
        AdditiveStatistics const& stats,
        AdditiveHyperparams const& hyper
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
