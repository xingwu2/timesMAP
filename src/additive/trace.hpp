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

#include "multiplicative/trace.hpp"

#include <string>
#include <vector>

// =================================================================================================
//     Additive Trace
// =================================================================================================

/**
 * @brief Helper class to handle the trace of a additive model chain, and compute statistics.
 *
 * As of now, the implementation is exactly the same as for the multiplicative model, so we just
 * re-use the class. If we need different behavior later on, we can implement this here.
 */
using AdditiveTrace = MultiplicativeTrace;

#endif // include guard
