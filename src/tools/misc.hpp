#ifndef TIMESMAP_TOOLS_MISC_H_
#define TIMESMAP_TOOLS_MISC_H_

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

#include "genesis/utils/text/string.hpp"

#include <iosfwd>
#include <string>
#include <stdexcept>
#include <vector>

// =================================================================================================
//      Formatting
// =================================================================================================

std::string format_columns(
    std::string const& left,
    std::string const& right,
    size_t left_width
);

void write_columns(
    std::ostream& out,
    std::string const& left,
    std::string const& right,
    size_t left_width,
    size_t right_width
);

// =================================================================================================
//      Misc
// =================================================================================================

/**
 * @brief Helper function to get the char representation for table separator chars,
 * given its textual description from the option.
 */
char translate_separator_char( std::string const& separator_char );

#endif // include guard
