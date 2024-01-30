#ifndef TIMESMAP_TOOLS_CLI_FORMATTER_H_
#define TIMESMAP_TOOLS_CLI_FORMATTER_H_

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

#include "tools/misc.hpp"

#include <string>

// =================================================================================================
//      CLI11 Formatter
// =================================================================================================

class CliFormatter : public CLI::Formatter
{
public:

    // -----------------------------------------------------------
    //     Overridden Formatter Functions
    // -----------------------------------------------------------

    virtual std::string make_subcommand( CLI::App const* sub) const override
    {
        auto const lcol = sub->get_name();
        auto const rcol = sub->get_description();
        return format_columns( "  " + lcol, rcol, get_column_width() );
    }

    virtual std::string make_option( CLI::Option const* opt, bool is_positional ) const override
    {
        auto const lcol = make_option_name(opt, is_positional) + make_option_opts(opt);
        auto const rcol = make_option_desc(opt);
        return format_columns( "  " + lcol, rcol, get_column_width() );
    }

};

#endif // include guard
