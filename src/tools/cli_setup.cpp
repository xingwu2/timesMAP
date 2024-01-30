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

#include "tools/cli_setup.hpp"

#include "tools/misc.hpp"

#include "genesis/utils/core/algorithm.hpp"
#include "genesis/utils/core/logging.hpp"
#include "genesis/utils/text/string.hpp"
#include "genesis/utils/tools/date_time.hpp"

#include <stdexcept>
#include <unordered_set>

// =================================================================================================
//      CLI11 Setup Internal Functions
// =================================================================================================

/**
 * @brief Magic constant for the left column width when printing information.
 */
static const size_t left_column_width_ = 35;

void print_header( CLI::App const* sub )
{
    // Print a nice header.
    LOG_BOLD << "timesMAP";
    LOG_BOLD;

    // Get the command usage line.
    std::string usage = sub->get_name();
    auto parent = sub->get_parent();
    while( parent ) {
        usage  = parent->get_name() + " " + usage;
        parent = parent->get_parent();
    }

    // Print basic command information.
    // LOG_BOLD << format_columns( "Invocation:", global_options.command_line(), left_column_width_ );
    // LOG_BOLD << format_columns( "Command:", usage, left_column_width_ );
    LOG_BOLD;
}

std::string get_option_value( CLI::Option const* option )
{
    std::string value;

    // Non-flags
    if( option->get_type_size() != 0 ) {

        // If the option was found on command line
        if( option->count() > 0 ) {
            value = CLI::detail::ini_join( option->results() );
            // value = genesus::utils::join( option->results(), " " );

        // Else use the default
        } else {
            value = option->get_default_str();
            // + " (default)";
        }

    // Flag, one passed
    } else if( option->count() == 1 ) {
        value = "true";

    // Flag, multiple passed
    } else if( option->count() > 1 ) {
        value = std::to_string(option->count());

    // Flag, not present
    } else if( option->count() == 0 ) {
        value = "false";
    }

    return value;
}

void print_option_values( CLI::App const* subcommand )
{
    // Store per-group output, so that it is properly sorted.
    // The vector keeps the order in which groups are added.
    // Its content are: group name and full output to be printed.
    std::vector<std::pair<std::string, std::string>> group_output;

    // Helper function to add text to an option group.
    // First check if the group was already used, and if not add it.
    auto get_group_content = [ &group_output ]( std::string const& name ) -> std::string& {
        for( auto& elem : group_output ) {
            if( elem.first == name ) {
                return elem.second;
            }
        }
        group_output.push_back({ name, "" });
        return group_output.back().second;
    };

    // Add output for each option.
    for( auto option : subcommand->get_options()) {

        // Do not add help option.
        if( option->get_name() == "-h,--help" || option->get_name() == "--help" ) {
            continue;
        }

        // Do not add options in the hidden group, using two ways to specify this.
        if(
            option->get_group().empty() ||
            genesis::utils::to_lower( option->get_group() ) == "hidden"
        ) {
            continue;
        }

        // Add the option to its group.
        auto const line = format_columns(
            "  " + option->get_name(),
            get_option_value( option ),
            left_column_width_
        );
        get_group_content( option->get_group() ) += line;
    }

    // Now we have a nicely sorted list of all options in their groups.
    // Print them!
    for( auto const& group : group_output ) {
        LOG_BOLD << group.first << ":";
        LOG_BOLD << group.second;
        LOG_BOLD;
    }
}

// =================================================================================================
//      CLI11 Setup
// =================================================================================================

std::function<void()> cli_callback_wrapper(
    CLI::App const*          subcommand,
    std::function<void()>    run_function
) {
    return [ subcommand, run_function ](){

        // Run the global options callback. Need to this before everything else,
        // so that the number of threads etc are properly set.
        // global_options.run_global();

        // Print out the full header, with all option values (including the number of threads
        // that was just set by the above global options callback).
        print_header( subcommand );
        print_option_values( subcommand );

        LOG_MSG << "Started "
                << genesis::utils::current_date() << " "
                << genesis::utils::current_time()
        ;
        LOG_BOLD;

        // Run the actual command callback function.
        run_function();

        LOG_BOLD;
        LOG_MSG << "Finished "
                << genesis::utils::current_date() << " "
                << genesis::utils::current_time()
        ;
    };
}

// =================================================================================================
//      Checks and Helpers
// =================================================================================================

void fix_cli_default_values( CLI::App& app )
{
    // Make all option capture their defaults now!
    for( auto option : app.get_options()) {
        // Stupid CLI only exposes const pointers... but the misuse here works,
        // as the function uses a non-const App.
        const_cast<CLI::Option*>( option )->capture_default_str();
    }

    // Recursively run this for subcommands.
    for( auto subcom : app.get_subcommands({}) ) {
        fix_cli_default_values( *subcom );
    }
}

void set_module_help_group( CLI::App& module, std::string const& group_name )
{
    for( auto subcom : module.get_subcommands({}) ) {

        // Get the current settings for the help flag.
        auto const name = subcom->get_help_ptr()->get_name(false, true);
        auto const desc = subcom->get_help_ptr()->get_description();

        // First remove it, then add it again. This way, it is the last one to be added,
        // which is nicer for the help message.
        subcom->set_help_flag();
        subcom->set_help_flag( name, desc );
        subcom->get_help_ptr()->group( group_name );
    }
}
