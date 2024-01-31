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

#include "additive/command.hpp"

#include "common/command.hpp"
#include "tools/cli_setup.hpp"

#include "genesis/utils/core/logging.hpp"
#include "genesis/utils/core/options.hpp"
#include "genesis/utils/core/thread_pool.hpp"

#include <memory>

// =================================================================================================
//      Command
// =================================================================================================

void setup_additive_command( CLI::App& app )
{
    // Make a sub-command
    auto sub = app.add_subcommand(
        "additive",
        "Run the additive model."
    );

    // Add all our options to it
    auto cli = std::make_shared<AdditiveCLI>();
    add_data_files_cli( *sub, cli->files );
    add_runparam_cli( *sub, cli->run );
    add_additive_hyper_params_cli( *sub, cli->hyper );

    // Set the function to be run when the command is issued
    sub->callback( cli_callback_wrapper(
        sub,
        [ cli ]() {
            run_additive_command( *cli );
        }
    ));
}

void run_additive_command( AdditiveCLI const& cli )
{
    run_multithreaded_command<AdditiveSampler>( cli );
}
