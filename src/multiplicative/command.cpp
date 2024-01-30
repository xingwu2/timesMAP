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

#include "multiplicative/command.hpp"

#include "tools/cli_setup.hpp"

#include "genesis/utils/core/logging.hpp"
#include "genesis/utils/core/options.hpp"
#include "genesis/utils/core/thread_pool.hpp"

#include <memory>

// =================================================================================================
//      Command
// =================================================================================================

void setup_multiplicative_command( CLI::App& app )
{
    // Make a sub-command
    auto sub = app.add_subcommand(
        "multiplicative",
        "Run the multiplicative model."
    );

    // Add all our options to it
    auto cli = std::make_shared<MultiplicativeCLI>();
    add_data_files_cli( *sub, cli->files );
    add_runparam_cli( *sub, cli->run );
    add_multiplicative_hyper_params_cli( *sub, cli->hyper );

    // Set the function to be run when the command is issued
    sub->callback( cli_callback_wrapper(
        sub,
        [ cli ]() {
            run_multiplicative_command( *cli );
        }
    ));
}

void run_multiplicative_command( MultiplicativeCLI const& cli )
{
    // Get the input
    auto const data = read_input_data( cli.files );
    initialize_hyperparams( data, const_cast<MultiplicativeHyperparams&>( cli.hyper ));

    // Create samplers for each chain.
    std::vector<MultiplicativeSampler> samplers;
    for( size_t i = 0; i < cli.run.num_chains; ++i ) {
        samplers.emplace_back(
            cli.run, i, cli.run.random_seed + i, data, cli.hyper
        );
    }

    // Set up the thread pool for processing.
    genesis::utils::Options::get().init_global_thread_pool( cli.run.num_threads );
    auto thread_pool = genesis::utils::Options::get().global_thread_pool();

    // Run all chains by adding them as tasks to the thread pool, and wait for them to finish.
    thread_pool->parallel_for(
        0, cli.run.num_chains,
        [&]( size_t i ){
            samplers[i].run_chain();
        }
    ).wait();

    // Summarize and write the trace
    for( size_t i = 0; i < cli.run.num_chains; ++i ) {
        samplers[i].get_trace().summarize_trace( data, cli.run.output_prefix + std::to_string( i ));
        samplers[i].get_trace().write_trace( cli.run.output_prefix + std::to_string( i ));
    }

    // TODO output median of means and other summaries
}
