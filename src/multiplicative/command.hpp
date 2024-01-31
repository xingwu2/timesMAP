#ifndef TIMESMAP_MULTIPLICATIVE_COMMAND_H_
#define TIMESMAP_MULTIPLICATIVE_COMMAND_H_

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

#include "common/input_data.hpp"
#include "common/run_params.hpp"
#include "multiplicative/sampler.hpp"

#include <string>
#include <vector>

// =================================================================================================
//     Multiplicative Command
// =================================================================================================

/**
 * @brief Collection of all input settings to run the model.
 *
 *  The values here typically are set by the user via the command line.
 *  The class contains the general data and run params, as well as the hyper params specific to the
 *  multiplicative model. Other models might adapt that as needed, using different hyper params.
 */
struct MultiplicativeCLI
{
    DataFiles files;
    Runparams run;
    MultiplicativeHyperparams hyper;
};

/**
 * @brief Add the "multiplicative" model command to a given app.
 */
void setup_multiplicative_command( CLI::App& app );

/**
 * @brief Callback function that is used by the above command to do the actual work.
 *
 * The command sets up the chains as needed, runs them, and writes out the results.
 */
void run_multiplicative_command( MultiplicativeCLI const& cli );

#endif // include guard
