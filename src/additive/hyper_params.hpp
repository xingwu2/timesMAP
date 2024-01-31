#ifndef TIMESMAP_ADDITIVE_HYPER_PARAMS_H_
#define TIMESMAP_ADDITIVE_HYPER_PARAMS_H_

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
#include "multiplicative/hyper_params.hpp"

#include <string>
#include <vector>

// =================================================================================================
//     Hyperparameters
// =================================================================================================

/**
 * @brief Collection of all hyper params for the additive model.
 *
 * This is currently _almost_ the same as the MultiplicativeHyperparams, except for one default
 * value. To avoid code duplication, we instead set the desired default value in the constructor
 * of the class here. That's a bit hacky, but well, for now probbaly better than duplicating the
 * whole class. Might need to refactor later.
 */
struct AdditiveHyperparams : public MultiplicativeHyperparams
{
    AdditiveHyperparams()
    {
        MultiplicativeHyperparams::pi_b_ratio = 50.0;
    }
};

// -------------------------------------------------------------------------
//     add_additive_hyper_params_cli
// -------------------------------------------------------------------------

/**
 * @brief Add all command line options for the additive hyper params to a given command.
 */
inline void add_additive_hyper_params_cli( CLI::App& app, AdditiveHyperparams& hyper )
{
    // The struct is currently derived from MultiplicativeHyperparams, and has the exact same
    // params, so we can re-use the function here.
    add_multiplicative_hyper_params_cli( app, hyper );
}

// -------------------------------------------------------------------------
//     initialize_hyperparams
// -------------------------------------------------------------------------

/**
 * @brief Init the hyper params using the command line options and data.
 *
 * This function differs from the multiplicative one, so we define it here fully.
 */
inline void initialize_hyperparams( Data const& data, AdditiveHyperparams& hyper )
{
    hyper.pi_b = data.num_snps / hyper.pi_b_ratio;
    if( hyper.pi_b < 1.0 ) {
        throw std::runtime_error( "Too few SNPs, chain would not converge." );
    }
}

#endif // include guard
