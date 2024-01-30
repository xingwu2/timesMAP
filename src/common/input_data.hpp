#ifndef TIMESMAP_COMMON_INPUT_DATA_H_
#define TIMESMAP_COMMON_INPUT_DATA_H_

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

#include "genesis/utils/containers/matrix.hpp"
#include "genesis/utils/containers/matrix/operators.hpp"

#include <string>
#include <vector>

// =================================================================================================
//     Data CLI
// =================================================================================================

/**
 * @brief Helper to collect the file paths for input.
 */
struct DataFiles
{
    std::string x_file;
    std::string y_file;
    std::string c_file;

    std::string separator_char = "tab";
};

void add_data_files_cli( CLI::App& app, DataFiles& files );

// =================================================================================================
//     Data Class
// =================================================================================================

/**
 * @brief All input data for the models, in concise matrices.
 *
 * We use a struct with all-public members here for simplicity, as we will only ever hand this
 * over to the sampler via const reference anyway.
 */
struct Data
{
    // Data
    genesis::utils::Matrix<unsigned char> x; // num_indiv x num_snps
    genesis::utils::Matrix<double> y;        // num_indiv
    genesis::utils::Matrix<double> c;        // num_indiv x num_covar

    // Helpers for intuitive naming
    size_t num_indiv;
    size_t num_snps;
    size_t num_covar;
};

// =================================================================================================
//     Input File Reading
// =================================================================================================

/**
 * @brief Read the genotype matrix (called `x`) from a simple tabular file.
 */
genesis::utils::Matrix<unsigned char> read_genotype_matrix(
    std::string const& x_file,
    char separator_char = '\t'
);

/**
 * @brief Read the phenotype matrix (called `y`) from a simple tabular file.
 */
genesis::utils::Matrix<double> read_phenotype_matrix(
    std::string const& y_file,
    char separator_char = '\t'
);

/**
 * @brief Read the covariate matrix (called `c`) from a simple tabular file.
 */
genesis::utils::Matrix<double> read_covariate_matrix(
    std::string const& c_file,
    char separator_char = '\t'
);

/**
 * @brief Wrapper function to read all three files from tabular input.
 */
Data read_input_data(
    std::string const& x_file,
    std::string const& y_file,
    std::string const& c_file,
    char separator_char = '\t'
);

/**
 * @brief Wrapper function to read all three files from tabular input using the CLI options.
 */
Data read_input_data( DataFiles const& files );

#endif // include guard
