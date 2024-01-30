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

#include "genesis/utils/containers/matrix/simple_reader.hpp"
#include "genesis/utils/text/string.hpp"

#include "tools/misc.hpp"

using namespace genesis;
using namespace genesis::utils;

// =================================================================================================
//     Data CLI
// =================================================================================================

void add_data_files_cli( CLI::App& app, DataFiles& files )
{
    app.add_option(
        "--genotype-matrix",
        files.x_file,
        "Input genotype matrix (called `x`), simple tabular file."
    )->group(
        "Input Data"
    )->check(
        CLI::ExistingPath
    )->required();

    app.add_option(
        "--phenotype-matrix",
        files.y_file,
        "Input phenotype matrix (called `y`), simple tabular file."
    )->group(
        "Input Data"
    )->check(
        CLI::ExistingPath
    )->required();

    app.add_option(
        "--covariate-matrix",
        files.c_file,
        "Input covariate matrix (called `c`), simple tabular file."
    )->group(
        "Input Data"
    )->check(
        CLI::ExistingPath
    )->required();

    app.add_option(
        "--separator-char",
        files.separator_char,
        "Separator char between fields of input tabular files."
    )->group(
        "Input Data"
    )->transform(
        CLI::IsMember({ "comma", "tab", "space", "semicolon" }, CLI::ignore_case )
    );
}

// =================================================================================================
//     Local Helper Functions
// =================================================================================================

template<typename T>
void check_matrix_range( Matrix<T> const& matrix, T min = 0, T max = 0 )
{
    for( auto e : matrix ) {
        if( ! std::isfinite( e )) {
            throw std::runtime_error( "Matrix has non-finite values" );
        }
        if(( min != 0 || max != 0 ) && ( e < min || e > max )) {
            throw std::runtime_error(
                "Matrix has values outside of [" + std::to_string( min ) + ", " +
                std::to_string( max ) + "]"
            );
        }
    }
}

// =================================================================================================
//     Input File Reading
// =================================================================================================

Matrix<unsigned char> read_genotype_matrix(
    std::string const& x_file,
    char separator_char
) {
    // Set up fast readers for simple matrices
    auto unsigned_char_reader = MatrixSimpleReader<unsigned char>( separator_char );
    unsigned_char_reader.parse_value_functor( parse_unsigned_integer<unsigned char> );

    // Read x (genotype matrix)
    LOG_INFO << "Reading genotype matrix X";
    auto x_mat = unsigned_char_reader.read( from_file( x_file ));
    LOG_INFO << "x[" << x_mat.rows() << "," << x_mat.cols() << "]";
    LOG_INFO << print( x_mat );
    check_matrix_range<unsigned char>( x_mat, 0, 2 );

    return x_mat;
}

Matrix<double> read_phenotype_matrix(
    std::string const& y_file,
    char separator_char
) {
    // Set up fast readers for simple matrices
    auto double_reader = MatrixSimpleReader<double>( separator_char );
    double_reader.parse_value_functor( parse_float<double> );

    // Read y (phenotype vector)
    LOG_INFO << "Reading phenotype vector Y";
    auto y_mat = double_reader.read( from_file( y_file ));
    // TODO we might want to allow multiple phenotypes later
    if( y_mat.cols() != 1 ) {
        throw std::invalid_argument( "Input y has more than one column" );
    }
    LOG_INFO << "y[" << y_mat.rows() << "]";
    LOG_INFO << print( y_mat );
    check_matrix_range( y_mat );

    return y_mat;
}

Matrix<double> read_covariate_matrix(
    std::string const& c_file,
    char separator_char
) {
    // Set up fast readers for simple matrices
    auto double_reader = MatrixSimpleReader<double>( separator_char );
    double_reader.parse_value_functor( parse_float<double> );

    // Read c (covariates)
    LOG_INFO << "Reading covariate matrix C";
    auto c_mat = double_reader.read( from_file( c_file ));
    LOG_INFO << "c[" << c_mat.rows() << "," << c_mat.cols() << "]";
    LOG_INFO << print( c_mat );
    check_matrix_range( c_mat );

    return c_mat;
}

Data read_input_data(
    std::string const& x_file,
    std::string const& y_file,
    std::string const& c_file,
    char separator_char
) {
    // Read the data
    Data data;
    data.x = read_genotype_matrix(  x_file, separator_char );
    data.y = read_phenotype_matrix( y_file, separator_char );
    data.c = read_covariate_matrix( c_file, separator_char );

    // Safety checks
    if( data.x.rows() != data.y.size() || data.x.rows() != data.c.rows() ) {
        throw std::invalid_argument( "Input data dimensions do not match" );
    }
    if( data.x.rows() == 0 || data.x.cols() == 0 || data.c.cols() == 0 ) {
        throw std::invalid_argument( "Input data is empty" );
    }

    // Set some nicely named helpers
    data.num_indiv = data.x.rows();
    data.num_snps  = data.x.cols();
    data.num_covar = data.c.cols();

    return data;
}

Data read_input_data( DataFiles const& files )
{
    auto const sep_char = translate_separator_char( files.separator_char );
    return read_input_data( files.x_file, files.y_file, files.c_file, sep_char );
}
