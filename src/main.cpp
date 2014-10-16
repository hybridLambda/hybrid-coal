/*
 * hybrid_coal is used to compute gene tree probabilities given species network under coalescent process.
 * 
 * Copyright (C) 2010, 2011 Sha (Joe) Zhu
 * 
 * This file is part of hybrid_coal.
 * 
 * hybrid_coal is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

/*! \file main.cpp
 *  \brief Main function of hybrid_coal
 */

#include "hybridcoal.hpp"
//#include"export.hpp"

//#include"all_gene_topo.hpp"
//#include"coal.hpp"

using namespace std;

//bool log_file_bool; /*! \todo fix log file */ 


//bool acc_bool;

//bool debug_bool;
//bool export_debug_bool;
//bool coal_debug_bool;
//bool utility_debug_bool;
//bool rm_debug_bool; 

/*! 
 * \fn int main(int argc, char *argv[])
 * \brief Main function for hybrid_coal 
 * */
int main(int argc, char *argv[]){

	if ( argc == 1 ) print_help(); 	//else, proceed
    
    try {
        HybridCoal run_hybridcoal ( argc, argv );
    }
    catch (const exception &e) {
      std::cerr << "Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }
}

