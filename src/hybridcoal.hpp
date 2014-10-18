/* 
 * hybrid-coal is used to compute gene tree probabilities given species network under coalescent process.
 * 
 * Copyright (C) 2010 -- 2014 Sha (Joe) Zhu
 * 
 * This file is part of hybrid-coal
 * 
 * hybrid-coal is free software: you can redistribute it and/or modify
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
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <net.hpp>
#include <sstream>      // std::stringstream

#ifndef HYBRDRIDCOAL_PARAM_INCLUDED
#define HYBRDRIDCOAL_PARAM_INCLUDED
using namespace std;

void print_example();
void print_help();
void print_option();

class HybridCoal{
    public:	
        /*! Constructors and Destructors */  
        HybridCoal(int argc, char *argv[]) : argc_(argc), argv_(argv) { this->init(); this->parse(); }        
        ~HybridCoal(){ dout << "~HybridCoal()" <<endl;};
        
        // ACTION
        void HybridCoal_core ( );
    private:
            /*! Members */ 
        int argc_;
        int argc_i;
        char * const* argv_;
        string tmp_input_str;
        string sp_str;

        bool print_tree_bool;
        bool plot_bool;
        bool symb_bool;        			
        bool latex_bool;
        bool print_gene_topo_bool;
        bool enumerate_gt_bool;
        bool maple_bool;
        bool all_gt_tree_bool ;
        string gt_file_name;
        vector <string> gt_tree_str_s;
        bool read_GENE_trees;
        string prefix;

        /*! Methods */              
        void check_gt_str_and_sign( string &gt_str );

        bool is_num ( const char *inchar );

        string read_input_line ( const char *inchar );
        //void read_gt();
        void init();
        void parse() ;
        void print();
        void read_input_lines(const char inchar[], vector <string> & out_vec);
        void read_sp_str( string & argv_i );
        void finalize();

        template < class T > T readNextInput() {
            ++argc_i;        
            if (argc_i >= argc_) throw std::invalid_argument( std::string( "Not enough parameters when parsing options: ") + argv_[argc_i-1]);
        
            char c;
            T input;
            std::stringstream ss(argv_[argc_i]);
            ss >> input;
            if (ss.fail() || ss.get(c)) throw std::invalid_argument( std::string( "Failed to parse option: ") + argv_[argc_i]); 
            return input;
        }

};
#endif
