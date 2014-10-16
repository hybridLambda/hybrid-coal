/* 
 * hybrid_coal is used to compute gene tree probabilities given species network under coalescent process.
 * 
 * Copyright (C) 2010 -- 2014 Sha (Joe) Zhu
 * 
 * This file is part of hybrid_coal
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
        ~HybridCoal();
        
        // ACTION
        void HybridCoal_core ( );
    private:
            /*! Members */ 
        int argc_;
        int argc_i;
        char * const* argv_;
        string tmp_input_str;

        //int num_sim_gt;
        bool print_tree_bool;
        bool plot_bool;
        //bool simulation_bool;
        //bool seg_bool;
        //bool tmrca_bool;
        //bool bl_bool;        
        //bool firstcoal_bool;        
        //bool freq_bool;
        //bool fst_bool;
        //size_t seed;	
        			
        string gt_file_name;
        //string mt_file_name;
        //vector <string> mt_tree_str_s;
        vector <string> gt_tree_str_s;
        //string seg_dir_name;
        //string extract_file_name;

       	//ofstream sim_gt_file_coal_unit;
        //ofstream sim_gt_file_mut_unit;
        //ofstream sim_gt_file_num_gener;
        //ofstream sim_gt_file_num_mut;
        //ofstream extract_file;
        bool read_GENE_trees;
        //bool read_mt_trees;        

		//vector <string> tax_name;
        //vector <double> monophyly;

        //action_board* simulation_jobs_;
        //SimulationParameters* parameters_;
        string prefix;

        /*! Methods */              
        //action_board* simulation_jobs() const { return this->simulation_jobs_; }
        //SimulationParameters* parameters() const { return this->parameters_;   }         
        
        //bool mono_fst_not_feasiable ( string flag );
        bool is_num ( const char *inchar );

        //string read_input_para ( const char *inchar, string in_str );
        string read_input_line ( const char *inchar );
        void init();
        void parse() ;
        void print();
        void read_input_lines(const char inchar[], vector <string> & out_vec);
        void read_sp_str( string & argv_i );
        //void read_sample_sizes();
        //void extract_mm_or_pop_param( string & mm_pop_string );
        //void create_new_site_data(string &gt_string_mut_num, int site_i);
        //void outtable_header( std::ofstream &output );
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
