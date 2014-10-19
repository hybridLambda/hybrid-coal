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

#include "node.hpp"
#include <valarray>
#include <fstream>

#ifndef NETWORK
#define NETWORK

/*! \brief Network class*/

class Net;

class Tree{
    friend class Net;
    friend class HybridCoal;
    friend class Frequency;
    friend class Figure;
    private:        
        string label_interior_node(string in_str);
        void enumerate_internal_branch( Node &node );
        
        size_t first_coal_rank();
        size_t current_enum_;
        void init(){
            this->current_enum_ = 0;
            this->is_Net = false;
            this->is_ultrametric = true;
            }

        bool start_of_tax_name(string in_str, size_t i);
        size_t Parenthesis_balance_index_backwards( string &in_str, size_t i );
        size_t Parenthesis_balance_index_forwards( string &in_str, size_t i );
        
        void check_Parenthesis(string &in_str);
        void check_labeled( string in_str );
        
        void check_isNet(); /*!< \brief To determin if a Net is network or not. \return is_Net */
        void check_isUltrametric(); /*!< \brief To determin if a Net is ultrametric or not. \return is_ultrametric */

        size_t first_coal_index ();
    
        string rewrite_internal_node_content( size_t i);       
        void connect_graph();
        void extract_tax_and_tip_names();
        
        void init_descendant();
        void init_node_clade();
        void rewrite_descendant();
        void rewrite_node_clade();
        
        string net_str; /*!< \brief species network string \todo this is new!!!*/
        size_t max_rank;
        vector< valarray <int> > descndnt;
        vector< valarray <int> > samples_below;
        vector<string> tip_name;
        vector < vector < int > > R_matrix;
        vector < vector < int > > M_matrix;
        vector < vector <size_t> > coal_hist_mat;
        void build_coal_hist_mat( Net & sp_net );
        void building_R_matrix();
        void building_M_matrix( Net & sp_net ) ;
        bool is_Net_() const { return this->is_Net ; }
        string extract_label(string &in_str, size_t i);
        void print_all_node();
        bool print_all_node_dout();
        bool print_matrix( vector < vector < int > > & mat );
        //Tree (){ this->init(); }
        Tree(string Tree_str);
        ~Tree(){};

        vector <string> tax_name;
        bool is_ultrametric; /*!< \brief true if the distances between tips and root are equal; false, otherwise */
        bool is_Net; /*!< \brief true if Net is a network; false if it's a tree */
    public:    
        vector < Node > NodeContainer;  /*!< \brief vector of nodes */
        void rewrite_node_content();
        string print_newick( Node * node );
};

class Net: public Tree {
    friend class HybridCoal;
        vector < double > brchlens_vec;
        vector < int > max_num_brch_vec;
        vector < vector < vector < double > > > gijoemat;
        vector < vector < int > > S_matrix;
        void assign_bl_to_vec();
        void build_gijoe();
        bool print_gijoemat();
        void building_S_matrix();
        
    public:
        Net (string Net_str) : Tree ( Net_str){};
        //Net () { this->init();}
        ~Net(){};
    };

string remove_interior_label(string in_str);
size_t end_of_label_or_bl( string &in_str, size_t i);
void readNextStringto( string &readto , int& argc_i, int argc_, char * const* argv_ );

string rm_one_child_root(string in_str);

#endif
