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


/*! \file rm_hybrid.hpp
 * \brief Header file for rm_hybrid.cpp */

#include"coal.hpp"
//#include<algorithm>

//extern bool maple_bool;
//extern bool symbolic;
//extern bool log_file_bool;
#ifndef RM_HYBRID
#define RM_HYBRID
//extern 

//bool debug_bool=rm_debug_bool;

class vec_Net_wiz_prior_p{
	public:
	vector <Net_wiz_prior_p> Net_vec;
	void clear(){Net_vec.clear();};
	
	vec_Net_wiz_prior_p(){
		vector <Net_wiz_prior_p> Net_vec;
	}
};

class rm_H_node{
	private:
		//bool maple_bool_local;
		string current_removing_net_string;
		string current_removing_net_string_enum;
		string left_hybrid_parameter;
		string right_hybrid_parameter;
		string left_hybrid_parameter_str;
		string right_hybrid_parameter_str;		
		string hybrid_parameter;
		string new_node_name; //remove this
		string prior_omega_string;
		string current_omega_string;
		class vec_Net_wiz_prior_p nchild_gt_one(int rm_node_index,bool maple_bool_local,int n_child);
		class vec_Net_wiz_prior_p nchild_eq_one(int rm_node_index,bool maple_bool_local);
		unsigned int new_node_name_i;
		//int rm_node_index; // remove this
		int current_root_enum;
		//int n_child; //remove this
		vector < vector <int> > current_lambda_sum;
		vector <int> current_prior_coal_hist;
		vector < valarray < int > > current_prior_coal_list;
		double left_hybrid_parameter_num;
		double right_hybrid_parameter_num;
		double current_omega;
		double prior_omega;		
		//vector < vector <int> > current_e_num_vec;
	public:
		class vec_Net_wiz_prior_p block_rm_H;
		rm_H_node(int rm_node_index_in,Net_wiz_prior_p new_Net_wiz_prior_p,bool maple_bool_local_in);
};

string nchild_eq_one_core(string in_str,int rm_node_index,size_t i_one_hybrid_child);






vector < valarray <int> > all_possible_comb(int n);
vector < vector < valarray < int > > > build_h_child(int n);
vector < vector < valarray < int > > > build_s_child(int n);
vector < valarray <int> > rearrange_A(vector < valarray <int> > A,int n);

string rm_one_child_interior_node(string old_net_string);
//vec_Net_wiz_prior_p simplify_Networks(string old_net_string, string gt_string);

void simplify_Networks_maple_block(const char* file_name,string old_net_string,string gt_tree_str,int gt_str_i);
void write_maple_block(const char* file_name,string sp_str,vector <string> maple_block_vec_str);
//void rm_zero_kids_hybrid_node(vector <Node*> sp_nodes_ptr_rm);
//void simplify_Networks_maple(string old_net_string, vector <string> gt_string);




Net initialize_enum_net(string in_str);
vector < vector <int> > initialize_lambda_sum(Net enum_net);
//Net_wiz_prior_p build_initial_Networks(string in_str);
vector <string> build_maple_block_gt_part(vector <string> maple_block_vec_str,vec_Net_wiz_prior_p my_rmd_networks,string gt_i,string gt_tree_str);
vector <string> build_maple_block_sp_part(vector <string> maple_block_vec_str, int old_my_rmd_network_size,int new_block_size,string gt_i,string sp_i);
vector <string> build_maple_block_omega_part(vector <string> maple_block_vec_str,int old_my_rmd_network_size,int new_block_size,vec_Net_wiz_prior_p my_rmd_networks);

vec_Net_wiz_prior_p simplify_Networks_one_block(bool hybrid_node,int rm_node_index,vec_Net_wiz_prior_p my_rmd_networks,size_t i,bool maple_bool_local,string gt_tree_str);
vec_Net_wiz_prior_p build_initial_Networks(string in_str);
vec_Net_wiz_prior_p simplify_Networks(string old_net_string, string gt_string);
vec_Net_wiz_prior_p rm_S_node(int rm_node_index,string gt_str,Net_wiz_prior_p new_Net_wiz_prior_p,bool maple_bool_local_in);
bool check_sp_coal_valid(bool coal_ed,Net my_gt_tree, vector < valarray <int> > new_coal_clade);

double compute_rm_S_omega(int n_child,int k_clade,vector < valarray < int > > brand_new_current_prior_coal_clades, Net gt_tree);

#endif
