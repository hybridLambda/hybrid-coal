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

/*! \file export.hpp
 * \brief Header file for export.cpp */

#include"utility.hpp"
#include"rm_hybrid.hpp"

#ifndef EXPORT
#define EXPORT
//extern 


void print_all_gt_topo(const char* file_name, vector <string> gt_tree_str_s);

string Fname_ext(const char* file_name, const char* ext);
string Fname_no_ext(const char* file_name, const char* ext);

void plot_in_latex(const char* file_name, Net net_dummy, int plot_option);
void plot_in_latex_file(const char* file_name, Net net_dummy, int plot_option);

void latex_header(const char* file_name);
void latex_nt_body1(const char* file_name,size_t topo_idx,string topo);
//void latex_nt_body2(const char* file_name,vec_Net_wiz_prior_p net_str_s,size_t i);
void latex_tre_body(const char* file_name,string gt_str,string sp_str);
void latex_tail(const char* file_name);

void plot_in_dot(const char* file_name, Net net_dummy, int plot_option);

void produce_maple(const char* file_name, string in_str, bool symb_bool, vector <string> gt_tree_str_s);
void maple_head(const char* file_name, string in_str,bool symb_bool);
void symb_maple(const char* file_name, string in_str);
void maple_tail(const char* file_name, size_t num_gt);

//void plot_enum_in_latex(const char* file_name, Net_wiz_prior_p net_dummy);

void prob_out_body(const char* file_name,size_t topo_idx,string topo,double prob,size_t num_tax);
void prob_out_tail(const char* file_name,size_t empty_space,double total_prob);

void list_sub_tree(string in_str,string gt_str);
void list_sub_net(string in_str,string gt_str);
void list_sub(string in_str,string gt_str);

int set_plot_option(bool plot_label,bool plot_branch);
valarray <int>  det_x_node (Net net_dummy);
void print_help();


#endif
