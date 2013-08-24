/* 
 * hybrid_coal is used to compute gene tree probabilities given species network under coalescent process.
 * 
 * hybrid_sim is used to simulate gene trees given species network under 
 * coalescent process.
 * 
 * Copyright (C) 2011 Sha (Joe) Zhu
 * 
 * This file is part of hybrid_sim and hybrid_coal
 * 
 * hybrid_sim is free software: you can redistribute it and/or modify
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


/*! \file utility.hpp
 * \brief Header file for network.cpp */


#include<iostream>
#include<string>
#include<sstream>
#include<fstream>
#include<vector>
#include<iomanip>
#include<valarray>
//#include <gmp.h>
#include"node.hpp"
//#include"net.hpp"
using namespace std;

#ifndef GLOBAL_H
#define GLOBAL_H	
//extern 
//extern bool log_file_bool;




bool start_of_tax_name(string in_str,size_t i);

void add_node(Node *parent_node, Node *child_node);
void find_tip(Node *current);
void find_hybrid_descndnt(Node *current);
bool find_descndnt(Node* current, string taxname);
bool find_descndnt2(Node* current, string taxname);
void rewrite_node_content(vector <Node*> Net_ptr);
//string construct_adding_new_Net_str(Net old_Net);
int ranking(Node *current);
//vector <string> all_n_tax_gene_tree(unsigned int tax_num);
double factorial(double a);	
double n_permu_a(double n, double a);
double n_choose_k(double n, double k);
double n_choose_k_old(double n, double k);

int factorial_int(int a);	
int n_permu_a_int(int n, int a);
int n_choose_k_int(int n, int k);
void appending_debug_file(string debug_file_input);
void appending_log_file(string log_file_input);
string remove_interior_label(string in_str);
string remove_brchlen(string in_str);
string tree_topo(string in_str);




string rm_and_hash_sign(string in_str);
string rm_and_sign(string in_str);
string rm_hash_sign(string in_str);
string rm_dash_sign(string in_str);


void check_and_remove(const char* file_name);

int my_exit();


size_t Parenthesis_balance_index_backwards(string in_str,size_t i);
size_t Parenthesis_balance_index_forwards(string in_str,size_t i);

//string extract_label(string in_str, size_t i);
size_t end_of_label_or_bl(string in_str, size_t i);

string write_para_into_tree(string sp_string, double para);


string extract_label(string in_str, size_t i);
size_t hybrid_hash_index(string in_str);
string extract_hybrid_label(string in_str);
string extract_hybrid_para_str(string in_str);
double extract_hybrid_para(string in_str);

string read_input_line(char inchar[]);
vector <string> read_input_lines(char inchar[]);
string read_input_para(char inchar[],string in_str);
bool is_num(char inchar[]);
string checking_labeled(string in_str);
string label_interior_node(string in_str);

string rewrite_net_str(string net_str, vector <string> tax_name, vector <string> tip_name);




void check_sp_str_dash_sign(string sp_str);
void check_gt_str_and_sign(string gt_str);


#endif
