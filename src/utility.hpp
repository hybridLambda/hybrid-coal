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

using namespace std;

#ifndef GLOBAL_H
#define GLOBAL_H	
//extern 
//extern bool log_file_bool;



/*! 
 */
 


//AC::AC(){
	//pAC=0.0;
	//}





/*! \brief Node of a tree or network, it also represent the branch between this node and its parent node
 */
class Node {
	public:
	//vector<int> descndnt;
	vector<int> des_num_descndnt_interior;
	vector<Node*> descndnt_interior_node; /*!< \brief list of pointers to its descndent interior nodes */
	vector<Node*> child; /*!< \brief list of pointers to its child nodes */	
	Node* parent1; /*!< \brief pointer to its parent node. */
	
	string clade; /*!< \brief clade at this node, \todo this should be modified to a vector <string> */
	
	string label; /*!< \brief String label of a node, each node has unique label */
	string node_content; /*!< \brief node content, the subtree string at this node */
	class Net * SubNetwork; /*!< \brief \todo pointer to the subtree of this node */
	unsigned int e_num; /*!< \brief numbering the branch */
	int rank; /*!< \brief rank of the node, tip node has rank one, the root has the highest rank */
	int num_child; /*!< \brief number of child \todo this can be replaced by child.size */
	int num_descndnt; /*!< \brief number of the tip nodes, that are descendant from this node */
	int num_descndnt_interior; /*!< \brief number of the interior nodes, that are descendant from this node \todo to be replaced by descndnt_interior_node.size()? */
	vector <double> path_time; 
	double absolute_time; /*!< \brief distance to the bottom of the tree */
	double brchlen1; /*!< \brief Branch length */
	bool visited;
	bool descndnt_of_hybrid; /*!< \brief Indicator of descendant of hybrid nodes. It's true, if it is a descendant of hybrid nodes; false, otherwise. */
	bool tip_bool; /*!< \brief Indicator of tip nodes. It's true, if it is a tip node, otherwise it is false. */
	
	unsigned int node_index; /*!< \brief node index in the array, \todo use this more often!!!*/
	
	
	/* These members apply to only hybrid nodes */
	bool hybrid; /*!< \brief Hybrid node only, indicator of a hybrid node */
	Node* parent2; /*!< \brief Hybrid node only, pointer to its second parent node. */
	unsigned int e_num2; /*!< \brief Hybrid node only, numbering the branch between the node and its second parent */
	//double prob_to_hybrid_left; /*!< \brief Hybrid node only, the probability that a lineage goes to the left */
	double brchlen2;/*!< \brief Hybrid node only, Branch length to the second parent*/

	string name; /*!< \brief Name of a node, this is not unique for nodes. e.g. if its label is A_1, name is A */
	
	//vector < class AC > AC_list;
		//vector < class CAC > CAC_list;

	//vector < unsigned int > AC_list_sizes;
	//vector < vector <unsigned int> > AC_list;
	//double pAC;
	
	//vector <unsigned int> Net_node_contains_gt_node1; /*!< Used while simulation, check if a Network node contains a gene tree node */
	//vector <unsigned int> Net_node_contains_gt_node2; /*!< Used while simulation, check if a Network node contains a gene tree node */
	
	Node(); /*!< \brief Initialize Node class*/
	void print_net_Node();
	void print_tree_Node();
	void clear();
	
		
};

/*! \brief Collection of Networks in a vector container*/
class Network_s{
	public:
	vector <class Net *> Net_vec;
	void clear(){Net_vec.clear();};
	
	Network_s(){
		vector <class Net *> Net_vec;
	}
	
};

/*! \brief Network class*/
class Net{	
	private:
	//string checking_labeled(string in_str);
	//string label_interior_node(string in_str);
	int enumerate_internal_branch(Node *current,int e_num_old);
	bool is_net_func(); /*!< \brief To determin if a Net is network or not. \return is_net */
	
	public:
	string net_str; /*!< \brief species network string \todo this is new!!!*/
	
	class Network_s SubNetworkS; /*!< \brief sub species networks \todo this is new!!!*/
	int max_rank;
	vector< valarray <int> > descndnt;
	vector< valarray <int> > descndnt2;
	vector<string> tax_name;
	vector<string> tip_name;
	//vector <Node*> Net_nodes_ptr; /*!< \brief pointers to the nodes of Net \todo use this instead of Net_nodes */
	vector <Node> Net_nodes;  /*!< \brief vector of nodes */
	bool is_net; /*!< \brief true if Net is a network; false if it's a tree */
	//bool is_ultrametric; /*!< \brief true if the distances between tips and root are equal; false, otherwise */
	void clear(); 
	void print_all_node();
	bool is_ultrametric_func(); /*!< \brief To determin if a Net is ultrametric or not. \return is_ultrametric */

	Net (){
		string net_str;
		vector <string> tax_name;
		vector <string> tip_name;
		vector <Node> Net_nodes;
		vector< valarray <int> > descndnt;
		vector< valarray <int> > descndnt2;
		//vector <Node*> Net_nodes_ptr;
	}
	
	//~Net(){
		//tax_name.clear();
		//tip_name.clear();
		//Net_nodes.clear();
	//}
	
	Net(string Net_str);
};

bool start_of_tax_name(string in_str,size_t i);

void add_node(Node *parent_node, Node *child_node);
void find_tip(Node *current);
void find_hybrid_descndnt(Node *current);
bool find_descndnt(Node* current, string taxname);
bool find_descndnt2(Node* current, string taxname);
void rewrite_node_content(vector <Node*> Net_ptr);
string construct_adding_new_Net_str(Net old_Net);
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
void checking_Parenthesis(string in_str);

string extract_label(string in_str, size_t i);
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
