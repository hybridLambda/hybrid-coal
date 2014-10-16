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

#include"node.hpp"
#include"utility.hpp"

#ifndef NET
#define NET


/*! \brief Network class*/
class Net{	
	private:
	//string checking_labeled(string in_str);
	//string label_interior_node(string in_str);
		//int enumerate_internal_branch(Node *current,int e_num_old);
		void is_net_func(); /*!< \brief To determin if a Net is network or not. \return is_net */
		void checking_Parenthesis(string in_str);
		void create_node_list();
		void merge_repeated_nodes_as_one();
		void linking_nodes();
		bool is_net_;
		int max_rank_;
		void find_tips();
		//void find_descndnt();
		void find_interior_descndnt();

		
		size_t num_tips_;
	public:
		size_t num_tips() const {return this->num_tips_;};
		void set_num_tips(size_t num){this->num_tips_ =  num;};
		
		void init();
		string net_str; /*!< \brief species network string \todo this is new!!!*/
		
		//class Network_s SubNetworkS; /*!< \brief sub species networks \todo this is new!!!*/
		
		vector< valarray <int> > descndnt;
		//vector< valarray <int> > descndnt2;
		vector<string> tax_name;
		vector<string> tip_name;
		
		//vector<string> extract_tip_names();
		void extract_tip_names();
		void extract_tax_names();
		//vector<string> extract_tax_names();
		
		//vector <Node*> Net_nodes_ptr; /*!< \brief pointers to the nodes of Net \todo use this instead of Net_nodes */
		vector <Node*> Net_nodes;  /*!< \brief vector of nodes */
		
		bool is_net() const {return this->is_net_;}; /*!< \brief true if Net is a network; false if it's a tree */
		void set_is_net(bool yesorno){this->is_net_ = yesorno;}
		
		int max_rank() const {return this->max_rank_;}
		void set_max_rank(int rank){this->max_rank_ = rank;}
		
		//bool is_ultrametric; /*!< \brief true if the distances between tips and root are equal; false, otherwise */
		void clear(); 
		bool print_all_node();
		
		//bool is_ultrametric_func(); /*!< \brief To determin if a Net is ultrametric or not. \return is_ultrametric */
	
		//Net (){
			//string net_str;
			//vector <string> tax_name;
			//vector <string> tip_name;
			//vector <Node> Net_nodes;
			////vector< valarray <int> > descndnt;
			////vector< valarray <int> > descndnt2;
			////vector <Node*> Net_nodes_ptr;
		//}
		Net();
		~Net(){
			net_str.clear();
			for (size_t i = 0; i < Net_nodes.size(); i++){delete Net_nodes[i];}
			//tax_name.clear();
			//tip_name.clear();
			//Net_nodes.clear();
		}
		
		Net(string Net_str);
};

string rm_one_child_root(string in_str);




//string write_para_into_tree(string sp_string, double para);

#endif
