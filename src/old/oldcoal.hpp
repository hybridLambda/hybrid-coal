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

/*! \file coal.hpp
 * \brief Header file for coal.cpp */


#include"utility.hpp"
#ifndef COAL
#define COAL
//extern bool latex_bool;
//extern bool maple_bool;
//extern bool debug_bool;
//extern 




/*! \brief species network with prior probabilities*/
class Net_wiz_prior_p{
	public:
	string s_net_string; /*!< species network string */
	string s_net_string_enum;
	int root_enum; //! \todo try remove this!, as enum of the root can be expressed in the expressions /
	vector < valarray < int > > prior_clade_list;
	vector <int> prior_coal_hist;
	vector < vector < int > > lambda_sum;
	//vector < vector <int> > e_num_vec;
//	vector <int> prior_number_coal_list;
	//string prior_prob_string;
	//double prior_prob_num;
	string omega_string;
	double omega;

	
	
	
	Net_wiz_prior_p(){
		vector < valarray < int > > prior_clade_list;
		vector <int> prior_coal_hist;
		vector < vector < int > > lambda_sum;
		//vector < vector <int> > e_num_vec;
		//prior_prob_num=1;
		//prior_prob_string="";
		omega=1;
		omega_string="";
		s_net_string="";
		root_enum=0;
	}
	
	void clear(){
		prior_clade_list.clear();
		prior_coal_hist.clear();
		lambda_sum.clear();
		//omega_string.clear();
		omega=1.0;
		omega_string="";
		s_net_string="";
		s_net_string_enum="";
		root_enum=0;
		//string s_net_string; /*!< species network string */
	//string s_net_string_enum;
	//int root_enum; //! \todo try remove this!, as enum of the root can be expressed in the expressions /
	//vector < valarray < int > > prior_clade_list;
	//vector <int> prior_coal_hist;
	//vector < vector < int > > lambda_sum;
	//string omega_string;
	//double omega;
		
	}
	
};




class gt_coal_in_st{
	private:
		vector < vector < vector < double > > > gijoe_matrix;
		vector < vector <int> > num_coal_in_branch;		
		void build_gijoe_matrix(vector < double > brchlens_vec, vector < int > max_num_brch_vec );
		void construct_prob_maple();
		void calc_prob_of_hist();
		Net my_sp_net;
		Net s_net_enum;
		//string prior_prob_str;
		string num_ranked_tree_str;
		string repeated_ranked_tree_str;
		double prior_prob;
		int num_ranked_tree;
		int repeated_ranked_tree;		
		string gt_str;
		string sp_str;
		
	public:
	
		vector < vector <int> > current_lambda_sum;
		vector < double > prob_of_hist;
		vector < vector <unsigned int> > coal_hist;
		vector < vector <int> > all_w;
		vector < vector <int> > all_d;
		vector < vector <int> > num_enter;
		vector < vector <int> > num_out;
		vector < string > prob_of_hist_maple;
		vector <int> prior_coal_hist;
		double probability;
		string prob_of_maple;
		bool sp_valid;
		void latex_print(const char* file_name);
		void print_coal_hist();
		//gt_coal_in_st new_gt_coal_in_st;
		//void print_S();
		//void print_M();
		//void print_R();
		~gt_coal_in_st();
		gt_coal_in_st();
		gt_coal_in_st(string gt_tree_str, Net_wiz_prior_p st_tree_wiz_prior_p);
		gt_coal_in_st(string gt_tree_str,string st_tree_str);			
		gt_coal_in_st(Net gt_tree, Net st_tree);
		gt_coal_in_st(Net gt_tree, Net_wiz_prior_p st_tree_wiz_prior_p);
		
		void clear(){
		current_lambda_sum.clear();
		prob_of_hist.clear();
		coal_hist.clear();
		all_w.clear();
		all_d.clear();
		num_enter.clear();
		num_out.clear();
		prob_of_hist_maple.clear();
		prior_coal_hist.clear();
		probability=0;
		prob_of_maple.clear();
		sp_valid=false;
		gijoe_matrix.clear();
		num_coal_in_branch.clear();
		my_sp_net.clear();
		s_net_enum.clear();
		
		num_ranked_tree_str="";
		repeated_ranked_tree_str="";
		prior_prob=1.0;
		num_ranked_tree=0;
		repeated_ranked_tree=0;		
		gt_str="";
		sp_str="";
		}
	
};






double gijoe(int u, int v, double T); /*!< \brief compute the probability of \a u lineages entering and \a v lineages exiting at a branch with length \a t. \todo move this function to rm_hybrid.hpp, in program coal, should be using gijoe matrix, would be faster.*/
void print_matrix(vector < vector < int > > mat);

bool check_multiple_sp_tip(Net pre_current_net);
bool check_gt_tip_already_coaled(Net my_gt_tree);


bool check_sp_valid(Net sp_net, Net gt_tree);

void print_prior_coal_hist(vector < valarray < int > > prior_clade_list);
#endif
