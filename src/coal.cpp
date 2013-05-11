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


/*! \file coal.cpp
 * \brief Remake of program COAL, to compute gene tree probabilities given species tree */


#include"coal.hpp"
bool coal_debug_bool=false;	






vector < vector <unsigned int> > build_coal_hist(
vector < vector <unsigned  int> > coal_hist,
unsigned int node_i,
vector < vector <unsigned int> > coal_hist_mat, 
vector < vector < int > > R_matrix)
{
	vector < vector <unsigned  int> > new_coal_hist;
	for (unsigned int coal_hist_i=0;coal_hist_i<coal_hist.size();coal_hist_i++){
		vector <unsigned int> coal_hist_dummy;
		coal_hist_dummy=coal_hist[coal_hist_i];
		unsigned int coal_hist_mat_node_i_i=0;
		for (unsigned int j_R_mat=0;j_R_mat<node_i;j_R_mat++){
			while (R_matrix[node_i][j_R_mat]==1 && coal_hist_mat[node_i][coal_hist_mat_node_i_i]< coal_hist_dummy[j_R_mat]){
				coal_hist_mat_node_i_i++;
			}
		}
		for (;coal_hist_mat_node_i_i<coal_hist_mat[node_i].size();coal_hist_mat_node_i_i++ ){
			vector <unsigned int> coal_hist_dummy_dummy;
			coal_hist_dummy_dummy=coal_hist_dummy;
			coal_hist_dummy_dummy.push_back(coal_hist_mat[node_i][coal_hist_mat_node_i_i]);
			new_coal_hist.push_back(coal_hist_dummy_dummy);
		}
	}
	
	if (coal_debug_bool){
		cout<<new_coal_hist.size()<<"new_coal_hist"<<endl;	
		for (unsigned int i=0;i<new_coal_hist.size();i++){
			for (unsigned int j=0;j<new_coal_hist[i].size();j++){
				cout<<new_coal_hist[i][j];
			}
			cout<<endl;
		}
	}
	
	if (node_i<coal_hist_mat.size()-1){
		node_i++;
		new_coal_hist=build_coal_hist(new_coal_hist,node_i,coal_hist_mat,R_matrix);		
	}

	return new_coal_hist;
}

void print_matrix(vector < vector < int > > mat){
	for (size_t i=0;i<mat.size();i++){
		for (size_t j=0;j<mat[i].size();j++){
			cout<<mat[i][j];
		}
		cout<<endl;
	}
	cout<<endl;		
}


vector < vector < int > > building_R_matrix(Net gt_tree){
	int gt_max_enum=gt_tree.Net_nodes.back().e_num;
	vector < vector < int > > R_matrix;
	for (int i_gt_enum=0;i_gt_enum<gt_max_enum;i_gt_enum++){
		vector <int> R_matrix_row;
		for (int j_gt_enum=0;j_gt_enum<i_gt_enum;j_gt_enum++){
			R_matrix_row.push_back(1);
		}
		R_matrix.push_back(R_matrix_row);
	}
	for (unsigned int i=0;i<gt_tree.Net_nodes.size()-1;i++){
		if (!gt_tree.Net_nodes[i].tip_bool){
			for (unsigned int j=0;j<i;j++){
				if (!gt_tree.Net_nodes[j].tip_bool){
					for (unsigned int i_des=0;i_des<gt_tree.descndnt[i].size();i_des++){
						int des_diff=gt_tree.descndnt[i][i_des]-gt_tree.descndnt[j][i_des];
						if (des_diff<0){
							R_matrix[gt_tree.Net_nodes[i].e_num-1][gt_tree.Net_nodes[j].e_num-1]=0;
							break;
						}
					}
				}
			}
		}
	}
	if (coal_debug_bool){
		cout<<"R matrix"<<endl;
		print_matrix(R_matrix);
	}
	return R_matrix;
}

vector < vector < int > > building_S_matrix(Net my_sp_net){
	int sp_max_enum=my_sp_net.Net_nodes.back().e_num;
	vector < vector < int > > S_matrix;
	
	for (int i_sp_enum=0;i_sp_enum<sp_max_enum;i_sp_enum++){
		vector < int > S_matrix_row(sp_max_enum,1);
		S_matrix_row[i_sp_enum]=0;
		S_matrix.push_back(S_matrix_row);
	}
	for (unsigned int i=0;i<my_sp_net.Net_nodes.size()-1;i++){
		if (!my_sp_net.Net_nodes[i].tip_bool){
			for (unsigned int j=0;j<my_sp_net.Net_nodes.size()-1;j++){
				if (!my_sp_net.Net_nodes[j].tip_bool){
					for (unsigned int i_des=0;i_des<my_sp_net.descndnt[i].size();i_des++){
						int des_diff=my_sp_net.descndnt[i][i_des]-my_sp_net.descndnt[j][i_des];
						if (des_diff<0){
							S_matrix[my_sp_net.Net_nodes[i].e_num-1][my_sp_net.Net_nodes[j].e_num-1]=0;
							break;
						}
					}
				}
			}
		}
		
	}
	if (coal_debug_bool){
		cout<<"S matrix"<<endl;
		print_matrix(S_matrix);
	}
	return S_matrix;
}
	
vector < vector < int > > building_M_matrix(
Net my_gt_tree,
Net my_sp_net)
{
	int sp_max_enum=my_sp_net.Net_nodes.back().e_num;
	int gt_max_enum=my_gt_tree.Net_nodes.back().e_num;
	vector < vector < int > > M_matrix;
	
	for (int i_sp_enum=0;i_sp_enum<sp_max_enum;i_sp_enum++){
		vector < int > M_matrix_row((gt_max_enum-1),1);
		M_matrix.push_back(M_matrix_row);
	
	}
	for (unsigned int i=0;i<my_sp_net.Net_nodes.size()-1;i++){
		if (!my_sp_net.Net_nodes[i].tip_bool){
			for (unsigned int j=0;j<my_gt_tree.Net_nodes.size()-1;j++){
				if (!my_gt_tree.Net_nodes[j].tip_bool){
					int des_diff_prod=1;
					for (unsigned int i_des=0;i_des<my_gt_tree.tax_name.size();i_des++){
						int des_diff=my_sp_net.descndnt[i][i_des]-my_gt_tree.descndnt[j][i_des]; /*! \todo this may not be right for multiple lineages per species*/
						if (des_diff<0){
							des_diff_prod=0;
							M_matrix[my_sp_net.Net_nodes[i].e_num-1][my_gt_tree.Net_nodes[j].e_num-1]=0;
							break;
						}
						else{  /*! \todo check this else!!! dont think this is needed */
							if (des_diff>0){
								des_diff_prod=0;}
						}
					}
					if (des_diff_prod==1){
						M_matrix[my_sp_net.Net_nodes[i].e_num-1][my_gt_tree.Net_nodes[j].e_num-1]=1;
					}
				}
			}
		}
	}
	
	if (coal_debug_bool){
		cout<<"M matrix"<<endl;
		print_matrix(M_matrix);
	}
	return M_matrix;
}

gt_coal_in_st::gt_coal_in_st(){	
	string num_ranked_tree_str;
	vector < vector <unsigned int> > coal_hist;
	vector < vector <int> > num_enter;
	vector < vector <int> > num_out;
	vector < vector <int> > num_coal_in_branch;		
	vector < vector <int> > all_w;
	vector < vector <int> > all_d;
	vector < double > prob_of_hist;
	vector < string > prob_of_hist_maple;
	string prob_of_maple;
	probability=0;
	sp_valid=false;
}
		
		
void gt_coal_in_st::construct_prob_maple(){
	for (unsigned int i_coal_hist=0;i_coal_hist<coal_hist.size();i_coal_hist++){
		string current_prob_maple;
		for (unsigned int i=0;i<all_w[i_coal_hist].size();i++){
			ostringstream current_branch,u_str,v_str,w_str,d_str;
			int lambda_subscript;
			for (unsigned int node_i=0; node_i<s_net_enum.Net_nodes.size();node_i++){
				if (s_net_enum.Net_nodes[node_i].e_num==(i+1)){			
					lambda_subscript=s_net_enum.Net_nodes[node_i].brchlen1;
					current_branch<<lambda_subscript;
					break;
				}
			}
			//ostringstream u_str;
			u_str<<num_enter[i_coal_hist][i];
			v_str<<num_out[i_coal_hist][i];
			w_str<<all_w[i_coal_hist][i];
			d_str<<all_d[i_coal_hist][i];
			if (i<all_w[i_coal_hist].size()-1){
				current_prob_maple=current_prob_maple+w_str.str()+"/"+d_str.str()+"*puvT("+u_str.str()+","+v_str.str()+",lambda["+current_branch.str()+"]";
				for (unsigned int lambda_sum_i=0; lambda_sum_i<current_lambda_sum[lambda_subscript-1].size();lambda_sum_i++){
					ostringstream current_branch_new;
					current_branch_new<<current_lambda_sum[lambda_subscript-1][lambda_sum_i];
					current_prob_maple=current_prob_maple+"+lambda["+current_branch_new.str()+"]";
				}
					//current_prob_maple=current_prob_maple+")";
				current_prob_maple=current_prob_maple+")*";
			}
			else{
				current_prob_maple=current_prob_maple+w_str.str()+"/"+d_str.str();
			}
		}
		prob_of_hist_maple.push_back(current_prob_maple);
		if (prob_of_maple.empty()){
			prob_of_maple=current_prob_maple;
		}
		else{
			prob_of_maple=prob_of_maple+"+"+current_prob_maple;
		}
	}
	prob_of_maple="("+prob_of_maple+");";
}


void gt_coal_in_st::calc_prob_of_hist(){
	probability=0.0;
	for (unsigned int i_coal_hist=0;i_coal_hist<coal_hist.size();i_coal_hist++){
		double current_prob_of_hist=1;
		
		for (size_t i=0;i<all_w[i_coal_hist].size();i++){
			double current_branch_lengths;
			size_t current_enum;
			for (size_t node_i=0; node_i<my_sp_net.Net_nodes.size();node_i++){
				if (my_sp_net.Net_nodes[node_i].e_num==(i+1)){
					current_branch_lengths=my_sp_net.Net_nodes[node_i].brchlen1;
					current_enum=my_sp_net.Net_nodes[node_i].e_num;
					break;
				}
			}	
			double current_gijoe;
			if ((i+1)==my_sp_net.Net_nodes.back().e_num){
				current_gijoe=1;
				}
				else{
				current_gijoe=gijoe_matrix[current_enum-1][num_enter[i_coal_hist][i]-1][num_out[i_coal_hist][i]-1];
//!!! this works!!! //current_gijoe=gijoe(num_enter[i_coal_hist][i], num_out[i_coal_hist][i], current_branch_lengths); 
					
			}
			//double current_gijoe=gijoe_matrix[num_enter[i_coal_hist][i]-1][num_out[i_coal_hist][i]-1];
			//cout<<"g["<<num_enter[i_coal_hist][i]<<","<< num_out[i_coal_hist][i]<<"]("<<current_branch_lengths<<")="<<current_gijoe<<"  ";
			//cout<<"gij_matrix:  "<<gijoe_matrix[current_enum-1][num_enter[i_coal_hist][i]-1][num_out[i_coal_hist][i]-1] <<endl;
			
			//if (current_gijoe==gijoe_matrix[current_enum-1][num_enter[i_coal_hist][i]-1][num_out[i_coal_hist][i]-1]){
				////cout<<"yes equal!!"<<endl;
			//}
			//else{
				//cout<<"No they are not equal!"<<endl;
				//cout<<"g["<<num_enter[i_coal_hist][i]<<","<< num_out[i_coal_hist][i]<<"]("<<current_branch_lengths<<")="<<current_gijoe<<"  ";
				//cout<<"gij_matrix:  "<<gijoe_matrix[current_enum-1][num_enter[i_coal_hist][i]-1][num_out[i_coal_hist][i]-1] <<endl;
			//}
			
			current_prob_of_hist=current_prob_of_hist*all_w[i_coal_hist][i]/all_d[i_coal_hist][i]*current_gijoe;
		}
		prob_of_hist.push_back(current_prob_of_hist);
		//cout<<endl<<current_prob_of_hist<<endl;
		probability=probability+current_prob_of_hist;
	}
	//probability=probability/num_ranked_tree*repeated_ranked_tree;
}



void gt_coal_in_st::latex_print(const char* file_name){
	if (sp_valid ){
		ofstream latex_file;
		latex_file.open (file_name, ios::out | ios::app | ios::binary);
			//cout<<"$(";
		latex_file<<"gt: {\\verb "<<gt_str<<" } \\\\"<<endl;
		latex_file<<"sp: {\\verb "<<sp_str<<" } \\\\"<<endl;
		
		latex_file<<"\\begin{tabular}{|ccc|}\\hline\n";
		
		for (unsigned int i_coal_hist=0;i_coal_hist<coal_hist.size();i_coal_hist++){	
			//cout<<"hea"<<endl;
			for (unsigned int prior_coal_hist_i=0;prior_coal_hist_i<prior_coal_hist.size();prior_coal_hist_i++){
				latex_file<<prior_coal_hist[prior_coal_hist_i]<<",";
			}
			for (unsigned int j_coal_hist=0;j_coal_hist<coal_hist[i_coal_hist].size();j_coal_hist++){
				if (coal_debug_bool){
					cout<<i_coal_hist<<" "<<j_coal_hist<<endl;
					cout<<coal_hist[i_coal_hist][j_coal_hist]<<endl;
				}
				latex_file<<coal_hist[i_coal_hist][j_coal_hist];
				if (j_coal_hist<coal_hist[i_coal_hist].size()-1){
					latex_file<<",";
					}
			}
			
			latex_file<<" & ";
			//cout<<"  ";
			latex_file<<prob_of_hist[i_coal_hist]<<" & "<<endl;
			//cout<<prob_of_hist[i_coal_hist]<<endl;
			latex_file<<"$(";
			for (size_t i=0;i<all_w[i_coal_hist].size();i++){
				//if (all_d[i_coal_hist][i]!=0){
				//if (all_w[i_coal_hist][i]/all_d[i_coal_hist][i]!=1){
				latex_file<<"\\frac{"<<all_w[i_coal_hist][i]<<"}{"<<all_d[i_coal_hist][i]<<"}";
					//cout<<all_w[i_coal_hist][i]<<"/"<<all_d[i_coal_hist][i];
				//}
				//}
				ostringstream current_branch;
				double lambda_subscript;
				if (coal_debug_bool){
					s_net_enum.print_all_node();
				}
				for (unsigned int node_i=0; node_i<s_net_enum.Net_nodes.size();node_i++){
					if (s_net_enum.Net_nodes[node_i].e_num==(i+1)){
						lambda_subscript=s_net_enum.Net_nodes[node_i].brchlen1;
						current_branch<<lambda_subscript;
		//cout<<"hea "<<i <<"  "<<s_net_enum.Net_nodes[node_i].brchlen1<<"  "<<current_branch.str() <<endl;	
						break;
		
					}
				}
		//cout<<"hea "<<i <<"  "<<"  "<<current_branch.str() <<endl;	
			
				if (i<all_w[i_coal_hist].size()-1){
					latex_file<<"p_{"<<num_enter[i_coal_hist][i]<<num_out[i_coal_hist][i]<<"}(\\lambda_{"<<current_branch.str()<<"}";
//if (){					
					//for (size_t lambda_sum_i=0; lambda_sum_i<current_lambda_sum[lambda_subscript-1].size();lambda_sum_i++){
						//ostringstream current_branch_new;
						//current_branch_new<<current_lambda_sum[lambda_subscript-1][lambda_sum_i];
						//latex_file<<"+\\lambda_{"+current_branch_new.str()+"}";
					//}
		//}
					latex_file<<")";//need to change, coalesced branch	
					//cout<<"p_{"<<num_enter[i_coal_hist][i]<<num_out[i_coal_hist][i]<<"}(lambda_{"<<i+1<<"})";	
				}
			
			}
			//cout<<")/"<< num_ranked_tree_str<<"*"<<repeated_ranked_tree_str<<endl;
			latex_file<<")$\\\\"<<endl;
			//cout<<")$"<<endl;
		}	
		latex_file<<"&"<<probability<<"& \\\\ \\hline\\end{tabular}\\\\ "<<endl;
		latex_file.close();
	}
	if (coal_debug_bool){
		cout<<"gt_coal_in_st::latex_print END"<<endl;
	}	
}

void gt_coal_in_st::print_coal_hist(){
	if (sp_valid){
		cout<<coal_hist.size()<<" Coalescent History"<<endl;
		for (unsigned int i_coal_hist=0;i_coal_hist<coal_hist.size();i_coal_hist++){	
			for (unsigned int j_coal_hist=0;j_coal_hist<coal_hist[i_coal_hist].size();j_coal_hist++){
				cout<<coal_hist[i_coal_hist][j_coal_hist];
				if (j_coal_hist<coal_hist[i_coal_hist].size()-1){
					cout<<",";
				}
			}
			cout<<endl;
		}
	}
}

gt_coal_in_st::~gt_coal_in_st(){
	coal_hist.clear();
	num_enter.clear();
	num_out.clear();
	num_coal_in_branch.clear();
	all_w.clear();
	all_d.clear();
}	





gt_coal_in_st::gt_coal_in_st(string gt_tree_str, Net_wiz_prior_p st_tree_wiz_prior_p){
	if (coal_debug_bool){
		cout<<"gt_coal_in_st::gt_coal_in_st(string gt_tree_str, Net_wiz_prior_p st_tree_wiz_prior_p) START"<<endl;
	}
	Net my_gt_tree(gt_tree_str);
	//gt_coal_in_st new_gt_coal_in_st(my_gt_tree, st_tree_wiz_prior_p);
	
		//Net my_gt_tree=gt_tree;
	Net load_Net(st_tree_wiz_prior_p.s_net_string);
	//cout<<st_tree_wiz_prior_p.s_net_string<<endl;
	my_sp_net=load_Net;
	//check if my_sp_net is valid
	sp_valid=true;
	prior_coal_hist=st_tree_wiz_prior_p.prior_coal_hist;
	
	//cout<<st_tree_wiz_prior_p.prior_clade_list.size()<<endl;
	//print_prior_coal_list(st_tree_wiz_prior_p.prior_clade_list);
	//gt_tree.print_all_node();
	
	for (size_t coal_clade_i=0;coal_clade_i<st_tree_wiz_prior_p.prior_clade_list.size();){
		bool coal_i_valid=false;
		for (size_t gt_node_i=0;gt_node_i<my_gt_tree.Net_nodes.size();gt_node_i++){
			valarray <bool> comp = (st_tree_wiz_prior_p.prior_clade_list[coal_clade_i] == my_gt_tree.descndnt[gt_node_i]);
			if ( comp.min() == true ){
				coal_i_valid=true;
				break;
			}
		}
		if (coal_i_valid){
			sp_valid=true;
			coal_clade_i++;
		}
		else{
			sp_valid=false;
			break;
		}
	}	
	
	
	if (sp_valid){
		gt_coal_in_st new_gt_coal_in_st(my_gt_tree, my_sp_net);
		probability=new_gt_coal_in_st.probability;//*st_tree_wiz_prior_p.omega;
		prob_of_maple=new_gt_coal_in_st.prob_of_maple;	
		all_w=new_gt_coal_in_st.all_w;
		all_d=new_gt_coal_in_st.all_d;
		coal_hist=new_gt_coal_in_st.coal_hist;
		prob_of_hist=new_gt_coal_in_st.prob_of_hist;
		num_enter=new_gt_coal_in_st.num_enter;
		num_out=new_gt_coal_in_st.num_out;
		gt_str=new_gt_coal_in_st.gt_str;
		sp_str=new_gt_coal_in_st.sp_str;
		//new_gt_coal_in_st.latex_print();
		//latex_print();
		
		Net  s_net_enum_dummy(st_tree_wiz_prior_p.s_net_string_enum);
		s_net_enum.clear();
		s_net_enum=s_net_enum_dummy;
		current_lambda_sum=st_tree_wiz_prior_p.lambda_sum;
		for (unsigned int coal_hist_i=0;coal_hist_i<coal_hist.size();coal_hist_i++){
			for (unsigned int coal_hist_ii=0;coal_hist_ii<coal_hist[coal_hist_i].size();coal_hist_ii++){
				for (unsigned int sp_node_i=0;sp_node_i<s_net_enum.Net_nodes.size();sp_node_i++){
					if (sp_node_i==s_net_enum.Net_nodes.size()-1){
						coal_hist[coal_hist_i][coal_hist_ii] = st_tree_wiz_prior_p.root_enum;
					}
					if (coal_hist[coal_hist_i][coal_hist_ii] == s_net_enum.Net_nodes[sp_node_i].e_num){
						coal_hist[coal_hist_i][coal_hist_ii] = s_net_enum.Net_nodes[sp_node_i].brchlen1;
						break;
					}
					
				}	
			}
		}
		construct_prob_maple();

	}
	else{
		probability=0;
		prob_of_hist_maple.push_back("0");
		prob_of_maple="0";		
	}
	
	
	
	
	
	//probability=new_gt_coal_in_st.probability;//*st_tree_wiz_prior_p.omega;
	//prob_of_maple=new_gt_coal_in_st.prob_of_maple;
	//all_d=new_gt_coal_in_st.all_d;
	//all_w=new_gt_coal_in_st.all_w;
	//sp_valid=new_gt_coal_in_st.sp_valid;
	//coal_hist=new_gt_coal_in_st.coal_hist;
	//prob_of_hist=new_gt_coal_in_st.prob_of_hist;
	//num_enter=new_gt_coal_in_st.num_enter;
	//num_out=new_gt_coal_in_st.num_out;
	//prior_coal_hist=new_gt_coal_in_st.prior_coal_hist;
	//gt_str=new_gt_coal_in_st.gt_str;
	//sp_str=new_gt_coal_in_st.sp_str;
	//Net  s_net_enum_dummy(st_tree_wiz_prior_p.s_net_string_enum);
	//s_net_enum.clear();
	//s_net_enum=s_net_enum_dummy;
	//current_lambda_sum=st_tree_wiz_prior_p.lambda_sum;
	//for (unsigned int coal_hist_i=0;coal_hist_i<coal_hist.size();coal_hist_i++){
		//for (unsigned int coal_hist_ii=0;coal_hist_ii<coal_hist[coal_hist_i].size();coal_hist_ii++){
			//for (unsigned int sp_node_i=0;sp_node_i<s_net_enum.Net_nodes.size();sp_node_i++){
				//if (sp_node_i==s_net_enum.Net_nodes.size()-1){
					//coal_hist[coal_hist_i][coal_hist_ii] = st_tree_wiz_prior_p.root_enum;
				//}
				//if (coal_hist[coal_hist_i][coal_hist_ii] == s_net_enum.Net_nodes[sp_node_i].e_num){
					//coal_hist[coal_hist_i][coal_hist_ii] = s_net_enum.Net_nodes[sp_node_i].brchlen1;
					//break;
				//}
				
			//}	
		//}
	//}
	
	//if (maple_bool){
		//construct_prob_maple();
	//}
	if (coal_debug_bool){
		cout<<probability<<endl;
		cout<<"gt_coal_in_st::gt_coal_in_st(string gt_tree_str, Net_wiz_prior_p st_tree_wiz_prior_p) END"<<endl;
		
	}	

}

void print_prior_coal_list(vector < valarray < int > > prior_clade_list){
	for (size_t i=0;i<prior_clade_list.size();i++){
		for (size_t j=0;prior_clade_list[i].size();j++){
			cout<<prior_clade_list[i][j];
		}
		cout<<" ";
	}	
	cout<<endl;
}	

	


gt_coal_in_st::gt_coal_in_st(Net gt_tree, Net st_tree){
	
	if (coal_debug_bool){
		cout<<"gt_coal_in_st::gt_coal_in_st(Net gt_tree, Net st_tree) start"<<endl;
		cout<<"sp: "<<construct_adding_new_Net_str(st_tree)<<endl;
	}
	Net pre_my_sp_net=st_tree;
	
	
	//my_sp_net=st_tree;
	Net my_gt_tree;
	bool multiple_sp_tip=false;
	for (unsigned int i=0;i<pre_my_sp_net.Net_nodes.size();i++){
		if (pre_my_sp_net.descndnt[i].sum()!=pre_my_sp_net.Net_nodes[i].num_descndnt){
			multiple_sp_tip=true;
			break;
		}
	}
	//st_tree.print_all_node();
	//string st_str=construct_adding_new_Net_str(st_tree);
	//size_t found=
	
	bool gt_tip_already_coaled=false;

	for (unsigned int i=0;i<gt_tree.Net_nodes.size();i++){
		if (gt_tree.descndnt[i].sum()!=gt_tree.Net_nodes[i].num_descndnt){
			gt_tip_already_coaled=true;
			break;
		}
	}
		
	sp_valid=true;
	//check if pre_my_sp_net is valid
	
	if (multiple_sp_tip){
		sp_valid= check_sp_valid(my_sp_net, gt_tree);
	}
	
	if (sp_valid){
		//double prior_omega_multiple_sp_tip=1;
		//checking stage2
		if (multiple_sp_tip && !gt_tip_already_coaled){
			string dummy_string=construct_adding_new_Net_str(gt_tree);
			Net load_gt_tree(dummy_string);
			//load_gt_tree.print_all_node();
			//my_gt_tree=load_gt_tree;
			vector <Node*> gt_nodes_ptr_rm;
			for (unsigned int i_link_ptr_rm=0;i_link_ptr_rm<load_gt_tree.Net_nodes.size();i_link_ptr_rm++){
				Node* new_node_ptr=NULL;
				gt_nodes_ptr_rm.push_back(new_node_ptr);
				gt_nodes_ptr_rm[i_link_ptr_rm]=&load_gt_tree.Net_nodes[i_link_ptr_rm];
			}
			for (unsigned int i_gt_node=0;i_gt_node<load_gt_tree.descndnt.size();i_gt_node++){
				for (unsigned int i_st_node=0;i_st_node<pre_my_sp_net.descndnt.size();i_st_node++){
					valarray <bool> comp = (load_gt_tree.descndnt[i_gt_node]==pre_my_sp_net.descndnt[i_st_node]);
					if (comp.min() && pre_my_sp_net.Net_nodes[i_st_node].tip_bool ){
						gt_nodes_ptr_rm[i_gt_node]->num_child=0;
						gt_nodes_ptr_rm[i_gt_node]->child.clear();
						//gt_nodes_ptr_rm[i_gt_node]->rank=1;
						gt_nodes_ptr_rm[i_gt_node]->tip_bool=true;
						//cout<<gt_nodes_ptr_rm[i_gt_node]->label<<&(gt_nodes_ptr_rm[i_gt_node]->label)<<endl;
						//cout<<gt_nodes_ptr_rm[i_gt_node]->clade<<&(gt_nodes_ptr_rm[i_gt_node]->clade)<<endl;
						gt_nodes_ptr_rm[i_gt_node]->label=gt_nodes_ptr_rm[i_gt_node]->clade;
						//cout<<gt_nodes_ptr_rm[i_gt_node]->label<<&(gt_nodes_ptr_rm[i_gt_node]->label)<<endl;
						gt_nodes_ptr_rm[i_gt_node]->node_content.clear();
					}
				}		
			}		
			//load_gt_tree.print_all_node();
			rewrite_node_content(gt_nodes_ptr_rm);
			//load_gt_tree.print_all_node();
			
			string new_gt_str=construct_adding_new_Net_str(load_gt_tree);
			//cout<<new_gt_str<<endl;
			new_gt_str=rm_and_sign(new_gt_str);
			//cout<<new_gt_str<<endl;
			string new_sp_str=construct_adding_new_Net_str(pre_my_sp_net);
			//cout<<new_sp_str<<endl;
			new_sp_str=rm_and_sign(new_sp_str);
			Net reprocess_my_sp_net(new_sp_str);
	
			Net new_gt_tree(new_gt_str);
			my_sp_net=reprocess_my_sp_net;
			my_gt_tree=new_gt_tree;
			//my_gt_tree.print_all_node();
			//cout<<new_gt_str<<"********************************************************* end"<<endl;
		}
		else{
			my_sp_net=st_tree;
			my_gt_tree=gt_tree;
			
		}
		gt_str=construct_adding_new_Net_str(my_gt_tree);
		
		sp_str=construct_adding_new_Net_str(my_sp_net);
		if (coal_debug_bool){
			cout<<"gt_str "<<gt_str<<endl;
			cout<<"sp_str "<<sp_str<<endl;
		}
		string enum_str=construct_adding_new_Net_str(my_sp_net);
		
		Net  s_net_enum_dummy(enum_str);
		//s_net_enum.clear();
		vector <Node*> old_Net_node_ptr;
		for (unsigned int i=0;i<s_net_enum_dummy.Net_nodes.size();i++){
			Node* new_node_ptr=NULL;
			old_Net_node_ptr.push_back(new_node_ptr);
			old_Net_node_ptr[i]=&s_net_enum_dummy.Net_nodes[i];	
			old_Net_node_ptr[i]->brchlen1=old_Net_node_ptr[i]->e_num;
		}
		
		s_net_enum=s_net_enum_dummy;
		
	
		//my_gt_tree.print_all_node();
		//my_sp_net.print_all_node();
		unsigned int sp_max_enum=my_sp_net.Net_nodes.back().e_num;
		int gt_max_enum=my_gt_tree.Net_nodes.back().e_num;
		
		//int M_matrix[sp_max_enum][gt_max_enum-1];
		//int S_matrix[sp_max_enum][sp_max_enum];
		vector < vector < int > > S_matrix=building_S_matrix(my_sp_net);
		vector < vector < int > > R_matrix=building_R_matrix(my_gt_tree);	
		vector < vector < int > > M_matrix=building_M_matrix(my_gt_tree,my_sp_net);
		vector < vector <unsigned int> > coal_hist_mat;
		int max_coal_hist_num=1;
		for (int i_gt_enum=0;i_gt_enum<gt_max_enum-1;i_gt_enum++){
			vector <unsigned int> coal_hist_vec;
			for (unsigned int i_sp_enum=0;i_sp_enum<sp_max_enum;i_sp_enum++){
				if (M_matrix[i_sp_enum][i_gt_enum]==1){
					coal_hist_vec.push_back(i_sp_enum+1);
				}
			}
			max_coal_hist_num=max_coal_hist_num*(coal_hist_vec.size());
			coal_hist_mat.push_back(coal_hist_vec);
		}
		
		for (unsigned int first_coal_mat_i=0;first_coal_mat_i<coal_hist_mat[0].size();first_coal_mat_i++){
			vector <unsigned int> coal_hist_dummy;
			coal_hist_dummy.push_back(coal_hist_mat[0][first_coal_mat_i]);
			coal_hist.push_back(coal_hist_dummy);
		}
		
		if (gt_max_enum-1>1){
			coal_hist=build_coal_hist(coal_hist,1,coal_hist_mat,R_matrix);
		}
		//for (unsigned int coal_hist_i=0;coal_hist_i<coal_hist.size();coal_hist_i++){
		//coal_hist[coal_hist_i].push_back(sp_max_enum);
		//}
		vector < double > brchlens_vec(sp_max_enum,0);
		vector < int > max_num_brch_vec(sp_max_enum,0);
		
		for (unsigned int node_i=0;node_i<my_sp_net.Net_nodes.size();node_i++){
			int index_enum=my_sp_net.Net_nodes[node_i].e_num;
			
			if (index_enum!=0){
				
				brchlens_vec[index_enum-1]=my_sp_net.Net_nodes[node_i].brchlen1;
				//cout<<my_sp_net.Net_nodes[node_i].label<<"  "<<my_sp_net.Net_nodes[node_i].brchlen1<<endl;
				max_num_brch_vec[index_enum-1]=my_sp_net.Net_nodes[node_i].num_descndnt;
				//cout<<index_enum<<" "<<brchlens_vec[index_enum-1]<<" "<<max_num_brch_vec[index_enum-1]<<endl;
			}
		}
	
		//for (unsigned int checking_i=0;checking_i<sp_max_enum;checking_i++){
			//cout<<checking_i+1<<" "<<max_num_brch_vec[checking_i]<<"  "<<brchlens_vec[checking_i]<<endl;
		//}
				//cout<<"HEA"<<endl;
		for (unsigned int i_coal_hist=0;i_coal_hist<coal_hist.size();i_coal_hist++){	
			//for (unsigned int j_coal_hist=0;j_coal_hist<coal_hist[i_coal_hist].size();j_coal_hist++){
				//cout<<coal_hist[i_coal_hist][j_coal_hist];
			//}
			vector <int> num_enter_i_coal;
			vector <int> num_out_i_coal;
			vector <int> num_coal_in_branch_i_coal;
			vector <int> w_i_coal;
			vector <int> d_i_coal;
			for (unsigned int i=0;i<sp_max_enum;i++){
				for (unsigned int ith_node=0;ith_node<my_sp_net.Net_nodes.size();ith_node++){
					if (!my_sp_net.Net_nodes[ith_node].tip_bool && my_sp_net.Net_nodes[ith_node].e_num==i+1){
						num_enter_i_coal.push_back(my_sp_net.Net_nodes[ith_node].num_descndnt);
					}
				}
			}	
			for (unsigned int i=0;i<sp_max_enum;i++){
				vector < int > clades_coal_in_branch;
				int num_coal_in_branch_i_coal_dummy=0;
				for (unsigned int coal_hist_position=0;coal_hist_position<coal_hist[i_coal_hist].size();coal_hist_position++){
					if (coal_hist[i_coal_hist][coal_hist_position]==i+1){
						num_coal_in_branch_i_coal_dummy++;
						clades_coal_in_branch.push_back(coal_hist_position);
						for (unsigned int i_S_mat=i;i_S_mat<(sp_max_enum);i_S_mat++){
							if (S_matrix[i_S_mat][i]==1){
								num_enter_i_coal[i_S_mat]--;
							}
						}	
					}			
				}
				num_coal_in_branch_i_coal.push_back(num_coal_in_branch_i_coal_dummy);	
				num_out_i_coal.push_back(num_enter_i_coal[i]-num_coal_in_branch_i_coal[i]);
				int d_dummy=1;
				int coal_dummy=num_enter_i_coal[i];
				for (int y=0;y<num_coal_in_branch_i_coal_dummy;y++){
					//d_dummy=d_dummy*n_choose_k(coal_dummy,2);
					d_dummy=d_dummy*coal_dummy*(coal_dummy-1)/2;
					coal_dummy--;
				}
				d_i_coal.push_back(d_dummy);
				
				int w_dummy=factorial(num_coal_in_branch_i_coal_dummy);
				vector < vector < int > > updated_R_mat;
				if (clades_coal_in_branch.size()>1){
					for (unsigned int updated_R_i=0;updated_R_i<clades_coal_in_branch.size();updated_R_i++){
						vector < int > updated_R_mat_row;
						for (unsigned int updated_R_j=0;updated_R_j<updated_R_i;updated_R_j++){
							updated_R_mat_row.push_back(R_matrix[clades_coal_in_branch[updated_R_i]][clades_coal_in_branch[updated_R_j]]);
						}
						updated_R_mat.push_back(updated_R_mat_row);
					}
					for (unsigned int updated_R_mat_i=0;updated_R_mat_i<updated_R_mat.size();updated_R_mat_i++){
						int sum_r=1;
						for (unsigned int updated_R_mat_j=0;updated_R_mat_j<updated_R_mat[updated_R_mat_i].size();updated_R_mat_j++){
							sum_r=sum_r+updated_R_mat[updated_R_mat_i][updated_R_mat_j];
						}
						w_dummy=w_dummy/sum_r;
					}
				}	
				
				w_i_coal.push_back(w_dummy);
				//if (w_dummy/d_dummy!=1){
					//cout<<w_dummy<<"/"<<d_dummy;
				//}
				//if (i<sp_max_enum-1){
					//cout<<"p_{"<<num_enter_i_coal[i]<<num_out_i_coal[i]<<"}(lambda_{"<<i+1<<"})";	
				//}
			}
			//cout<<endl;
			num_enter.push_back(num_enter_i_coal);
			num_out.push_back(num_out_i_coal);
			num_coal_in_branch.push_back(num_coal_in_branch_i_coal);
			all_w.push_back(w_i_coal);
			all_d.push_back(d_i_coal);	
			coal_hist[i_coal_hist].push_back(sp_max_enum);
		}
		
		//vector < vector < vector < double > > > gijoe_matrix=
		build_gijoe_matrix(brchlens_vec, max_num_brch_vec);
		calc_prob_of_hist();
		
	}
	else{
		probability=0;
		prob_of_hist_maple.push_back("0");
		prob_of_maple="0";
	}
	if (coal_debug_bool){
		cout<<probability<<endl;
		cout<<"gt_coal_in_st::gt_coal_in_st(Net gt_tree, Net st_tree) end"<<endl;
		
	}
}


double gijoe(
int u, /*!< number of branch in */
int v, /*!< number of branch out */
double T) /*!< branch length*/
{
	//if (T==0){
		//return 1;
	//}
	//else{
	double sums=0;
	for (int k=v;k<=u;k++){
		double prods=exp(-k*(k-1)*T/2)*(2*k-1)*pow(-1.0,1.0*(k-v))/factorial(v*1.0)/factorial((k-v)*1.0)/(v+k-1);
		for (int y=0;y<k;y++){
			prods=prods*(v+y)*(u-y)/(u+y);
		}
		sums=sums+prods;
	}
	return sums;//}
}




void gt_coal_in_st::build_gijoe_matrix(
vector < double > brchlens_vec,
 vector < int > max_num_brch_vec )
 {
	 gijoe_matrix.clear();
	double v1,v2,v3,u1,u2;
	int u,v,h,k;
	for (size_t b=0;b<brchlens_vec.size();b++){
		vector < vector < double > > gijoe_matrix_b;
		vector < double > empty_gijoe_matrix_b_1;
		gijoe_matrix_b.push_back(empty_gijoe_matrix_b_1);
		for (u=2;u<=max_num_brch_vec[b];u++){
			
			vector < double > gijoe_matrix_b_u;
				for (v=1;v<=u;v++){
					double gi_joe_buv;
					//if (brchlens_vec[b]==0){
					if (b == brchlens_vec.size()-1){
						gi_joe_buv=1;
						}
					else{
						gi_joe_buv=0;
						for(k=v; k <= u; k++) {
							v1 = 1.0;
							for(h = v; h <= v+k-2; h++){
								v1 *= h;
							}
							v2 = 1.0; 
							for(h = v; h > 1; h--){         
								v2 *= h;
							}
							v3 = 1.0;
							for(h = k-v; h > 1; h--){
								v3 *= h;
							}
							u1 = 1.0;
							for(h = u; h >= u - k + 1; h--){
								u1 *= h;
							}
							u2 = 1.0;
							for(h = u; h <= u+k-1; h++){
								u2 *= h;
							}
							gi_joe_buv += exp(.5*k*(1.0-k)*brchlens_vec[b])*(2.0*k-1.0)*pow(-1.0,(k-v)*1.0)*v1*u1/(v2*v3*u2);
						}
						
						
					}
					
					//cout<<"b:"<<b+1<<" u:"<<u<<" v:"<<v<<endl;
					//cout<<gi_joe_buv<<"  "<<gijoe(u,v,brchlens_vec[b])<<endl;
					//if (gi_joe_buv==gijoe(u,v,brchlens_vec[b])){cout<<"yes"<<endl;}
					//else{cout<<"no"<<endl;}
					//cout<<"g_["<<u<<","<<v<<"]("<<brchlens_vec[b]<<")="<<gi_joe_buv<<"  ";
					gijoe_matrix_b_u.push_back(gi_joe_buv);
					//cout<<"      "<<gijoe_matrix_b_u.size()<<endl;
				}
				//cout<<endl;
			gijoe_matrix_b.push_back(gijoe_matrix_b_u);
			//cout<<"  "<<gijoe_matrix_b.size()<<endl;
		}
		//cout<<endl<<endl;
		gijoe_matrix.push_back(gijoe_matrix_b);
		//cout<<gijoe_matrix.size()<<endl;
	}
		//for (unsigned int mat_b=0;mat_b<gijoe_matrix.size();mat_b++){
			//for (unsigned int mat_i=0;mat_i<gijoe_matrix[mat_b].size();mat_i++){
				//for (unsigned int mat_j=0;mat_j<gijoe_matrix[mat_b][mat_i].size();mat_j++){
					//cout<<"g_["<<mat_i<<","<<mat_j<<"]("<<brchlens_vec[mat_b]<<")="<<gijoe_matrix[mat_b][mat_i][mat_j]<<"  ";
				
				//}
			//cout<<endl;
				
			//}
			//cout<<endl;
		//}
	//return gijoe_matrix;
}

bool check_sp_valid(Net sp_net, Net gt_tree){
	bool sp_valid_out=true;
	for (unsigned int coal_clade_i=0;coal_clade_i<sp_net.descndnt.size();){
		if (sp_net.Net_nodes[coal_clade_i].tip_bool){
			bool coal_i_valid=false;
			for (unsigned int gt_node_i=0;gt_node_i<gt_tree.Net_nodes.size();gt_node_i++){
				valarray <bool> comp = (sp_net.descndnt[coal_clade_i] == gt_tree.descndnt[gt_node_i]);
				if ( comp.min() == true ){
					coal_i_valid=true;
					break;
				}
			}
			if (coal_i_valid){
				sp_valid_out=true;
				coal_clade_i++;
			}
			else{
				sp_valid_out=false;
				break;
			}
		}
		else{
			coal_clade_i++;}
	}	
	return sp_valid_out;
}



bool check_multiple_sp_tip(Net pre_current_net){
	bool multiple_sp_tip=false;
	for (unsigned int i=0;i<pre_current_net.Net_nodes.size();i++){
		if (pre_current_net.descndnt[i].sum()!=pre_current_net.Net_nodes[i].num_descndnt){
			multiple_sp_tip=true;
			break;
		}
	}
	return 	multiple_sp_tip;
}

bool check_gt_tip_already_coaled(Net my_gt_tree){
	bool gt_tip_already_coaled=false;
	for (unsigned int i=0;i<my_gt_tree.Net_nodes.size();i++){
		if (my_gt_tree.descndnt[i].sum()!=my_gt_tree.Net_nodes[i].num_descndnt){
			gt_tip_already_coaled=true;
			break;
		}
	}
	
	return gt_tip_already_coaled;
}
