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

/*! \file main.cpp
 *  \brief Main function of hybrid_coal
 */

#include"export.hpp"

#include"all_gene_topo.hpp"
#include"coal.hpp"

//using namespace std;

//bool log_file_bool; /*! \todo fix log file */ 


//bool acc_bool;

//bool debug_bool;
//bool export_debug_bool;
//bool coal_debug_bool;
//bool utility_debug_bool;
//bool rm_debug_bool; 

/*! 
 * \fn int main(int argc, char *argv[])
 * \brief Main function for hybrid_coal 
 * */
int main(int argc, char *argv[]){
bool maple_bool;/*! \todo get rid of this */  //=false; 

	//system("reset");
	//check_and_remove("debug_file");
	check_and_remove("log_file");
	
	//debug_bool=false;
	
	bool main_debug_bool=false; 	
	bool export_debug_bool=false;
	bool coal_debug_bool=false;
	bool utility_debug_bool=false;
	bool rm_debug_bool=false;
	
	bool acc_bool=false;
	
//bool latex_bool;
	bool latex_bool=false;
	bool pdf_bool=false;	
	string latex_F_name="latex_prob.tex";
	
	maple_bool=false;
	string maple_F_name="maple_prob.mw";
	
	bool print_gene_topo_bool=false;
	string gtopo_F_name="GENE_topo";
	
	bool help=false;
	string net_str;
	vector <string> gt_tree_str_s;
	
	bool all_gt_tree_bool=true;
	bool symb_bool=false;
	bool print_tree=false;
		
	bool list_sub_network_bool=false;

	int plot_option;//=0;
	bool plot_label=false;
	bool plot_branch=false;
	bool plot_bool=false;
	string tex_fig_name="texfigure.tex";
	bool dot_bool=false;
	string dot_fig_name="figure.dot";
	bool out_bool=false;
	string out_file="out_coal";

	for (int argc_i=0;argc_i<argc;argc_i++){
		string argv_i(argv[argc_i]);
		
		if (argv_i=="-gt"){
			gt_tree_str_s=read_input_lines(argv[argc_i+1]);
			all_gt_tree_bool=false;
			string gt_log=argv[argc_i+1];
			gt_log="Gene trees Input: "+gt_log;
			appending_log_file(gt_log);
		}
		if (argv_i=="-sp"){
			net_str=read_input_line(argv[argc_i+1]);
			string sp_log=argv[argc_i+1];
			sp_log="Species Input: "+sp_log;
			appending_log_file(sp_log);
			sp_log="Species structure: "+net_str;
			appending_log_file(sp_log);

		}
		if (argv_i=="-latex"){
			latex_bool=true;			
		}
		
		if (argv_i=="-latexF"){
			latex_bool=true;
			latex_F_name=Fname_ext(argv[argc_i+1],".tex");
		}
		if (argv_i=="-out"){
			out_bool=true;
			out_file=argv[argc_i+1];
		}
		
		if (argv_i=="-pdf"){
			latex_bool=true;
			pdf_bool=true;
		}
		//if (argv_i=="-history"){
			//coal_history=true;
		//}
		if (argv_i=="-maple"){
			maple_bool=true;
		}
		if (argv_i=="-mapleF"){
			maple_bool=true;
			maple_F_name=Fname_ext(argv[argc_i+1],".mw");
		}
		
		if (argv_i=="-h" || argv_i=="-help" ){
			help=true;
		}

		//if (argv_i=="-debug"){
			//debug_bool=true;
		//}
		if (argv_i=="-debug"){
			//samples_bool=true;
			for (int argc_j=argc_i+1;argc_j<argc;argc_j++){
				string s(argv[argc_j]);
				if (s[0]=='-'){
					break;
				}
				if (s=="coal"){
					coal_debug_bool=true;
				}
				if (s=="utility"){
					utility_debug_bool=true;
				}
				if (s=="rm"){
					rm_debug_bool=true;
				}
				if (s=="export"){
					export_debug_bool=true;
				}
				if (s=="main"){
					main_debug_bool=true;
				}
			}
		}
		
		if (argv_i=="-symb"){
			symb_bool=true;
		}
		if (argv_i=="-print"){
			print_tree=true;
		}
		if (argv_i=="-plot"){
			plot_bool=true;
		}
		if (argv_i=="-plotF"){
			plot_bool=true;
			tex_fig_name=argv[argc_i+1];
		}
		
		
		if (argv_i=="-label"){
			plot_label=true;
		}

		if (argv_i=="-branch"){
			plot_branch=true;
		}		

		if (argv_i=="-dot"){
			dot_bool=true;
		}
		if (argv_i=="-dotF"){
			dot_bool=true;
			dot_fig_name=argv[argc_i+1];
		}
		

		if (argv_i=="-gtopo"){
			print_gene_topo_bool=true;
		}
		if (argv_i=="-gtopoF"){
			print_gene_topo_bool=true;
			gtopo_F_name=argv[argc_i+1];
		}
		
		
		if (argv_i=="-sub"){
			list_sub_network_bool=true;
		}
		if (argv_i=="-acc"){
			acc_bool=true;
		}
		
	}
		
	if (help || argc==1){
		print_help();
	}
	else{
		
		if (net_str.size()<0){
			exit(1);
		}
		
		//process net_str, check if gt have multiple lineages
		if (!all_gt_tree_bool){
			Net test_gt(gt_tree_str_s[0]);
			//Net test_sp(net_str);
			if (test_gt.tip_name.size()>test_gt.tax_name.size()){
				net_str=rewrite_net_str(net_str, test_gt.tax_name,test_gt.tip_name );
			}
			for (size_t i=0;i<gt_tree_str_s.size();i++){
				gt_tree_str_s[i]=rm_dash_sign(gt_tree_str_s[i]);
				//cout<<gt_tree_str_s[i]<<endl;
			}
			cout<<net_str<<endl;
		}
		
		
		
		
		Net net_dummy(net_str);
		if (print_tree){
			net_dummy.print_all_node();
			appending_log_file("Tree printed");
		}
		plot_option=set_plot_option(plot_label,plot_branch);

		if (plot_bool){
			plot_in_latex_file(tex_fig_name.c_str(), net_dummy,plot_option);	
		}
		
		if (dot_bool){
			plot_in_dot(dot_fig_name.c_str(), net_dummy,plot_option);			
		}
			
		if (all_gt_tree_bool){
			gt_tree_str_s=all_n_tax_gene_tree(net_dummy.tip_name);
			if (print_gene_topo_bool){
				print_all_gt_topo(gtopo_F_name.c_str(),gt_tree_str_s);
			}
		}

		if (list_sub_network_bool && gt_tree_str_s.size()==1){
			list_sub(net_str,gt_tree_str_s[0]);
		}

		if (list_sub_network_bool && all_gt_tree_bool){
			string 	empty_str="";
			list_sub(net_str,empty_str);
		}
				
		if (print_tree || plot_bool || dot_bool || print_gene_topo_bool || list_sub_network_bool && gt_tree_str_s.size()>1){
			//return 0;
			return my_exit();
		}		
					
				
		if (latex_bool){
			
			latex_header(latex_F_name.c_str());	
			plot_in_latex(latex_F_name.c_str(), net_str,1);		
		}
		
		check_and_remove(out_file.c_str());
		
		if (net_dummy.is_net){		
		
			double total_prob=0;
			for (size_t topo_i=0;topo_i<gt_tree_str_s.size();topo_i++){
				string current_gene_topo=gt_tree_str_s[topo_i];
				//check_gt_str(current_gene_topo);

				if (main_debug_bool){
					cout<<current_gene_topo<<endl;
				}
				
				vec_Net_wiz_prior_p net_str_s=simplify_Networks(net_str,current_gene_topo);
				//vec_Net_wiz_prior_p * net_str_s_pt;
				//net_str_s_pt= new vec_Net_wiz_prior_p;
				//vec_Net_wiz_prior_p net_str_s=simplify_Networks(net_str,current_gene_topo);
				//net_str_s_pt=&net_str_s;
				
				if (latex_bool){
					latex_nt_body1(latex_F_name.c_str(),topo_i+1,current_gene_topo);
				}
				double total_prob_topo_i=0;
				for (size_t i=0;i<net_str_s.Net_vec.size();i++){
					gt_coal_in_st gt_in_sp(current_gene_topo,net_str_s.Net_vec[i]);
					//gt_coal_in_st gt_in_sp(current_gene_topo,net_str_s_pt->Net_vec[i]);
					if (main_debug_bool){
						cout<<current_gene_topo<<endl;
						//cout<<net_str_s_pt->Net_vec[i].s_net_string<<endl;
						cout<<topo_i<<"  "<<i<<"  "<<net_str_s.Net_vec[i].omega <<"  "<<net_str_s.Net_vec[i].s_net_string<<endl;
					}
					total_prob_topo_i=total_prob_topo_i+net_str_s.Net_vec[i].omega*gt_in_sp.probability;
					//total_prob_topo_i=total_prob_topo_i+net_str_s_pt->Net_vec[i].omega*gt_in_sp.probability;
					if (latex_bool){
						latex_nt_body2(latex_F_name.c_str(),net_str_s, i);
						gt_in_sp.latex_print(latex_F_name.c_str());
					}
					gt_in_sp.clear();
					net_str_s.Net_vec[i].clear();
				}
				
				prob_out_body(out_file.c_str(),topo_i+1,gt_tree_str_s[topo_i],total_prob_topo_i,net_dummy.tax_name.size());
				
				total_prob=total_prob+total_prob_topo_i;
				
				net_str_s.clear();
			}
			prob_out_tail(out_file.c_str(),gt_tree_str_s[0].size(), total_prob);
			if (main_debug_bool){
				cout<<"End of gene tree probability given species network"<<endl;
			}
		}
		else{
			if (main_debug_bool){
				cout<<"Start of gene tree probability given species tree"<<endl;
			}
			net_dummy.clear();
			net_str=rm_one_child_root(net_str);
			
			Net new_net_dummy(net_str);
		
			net_dummy=new_net_dummy;
			double total_prob=0;
			vector < valarray <int> > prior_coal_clades_dummy;
			vector <int> prior_coal_hist;

			
			
			double total_prob_topo_i=0;

			if (prior_coal_clades_dummy.size()>0){
				Net_wiz_prior_p coaled_sp_wiz_prior;
	Net_wiz_prior_p initial_Networks;	
	Net old_net_enum=initialize_enum_net(net_str);
	coaled_sp_wiz_prior.s_net_string=net_str;
	coaled_sp_wiz_prior.omega=1; // fix this!!!!!!!!!!!!!!!
	coaled_sp_wiz_prior.lambda_sum=initialize_lambda_sum(old_net_enum);
	coaled_sp_wiz_prior.s_net_string_enum=construct_adding_new_Net_str(old_net_enum);
	coaled_sp_wiz_prior.root_enum=old_net_enum.Net_nodes.back().e_num;
	coaled_sp_wiz_prior.prior_clade_list=prior_coal_clades_dummy;//=initial_prior_clade_list;
	coaled_sp_wiz_prior.prior_coal_hist=prior_coal_hist;//=initial_prior_coal_hist;
				
				for (size_t topo_i=0;topo_i<gt_tree_str_s.size();topo_i++){
					//gt_coal_in_st gt_in_sp(current_gene_topo,net_str_s.Net_vec[i]);
					string current_gene_topo=gt_tree_str_s[topo_i];
					gt_coal_in_st gt_in_sp(current_gene_topo, coaled_sp_wiz_prior);
					//total_prob_topo_i=total_prob_topo_i+net_str_s.Net_vec[i].omega*gt_in_sp.probability;
					
					if (latex_bool){
						//latex_nt_body2(latex_F_name.c_str(),net_str_s, 0);
						gt_in_sp.latex_print(latex_F_name.c_str());
					}
					total_prob_topo_i=total_prob_topo_i+coaled_sp_wiz_prior.omega*gt_in_sp.probability;
					
					prob_out_body(out_file.c_str(),topo_i+1,gt_tree_str_s[topo_i],gt_in_sp.probability,net_dummy.tax_name.size());
					total_prob=total_prob+total_prob_topo_i;
					gt_in_sp.clear();
				}
			
				
			}
			else{
				for (size_t topo_i=0;topo_i<gt_tree_str_s.size();topo_i++){
					Net current_gt(gt_tree_str_s[topo_i]);
					gt_coal_in_st gt_in_sp(current_gt,net_dummy);
					if (latex_bool){
						latex_tre_body(latex_F_name.c_str(),gt_tree_str_s[topo_i],net_str);
						gt_in_sp.latex_print(latex_F_name.c_str());
					}
					
					prob_out_body(out_file.c_str(),topo_i+1,gt_tree_str_s[topo_i],gt_in_sp.probability,net_dummy.tax_name.size());
				
					total_prob=total_prob+gt_in_sp.probability;
					gt_in_sp.clear();
				}
				
			}
			prob_out_tail(out_file.c_str(),gt_tree_str_s[0].size(), total_prob);
			
			if (main_debug_bool){
				cout<<"End of gene tree probability given species tree"<<endl;
			}
		}

		if (latex_bool){
			latex_tail(latex_F_name.c_str());
		}
		if (maple_bool){
			produce_maple(maple_F_name.c_str(),net_str, symb_bool,  gt_tree_str_s);
		}
		
	}
	
	if (pdf_bool){
		string command="pdflatex "+latex_F_name;
		system(command.c_str());
		//system("pdflatex latex_prob.tex");
	}
		return my_exit();
	//return 0;
}

