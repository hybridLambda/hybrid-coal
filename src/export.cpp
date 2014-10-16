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

/*! \file export.cpp
 *  \brief Functions to write output into maple, latex and other files
 */
 
 
#include"export.hpp"
bool export_debug_bool=false;

void print_all_gt_topo(const char* file_name,vector <string> gt_tree_str_s){
	check_and_remove(file_name);
	ofstream out_file;
	out_file.open(file_name, ios::out | ios::app | ios::binary);
	for (size_t i=0;i<gt_tree_str_s.size();i++){			
		out_file<<gt_tree_str_s[i]<<";"<<endl;		
	}
	//out_file<<gt_tree_str_s.size()<<endl;
	out_file.close();
	string appending_str(file_name);
	ostringstream num_topo_strm;
	num_topo_strm<< gt_tree_str_s.size();
	appending_str=num_topo_strm.str()+" gene tree topologies listed in file: "+appending_str;
	appending_log_file(appending_str);
}



void latex_header(const char* file_name){
	check_and_remove(file_name);
	ofstream latex_file;
	latex_file.open (file_name, ios::out | ios::app | ios::binary); 
	latex_file <<"\\documentclass[9pt]{article}\n";
	latex_file <<"\\usepackage{tikz,graphics,graphicx,lscape,fullpage,multicol,setspace}\n \\singlespacing\n \\begin{document}\n ";	
	latex_file<<"\\ifx\\du\\undefined\\newlength{\\du}\\fi\\setlength{\\du}{35\\unitlength}\n";
	latex_file.close();
	if (export_debug_bool){
		cout<<"latex_header END"<<endl;
	}
}

void latex_nt_body1(const char* file_name,size_t topo_idx,string gt_topo){
	ofstream latex_file;
	latex_file.open (file_name, ios::out | ios::app | ios::binary);
	latex_file<<"\\noindent ********************   "<<topo_idx<<"   **********************\\\\"<<endl;
	latex_file<<"\\begin{verbatim}GT: "<<gt_topo<<"\\end{verbatim}"<<endl;
	latex_file.close();	
	if (export_debug_bool){
		cout<<"latex_nt_body1 END"<<endl;
	}
}

//void latex_nt_body2(const char* file_name,vec_Net_wiz_prior_p net_str_s,size_t i){
	//ofstream latex_file;
	//latex_file.open (file_name, ios::out | ios::app | ios::binary);
	
	//latex_file<<"\\begin{verbatim}ST: "<<tree_topo(net_str_s.Net_vec[i].s_net_string)<<"\\end{verbatim} "<<endl;
	//latex_file<<"$"<<net_str_s.Net_vec[i].omega_string<<"$\\\\"<<endl;
	//latex_file<<"$\\omega="<<net_str_s.Net_vec[i].omega<<"$\\\\"<<endl;
	//latex_file.close();
	//////plot_enum_in_latex("latex_prob.tex",net_str_s.Net_vec[i]);
	//plot_in_latex(file_name, net_str_s.Net_vec[i].s_net_string_enum,2);		
	//////plot_in_latex("latex_prob.tex", net_str_s.Net_vec[i].s_net_string,1);		

	//if (export_debug_bool){
		//cout<<"latex_nt_body2 END"<<endl;
	//}
//}

void latex_tre_body(const char* file_name,string gt_str,string sp_str){
	ofstream latex_file;
	latex_file.open (file_name, ios::out | ios::app | ios::binary);
	latex_file<<"\\noindent ****************************************************\\\\"<<endl;
	latex_file<<"\\begin{verbatim}Load in GT: "<<gt_str<<"\\end{verbatim}"<<endl;
	latex_file<<"\\begin{verbatim}ST: "<<tree_topo(sp_str)<<"\\end{verbatim}"<<endl;
	latex_file.close();
	if (export_debug_bool){
		cout<<"latex_tre_body END"<<endl;
	}
}

void latex_tail(const char* file_name){
	ofstream latex_file;
	latex_file.open (file_name, ios::out | ios::app | ios::binary);
	latex_file <<"\\end{document}\n";
	latex_file.close();
	string file_str(file_name);
	string appending_log_str="Coalescent history generated in file: "+file_str;
	appending_log_file(appending_log_str);
	if (export_debug_bool){
		cout<<"latex_tail END"<<endl;
	}	
}


void list_sub_tree(string in_str,string gt_str){
	//check_and_remove("sub_trees");
	//ofstream sub_trees_file;
	//sub_trees_file.open("sub_trees", ios::out | ios::app | ios::binary);
	//////string empty_str="";
	//////vec_Net_wiz_prior_p net_str_s=simplify_Networks(in_str,empty_str);
	//vec_Net_wiz_prior_p net_str_s=simplify_Networks(in_str,gt_str);
	
	//for (size_t i=0;i<net_str_s.Net_vec.size();i++){
		//sub_trees_file<<net_str_s.Net_vec[i].s_net_string<<"\t"<<net_str_s.Net_vec[i].omega <<endl;
		
	//}
	//////sub_trees_file<<net_str_s.Net_vec.size()<<" "<<" "<<endl;
	//sub_trees_file.close();
	//appending_log_file("Simplified species tree are listed in: sub_trees");
	
}



void list_sub_net(string in_str,string gt_str){
	//check_and_remove("sub_networks");
	//check_and_remove("reverse_sub_net");
	////string 	empty_str="";
	////simplify_Networks_maple_block("reverse_sub_net",in_str, empty_str,1);
	//simplify_Networks_maple_block("reverse_sub_net",in_str, gt_str,1);

	//ifstream maple_file;
	//maple_file.open("reverse_sub_net");
	//ofstream sub_network_file;
	//sub_network_file.open("sub_networks", ios::out | ios::app | ios::binary);
	//string maple_file_content;
	//getline (maple_file,maple_file_content);
	//int count_subnetwork=0;
	//while (maple_file_content.size()>0){   
		////getline (gt_tree_file,gt_tree_str);
		//if (maple_file_content.find(";#s[")!=string::npos){
			//count_subnetwork++;
			//sub_network_file<<maple_file_content.substr(maple_file_content.find("]=")+2,maple_file_content.size())<<endl;
		//}
		////gt_tree_str_s.push_back(gt_tree_str);
		//getline (maple_file,maple_file_content);
	//}
	//maple_file.close();
	//check_and_remove("reverse_sub_net");		
	//sub_network_file<<count_subnetwork<<endl;
	//sub_network_file.close();
	//appending_log_file("Simplified species tree are listed in: sub_networks");	
}



void list_sub(string in_str,string gt_str){
	list_sub_tree(in_str,gt_str);
	list_sub_net(in_str,gt_str);
}


//void produce_maple(const char* file_name, string in_str, bool symb_bool, vector <string> gt_tree_str_s){
	//if (export_debug_bool){
		//cout<<"produce_maple START"<<endl;
	//}
	//maple_head(file_name,in_str,symb_bool);
	
	//for (size_t i=0;i<gt_tree_str_s.size();i++){
		////vector <string> maple_block=simplify_Networks_maple_block(old_net_string, gt_tree_str_s[i],i);
		//simplify_Networks_maple_block(file_name,in_str, gt_tree_str_s[i],i);
		////for (size_t maple_block_i=0;maple_block_i<maple_block.size();maple_block_i++){
			////maple_file_rv_vec.push_back(maple_block[maple_block_i]);
			
		////}
	//}			
	//maple_tail(file_name,gt_tree_str_s.size()-1);
//}

void maple_head(const char* file_name, string in_str,bool symb_bool){
	check_and_remove(file_name);
	ofstream maple_file;
	maple_file.open(file_name, ios::out | ios::app | ios::binary);	
	maple_file<<"restart;"<<endl;
	maple_file<<"puvT:=(u,v,T)->sum(exp(-k*(k-1)*T/2) *(2*k-1)*(-1)^(k-v) / (v! * (k-v)!*(v+k-1)) *product((v+y)*(u-y)/(u+y),y=0..k-1),k=v..u);"<<endl;
	maple_file<<"#s="<<in_str<<endl;

	if (!symb_bool){
		symb_maple(file_name,in_str);
	}
	maple_file.close();	
}

void symb_maple(const char* file_name,string in_str){
	//Net * net_dummy = new Net(in_str);
	//ofstream maple_file;
	//maple_file.open(file_name, ios::out | ios::app | ios::binary);	
	//maple_file<<"unprotect(gamma); "<<endl;
	
	//for (size_t i=0;i<net_dummy->Net_nodes.size()-1;i++){
		//if (!net_dummy->Net_nodes[i]->tip()){
			//maple_file<<"lambda["<<net_dummy->Net_nodes[i].e_num<<"]:="<<net_dummy->Net_nodes[i].brchlen1<<";"<<endl;
			//if (net_dummy->Net_nodes[i].hybrid){
				//maple_file<<"lambda["<<net_dummy->Net_nodes[i].e_num2<<"]:="<<net_dummy->Net_nodes[i].brchlen2<<";"<<endl;		
				//size_t new_node_name_i=net_dummy->Net_nodes[i].label.find("#");
				////for (new_node_name_i=0;new_node_name_i<net_dummy->Net_nodes[i].label.size();new_node_name_i++){
					////if (net_dummy->Net_nodes[i].label[new_node_name_i]=='#'){
						////break;
					////}
				////}
				//maple_file<<"gamma["<<net_dummy->Net_nodes[i].label.substr(1,new_node_name_i-1) <<"]:="<<net_dummy->Net_nodes[i].label.substr(new_node_name_i+1,net_dummy->Net_nodes[i].label.size()-1) <<";"<<endl;		
			//}
		//}
	//}	
	//maple_file.close();	
}



void maple_tail(const char* file_name, size_t num_gt){
	ofstream maple_file;
	maple_file.open(file_name, ios::out | ios::app | ios::binary);	
	maple_file<<"numgt:="<<num_gt+1<<";"<<endl;
	maple_file<<"total:=0;for i from 0 by 1 to "<<num_gt<<" do total:=total+p(g[i],s[0]) end do;"<<endl;
	maple_file<<"writeto(\"maple_prob_out.txt\");"<<endl;
	maple_file<<"simplify(total);"<<endl;
	maple_file<<"writeto(terminal);"<<endl;
	maple_file.close();
	string appending_log_str(file_name);
	appending_log_str="Maple file produced in: "+appending_log_str;
	appending_log_file(appending_log_str);
	
}


void prob_out_body(const char* file_name,size_t topo_idx,string topo,double prob,size_t num_tax){
	if (num_tax<5){
		cout<<setw (4)<<topo_idx<<"  "<<topo<<"   "<<prob<<endl;
	}
	ofstream coal_out_file;
	coal_out_file.open (file_name, ios::out | ios::app | ios::binary); 
	coal_out_file<<setw (4)<<topo_idx<<"  "<<topo<<"   "<<prob<<"\n";
	coal_out_file.close();
}

void prob_out_tail(const char* file_name,size_t empty_space,double total_prob){
	ofstream coal_out_file;
	coal_out_file.open (file_name, ios::out | ios::app | ios::binary); 
	//coal_out_file<<setw(12)<<"Total"<<setw(empty_space)<<total_prob<<endl;	
	coal_out_file.close();
	cout<<setw(12)<<"Total"<<setw(empty_space)<<total_prob<<endl;	
	ostringstream total_prob_strm;
	total_prob_strm<<total_prob;
	string appending_log_str="Total probability: "+total_prob_strm.str();
	appending_log_file(appending_log_str);
	appending_log_str=file_name;
	appending_log_str="Gene tree probabilities produced in file: "+appending_log_str;
	appending_log_file(appending_log_str);

}

