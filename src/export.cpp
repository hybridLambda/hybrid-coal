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



/*! \brief Core function of drawing a network in .tex files. 
 */
void plot_in_latex(const char* file_name /*! Name for the figure file */ , 
	Net* net_dummy,
// string net_str /*! Input network written in extended newick form */,
	int plot_option /*! '0' -- do not label the branches; '1' -- enumerate the branches by  postorder tree traversal; '2' -- label the branch lengths of the network*/
	){
	//Net net_dummy(net_str);
	ofstream latex_file;
	latex_file.open (file_name, ios::out | ios::app | ios::binary); 
	latex_file <<"\\begin{center}\n";	
	latex_file <<"\\begin{tikzpicture}[thick]\n";
	valarray <int>  x_node=det_x_node (net_dummy);
	for (size_t node_i=0;node_i<net_dummy->Net_nodes.size();node_i++){
		Node * current_node =net_dummy->Net_nodes[node_i];
		string sp_node_label=net_dummy->net_str.substr(current_node->nodeName_start(), current_node->nodeName_length());
		sp_node_label=rm_and_hash_sign(sp_node_label);
		if (current_node->tip()){
			latex_file<<"\\node at ("<<x_node[node_i]<<"\\du,"<<current_node->rank()<<"\\du) [circle,draw] ("<<sp_node_label <<") {$"<<sp_node_label <<"$};\n";
		//latex_file<<"\\node at ("<<x_node[node_i]<<"\\du,"<<y_node[node_i]<<"\\du) [circle,draw] ("<<sp_node_label <<") {$"<<sp_node_label <<"$};\n";
		}
		else{
			latex_file<<"\\node at ("<<x_node[node_i]<<"\\du,"<<current_node->rank()<<"\\du) [circle,fill=orange,draw] ("<<sp_node_label <<") {$"<<sp_node_label <<"$};\n";
		}
	}

	for (size_t node_i=0;node_i<net_dummy->Net_nodes.size()-1;node_i++){
		Node * current_node =net_dummy->Net_nodes[node_i];
		string sp_node_label=net_dummy->net_str.substr(current_node->nodeName_start(), current_node->nodeName_length());
		//string sp_node_label=net_dummy->Net_nodes[node_i].label;
		sp_node_label=rm_and_hash_sign(sp_node_label);
		Node * parent1 = current_node->parent1();
		string sp_node_parent1_label=net_dummy->net_str.substr(parent1->nodeName_start(), parent1->nodeName_length());
		sp_node_parent1_label=rm_and_hash_sign(sp_node_parent1_label);
		if (!net_dummy->Net_nodes[node_i]->tip()){
			if (plot_option==1){
				latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent1_label<<") node [midway,left]{"<< current_node->e_num() <<"};\n";
			}
			else{
				if (plot_option==2){
					latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent1_label<<") node [midway,left]{"<< current_node->brchlen1() <<"};\n";
				}
				else{
					latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent1_label<<");\n";	
				}	
			}
			if (current_node->parent2()){
				Node * parent2 = current_node->parent2();
				string sp_node_parent2_label=net_dummy->net_str.substr(parent2->nodeName_start(), parent2->nodeName_length());
				//string sp_node_parent2_label=net_dummy->Net_nodes[node_i].parent2->label;	
				sp_node_parent2_label=rm_and_hash_sign(sp_node_parent2_label);
				if (plot_option==1){
					latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent2_label<<") node [midway,left]{"<< current_node->e_num2() <<"};\n";
				}
				else{
					if (plot_option==2){
						latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent2_label<<") node [midway,left]{"<< current_node->brchlen2() <<"};\n";
					}
					else{
						latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent2_label<<");\n";
					}		
				}
			}
		}
		else{
			if (plot_option==2){
				latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent1_label<<") node [midway,left]{"<< current_node->brchlen1() <<"};\n";
			}
			else{
				latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent1_label<<");\n";	
			}	
			if (current_node->parent2()){
				Node * parent2 = current_node->parent2();
				string sp_node_parent2_label=net_dummy->net_str.substr(parent2->nodeName_start(), parent2->nodeName_length());
				//string sp_node_parent2_label=net_dummy->Net_nodes[node_i].parent2->label;
				sp_node_parent2_label=rm_and_hash_sign(sp_node_parent2_label);
				if (plot_option==2){
					latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent2_label<<") node [midway,left]{"<< current_node->brchlen2() <<"};\n";
				}
				else{
					latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent2_label<<");\n";
				}		
			}
		}
	}
	latex_file <<"\\end{tikzpicture}\n\n";
	latex_file <<"\\end{center}\n";
	latex_file.close();	
}


/*! \brief Produce a dot file, which is used to draw the network, and compile the dot file to a pdf file.
 */
void plot_in_dot(const char* file_name /*! Name for the figure file */,
	Net* net_dummy,
// string net_str /*! Input network written in extended newick form */,
	int plot_option /*! '0' -- do not label the branches; '1' -- enumerate the branches by  postorder tree traversal; '2' -- label the branch lengths of the network*/){
	string file_name_no_dot=Fname_no_ext(file_name,".dot");
	string file_name_with_dot=Fname_ext(file_name,".dot");
	ofstream dot_file;
	check_and_remove(file_name_with_dot.c_str());
	dot_file.open (file_name_with_dot.c_str(), ios::out | ios::app | ios::binary); 
			
	dot_file <<"graph G {\n rankdir=BT; ratio=compress;\n";//page="14,14"; determines the size of the ps output
	//valarray <int>  x_node=det_x_node (net_dummy);

	for (size_t node_i=0;node_i<net_dummy->Net_nodes.size()-1;node_i++){
		Node * current_node =net_dummy->Net_nodes[node_i];
		string sp_node_label=net_dummy->net_str.substr(current_node->nodeName_start(), current_node->nodeName_length());
		//string sp_node_label=net_dummy->Net_nodes[node_i].label;
		sp_node_label=rm_and_hash_sign(sp_node_label);
		//string sp_node_parent1_label=net_dummy->Net_nodes[node_i].parent1->label;
		Node * parent1 = current_node->parent1();
		string sp_node_parent1_label=net_dummy->net_str.substr(parent1->nodeName_start(), parent1->nodeName_length());
		sp_node_parent1_label=rm_and_hash_sign(sp_node_parent1_label);
		if (!current_node->tip()){		
			if (plot_option==1){
				dot_file<<sp_node_label<<" -- "<<sp_node_parent1_label<<"[label=\""<< current_node->e_num() <<"\"];\n";
			}
			else{
				if (plot_option==2){
					dot_file<<sp_node_label<<" -- "<<sp_node_parent1_label<<"[label=\""<< current_node->brchlen1() <<"\"];\n";	
				}
				else{
					dot_file<<sp_node_label<<" -- "<<sp_node_parent1_label<<";\n";//
				}	
			}
			if (current_node->parent2()){
				Node * parent2 = current_node->parent2();
				string sp_node_parent2_label=net_dummy->net_str.substr(parent2->nodeName_start(), parent2->nodeName_length());
				//string sp_node_parent2_label=net_dummy->Net_nodes[node_i].parent2->label;
				sp_node_parent2_label=rm_and_hash_sign(sp_node_parent2_label);
				if (plot_option==1){
					dot_file<<sp_node_label<<" -- "<<sp_node_parent2_label<<"[label=\""<< current_node->e_num2() <<"\"];\n";
				}
				else{
					if (plot_option==2){
						dot_file<<sp_node_label<<" -- "<<sp_node_parent2_label<<"[label=\""<< current_node->brchlen2() <<"\"];\n";	
					}
					else{
						dot_file<<sp_node_label<<" -- "<<sp_node_parent2_label<<";\n";//<<"[label=\""<< net_dummy->Net_nodes[node_i].e_num2 <<"\"];\n";
					}	
				}	
				
			}
		}
		else{
			if (plot_option==2){
				dot_file<<sp_node_label<<" -- "<<sp_node_parent1_label<<"[label=\""<< current_node->brchlen1() <<"\"];\n";	
			}
			else{
				dot_file<<sp_node_label<<" -- "<<sp_node_parent1_label<<";\n";//
			}	
			
			if (current_node->parent2()){
				Node * parent2 = current_node->parent2();
				string sp_node_parent2_label=net_dummy->net_str.substr(parent2->nodeName_start(), parent2->nodeName_length());
				//string sp_node_parent2_label=net_dummy->Net_nodes[node_i].parent2->label;
				sp_node_parent2_label=rm_and_hash_sign(sp_node_parent2_label);					
				if (plot_option==2){
					dot_file<<sp_node_label<<" -- "<<sp_node_parent2_label<<"[label=\""<< current_node->brchlen2() <<"\"];\n";	
				}
				else{
					dot_file<<sp_node_label<<" -- "<<sp_node_parent2_label<<";\n";//<<"[label=\""<< net_dummy->Net_nodes[node_i].e_num2 <<"\"];\n";
				}
			}
		}
	}
	//cout<<"hea"<<endl;
	//if (net_dummy->is_ultrametric_func()){
		for (int rank_i=net_dummy->Net_nodes.back()->rank();rank_i>0;rank_i--){
			dot_file<<"{ rank=same; ";
			vector <int> x_node_dummy;
			vector <size_t> x_node_dummy_index;
			for (size_t node_i=0;node_i<net_dummy->Net_nodes.size();node_i++){
				Node * current_node =net_dummy->Net_nodes[node_i];
				if (current_node->rank()==rank_i){
					
					string sp_node_label=net_dummy->net_str.substr(current_node->nodeName_start(), current_node->nodeName_length());
					//string sp_node_label=net_dummy->Net_nodes[node_i].label;
					sp_node_label=rm_and_hash_sign(sp_node_label);
					dot_file<<sp_node_label<<" ";
				}	
			}
			dot_file<<"} ;\n";
		}
	//}
	
	dot_file <<"}\n";
	dot_file.close();

	string command="dot -Tps "+file_name_no_dot+".dot -o "+file_name_no_dot+".ps";
	int sys=system(command.c_str());
	command="convert "+file_name_no_dot+".ps -resize 100\% "+file_name_no_dot+".pdf";
	//cout<<command<<endl;
	sys=system(command.c_str());

	string appending_log_str="Dot figure generated in file: "+file_name_no_dot+".pdf";
	appending_log_file(appending_log_str);
}


string Fname_ext(const char* file_name,const char* ext){
	string in_str(file_name);
	string ext_str(ext);
	size_t found=in_str.find(ext_str);
	if (found!=string::npos){
		return in_str;
	}
	else{
		return in_str+ext_str;
	}
}

string Fname_no_ext(const char* file_name,const char* ext){
	string in_str(file_name);
	string ext_str(ext);
	size_t found=in_str.find(ext_str);
	if (found!=string::npos){
		return in_str.substr(0,found);
	}
	else{
		return in_str;
	}
}


/*! \brief Produce a tex file, which is used to draw the network 
 */
void plot_in_latex_file(const char* file_name /*! Name for the figure file */ , 
	Net* net_dummy,
// string net_str /*! Input network written in extended newick form */,
	int plot_option /*! '0' -- do not label the branches; '1' -- enumerate the branches by  postorder tree traversal; '2' -- label the branch lengths of the network*/ 
	){
	ofstream latex_file;
	string file_name_no_dot=Fname_no_ext(file_name,".tex");
	string file_name_with_dot=Fname_ext(file_name,".tex");
	latex_header(file_name_with_dot.c_str());
	plot_in_latex(file_name_with_dot.c_str(), net_dummy,plot_option);	
	latex_tail(file_name_with_dot.c_str());
	string command="pdflatex "+file_name_with_dot;
	int sys=system(command.c_str());
	string appending_log_str="Network figure generated in file: "+file_name_no_dot+".pdf";
	appending_log_file(appending_log_str);

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


/*! \brief Set plot options */
/*! \return int '0' -- do not label the branches; '1' -- enumerate the branches by  postorder tree traversal; '2' -- label the branch lengths of the network*/ 
int set_plot_option(bool plot_label,bool plot_branch){
	int plot_option=0;
	if (plot_label){
		plot_option=1;
		appending_log_file("Internal branches are labelled by post-order tree traversal.");
	}
	else{
		if (plot_branch){
			plot_option=2;
			appending_log_file("Branch lengths are labelled.");
		}
	}
	
	//int plot_option;
	//if (plot_label){
		//plot_option=1;
		//appending_log_file("Internal branches are labelled by post-order tree traversal.");
	//}
	//else if (plot_branch){
			//plot_option=2;
			//appending_log_file("Branch lengths are labelled.");
		//}
		
	//else{
		//plot_option=0;
	//}
	return plot_option;
}

/*! \brief When drawing network in .tex files, detemine the x coordinates of nodes \todo improve this!
 */
valarray <int>  det_x_node (Net* net_dummy){
	valarray <int>  x_node (net_dummy->Net_nodes.size());
	x_node[x_node.size()-1]=0;
	
	//valarray <int>  y_node (net_dummy->Net_nodes.size());
	//y_node[x_node.size()-1]=net_dummy->Net_nodes.back().e_num;
	
	for (int rank_i=net_dummy->Net_nodes.back()->rank();rank_i>0;rank_i--){
		vector <int> x_node_dummy;
		vector <size_t> x_node_dummy_index;
		for (size_t node_i=0;node_i<net_dummy->Net_nodes.size();node_i++){
			Node * current_node = net_dummy->Net_nodes[node_i];
			if (current_node->rank()==rank_i){
				size_t n_child=current_node->child.size();
				int parent_x=x_node[node_i];
				//int parent_y=y_node[node_i];
				int start_child_x=parent_x-floor(n_child/2);
				//int start_child_x=0;
				bool odd_num_child=false;
				if ((n_child % 2) == 1){
					odd_num_child=true;
				}
				if (odd_num_child){
					for (size_t child_i=0;child_i<n_child;child_i++){
						for (size_t node_j=0;node_j<net_dummy->Net_nodes.size();node_j++){
							if (net_dummy->Net_nodes[node_j]==current_node->child[child_i]){			
								//int child_y=y_node[node_j];
								if (start_child_x==parent_x){
									x_node[node_j]=parent_x;
									//y_node[node_j]=parent_y-1;
									start_child_x++;
								}
								else{
									//x_node[node_j]=start_child_x*(parent_y-child_y)+parent_x;
									x_node[node_j]=start_child_x;
									//y_node[node_j]=parent_y-1;
								}
								start_child_x++;
							}
						}
					}
				}
				else{
					for (size_t child_i=0;child_i<n_child;child_i++){
						for (size_t node_j=0;node_j<net_dummy->Net_nodes.size();node_j++){
							if (net_dummy->Net_nodes[node_j] == current_node->child[child_i]){
								if (start_child_x==parent_x){										
									start_child_x++;
								}
								//	x_node[node_j]=start_child_x*(parent_y-child_y)+parent_x;
								x_node[node_j]=start_child_x;
								//y_node[node_j]=parent_y-1;
								start_child_x++;
							}
						}
					}
					
				}
				x_node_dummy.push_back(x_node[node_i]);
				x_node_dummy_index.push_back(node_i);
			}
		}
		if (x_node_dummy.size()>1){
			bool need_to_shift=true;
			while (need_to_shift){
				for (size_t x_node_dummy_i=0;x_node_dummy_i<x_node_dummy.size();x_node_dummy_i++){
					int current_x_node_dummy=x_node_dummy[x_node_dummy_i];
					for (size_t x_node_dummy_j=x_node_dummy_i+1;x_node_dummy_j<x_node_dummy.size();x_node_dummy_j++){
						if (current_x_node_dummy==x_node_dummy[x_node_dummy_j]){
							if (x_node_dummy[x_node_dummy_j]>0){
								x_node_dummy[x_node_dummy_j]++;
								x_node[x_node_dummy_index[x_node_dummy_j]]++;
							}
							else{
								x_node_dummy[x_node_dummy_j]--;
								x_node[x_node_dummy_index[x_node_dummy_j]]--;
							}
						}
					}
				}
				need_to_shift=false;
				for (size_t x_node_dummy_i=0;x_node_dummy_i<x_node_dummy.size();x_node_dummy_i++){
					int current_x_node_dummy=x_node_dummy[x_node_dummy_i];
					for (size_t x_node_dummy_j=x_node_dummy_i+1;x_node_dummy_j<x_node_dummy.size();x_node_dummy_j++){
						if (current_x_node_dummy==x_node_dummy[x_node_dummy_j]){
							need_to_shift=true;
							break;
						}		
					}
					if (need_to_shift){
						break;
					}
				}
			}
		}
	}
	return x_node;
	
}



void print_help(){
	cout<<"*****************************************************************"<<endl;
	cout<<"*			Hybrid_coal beta 0.1			*"<<endl;
	cout<<"*			  Author: Joe ZHu			*"<<endl;
	cout<<"*****************************************************************"<<endl;
	cout<<endl<<endl;

	cout<<"-h -- Help. List the following content."<<endl;
	cout<<"-gt INPUT -- Input the gene tree string through command line or a file."<<endl;
	cout<<"-sp INPUT -- Input the species network/tree string through command line or a file."<<endl;
	cout<<"-maple/-mapleF -- Generate a Maple executeable script file to calculate the gene tree probabilities given species networks."<<endl;
	cout<<"-latex/-latexF -- Generate the coalescent history of a gene tree within a species network."<<endl;
	cout<<"-symb -- To enable the Maple script calculate the symbolic gene tree probabilities."<<endl;
	cout<<"-plot/-dot [option] Use LaTEX(-plot) or Dot (-dot) to draw the input (defined by -sp) network(tree)."<<endl;
	cout<<"-plotF/-dotF FILE -- Generated figure will be saved in FILE"<<endl;
	//cout<<"-print -- Prints out the node content "<<endl;
	cout<<"-gtopo -- To generate the gene tree topologies of a given set of taxa."<<endl;
	//cout<<"-sub -- Produce file sub_networks, which contanis all sub trees of a species network"<<endl;
	cout<<"-out -- Specify the file name of the output probability"<<endl;

	cout<<endl;
	cout<<"Examples:"<<endl;
	cout<<"hybrid_coal -sp '((((B:1,C:1)s1:1)h1#.5:1,A:3)s2:1,(h1#.5:1,D:3)s3:1)r;'"<<endl;	
	cout<<"hybrid_coal -sp trees/4_tax_sp_nt1_para -gt '(((A,D),C),B)'"<<endl;	
	cout<<"hybrid_coal -sp trees/4_tax_sp_nt1_para -gt trees/4_tax_gt4.tre -maple"<<endl;	
	cout<<"hybrid_coal -sp trees/4_tax_sp_nt1_para -gt trees/4_tax_gt4.tre -latex"<<endl;
	cout<<"hybrid_coal -sp trees/4_tax_sp_nt1_para -maple -latex"<<endl;	
	cout<<"hybrid_coal -sp trees/4_tax_sp_nt1_para -plot"<<endl;
	cout<<"hybrid_coal -sp trees/4_tax_sp_nt1_para -plot -branch"<<endl;
	cout<<"hybrid_coal -sp trees/4_tax_sp_nt1_para -plot -label"<<endl;
	cout<<"hybrid_coal -sp trees/4_tax_sp_nt1_para -dot"<<endl;
	cout<<"hybrid_coal -sp trees/4_tax_sp_nt1_para -print"<<endl;
	cout<<"hybrid_coal -sp trees/4_tax_sp_nt1_para -gtopo"<<endl;
	cout<<"hybrid_coal -sp trees/4_tax_sp_nt1_para -sub"<<endl;
	cout<<endl;	
	
}
