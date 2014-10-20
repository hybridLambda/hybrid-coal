/* 
 * hybrid-coal is used to compute gene tree probabilities given species network under coalescent process.
 * 
 * Copyright (C) 2010 -- 2014 Sha (Joe) Zhu
 * 
 * This file is part of hybrid-coal
 * 
 * hybrid-coal is free software: you can redistribute it and/or modify
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



