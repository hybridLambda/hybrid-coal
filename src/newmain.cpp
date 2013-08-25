#include"net.hpp"
#include"node.hpp"
#include"utility.hpp"

int main(int argc, char *argv[]){
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
	
	//maple_bool=false;
	//string maple_F_name="maple_prob.mw";
	
	bool print_gene_topo_bool=false;
	string gtopo_F_name="GENE_topo";
	
	bool help=false;
	string net_str;
	vector <string> gt_tree_str_s;
	
	bool all_gt_tree_bool=true;
	bool symb_bool=false;
	bool print_tree=false;
		
	bool list_sub_network_bool=false;

	//int plot_option;//=0;
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
		//if (argv_i=="-latex"){
			//latex_bool=true;			
		//}
		
		//if (argv_i=="-latexF"){
			//latex_bool=true;
			//latex_F_name=Fname_ext(argv[argc_i+1],".tex");
		//}
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
		//if (argv_i=="-maple"){
			//maple_bool=true;
		//}
		//if (argv_i=="-mapleF"){
			//maple_bool=true;
			//maple_F_name=Fname_ext(argv[argc_i+1],".mw");
		//}
		
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
	if (argc==1){net_str="((((B:1,C:1)s1:1)h1#.5:1,A:3)s2:1,(h1#.5:1,D:3)s3:1)r;";}
	
	
	Net* net_dummy = new Net(net_str);
		if (print_tree){
			net_dummy->print_all_node();
			appending_log_file("Tree printed");
		}
		//plot_option=set_plot_option(plot_label,plot_branch);

		//if (plot_bool){
			//plot_in_latex_file(tex_fig_name.c_str(), net_dummy,plot_option);	
		//}
		
		//if (dot_bool){
			//plot_in_dot(dot_fig_name.c_str(), net_dummy,plot_option);			
		//}
			
		//if (all_gt_tree_bool){
			//gt_tree_str_s=all_n_tax_gene_tree(net_dummy.tip_name);
			//if (print_gene_topo_bool){
				//print_all_gt_topo(gtopo_F_name.c_str(),gt_tree_str_s);
			//}
		//}

		//if (list_sub_network_bool && gt_tree_str_s.size()==1){
			//list_sub(net_str,gt_tree_str_s[0]);
		//}

		//if (list_sub_network_bool && all_gt_tree_bool){
			//string 	empty_str="";
			//list_sub(net_str,empty_str);
		//}
				
		//if (print_tree || plot_bool || dot_bool || print_gene_topo_bool || list_sub_network_bool && gt_tree_str_s.size()>1){
			////return 0;
			//return my_exit();
		//}		

}
