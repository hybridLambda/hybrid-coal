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

#include "hybridcoal.hpp"
#include "plot/figure.hpp"

HybridCoal::init(){}

void HybridCoal::parse(){
    
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


    
}


void HybridCoal::finalize(){
    
    this->parameters()->finalize( );    
   
    if ( this->print_tree_bool ) 
        this->print();
   
    if ( this->plot_bool ){
        Figure figure_para ( this->argc_, this->argv_ );
        figure_para.figure_file_prefix = this->prefix;
        figure_para.finalize();
        figure_para.plot( this->parameters()->net_str );
        exit(EXIT_SUCCESS);
    }
}

void HybridCoal::print(){
    Net net( this->parameters_->net_str );
    net.print_all_node();
    exit(EXIT_SUCCESS);
}

void print_example(){
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


void print_option(){
	cout<<setw(20)<<"-h or -help"<<"  --  "<<"Help. List the following content."<<endl;
	cout<<setw(20)<<"-gt STR"<<"  --  "<<"Input the gene tree string string through command line or a file."<<endl;
	cout<<setw(20)<<"-sp STR"<<"  --  "<<"Input the species network/tree string through command line or a file."<<endl;
    cout<<setw(20)<<"-maple [option]"<<"  --  "<<"Generate a Maple executeable script file to calculate the gene tree probabilities given species networks."<<endl;
    cout<<setw(20)<<"[-symb]"<<"  --  "<<"To enable the Maple script calculate the symbolic gene tree probabilities."<<endl;
	cout<<setw(20)<<"-latex"<<"  -- Generate the coalescent history of a gene tree within a species network."<<endl;
	cout<<setw(20)<<"-gtopo"<<"  -- To generate the gene tree topologies of a given set of taxa."<<endl;
	//cout<<"-sub -- Produce file sub_networks, which contanis all sub trees of a species network"<<endl;

	cout<<setw(20)<<"-plot/-dot [option]"<<"  --  "<<"Use LaTEX(-plot) or Dot (-dot) to draw the input (defined by -spcu) network(tree)."<<endl;
	cout<<setw(20)<<"    [-branch]"<<"  --  "<<"Branch lengths will be labelled in the figure."<<endl;
	cout<<setw(20)<<"-o STR "<<"  --  "<<"Specify the file name prefix of the output."<<endl;
}


/*! \brief hybrid-coal help file*/
void print_help(){
	cout<<endl;
	cout<<endl;
	cout<<"*****************************************************************"<<endl;
	cout<<"*                      hybrid-coal beta 0.1                     *"<<endl;
	cout<<"*                         Author: Joe ZHU                       *"<<endl;
	cout<<"*****************************************************************"<<endl;
	cout<<endl<<endl;
	print_option();
	print_example();
    exit (EXIT_SUCCESS);
}

