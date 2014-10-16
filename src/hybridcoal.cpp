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
HybridCoal::init(){}


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
    cout<<"-h -- Help. List the following content."<<endl;
	cout<<"-gt INPUT -- Input the gene tree string through command line or a file."<<endl;
	cout<<"-sp INPUT -- Input the species network/tree string through command line or a file."<<endl;
	cout<<"-maple/-mapleF -- Generate a Maple executeable script file to calculate the gene tree probabilities given species networks."<<endl;
	cout<<"-latex/-latexF -- Generate the coalescent history of a gene tree within a species network."<<endl;
	cout<<"-symb -- To enable the Maple script calculate the symbolic gene tree probabilities."<<endl;
	//cout<<"-print -- Prints out the node content "<<endl;
	cout<<"-gtopo -- To generate the gene tree topologies of a given set of taxa."<<endl;
	//cout<<"-sub -- Produce file sub_networks, which contanis all sub trees of a species network"<<endl;

	cout<<setw(20)<<"-plot/-dot [option]"<<"  --  "<<"Use LaTEX(-plot) or Dot (-dot) to draw the input (defined by -spcu) network(tree)."<<endl;
	cout<<setw(20)<<"      -branch"<<"  --  "<<"Branch lengths will be labelled in the figure."<<endl;

	cout<<setw(20)<<"-o STR "<<"  --  "<<"Specify the file name prefix of the output."<<endl;

	//cout<<setw(20)<<"-h or -help"<<"  --  "<<"Help. List the following content."<<endl;
	//cout<<setw(20)<<"-spcu STR"<<"  --  "<<"Input the species network/tree string through command line or a file."<<endl;
	//cout<<setw(26)<<" "<<"Branch lengths of the INPUT are in coalescent unit."<<endl;
	//cout<<setw(20)<<"-spng STR"<<"  --  "<<"Input the species network/tree string through command line or a file. "<<endl;
	//cout<<setw(26)<<" "<<"Branch lengths of the INPUT are in number of generation."<<endl;
	//cout<<setw(20)<<"-pop STR/FLT"<<"  --  "<<"Population sizes are defined by a single numerical constant, "<<endl;
	//cout<<setw(26)<<" "<<"or a string which specifies the population size on each branch. "<<endl;
	//cout<<setw(26)<<" "<<"The string can be input through command line or a file. "<<endl;
	//cout<<setw(26)<<" "<<"By default, population size 10,000 is used."<<endl;
}


/*! \brief hybrid-Lambda help file*/
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

