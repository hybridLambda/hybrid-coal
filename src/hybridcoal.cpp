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

#include "hybridcoal.hpp"
#include "plot/figure.hpp"

void HybridCoal::init(){
    this->prefix    = "OUT";    
    this->symb_bool = false;
    this->print_tree_bool = false;
    this->print_gene_topo_bool = false; 
    this->plot_bool       = false;
    this->latex_bool = false;
    this->argc_i = 1;
    this->maple_bool = false;
    this->all_gt_tree_bool = true;
}

void HybridCoal::read_sp_str( string & argv_i ){
    readNextStringto( this->tmp_input_str , this->argc_i, this->argc_,  this->argv_ );
    this->sp_str = read_input_line( tmp_input_str.c_str() );
}

void HybridCoal::read_input_lines(const char inchar[], vector <string> & out_vec){
    ifstream in_file( inchar );
    string out_str;
    if ( in_file.good() ){
        getline ( in_file, out_str );
        while ( out_str.size() > 0 ){   
            out_vec.push_back( out_str );
            getline ( in_file, out_str );
        }
    } else {
        string dummy_str(inchar);
        if (dummy_str.find('(')!=string::npos && dummy_str.find(')')!=string::npos){
            out_str=dummy_str;
            out_vec.push_back(out_str);
        } else{
            throw std::invalid_argument("Invalid input file. " + string (inchar) );
        }
    }
    in_file.close();
}

void HybridCoal::check_gt_str_and_sign( string &gt_str ){
	if ( gt_str.find('&') != string::npos ){
        throw std::invalid_argument ( "Error: gene tree " + gt_str + " can not be a gene tree string" );
	}	
}



void HybridCoal::parse(){
    while (argc_i < argc_){	
        std::string argv_i(argv_[argc_i]);
        if ( argv_i == "-h" || argv_i == "-help" ){ print_help(); }		
		else if (argv_i=="-gt"){
            readNextStringto( this->gt_file_name , this->argc_i, this->argc_,  this->argv_ );
		}
		else if (argv_i=="-latex"){  this->latex_bool=true;	}
        if ( argv_i == "-sp"){ this->read_sp_str(argv_i); }

        else if ( argv_i =="-o" ) readNextStringto( this->prefix , this->argc_i, this->argc_,  this->argv_ );

		// else if (argv_i=="-pdf"){
			//latex_bool=true;
			//pdf_bool=true;
		//}
		else if ( argv_i == "-maple" ){ this->maple_bool=true; }
		//if (argv_i=="-mapleF"){
			//maple_bool=true;
			//maple_F_name=Fname_ext(argv[argc_i+1],".mw");
		//}
		else if ( argv_i == "-symb" ){ this->symb_bool=true; }
        else if ( argv_i == "-plot" || argv_i == "-dot" ){ this->plot_bool = true; }
        else if ( argv_i == "-label" || argv_i == "-branch" ){ argc_i++; continue; }        
        else if ( argv_i == "-print" ){ this->print_tree_bool = true; }
		else if ( argv_i == "-gtopo" ){ this->print_gene_topo_bool = true; }
		//else if (argv_i=="-sub"){ list_sub_network_bool=true;}
		//else if (argv_i=="-acc"){ acc_bool=true;	}
        else { throw std::invalid_argument ( "Unknown flag:" + argv_i); }
        //else { cout <<"  need to change this !!!" << argv_i<<endl; argc_i++;continue; } // need to change this !!!
        argc_i++;	
	}
    this->finalize();
}


string HybridCoal::read_input_line(const char *inchar){
	string out_str;
    ifstream in_file( inchar );
	if (in_file.good())	getline ( in_file, out_str); 
	else{
		string dummy_str(inchar);
		if (dummy_str.find('(')!=string::npos && dummy_str.find(')')!=string::npos) out_str=dummy_str;
		else  throw std::invalid_argument("Invalid input file. " + string (inchar) );
	}
	in_file.close();			
    return 	out_str;
}


void HybridCoal::finalize(){
    if ( this->print_tree_bool ) {
        this->print();
        return;
    }
   
    if ( this->plot_bool ){
        Figure figure_para ( this->argc_, this->argv_ );
        figure_para.figure_file_prefix = this->prefix;
        figure_para.finalize();
        figure_para.plot( this->sp_str );
        //exit(EXIT_SUCCESS);
        return;
    }
    
    if ( this->gt_file_name.size() > 0 ) {
        this->read_input_lines( this->gt_file_name.c_str(), this->gt_tree_str_s);
    }
}

void HybridCoal_core(){
    
}

void HybridCoal::print(){
    Net net( this->sp_str );
    net.print_all_node();
    //exit(EXIT_SUCCESS);
}

void print_example(){
	cout<<"Examples:"<<endl;
	cout<<"hybrid-coal -sp '((((B:1,C:1)s1:1)h1#.5:1,A:3)s2:1,(h1#.5:1,D:3)s3:1)r;'"<<endl;	
	cout<<"hybrid-coal -sp trees/4_tax_sp_nt1_para -gt '(((A,D),C),B)'"<<endl;	
	cout<<"hybrid-coal -sp trees/4_tax_sp_nt1_para -gt trees/4_tax_gt4.tre -maple"<<endl;	
	cout<<"hybrid-coal -sp trees/4_tax_sp_nt1_para -gt trees/4_tax_gt4.tre -latex"<<endl;
	cout<<"hybrid-coal -sp trees/4_tax_sp_nt1_para -maple -latex"<<endl;	
	cout<<"hybrid-coal -sp trees/4_tax_sp_nt1_para -plot"<<endl;
	cout<<"hybrid-coal -sp trees/4_tax_sp_nt1_para -plot -branch"<<endl;
	cout<<"hybrid-coal -sp trees/4_tax_sp_nt1_para -plot -label"<<endl;
	cout<<"hybrid-coal -sp trees/4_tax_sp_nt1_para -dot"<<endl;
	cout<<"hybrid-coal -sp trees/4_tax_sp_nt1_para -print"<<endl;
	cout<<"hybrid-coal -sp trees/4_tax_sp_nt1_para -gtopo"<<endl;
	cout<<"hybrid-coal -sp trees/4_tax_sp_nt1_para -sub"<<endl;
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

