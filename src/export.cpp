/* 
 * hybrid-coal is used to compute gene tree probabilities given species network under coalescent process.
 * 
 * Copyright (C) 2010 -- 2015 Sha (Joe) Zhu
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
    ////string     empty_str="";
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
