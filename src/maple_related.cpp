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
            //maple_file<<"lambda["<<net_dummy->Net_nodes[i].edge<<"]:="<<net_dummy->Net_nodes[i].brchlen1<<";"<<endl;
            //if (net_dummy->Net_nodes[i].hybrid){
                //maple_file<<"lambda["<<net_dummy->Net_nodes[i].edge2<<"]:="<<net_dummy->Net_nodes[i].brchlen2<<";"<<endl;        
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



void simplify_Networks_maple_block(const char* file_name,string sp_str,string gt_tree_str,int gt_str_i){
    bool maple_bool_local=true;
    if (rm_debug_bool){
        string bug_statement="Start maple block for  "+sp_str;
    }
    vec_Net_wiz_prior_p my_rmd_networks=build_initial_Networks(sp_str);

    ostringstream gt_i_strm;
    gt_i_strm<<gt_str_i;
    vector <string> maple_block_vec_str;
    for (size_t i=0;i<my_rmd_networks.Net_vec.size();i++){
        string maple_file_rv_vec_dummy;
        int old_my_rmd_network_size=my_rmd_networks.Net_vec.size(); 
        ostringstream i_str;
        i_str<<i;
        Net current_net(my_rmd_networks.Net_vec[i].s_net_string);    
        int rm_node_index=choose_rm_node(current_net);
        
        if (rm_node_index>=0){
            my_rmd_networks=simplify_Networks_one_block(current_net.Net_nodes[rm_node_index].hybrid,rm_node_index,my_rmd_networks,i,maple_bool_local,gt_tree_str);
            int new_block_size=my_rmd_networks.Net_vec.size()- old_my_rmd_network_size;
            maple_block_vec_str=build_maple_block_sp_part(maple_block_vec_str, old_my_rmd_network_size,new_block_size,gt_i_strm.str(),i_str.str());
            maple_block_vec_str=build_maple_block_omega_part(maple_block_vec_str,old_my_rmd_network_size,new_block_size,my_rmd_networks);

        }
    }
    
    if (gt_tree_str.size()>0){
        maple_block_vec_str=build_maple_block_gt_part(maple_block_vec_str,my_rmd_networks,gt_i_strm.str(),gt_tree_str);
    }

    write_maple_block(file_name,sp_str,maple_block_vec_str);
    
    maple_block_vec_str.clear();
}


vector <string> build_maple_block_sp_part(vector <string> maple_block_vec_str, int old_my_rmd_network_size,int new_block_size,string gt_i,string sp_i){
    string maple_file_rv_vec_dummy="p(g["+gt_i+"],s["+sp_i+"]):=";
    for (int i_sp_vec=0;i_sp_vec<new_block_size;i_sp_vec++ ){
        ostringstream net_counter_str;
        net_counter_str<<old_my_rmd_network_size+i_sp_vec;
        maple_file_rv_vec_dummy=maple_file_rv_vec_dummy+"omega["+net_counter_str.str()+"]*p(g["+gt_i+"],s["+net_counter_str.str()+"])";
        if (i_sp_vec<(new_block_size-1)){
            maple_file_rv_vec_dummy=maple_file_rv_vec_dummy+"+";
        }
    }
    if ( maple_file_rv_vec_dummy[maple_file_rv_vec_dummy.size()-1]== '='){
        maple_file_rv_vec_dummy=maple_file_rv_vec_dummy+"0;";
    }
    else{
        maple_file_rv_vec_dummy=maple_file_rv_vec_dummy+";";
    }
    maple_block_vec_str.push_back(maple_file_rv_vec_dummy);    
    if (rm_debug_bool){
        cout<<"End of build_maple_block_sp_part"<<endl;
    }
    return maple_block_vec_str;
    
}

vector <string> build_maple_block_omega_part(vector <string> maple_block_vec_str,int old_my_rmd_network_size,int new_block_size,vec_Net_wiz_prior_p my_rmd_networks){
    for (int i_sp_vec=0;i_sp_vec<new_block_size;i_sp_vec++){
        int net_counter=old_my_rmd_network_size+i_sp_vec;
        ostringstream net_counter_str;
        net_counter_str<<net_counter;
        string maple_file_rv_vec_dummy="omega["+net_counter_str.str()+"]:="+my_rmd_networks.Net_vec[net_counter].omega_string+";"+"#s["+net_counter_str.str()+"]="+my_rmd_networks.Net_vec[net_counter].s_net_string;
        maple_block_vec_str.push_back(maple_file_rv_vec_dummy);
    }
    if (rm_debug_bool){
        cout<<"End of build_maple_block_omega_part"<<endl;
    }
    return maple_block_vec_str;
}

vector <string> build_maple_block_gt_part(vector <string> maple_block_vec_str,vec_Net_wiz_prior_p my_rmd_networks,string gt_i,string gt_tree_str){
    for (size_t i=0;i<my_rmd_networks.Net_vec.size();i++){
        Net current_net(my_rmd_networks.Net_vec[i].s_net_string);
        if (!current_net.is_net){
            ostringstream i_str;
            i_str<<i;
            gt_coal_in_st trial1(gt_tree_str,my_rmd_networks.Net_vec[i]);
            string current_p_gs="p(g["+gt_i+"],s["+i_str.str()+"]):="+trial1.prob_of_maple;
            maple_block_vec_str.push_back(current_p_gs);
        }
    }    
    if (rm_debug_bool){    
        cout<<"End of build_maple_block_gt_part"<<endl;
    }
    return maple_block_vec_str;
}


void write_maple_block(const char* file_name,string sp_str,vector <string> maple_block_vec_str){
    ofstream maple_file;
    maple_file.open(file_name, ios::out | ios::app | ios::binary);    
    for (int rit=int(maple_block_vec_str.size()-1) ; rit >=0; rit-- ){
        maple_file<<maple_block_vec_str[rit]<<endl;        
    }
    maple_file<<"#s[0]="<< sp_str<<endl;
    maple_file.close();
    if (rm_debug_bool){
        cout<<"End maple block for  "<<sp_str<<endl;
    }        
}
