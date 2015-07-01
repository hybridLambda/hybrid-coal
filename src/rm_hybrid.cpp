/*
 * hybrid-coal is used to compute gene tree probabilities given species network under coalescent process.
 * 
 * Copyright (C) 2010, 2011 Sha (Joe) Zhu
 * 
 * This file is part of hybrid-coal.
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

/*! \file rm_hybrid.cpp
 * \brief Remove hybrid node and descndent node of a hybrid node */

#include"rm_hybrid.hpp"
bool rm_debug_bool=false;

Net sub_brchlen_by_enum(string old_net_string, vector < vector <int> > e_num_vec_dummy){
//string sub_brchlen_by_enum(string old_net_string, vector < vector <int> > e_num_vec_dummy){
    Net old_Net(old_net_string);
    if (e_num_vec_dummy.size()!=old_Net.Net_nodes.size()){
        cout<<"sub_brchlen_by_enum!!!!"<<old_Net.Net_nodes.size()<<"  "<<e_num_vec_dummy.size() <<endl;
    }
    vector <Node*> old_Net_node_ptr;
    for (size_t i=0;i<old_Net.Net_nodes.size();i++){
        Node* new_node_ptr=NULL;
        old_Net_node_ptr.push_back(new_node_ptr);
        old_Net_node_ptr[i]=&old_Net.Net_nodes[i];    
        old_Net_node_ptr[i]->brchlen1=e_num_vec_dummy[i][0];
        if (old_Net_node_ptr[i]->parent2){
            old_Net_node_ptr[i]->brchlen2=e_num_vec_dummy[i][1];
        }
    
    }    
    rewrite_node_content(old_Net_node_ptr);
    return old_Net;
}


void print_H_S_A_matrix(vector < vector < valarray < int > > > A_matrix){
    for (size_t A_matrix_i=0;A_matrix_i<A_matrix.size();A_matrix_i++){
        for (size_t A_matrix_i_i=0;A_matrix_i_i<A_matrix[A_matrix_i].size();A_matrix_i_i++){
            for (size_t A_matrix_i_i_i=0;A_matrix_i_i_i<A_matrix[A_matrix_i][A_matrix_i_i].size();A_matrix_i_i_i++){
                cout<<A_matrix[A_matrix_i][A_matrix_i_i][A_matrix_i_i_i];
            }
            cout<<" ";
        }
        cout<<endl;
    }
}


int choose_rm_node(Net current_Net){
    size_t first_rm=current_Net.Net_nodes.size()-1;
    for (size_t i=0;i<current_Net.Net_nodes.size();i++){
        if (((current_Net.Net_nodes[i].hybrid+current_Net.Net_nodes[i].descndnt_of_hybrid)*(1-current_Net.Net_nodes[i].tip_bool))>=1 && current_Net.Net_nodes[i].rank < current_Net.Net_nodes[first_rm].rank){ /*! \todo Check if a tip node should be removed or not */
            first_rm=i;
        }
    }
    if (first_rm==current_Net.Net_nodes.size()-1){
        return -1;}
    else{
        return first_rm;
    }
}


vector <int> disjoint_list_s(int n, valarray <int> A_i,int i,vector <valarray <int> >A){
    vector <int> disjoint_list_s_return;
    for (size_t j=pow(2.0,1.0*(n-1))-1;j<A.size();j++){
        valarray <int> A_j=A[j];
        int a_inter_b=0;
        for (int A_j_i=0;A_j_i<n;A_j_i++){
             if (A_i[A_j_i]!=0 && A_i[A_j_i]==A_j[A_j_i]){
                a_inter_b=1;
                break;
            }    
        }
        if (a_inter_b==0){
            disjoint_list_s_return.push_back(j);
        }
    }
    return  disjoint_list_s_return;
}

int disjoint_list_h(int n,int i,vector <valarray <int> >A){
     valarray <int> A_i=A[i];
    int disjoint_list_h_return;
    for (size_t j=0;j<pow(2.0,1.0*n);j++){
        valarray <int> A_j=A[j];
        valarray <int> A_ij_sum=A_i+A_j;
        int a_inter_b=0;
        for (int A_j_i=0;A_j_i<n;A_j_i++){
             if (A_i[A_j_i]!=0 && A_i[A_j_i]==A_j[A_j_i]){
                a_inter_b=1;
                break;
            }    
        }
        if (a_inter_b==0 && A_ij_sum.sum()==n){
            disjoint_list_h_return=j;
        }
    }
    return  disjoint_list_h_return;
}


vector < valarray <int> > all_possible_comb(int n){
    vector <int> Avec;
    for (int ii=0;ii<n;ii++){
        int i=pow(1.0*2,1.0*ii);
        while (i<=pow(1.0*2,1.0*n)){
        vector <int> ones(pow(1.0*2,1.0*ii),1);
        vector <int> zeros(pow(1.0*2,1.0*ii),0);
        Avec.insert(Avec.begin(),zeros.begin(),zeros.end());
        Avec.insert(Avec.begin(),ones.begin(),ones.end());
        i=i+pow(1.0*2,1.0*(ii+1));
        }
    }
    vector < valarray <int> > A;
    for (int ii=0;ii<(pow(1.0*2,1.0*n));ii++){
        int    A_row_array[n];
        int A_row_array_ind=0;
        for (size_t i=0;i<Avec.size();i++){
            if (fmod(i+1,pow(1.0*2,1.0*n))==ii){
                A_row_array[A_row_array_ind]=Avec[i];
                A_row_array_ind++;
            }
        }
        valarray <int> A_row(A_row_array,n);
        A.push_back(A_row);
    }
    return A;
}

vector < valarray <int> > rearrange_A(vector < valarray <int> > A,int n){
    vector <int> A_sum;
    for (int i=0; i<pow(1.0*2,1.0*(n-1));i++){
    A_sum.push_back(A[i].sum());
    }
    int numLength =  pow(1.0*2,1.0*(n-1)); 
    for (int i=0; i<numLength-1;i++){
        for(int j = (i+1); j < numLength; j++){
            if (A_sum[i] < A_sum[j]){
                int temp= A_sum[i];          // swap
                A_sum[i] = A_sum[j];
                A_sum[j] = temp;
                valarray<int> temp_valarray= A[i];          // swap
                A[i] = A[j];
                A[j] = temp_valarray;
            }
        }
    }    
    return A;
}




vector < vector < valarray <int> > > build_s_child(int n){
    vector < valarray <int> > A_old=all_possible_comb(n);
    vector < valarray <int> > A=rearrange_A(A_old,n);
    vector < vector <vector <int> > > all_list;
    vector <int> empty_vec_dim1;
    vector < vector <int> > empty_vec_dim2;
    empty_vec_dim2.push_back(empty_vec_dim1);
    all_list.push_back(empty_vec_dim2);
    for (size_t i=1; i<pow(1.0*2,1.0*(n-1))+1;i++){
        valarray <int> A_i=A[i];
        vector <int> disjoint_set_list=disjoint_list_s(n,A_i,i,A);
        vector < vector <int> > new_list;
        for (size_t disjoint_i=0;disjoint_i<disjoint_set_list.size();disjoint_i++){
            valarray <int> new_A=A_i+A[disjoint_set_list[disjoint_i]];
            for (size_t upto_i=0;upto_i<i;upto_i++){
                valarray<bool> comp = (new_A==A[upto_i]);
                if (comp.min() == true){
                    vector < vector <int> > current_list=all_list[upto_i];
                    for (size_t current_list_i=0;current_list_i<current_list.size();current_list_i++){
                        current_list[current_list_i].push_back(disjoint_set_list[disjoint_i]);
                        }
                    new_list.insert(new_list.begin(),current_list.begin(),current_list.end());
                }
            }
        }
        all_list.push_back(new_list);
    }
    vector < vector <int> > all_list_in_vec;
    for (size_t print_all_list_i=0;print_all_list_i<all_list.size();print_all_list_i++){
        for (size_t print_all_list_i_i=0;print_all_list_i_i<all_list[print_all_list_i].size();print_all_list_i_i++){
            all_list[print_all_list_i][print_all_list_i_i].push_back(print_all_list_i);
            for (size_t print_all_list_i_i_i=0;print_all_list_i_i_i<all_list[print_all_list_i][print_all_list_i_i].size();print_all_list_i_i_i++){
                all_list_in_vec.push_back(all_list[print_all_list_i][print_all_list_i_i]);
            }
        }
    }

    for (size_t i=0;i<all_list_in_vec.size();i++){
        sort(all_list_in_vec[i].begin(),all_list_in_vec[i].end());
    }

    vector < vector <int> > new_all_list;
    new_all_list.push_back(all_list_in_vec[0]);
    for (size_t i=0;i<all_list_in_vec.size();i++){
        int all_compare=0;
        for (size_t ii=0;ii<new_all_list.size();ii++){
            if (all_list_in_vec[i].size()==new_all_list[ii].size()){
                int compare=1;
                for (size_t iii=0;iii<all_list_in_vec[i].size();iii++){
                    if (all_list_in_vec[i][iii]!=new_all_list[ii][iii]){
                        compare=0;
                        break;
                    }
                }
                all_compare=all_compare+compare;
            }
        }
        if (all_compare==0){
            new_all_list.push_back(all_list_in_vec[i]);

        }

    }
    
    vector < vector < valarray <int> > > descdent_all_list;
    for (size_t i=0;i<new_all_list.size();i++){
        vector < valarray <int> > current_descent;
        for (size_t ii=0;ii<new_all_list[i].size();ii++){
            valarray <int> descent(0,n);
            for (int iii=0;iii<n;iii++){
                descent[iii]=A[new_all_list[i][ii]][iii];
            }
            current_descent.push_back(descent);
        }
        descdent_all_list.push_back(current_descent);
    }
        
    return descdent_all_list;
    
}



vector < vector < valarray < int > > > build_h_child(int n){
    vector < vector < valarray < int > > > descdent_all_list;
    vector < valarray < int > > A=all_possible_comb(n);
    for (size_t i=0; i<pow(1.0*2,1.0*n);i++){
        int i_comp=disjoint_list_h(n,i,A);
        vector < valarray < int > > current_descdent_list;
        valarray <int> current_descdent_list_i=A[i];
        valarray <int> current_descdent_list_i_comp=A[i_comp];
        current_descdent_list.push_back(current_descdent_list_i);
        current_descdent_list.push_back(current_descdent_list_i_comp);
        descdent_all_list.push_back(current_descdent_list);
    }

    return descdent_all_list;
}


void rm_zero_kids_hybrid_node(vector <Node*> sp_nodes_ptr_rm){
    //check hybrid node has zero kids
    for (size_t sp_nodes_ptr_rm_i=0;sp_nodes_ptr_rm_i<sp_nodes_ptr_rm.size();sp_nodes_ptr_rm_i++){
        if ( ((sp_nodes_ptr_rm[sp_nodes_ptr_rm_i]->num_child==0) && sp_nodes_ptr_rm[sp_nodes_ptr_rm_i]->hybrid==1 ) ){
            Node * rm_parent_dummy_hybrid[2];
            rm_parent_dummy_hybrid[0]=NULL;
            rm_parent_dummy_hybrid[1]=NULL;
            rm_parent_dummy_hybrid[0]=sp_nodes_ptr_rm[sp_nodes_ptr_rm_i]->parent1;
            rm_parent_dummy_hybrid[1]=sp_nodes_ptr_rm[sp_nodes_ptr_rm_i]->parent2;
            int which_child_node_hybrid[2];
            for (int parent_i_hybrid=0;parent_i_hybrid<2;parent_i_hybrid++){                
                for (int child_i_hybrid=0;rm_parent_dummy_hybrid[parent_i_hybrid]->child.size();child_i_hybrid++){
                    if (rm_parent_dummy_hybrid[parent_i_hybrid]->child[child_i_hybrid]->label==sp_nodes_ptr_rm[sp_nodes_ptr_rm_i]->label){
                        which_child_node_hybrid[parent_i_hybrid]=child_i_hybrid;
                        break;
                    }    
                }
                rm_parent_dummy_hybrid[parent_i_hybrid]->num_child--;
                rm_parent_dummy_hybrid[parent_i_hybrid]->child.erase(rm_parent_dummy_hybrid[parent_i_hybrid]->child.begin()+which_child_node_hybrid[parent_i_hybrid]);
            }
            sp_nodes_ptr_rm[sp_nodes_ptr_rm_i]->clear();
        }
    }
}

vector < vector <int> > rm_one_child_interior_lambda_sum(string old_net_string, vector < vector < int > > old_lambda_sum){
    Net old_Net(old_net_string);
    //old_Net.print_all_node();
    vector <Node*> sp_nodes_ptr;
    for (size_t i=0;i<old_Net.Net_nodes.size();i++){
        Node* new_node_ptr=NULL;
        sp_nodes_ptr.push_back(new_node_ptr);
        sp_nodes_ptr[i]=&old_Net.Net_nodes[i];
    }
        
    for (size_t sp_nodes_ptr_i=0;sp_nodes_ptr_i<sp_nodes_ptr.size();sp_nodes_ptr_i++){
        if (sp_nodes_ptr[sp_nodes_ptr_i]->num_child==1 && sp_nodes_ptr[sp_nodes_ptr_i]->parent2==NULL && !sp_nodes_ptr[sp_nodes_ptr_i]->child[0]->tip_bool){
            old_lambda_sum[sp_nodes_ptr[sp_nodes_ptr_i]->child[0]->brchlen1-1].push_back(sp_nodes_ptr[sp_nodes_ptr_i]->brchlen1);
        }
    }
    old_Net.clear();
    return old_lambda_sum;    
}



string rm_one_child_interior_node(string old_net_string,bool e_num_bool){
    Net old_Net(old_net_string);
    vector <Node*> sp_nodes_ptr;
    for (size_t i=0;i<old_Net.Net_nodes.size();i++){
        Node* new_node_ptr=NULL;
        sp_nodes_ptr.push_back(new_node_ptr);
        sp_nodes_ptr[i]=&old_Net.Net_nodes[i];
    }
    for (size_t sp_nodes_ptr_i=0;sp_nodes_ptr_i<sp_nodes_ptr.size();sp_nodes_ptr_i++){
        if (sp_nodes_ptr[sp_nodes_ptr_i]->num_child==1 && sp_nodes_ptr[sp_nodes_ptr_i]->parent2==NULL){
            Node * rm_one_parent_dummy=sp_nodes_ptr[sp_nodes_ptr_i]->parent1;
            int which_one_child_parent;
            for (int child_i=0;rm_one_parent_dummy->child.size();child_i++){
                if (rm_one_parent_dummy->child[child_i]->label==sp_nodes_ptr[sp_nodes_ptr_i]->label){
                which_one_child_parent=child_i;
                break;
                }                
            }
            rm_one_parent_dummy->num_child--;
            rm_one_parent_dummy->child.erase(rm_one_parent_dummy->child.begin()+which_one_child_parent);
            sp_nodes_ptr[sp_nodes_ptr_i]->child[0]->parent1=NULL;
            add_node(rm_one_parent_dummy,sp_nodes_ptr[sp_nodes_ptr_i]->child[0]);
            if (!e_num_bool){
                sp_nodes_ptr[sp_nodes_ptr_i]->child[0]->brchlen1=sp_nodes_ptr[sp_nodes_ptr_i]->child[0]->brchlen1+sp_nodes_ptr[sp_nodes_ptr_i]->brchlen1;
            }
            sp_nodes_ptr[sp_nodes_ptr_i]->clear();
        }
    }
    rewrite_node_content(sp_nodes_ptr);
    return construct_adding_new_Net_str(old_Net);
    
}

vector <int> build_left_or_right_vec(int n_child,vector < valarray <int> > A_matrix_ith_row,size_t A_matrix_i){
    vector <int> left_or_right(2,0);
    for (size_t parent_i=0;parent_i<2;parent_i++){
        if (A_matrix_i>1){                
            for (int A_matrix_i_i_i=0;A_matrix_i_i_i<n_child;A_matrix_i_i_i++){
                if (A_matrix_ith_row[parent_i][A_matrix_i_i_i]==1){
                    left_or_right[parent_i]++;
                }
            }        
        }
        else{
            if (A_matrix_i==parent_i){
                for (int A_matrix_i_i_i=0;A_matrix_i_i_i<n_child;A_matrix_i_i_i++){
                    left_or_right[parent_i]++;
                }
            }
        }
    }
    return left_or_right;
}
        

vector <string> build_left_or_right_string(vector <int> left_or_right){
    vector <string> left_or_right_string;
    for (size_t i=0;i<left_or_right.size();i++){
        ostringstream left_or_right_string_dummy;
        left_or_right_string_dummy<< left_or_right[i];
        left_or_right_string.push_back(left_or_right_string_dummy.str());
    }
    return left_or_right_string;
}


string nchild_gt_one_core(string current_removing_net_string,string new_node_name,size_t new_node_name_i,size_t rm_node_index,int n_child, vector < valarray <int> > A_matrix_ith_row,size_t A_matrix_i){
    Net current_removing_net(current_removing_net_string);    
    vector <Node*> sp_nodes_ptr_rm;
    Node new_parent_dummy_L;
    Node new_parent_dummy_R;
    new_parent_dummy_L.label=new_node_name.substr(0,new_node_name_i)+"L";
    new_parent_dummy_R.label=new_node_name.substr(0,new_node_name_i)+"R";
    new_parent_dummy_L.brchlen1=current_removing_net.Net_nodes[rm_node_index].brchlen1;
    new_parent_dummy_R.brchlen1=current_removing_net.Net_nodes[rm_node_index].brchlen2;    
    valarray <int> new_parrent_dummy_descndent(0,current_removing_net.tax_name.size());
    
    current_removing_net.Net_nodes.push_back(new_parent_dummy_L);
    current_removing_net.Net_nodes.push_back(new_parent_dummy_R);
    current_removing_net.descndnt.push_back(new_parrent_dummy_descndent);
    current_removing_net.descndnt.push_back(new_parrent_dummy_descndent);    

    for (size_t i_link_ptr_rm=0;i_link_ptr_rm<current_removing_net.Net_nodes.size();i_link_ptr_rm++){
        Node* new_node_ptr=NULL;
        sp_nodes_ptr_rm.push_back(new_node_ptr);
        sp_nodes_ptr_rm[i_link_ptr_rm]=&current_removing_net.Net_nodes[i_link_ptr_rm];        
    }
    
    Node * rm_parent_dummy[2];
    rm_parent_dummy[0]=NULL;
    rm_parent_dummy[1]=NULL;
    rm_parent_dummy[0]=sp_nodes_ptr_rm[rm_node_index]->parent1;
    rm_parent_dummy[1]=sp_nodes_ptr_rm[rm_node_index]->parent2;
    int which_child_node[2];

    for (size_t parent_i=0;parent_i<2;parent_i++){
        for (int child_i=0;rm_parent_dummy[parent_i]->child.size();child_i++){
            if (rm_parent_dummy[parent_i]->child[child_i]->label==sp_nodes_ptr_rm[rm_node_index]->label){                                
                which_child_node[parent_i]=child_i;
                break;
            }    
        }

        rm_parent_dummy[parent_i]->num_child--;
        rm_parent_dummy[parent_i]->child.erase(rm_parent_dummy[parent_i]->child.begin()+which_child_node[parent_i]);
        
        if (A_matrix_i>1){                
            for (int A_matrix_i_i_i=0;A_matrix_i_i_i<n_child;A_matrix_i_i_i++){
                sp_nodes_ptr_rm[rm_node_index]->child[A_matrix_i_i_i]->parent1=NULL;
                if (A_matrix_ith_row[parent_i][A_matrix_i_i_i]==1){
                    add_node(sp_nodes_ptr_rm[sp_nodes_ptr_rm.size()+parent_i-2],sp_nodes_ptr_rm[rm_node_index]->child[A_matrix_i_i_i]);
                    sp_nodes_ptr_rm[rm_node_index]->child[A_matrix_i_i_i]->descndnt_of_hybrid=0;        
                }
            }        
            add_node(rm_parent_dummy[parent_i],sp_nodes_ptr_rm[sp_nodes_ptr_rm.size()+parent_i-2]);            
        }
        else{
            if (A_matrix_i==parent_i){
                for (int A_matrix_i_i_i=0;A_matrix_i_i_i<n_child;A_matrix_i_i_i++){
                    sp_nodes_ptr_rm[rm_node_index]->child[A_matrix_i_i_i]->parent1=NULL;
                    add_node(sp_nodes_ptr_rm[sp_nodes_ptr_rm.size()+parent_i-2],sp_nodes_ptr_rm[rm_node_index]->child[A_matrix_i_i_i]);
                    sp_nodes_ptr_rm[rm_node_index]->child[A_matrix_i_i_i]->descndnt_of_hybrid=0;
                }
                add_node(rm_parent_dummy[parent_i],sp_nodes_ptr_rm[sp_nodes_ptr_rm.size()+parent_i-2]);
            }
        }
    }

    sp_nodes_ptr_rm[rm_node_index]->parent1=NULL;
    sp_nodes_ptr_rm[rm_node_index]->parent2=NULL;

    //check hybrid node has zero kids
    rm_zero_kids_hybrid_node(sp_nodes_ptr_rm);
    
    for (size_t sp_nodes_ptr_rm_i=0;sp_nodes_ptr_rm_i<sp_nodes_ptr_rm.size();sp_nodes_ptr_rm_i++){
        sp_nodes_ptr_rm[sp_nodes_ptr_rm_i]->num_descndnt=0;
        for (size_t tax_name_i=0;tax_name_i<current_removing_net.tax_name.size();tax_name_i++){
            if (find_descndnt(sp_nodes_ptr_rm[sp_nodes_ptr_rm_i],current_removing_net.tax_name[tax_name_i])){
                current_removing_net.descndnt[sp_nodes_ptr_rm_i][tax_name_i]=1;
                }
            sp_nodes_ptr_rm[sp_nodes_ptr_rm_i]->num_descndnt=current_removing_net.descndnt[sp_nodes_ptr_rm_i].sum();
        }
    }
    
    current_removing_net.Net_nodes[rm_node_index].clear();
    Node * swap_dummy;
    swap_dummy=sp_nodes_ptr_rm.back();
    find_hybrid_descndnt(sp_nodes_ptr_rm.back());
    sp_nodes_ptr_rm.back()=sp_nodes_ptr_rm[sp_nodes_ptr_rm.size()-3];
    sp_nodes_ptr_rm[sp_nodes_ptr_rm.size()-3]=swap_dummy;

    rewrite_node_content(sp_nodes_ptr_rm);
    string adding_new_Net_string;    
    adding_new_Net_string=current_removing_net.Net_nodes[current_removing_net.Net_nodes.size()-3].node_content;
    adding_new_Net_string=adding_new_Net_string+current_removing_net.Net_nodes[current_removing_net.Net_nodes.size()-3].label;
    adding_new_Net_string.push_back(';');

    return adding_new_Net_string;
}


rm_H_node::rm_H_node(int rm_node_index,Net_wiz_prior_p new_Net_wiz_prior_p,bool maple_bool_local){
    //maple_bool_local=maple_bool_local_in;
    if (rm_debug_bool){
        string bug_statement="        Start of rm_H_node from " + new_Net_wiz_prior_p.s_net_string;
    }
    //rm_node_index=rm_node_index_in;
    Net current_net(new_Net_wiz_prior_p.s_net_string);
    
    current_prior_coal_list=new_Net_wiz_prior_p.prior_clade_list;
    prior_omega=new_Net_wiz_prior_p.omega;
    prior_omega_string=new_Net_wiz_prior_p.omega_string;    
    current_removing_net_string=new_Net_wiz_prior_p.s_net_string;
    current_removing_net_string_enum=new_Net_wiz_prior_p.s_net_string_enum;
    current_prior_coal_hist=new_Net_wiz_prior_p.prior_coal_hist;
    
    new_node_name=current_net.Net_nodes[rm_node_index].label;
    //n_child=current_net.Net_nodes[rm_node_index].num_child;
    
    
    current_root_enum=new_Net_wiz_prior_p.root_enum;
    current_lambda_sum=new_Net_wiz_prior_p.lambda_sum;
    

    new_node_name_i=hybrid_hash_index(new_node_name);
    //new_node_name_i=net_dummy.Net_nodes[i].label.find("#");
    if (new_node_name_i<new_node_name.size()){// check if this is neccessary
        //hybrid_parameter=new_node_name.substr(new_node_name_i+1,new_node_name.size()-1);
        hybrid_parameter=extract_hybrid_para_str(new_node_name);
        left_hybrid_parameter=hybrid_parameter;
        if (maple_bool_local){
            //left_hybrid_parameter_str="gamma["+new_node_name.substr(2,new_node_name_i-2)+"]";
            left_hybrid_parameter_str="gamma["+new_node_name.substr(1,new_node_name_i-1)+"]";
        }
        else{
            left_hybrid_parameter_str="\\gamma"+new_node_name.substr(1,new_node_name_i-1);
        }
        right_hybrid_parameter="(1-"+hybrid_parameter+")";
        right_hybrid_parameter_str="(1-"+left_hybrid_parameter_str+")";
        
    }
    
    
    left_hybrid_parameter_num=extract_hybrid_para(new_node_name);
    right_hybrid_parameter_num=1-left_hybrid_parameter_num;
    
    if (current_net.Net_nodes[rm_node_index].num_child>1){
        block_rm_H=nchild_gt_one(rm_node_index,maple_bool_local,current_net.Net_nodes[rm_node_index].num_child);
    }
    else{
        block_rm_H=nchild_eq_one(rm_node_index,maple_bool_local);
    }
    if (rm_debug_bool){
        cout<<"        End of rm_H_node "<< current_net.Net_nodes[rm_node_index].label<< " from " << new_Net_wiz_prior_p.s_net_string<<endl;
    }
}

vec_Net_wiz_prior_p rm_H_node::nchild_gt_one(int rm_node_index,bool maple_bool_local,int n_child){
    vec_Net_wiz_prior_p my_rmd_networks;
    vector < vector < valarray <int> > > A_matrix=build_h_child(n_child);
    for (size_t A_matrix_i=0;A_matrix_i<A_matrix.size();A_matrix_i++){        
        string adding_new_Net_string=nchild_gt_one_core(current_removing_net_string,new_node_name, new_node_name_i,rm_node_index,n_child, A_matrix[A_matrix_i],A_matrix_i);
        string adding_new_Net_string_enum=nchild_gt_one_core(current_removing_net_string_enum,new_node_name, new_node_name_i,rm_node_index,n_child, A_matrix[A_matrix_i],A_matrix_i);
        //cout<<adding_new_Net_string<<endl;
        //cout<<adding_new_Net_string_enum<<endl;
        vector <int> left_or_right=build_left_or_right_vec(n_child,A_matrix[A_matrix_i],A_matrix_i);
        vector <string> left_or_right_string=build_left_or_right_string(left_or_right);
        vector < vector < int > > new_lambda_sum=rm_one_child_interior_lambda_sum(adding_new_Net_string_enum, current_lambda_sum);
                
        bool enum_false=false;
        string adding_new_Net_string_updated=rm_one_child_interior_node(adding_new_Net_string,enum_false);    
        
        bool enum_true=true;
        string adding_new_Net_string_updated_enum=rm_one_child_interior_node(adding_new_Net_string_enum,enum_true);                        
    
        string new_left_hybrid="("+left_hybrid_parameter_str+"^"+left_or_right_string[0]+")";
        string new_right_hybrid="("+right_hybrid_parameter_str+"^"+left_or_right_string[1]+")";                            
        
        if (maple_bool_local){
            //current_omega=prior_omega*pow(1.0*left_hybrid_parameter_num,1.0*left_or_right[0])*pow(1.0*right_hybrid_parameter_num,1.0*left_or_right[1]);
            current_omega_string=new_left_hybrid+"*"+new_right_hybrid;            
        }
        else{
            current_omega=prior_omega*pow(1.0*left_hybrid_parameter_num,1.0*left_or_right[0])*pow(1.0*right_hybrid_parameter_num,1.0*left_or_right[1]);                            
            current_omega_string=prior_omega_string+new_left_hybrid+new_right_hybrid;            
        }

        Net_wiz_prior_p current_my_rmd_networks;
        current_my_rmd_networks.s_net_string=adding_new_Net_string_updated;
        current_my_rmd_networks.omega_string=current_omega_string;
        current_my_rmd_networks.omega=current_omega;
        current_my_rmd_networks.prior_clade_list=current_prior_coal_list;
        current_my_rmd_networks.s_net_string_enum=adding_new_Net_string_updated_enum;
        current_my_rmd_networks.root_enum=current_root_enum;
        current_my_rmd_networks.prior_coal_hist=current_prior_coal_hist;
        current_my_rmd_networks.lambda_sum=new_lambda_sum;
        my_rmd_networks.Net_vec.push_back(current_my_rmd_networks);
    }
    return my_rmd_networks;
}


string nchild_eq_one_core(string in_str,int rm_node_index,size_t i_one_hybrid_child){        
    Net current_removing_net(in_str);
    vector <Node*> sp_nodes_ptr_rm;
    for (size_t i_link_ptr_rm=0;i_link_ptr_rm<current_removing_net.Net_nodes.size();i_link_ptr_rm++){
        Node* new_node_ptr=NULL;
        sp_nodes_ptr_rm.push_back(new_node_ptr);
        sp_nodes_ptr_rm[i_link_ptr_rm]=&current_removing_net.Net_nodes[i_link_ptr_rm];
    }
    
    Node * rm_parent_dummy[2];
    rm_parent_dummy[0]=NULL;
    rm_parent_dummy[1]=NULL;
    rm_parent_dummy[0]=sp_nodes_ptr_rm[rm_node_index]->parent1;
    rm_parent_dummy[1]=sp_nodes_ptr_rm[rm_node_index]->parent2;
    
    int which_child_node[2];
    int parent_i=0;
    for (;parent_i<2;parent_i++){
        for (int child_i=0;rm_parent_dummy[parent_i]->child.size();child_i++){
        if (rm_debug_bool){
            cout<<rm_parent_dummy[parent_i]->label<<endl;
            cout<<rm_parent_dummy[parent_i]->child.size()<<endl;
            cout<<rm_parent_dummy[parent_i]->child[0]->label<<endl;
            cout<<sp_nodes_ptr_rm[rm_node_index]->label<<endl;
        }
            if (rm_parent_dummy[parent_i]->child[child_i]->label==sp_nodes_ptr_rm[rm_node_index]->label){
                which_child_node[parent_i]=child_i;
                break;
            }    
            
        }
        rm_parent_dummy[parent_i]->num_child--;
        rm_parent_dummy[parent_i]->child.erase(rm_parent_dummy[parent_i]->child.begin()+which_child_node[parent_i]);    
    }
    sp_nodes_ptr_rm[rm_node_index]->child[0]->descndnt_of_hybrid=0;
    sp_nodes_ptr_rm[rm_node_index]->child[0]->parent1=NULL;
    sp_nodes_ptr_rm[rm_node_index]->child[0]->parent2=NULL;        
    add_node(rm_parent_dummy[i_one_hybrid_child],sp_nodes_ptr_rm[rm_node_index]->child[0]);
    sp_nodes_ptr_rm[rm_node_index]->child.clear();

    //check hybrid node has zero kids
    rm_zero_kids_hybrid_node(sp_nodes_ptr_rm);

    sp_nodes_ptr_rm[rm_node_index]->parent1=NULL;
    sp_nodes_ptr_rm[rm_node_index]->parent2=NULL;
    rewrite_node_content(sp_nodes_ptr_rm);
    string adding_new_Net_string=construct_adding_new_Net_str(current_removing_net);
    return adding_new_Net_string;
}




vec_Net_wiz_prior_p rm_H_node::nchild_eq_one(int rm_node_index,bool maple_bool_local){
    vec_Net_wiz_prior_p my_rmd_networks;
    for (size_t i_one_hybrid_child=0;i_one_hybrid_child<2;i_one_hybrid_child++){
        string adding_new_Net_string=nchild_eq_one_core(current_removing_net_string, rm_node_index, i_one_hybrid_child);
        string adding_new_Net_string_enum=nchild_eq_one_core(current_removing_net_string_enum, rm_node_index, i_one_hybrid_child);

        vector < vector < int > > new_lambda_sum=rm_one_child_interior_lambda_sum(adding_new_Net_string_enum, current_lambda_sum);

        bool enum_false=false;
        string adding_new_Net_string_updated=rm_one_child_interior_node(adding_new_Net_string,enum_false);                
        bool enum_true=true;
        string adding_new_Net_string_updated_enum=rm_one_child_interior_node(adding_new_Net_string_enum,enum_true);                    
    
        double uni_hybrid_paramter_num;
        string uni_hybrid_paramter_string;    
        
        if (i_one_hybrid_child==0){
            uni_hybrid_paramter_string=left_hybrid_parameter_str;
            uni_hybrid_paramter_num=left_hybrid_parameter_num;
        }
        else{
            uni_hybrid_paramter_string=right_hybrid_parameter_str;
            uni_hybrid_paramter_num=right_hybrid_parameter_num;
        }
    
        if (maple_bool_local){
            current_omega_string="("+uni_hybrid_paramter_string+")";
            //current_omega=prior_omega*uni_hybrid_paramter_num;
        }
        else{
            current_omega_string=prior_omega_string+"("+uni_hybrid_paramter_string+")";
            current_omega=prior_omega*uni_hybrid_paramter_num;
        }    
                    
        Net_wiz_prior_p current_my_rmd_networks;
        current_my_rmd_networks.s_net_string=adding_new_Net_string_updated;
        current_my_rmd_networks.omega_string=current_omega_string;
        current_my_rmd_networks.omega=current_omega;
        current_my_rmd_networks.prior_clade_list=current_prior_coal_list;
        current_my_rmd_networks.s_net_string_enum=adding_new_Net_string_updated_enum;
        current_my_rmd_networks.root_enum=current_root_enum;
        current_my_rmd_networks.prior_coal_hist=current_prior_coal_hist;
        current_my_rmd_networks.lambda_sum=new_lambda_sum;
        
        my_rmd_networks.Net_vec.push_back(current_my_rmd_networks);
        
    }
    return my_rmd_networks;
}

bool check_sp_coal_valid(bool coal_ed,Net my_gt_tree, vector < valarray <int> > new_coal_clade){
    bool sp_coal_valid=false;
    if (coal_ed){
        for (size_t coal_clade_i=0;coal_clade_i<new_coal_clade.size();){
            bool coal_i_valid=false;
            for (size_t gt_node_i=0;gt_node_i<my_gt_tree.Net_nodes.size();gt_node_i++){
                valarray <bool> comp = (new_coal_clade[coal_clade_i] == my_gt_tree.descndnt[gt_node_i]);
                if ( comp.min() == true ){
                    coal_i_valid=true;
                    break;
                }
            }
            
            if (coal_i_valid){
                sp_coal_valid=true;
                coal_clade_i++;
                }
            else{
                sp_coal_valid=false;
                break;
            }
        }    
    }
    else{
        sp_coal_valid=true;
    }    
    return sp_coal_valid;
}


vec_Net_wiz_prior_p rm_S_node(int rm_node_index,string gt_str,Net_wiz_prior_p new_Net_wiz_prior_p,bool maple_bool_local_in){
    bool maple_bool_local=maple_bool_local_in;
    
    if (rm_debug_bool){
        cout<<"        Start of rm_S_node from " <<new_Net_wiz_prior_p.s_net_string<<endl;
    }
    
    vec_Net_wiz_prior_p my_rmd_networks;
    Net current_net;
    
    Net pre_current_net(new_Net_wiz_prior_p.s_net_string);

    //Net gt_tree;
    
    Net my_gt_tree(gt_str);
    string updated_gt_str;
    
    //string net_str=new_Net_wiz_prior_p.s_net_string;
    vector <string> original_tax_name=pre_current_net.tax_name;
    //for (size_t tax_i=0;tax_i<original_tax_name.size();tax_i++){
        //cout<<original_tax_name[tax_i]<<endl;
    //}
    //cout<<original_tax_name.size()<<endl;
    
    bool multiple_sp_tip=false;
    bool gt_tip_already_coaled=false;
    if (gt_str.size()>0){
        multiple_sp_tip=check_multiple_sp_tip(pre_current_net);
        gt_tip_already_coaled=check_gt_tip_already_coaled(my_gt_tree);
    }
    
    bool sp_valid_out=true;
    //check if my_Net is valid
    
    if (multiple_sp_tip && gt_str.size()>0){
        sp_valid_out=check_sp_valid(pre_current_net,my_gt_tree);    
    }

    //pre_current_net.print_all_node();

    
    if (sp_valid_out){

        if (multiple_sp_tip && !gt_tip_already_coaled){
            
            string dummy_string=construct_adding_new_Net_str(my_gt_tree);
            Net load_gt_tree(dummy_string);
            vector <Node*> gt_nodes_ptr_rm;
            for (size_t i_link_ptr_rm=0;i_link_ptr_rm<load_gt_tree.Net_nodes.size();i_link_ptr_rm++){
                Node* new_node_ptr=NULL;
                gt_nodes_ptr_rm.push_back(new_node_ptr);
                gt_nodes_ptr_rm[i_link_ptr_rm]=&load_gt_tree.Net_nodes[i_link_ptr_rm];
            }
            for (size_t i_gt_node=0;i_gt_node<load_gt_tree.descndnt.size();i_gt_node++){
                for (size_t i_st_node=0;i_st_node<pre_current_net.descndnt.size();i_st_node++){
                    valarray <bool> comp = (load_gt_tree.descndnt[i_gt_node]==pre_current_net.descndnt[i_st_node]);
                    if (comp.min() && pre_current_net.Net_nodes[i_st_node].tip_bool ){
                        gt_nodes_ptr_rm[i_gt_node]->num_child=0;
                        gt_nodes_ptr_rm[i_gt_node]->child.clear();
                        //gt_nodes_ptr_rm[i_gt_node]->rank=1;
                        gt_nodes_ptr_rm[i_gt_node]->tip_bool=true;
                        //cout<<gt_nodes_ptr_rm[i_gt_node]->label<<&(gt_nodes_ptr_rm[i_gt_node]->label)<<endl;
                        //cout<<gt_nodes_ptr_rm[i_gt_node]->clade<<&(gt_nodes_ptr_rm[i_gt_node]->clade)<<endl;
                        gt_nodes_ptr_rm[i_gt_node]->label=gt_nodes_ptr_rm[i_gt_node]->clade;
                        //cout<<gt_nodes_ptr_rm[i_gt_node]->label<<&(gt_nodes_ptr_rm[i_gt_node]->label)<<endl;
                        gt_nodes_ptr_rm[i_gt_node]->node_content.clear();
                    }
                }        
            }        
            //load_gt_tree.print_all_node();
            rewrite_node_content(gt_nodes_ptr_rm);
            //load_gt_tree.print_all_node();
            
            
            //string new_gt_str=construct_adding_new_Net_str(load_gt_tree);

            
            //Net new_gt_tree(new_gt_str);
            //gt_tree=new_gt_tree;
            updated_gt_str=construct_adding_new_Net_str(load_gt_tree);            

            string new_nt_str=construct_adding_new_Net_str(pre_current_net);
//pre_current_net.print_all_node();
//cout<<new_nt_str<<endl;
//!!!!!!!!!
            
//new_nt_str=rm_and_sign(new_nt_str);
            Net new_current_net(new_nt_str);
            //new_current_net.print_all_node();
            current_net=new_current_net;
            //my_gt_tree.print_all_node();
            //cout<<new_gt_str<<"********************************************************* end"<<endl;
        }
        else{
            //gt_tree=my_gt_tree;
            updated_gt_str=gt_str;
            current_net=pre_current_net;
            
        }
                //cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        //gt_tree.print_all_node();
        //cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
            
        
        //current_net.print_all_node();
            
        vector < valarray < int > > current_prior_coal_list=new_Net_wiz_prior_p.prior_clade_list;
        vector <int> current_prior_coal_hist=new_Net_wiz_prior_p.prior_coal_hist;
        vector < vector <int> > current_lambda_sum=new_Net_wiz_prior_p.lambda_sum;
    
    
        double prior_prior_prob_num=new_Net_wiz_prior_p.omega;
        string prior_prior_prob_string=new_Net_wiz_prior_p.omega_string;
        int n_child=current_net.Net_nodes[rm_node_index].num_child;
        
        //*** section of removing the descdent of the hybrid node.
        string current_removing_net_string=construct_adding_new_Net_str(current_net);
        vector < vector < valarray <int> > > A_matrix=build_s_child(n_child);
        string current_s_net_string_enum=new_Net_wiz_prior_p.s_net_string_enum;
        int current_root_enum=new_Net_wiz_prior_p.root_enum;
        
            //gt_tree.print_all_node();
                //cout<<"gt_tree.print_all_node(); not changed?"<<endl;
        for (size_t A_matrix_i=0;A_matrix_i<A_matrix.size();A_matrix_i++){
            Net current_removing_net(current_removing_net_string);
            Net current_removing_net_enum(current_s_net_string_enum);
            //cout<<current_s_net_string_enum<<endl;
            //vector < vector <int> > new_lambda_sum=current_lambda_sum;
            vector <Node*> sp_nodes_ptr_rm;
            vector <Node*> sp_nodes_ptr_rm_enum;
    
            for (size_t i_link_ptr_rm=0;i_link_ptr_rm<current_removing_net.Net_nodes.size();i_link_ptr_rm++){
                Node* new_node_ptr=NULL;
                sp_nodes_ptr_rm.push_back(new_node_ptr);
                sp_nodes_ptr_rm[i_link_ptr_rm]=&current_removing_net.Net_nodes[i_link_ptr_rm];
                Node* new_node_ptr_enum=NULL;
                sp_nodes_ptr_rm_enum.push_back(new_node_ptr_enum);
                sp_nodes_ptr_rm_enum[i_link_ptr_rm]=&current_removing_net_enum.Net_nodes[i_link_ptr_rm];
            }
            int which_child_node;
            
            Node *rm_parent_dummy=sp_nodes_ptr_rm[rm_node_index]->parent1;
            Node *rm_parent_dummy_enum=sp_nodes_ptr_rm_enum[rm_node_index]->parent1;
    
            for (int child_i=0;sp_nodes_ptr_rm[rm_node_index]->parent1->child.size();child_i++){
                if (rm_parent_dummy->child[child_i]->label==sp_nodes_ptr_rm[rm_node_index]->label){
                which_child_node=child_i;
                break;
                }    
            }
            rm_parent_dummy->num_child--;
            rm_parent_dummy_enum->num_child--;
            
            rm_parent_dummy->child.erase(rm_parent_dummy->child.begin()+which_child_node);
            rm_parent_dummy_enum->child.erase(rm_parent_dummy_enum->child.begin()+which_child_node);
    
            vector < Node* > rm_child_dummy;
            vector < Node* > rm_child_dummy_enum;
    
            bool coal_ed=false;
            
            int brch_lambda=sp_nodes_ptr_rm_enum[rm_node_index]->brchlen1;
    
            vector <int> new_prior_coal_hist=current_prior_coal_hist;
            vector < valarray < int > > new_current_prior_coal_clades=current_prior_coal_list;
            vector < valarray < int > > brand_new_current_prior_coal_clades;
            
            
            for (size_t A_matrix_i_i=0;A_matrix_i_i<A_matrix[A_matrix_i].size();A_matrix_i_i++){
                Node* current_rm_child_dummy=NULL;    
                Node* current_rm_child_dummy_enum=NULL;    
                int howmany_coaled=0;
                bool A_matrix_i_i_coaled=false;        
                
                for (int A_matrix_i_i_i=0;A_matrix_i_i_i<n_child;A_matrix_i_i_i++){
                    if (A_matrix[A_matrix_i][A_matrix_i_i][A_matrix_i_i_i]==1){
                        if (current_rm_child_dummy){
                            //current_rm_child_dummy->label.push_back('&');
                            current_rm_child_dummy->label=current_rm_child_dummy->label+ "&" + sp_nodes_ptr_rm[rm_node_index]->child[A_matrix_i_i_i]->label;
                            current_rm_child_dummy_enum->label=current_rm_child_dummy->label;
                            coal_ed=true;
                            howmany_coaled++;
                            A_matrix_i_i_coaled=true;
                            sp_nodes_ptr_rm[rm_node_index]->child[A_matrix_i_i_i]->clear();
                            sp_nodes_ptr_rm_enum[rm_node_index]->child[A_matrix_i_i_i]->clear();
    
                        }
                        else{
                            current_rm_child_dummy=sp_nodes_ptr_rm[rm_node_index]->child[A_matrix_i_i_i];
                            current_rm_child_dummy_enum=sp_nodes_ptr_rm_enum[rm_node_index]->child[A_matrix_i_i_i];
                        }
                    }
                }
                //gt_tree.print_all_node();
                //cout<<"gt_tree.print_all_node(); changed?"<<endl;
                if (rm_debug_bool){
                    cout<<"HEA"<<current_rm_child_dummy->label<<endl;
                }
                rm_child_dummy.push_back(current_rm_child_dummy);    
                rm_child_dummy_enum.push_back(current_rm_child_dummy_enum);
                
                if (A_matrix_i_i_coaled){
                    valarray <int> A_matrix_i_i_valarray (0,current_removing_net.tax_name.size());
                    vector <string> contained_tax_s;
                    for (size_t current_rm_child_dummy_label_i=0;current_rm_child_dummy_label_i<current_rm_child_dummy->label.size();){
                            string contained_tax;
                            size_t current_rm_child_dummy_label_j;
                            for (current_rm_child_dummy_label_j=current_rm_child_dummy_label_i;current_rm_child_dummy_label_j<current_rm_child_dummy->label.size();current_rm_child_dummy_label_j++){
                                if (current_rm_child_dummy->label[current_rm_child_dummy_label_j+1]=='&'){
                                    break;
                                }
                            }
                            contained_tax=current_rm_child_dummy->label.substr(current_rm_child_dummy_label_i,current_rm_child_dummy_label_j-current_rm_child_dummy_label_i+1);
                            contained_tax_s.push_back(contained_tax);
                            current_rm_child_dummy_label_i=current_rm_child_dummy_label_j+2;
                    }
                    for (size_t contained_tax_i=0;contained_tax_i<contained_tax_s.size();contained_tax_i++){
                        for (size_t tax_name_i=0;tax_name_i<current_removing_net.tax_name.size();tax_name_i++){
                            if (contained_tax_s[contained_tax_i]==current_removing_net.tax_name[tax_name_i]){
                                A_matrix_i_i_valarray[tax_name_i]=1;
                            }
                        }
                    }
                    
                    for (int howmany_i=0;howmany_i<howmany_coaled;howmany_i++){
                        new_prior_coal_hist.push_back(brch_lambda);
                    }
                    string A_matrix_i_i_clade=current_rm_child_dummy->label;//="(";

                    valarray <int> A_matrix_i_i_valarray_new(original_tax_name.size());

                    for (size_t tax_i=0;tax_i<original_tax_name.size();tax_i++){
                        size_t found;
                        //cout<<original_tax_name[tax_i]<<endl;
                        found=A_matrix_i_i_clade.find(original_tax_name[tax_i]);
                        if (found!=string::npos){
                            A_matrix_i_i_valarray_new[tax_i]++;
                        }
                    }
                    //cout<<A_matrix_i_i_clade<<endl;
                    new_current_prior_coal_clades.push_back(A_matrix_i_i_valarray_new);    
                    
                    //current_prior_coal_list.push_back(A_matrix_i_i_valarray);    
                    brand_new_current_prior_coal_clades.push_back(A_matrix_i_i_valarray);    
                }                    
            }
            
        
    
            bool sp_coal_valid=false;
            if (gt_str.size()>0){
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                
                //sp_coal_valid=check_sp_coal_valid( coal_ed, gt_tree, new_current_prior_coal_clades);
                sp_coal_valid=check_sp_coal_valid( coal_ed, my_gt_tree, new_current_prior_coal_clades);
            }
            else{
                sp_coal_valid=true;
            }
            
            
            //check!!!! i dont think this is needed!!!
            //if (gt_tree.Net_nodes.size()==0 && coal_ed){
                //sp_coal_valid=true;
            //}
                        
            if (sp_coal_valid){
                

                for (size_t dummy_i=0;dummy_i<rm_child_dummy.size();dummy_i++){
                    rm_child_dummy[dummy_i]->brchlen1=rm_child_dummy[dummy_i]->brchlen1+sp_nodes_ptr_rm[rm_node_index]->brchlen1;
                    rm_child_dummy[dummy_i]->parent1=NULL;
                    rm_child_dummy[dummy_i]->parent2=NULL;
                    add_node(rm_parent_dummy,rm_child_dummy[dummy_i]);
                    rm_child_dummy_enum[dummy_i]->parent1=NULL;
                    rm_child_dummy_enum[dummy_i]->parent2=NULL;
                    add_node(rm_parent_dummy_enum,rm_child_dummy_enum[dummy_i]);    
                }

                
                double current_brchlen=current_removing_net.Net_nodes[rm_node_index].brchlen1;
                current_removing_net.Net_nodes[rm_node_index].clear();
                current_removing_net_enum.Net_nodes[rm_node_index].clear();
                //new_current_enum_vec.erase(new_current_enum_vec.begin()+rm_node_index);//[rm_node_index].clear();
                rewrite_node_content(sp_nodes_ptr_rm);
                rewrite_node_content(sp_nodes_ptr_rm_enum);
                string adding_new_Net_string=construct_adding_new_Net_str(current_removing_net);
                string adding_new_Net_string_enum=construct_adding_new_Net_str(current_removing_net_enum);

                //Net adding_tree_dummy_enum(adding_new_Net_string_enum);

                
                //int k_clade=A_matrix[A_matrix_i].size();
                ////int num_w=factorial(k_clade);
                //double new_omega_prob=1.0*factorial(n_child-k_clade);
                ////cout<<new_omega_prob<<endl;
                ////debug_file<<"check here"<<endl;
                //if (rm_debug_bool){
                    //cout<<"*********"<<endl;
                    //cout<<current_net.Net_nodes[rm_node_index].node_content<<endl;
                    //cout<< n_child-k_clade<< " coal "<<endl;
                //}
                ////gt_tree.print_all_node();
                ////cout<<"gt_tree.print_all_node(); ova"<<endl;
                //for (size_t coal_clade_i=0;coal_clade_i<brand_new_current_prior_coal_clades.size();coal_clade_i++){
                    //if (brand_new_current_prior_coal_clades[coal_clade_i].sum()>1){
                    ////for (size_t jj=0;jj<brand_new_current_prior_coal_clades[coal_clade_i].size();jj++){
                        ////cout<<brand_new_current_prior_coal_clades[coal_clade_i][jj];
                    ////}
                    ////cout<<endl;
                        //for (size_t gt_node_i=0;gt_node_i<gt_tree.Net_nodes.size();gt_node_i++){
                                ////cout<<gt_node_i<<" gt_tree.Net_nodes[gt_node_i].label "<<gt_tree.Net_nodes[gt_node_i].label<<endl;
                            //valarray <bool> comp = (brand_new_current_prior_coal_clades[coal_clade_i] == gt_tree.descndnt[gt_node_i]);
                            //if ( comp.min() == true ){
                                //new_omega_prob=new_omega_prob/(1.0*(1+gt_tree.Net_nodes[gt_node_i].num_descndnt_interior));
            ////!!!!!!!!!!!!!!!!!!! error is here!!!!!!                

            /////////                
                                ////cout<<gt_node_i<<" gt_tree.Net_nodes[gt_node_i].label "<<gt_tree.Net_nodes[gt_node_i].label<<endl;
                                ////cout<<"gt_tree.Net_nodes[gt_node_i].interior_nodes_below.size()"<<gt_tree.Net_nodes[gt_node_i].interior_nodes_below.size()<<endl;    
                                ////cout<<"gt_tree.Net_nodes[gt_node_i].num_descndnt_interior"<<gt_tree.Net_nodes[gt_node_i].num_descndnt_interior<<endl;    
                                //if (gt_tree.Net_nodes[gt_node_i].interior_nodes_below.size()>1){
                                    //for (size_t interior_i=0;interior_i<gt_tree.Net_nodes[gt_node_i].interior_nodes_below.size();interior_i++){
                                ////if (gt_tree.Net_nodes[gt_node_i].num_descndnt_interior>0){
                                    ////for (size_t interior_i=0;interior_i<gt_tree.Net_nodes[gt_node_i].num_descndnt_interior;interior_i++){
                                    
                                        ////cout<<gt_tree.Net_nodes[gt_node_i].interior_nodes_below[interior_i]->num_descndnt_interior<<endl;
                                        //Node * dummy=gt_tree.Net_nodes[gt_node_i].interior_nodes_below[interior_i];
                                        //double num_int_des=double(dummy->num_descndnt_interior);
                                        ////cout<<"hea"<<endl;
                                        //cout<<gt_tree.Net_nodes[gt_node_i].label<<"  "<<gt_tree.Net_nodes[gt_node_i].interior_nodes_below[interior_i]->label<<"  "  <<num_int_des<<endl;
                                        //num_int_des=num_int_des+1.0;
                                        //new_omega_prob=new_omega_prob/(1.0*(num_int_des));
                                    //}
                                //}
                                
                                

                            //}
                            ////for (size_t interior_i=0;interior_i<gt_tree.Net_nodes[gt_node_i].interior_nodes_below.size();interior_i++){
                                ////new_omega_prob=new_omega_prob/(1.0*(1+(gt_tree.Net_nodes[gt_node_i].interior_nodes_below[interior_i]->num_descndnt_interior)));
                                
                            ////}
                            
                        //}
                    //}
                //}

                //cout<<new_omega_prob<<endl;
                Net gt_tree(updated_gt_str);
                int k_clade=A_matrix[A_matrix_i].size();
                double new_omega_prob=compute_rm_S_omega(n_child,k_clade, brand_new_current_prior_coal_clades, gt_tree);
                if (rm_debug_bool){
                    cout<<new_omega_prob<<"/";
                }
                ostringstream omega_w_str;
                omega_w_str<<new_omega_prob;
                ostringstream omega_d_str;
                int checking_n_minus_k=0;
                int omega_d=1;
                for (double coal_i=n_child;coal_i>double(k_clade);coal_i--){
                    //for (int coal_i=n_child;coal_i>k_clade;coal_i--){
                        //new_trial_prob=new_trial_prob/(coal_i*1.0)/((coal_i-1)*1.0)*2.0;
                        //new_trial_prob=new_trial_prob/coal_i/(coal_i-1)*2.0;
                        new_omega_prob=new_omega_prob/coal_i/(coal_i-1)*2.0;
                        checking_n_minus_k++;
                        omega_d=omega_d*coal_i*(coal_i-1)/2;
                }
                omega_d_str<<omega_d;
                if (rm_debug_bool){
                    cout<<omega_d<<endl;
                    cout<<"checking #coal"<<checking_n_minus_k<<endl;
                }
                
                
                ostringstream num_in;
                num_in<<n_child;
                ostringstream num_out;
                num_out<<k_clade;
                string current_omega_string;
                double current_omega;
                ostringstream brch_num;
                brch_num<<brch_lambda;
//!!!!!!!!!!!!!!!!!!! error is here!!!!!!                
                if (maple_bool_local){
                    //if (current_lambda_sum[brch_lambda-1].size()>0){
                        current_omega_string=omega_w_str.str()+"/"+omega_d_str.str()+"*puvT("+num_in.str()+","+num_out.str()+","+"lambda["+brch_num.str()+"]";
                        for (size_t lambda_sum_i=0; lambda_sum_i<current_lambda_sum[brch_lambda-1].size();lambda_sum_i++){
                            ostringstream brch_num_new;
                            brch_num_new<<current_lambda_sum[brch_lambda-1][lambda_sum_i];
                            current_omega_string=current_omega_string+"+lambda["+brch_num_new.str()+"]";
                        }
                        current_omega_string=current_omega_string+")";
                    //}
                    //else{
                        //current_omega_string="puvT("+num_in.str()+","+num_out.str()+","+"lambda["+brch_num.str()+"])";
                    //}
                    //current_omega=prior_prior_prob_num*gijoe(n_child,k_clade,current_brchlen);
                
                }
                else{                
            
                    current_omega=prior_prior_prob_num*gijoe(n_child,A_matrix[A_matrix_i].size(),current_brchlen)*new_omega_prob;
                    current_omega_string=prior_prior_prob_string+"\\frac{"+omega_w_str.str()+ "}{"+omega_d_str.str()+"}p_{"+num_in.str()+num_out.str()+"}(\\lambda_{"+brch_num.str()+"}";
                    for (size_t lambda_sum_i=0; lambda_sum_i<current_lambda_sum[brch_lambda-1].size();lambda_sum_i++){
                        ostringstream brch_num_new;
                        brch_num_new<<current_lambda_sum[brch_lambda-1][lambda_sum_i];
                        current_omega_string=current_omega_string+"+\\lambda_{"+brch_num_new.str()+"}";
                    }
                    current_omega_string=current_omega_string+")";

                }


                Net_wiz_prior_p current_my_rmd_networks;
                if (rm_debug_bool){
                    cout<<adding_new_Net_string<<"old"<<endl;
                    cout<<adding_new_Net_string_enum<<endl;
                    cout<<"current_omega  "<<current_omega<<endl;
                    cout<<current_omega_string<<endl;
                    cout<<current_root_enum<<endl;
                    cout<<adding_new_Net_string_enum<<endl;
                }

                current_my_rmd_networks.s_net_string=adding_new_Net_string;
                
                current_my_rmd_networks.omega_string=current_omega_string;
                current_my_rmd_networks.omega=current_omega;
                
                current_my_rmd_networks.prior_clade_list=new_current_prior_coal_clades;
                current_my_rmd_networks.s_net_string_enum=adding_new_Net_string_enum;

                current_my_rmd_networks.root_enum=current_root_enum;
                current_my_rmd_networks.prior_coal_hist=new_prior_coal_hist;
                current_my_rmd_networks.lambda_sum=current_lambda_sum;
                my_rmd_networks.Net_vec.push_back(current_my_rmd_networks);
                
            }    
            current_removing_net.clear();
            current_removing_net_enum.clear();

        }
        A_matrix.clear();
        //*** section of removing the descdent of the hybrid node.
        if (rm_debug_bool){
            cout<<"        End of rm_S_node "<< current_net.Net_nodes[rm_node_index].label<< " from " << new_Net_wiz_prior_p.s_net_string<<endl;
        }
        current_net.clear();
    }
    return my_rmd_networks;
}

Net initialize_enum_net(string in_str){
    Net enum_net(in_str);
    vector <Node*> enum_net_ptr;
    for (size_t node_i=0;node_i<enum_net.Net_nodes.size();node_i++){
        Node* new_node_ptr=NULL;
        enum_net_ptr.push_back(new_node_ptr);
        enum_net_ptr[node_i]=&enum_net.Net_nodes[node_i];    
        enum_net_ptr[node_i]->brchlen1=enum_net_ptr[node_i]->e_num;
        if (enum_net_ptr[node_i]->parent2){
            enum_net_ptr[node_i]->brchlen2=enum_net_ptr[node_i]->e_num2;
        }
    }
    rewrite_node_content(enum_net_ptr);
    return enum_net;
}

vector < vector <int> > initialize_lambda_sum(Net enum_net){
    vector < vector <int> > lambda_sum_initial;
    for (size_t enum_i=0;enum_i<enum_net.Net_nodes.back().e_num;enum_i++){
        vector <int> lambda_sum_dummy;
        lambda_sum_initial.push_back(lambda_sum_dummy);
    }
    return lambda_sum_initial;
}

vec_Net_wiz_prior_p build_initial_Networks(string in_str){
    Net_wiz_prior_p initial_Networks;    
    vec_Net_wiz_prior_p my_rmd_networks;    
    Net old_net_enum=initialize_enum_net(in_str);
    initial_Networks.s_net_string=in_str;
    initial_Networks.omega=1;
    initial_Networks.lambda_sum=initialize_lambda_sum(old_net_enum);
    initial_Networks.s_net_string_enum=construct_adding_new_Net_str(old_net_enum);
    initial_Networks.root_enum=old_net_enum.Net_nodes.back().e_num;
    //vector < valarray < int > > initial_prior_clade_list;
    //vector <int> initial_prior_coal_hist;
    //vector < vector < int > > initial_lambda_sum;
    initial_Networks.prior_clade_list.clear();//=initial_prior_clade_list;
    initial_Networks.prior_coal_hist.clear();//=initial_prior_coal_hist;
    //initial_Networks.lambda_sum.clear();//=initial_lambda_sum;
    
    my_rmd_networks.Net_vec.push_back(initial_Networks);
    return my_rmd_networks;
}    



vec_Net_wiz_prior_p simplify_Networks_one_block(bool hybrid_node,int rm_node_index,vec_Net_wiz_prior_p my_rmd_networks,size_t i,bool maple_bool_local,string gt_tree_str){
    string rm_string=my_rmd_networks.Net_vec[i].s_net_string;
    if (rm_debug_bool){
        cout<<"    Start rm_node for "<<rm_string<<endl;
    }
    if (hybrid_node){
        rm_H_node new_block(rm_node_index,my_rmd_networks.Net_vec[i],maple_bool_local);
        for (int block_network_i=0;block_network_i<new_block.block_rm_H.Net_vec.size();block_network_i++){
            my_rmd_networks.Net_vec.push_back(new_block.block_rm_H.Net_vec[block_network_i]);                
        }
    }
    else {
        vec_Net_wiz_prior_p new_block=rm_S_node(rm_node_index,gt_tree_str,my_rmd_networks.Net_vec[i],maple_bool_local);
        for (int block_network_i=0;block_network_i<new_block.Net_vec.size();block_network_i++){
            my_rmd_networks.Net_vec.push_back(new_block.Net_vec[block_network_i]);
        }
    }
    if (rm_debug_bool){
        cout<<"    End rm_node for "<<rm_string<<endl;
    }
    return my_rmd_networks;
}




vec_Net_wiz_prior_p simplify_Networks(string sp_str, string gt_string){
    bool maple_bool_local=false;
    if (rm_debug_bool){
        cout<<"Start simplify_Networks for "<<sp_str<<endl;        
    }
    vec_Net_wiz_prior_p my_rmd_networks=build_initial_Networks(sp_str);

    string adding_new_Net_str;
    bool needs_to_rm=false;
    while (!needs_to_rm){
        for (size_t i=0;i<my_rmd_networks.Net_vec.size();){
            Net current_net(my_rmd_networks.Net_vec[i].s_net_string);    
            int rm_node_index=choose_rm_node(current_net);
            if (rm_node_index>=0){
                my_rmd_networks=simplify_Networks_one_block(current_net.Net_nodes[rm_node_index].hybrid,rm_node_index,my_rmd_networks,i,maple_bool_local,gt_string);
                my_rmd_networks.Net_vec.erase(my_rmd_networks.Net_vec.begin()+i);
            }
            else{
                i++;
            }
        }
        needs_to_rm=true;
        for (size_t j=0;j<my_rmd_networks.Net_vec.size();){
            Net checking_rm_net(my_rmd_networks.Net_vec[j].s_net_string);
            int rm_node_index=choose_rm_node(checking_rm_net);
            if (rm_node_index>=0){
            needs_to_rm=false;
            break;}
            else{
                j++;
            }
        }        
    }
    
        
    if (rm_debug_bool){
        string bug_statement="End simplify_Networks for "+sp_str;
    }
    return my_rmd_networks;
}





double compute_rm_S_omega(int n_child,int k_clade,vector < valarray < int > > brand_new_current_prior_coal_clades, Net gt_tree){
    //int k_clade=A_matrix[A_matrix_i].size();
    
    double new_omega_prob=1.0*factorial(n_child-k_clade);
    //cout<<new_omega_prob<<endl;
    if (rm_debug_bool){
        cout<<"*********"<<endl;
        //cout<<current_net.Net_nodes[rm_node_index].node_content<<endl;
        cout<< n_child-k_clade<< " coal "<<endl;
    }
    //gt_tree.print_all_node();
    //cout<<"gt_tree.print_all_node(); ova"<<endl;
    for (size_t coal_clade_i=0;coal_clade_i<brand_new_current_prior_coal_clades.size();coal_clade_i++){
        if (brand_new_current_prior_coal_clades[coal_clade_i].sum()>1){
        //for (size_t jj=0;jj<brand_new_current_prior_coal_clades[coal_clade_i].size();jj++){
            //cout<<brand_new_current_prior_coal_clades[coal_clade_i][jj];
        //}
        //cout<<endl;
            for (size_t gt_node_i=0;gt_node_i<gt_tree.Net_nodes.size();gt_node_i++){
                    //cout<<gt_node_i<<" gt_tree.Net_nodes[gt_node_i].label "<<gt_tree.Net_nodes[gt_node_i].label<<endl;
                valarray <bool> comp = (brand_new_current_prior_coal_clades[coal_clade_i] == gt_tree.descndnt[gt_node_i]);
                if ( comp.min() == true ){
                    new_omega_prob=new_omega_prob/(1.0*(1+gt_tree.Net_nodes[gt_node_i].num_descndnt_interior));
//!!!!!!!!!!!!!!!!!!! error is here!!!!!!                

///////                
                    //cout<<gt_node_i<<" gt_tree.Net_nodes[gt_node_i].label "<<gt_tree.Net_nodes[gt_node_i].label<<endl;
                    //cout<<"gt_tree.Net_nodes[gt_node_i].interior_nodes_below.size()"<<gt_tree.Net_nodes[gt_node_i].interior_nodes_below.size()<<endl;    
                    //cout<<"gt_tree.Net_nodes[gt_node_i].num_descndnt_interior"<<gt_tree.Net_nodes[gt_node_i].num_descndnt_interior<<endl;    
                    if (gt_tree.Net_nodes[gt_node_i].interior_nodes_below.size()>1){
                        for (size_t interior_i=0;interior_i<gt_tree.Net_nodes[gt_node_i].interior_nodes_below.size();interior_i++){
                    //if (gt_tree.Net_nodes[gt_node_i].num_descndnt_interior>0){
                        //for (size_t interior_i=0;interior_i<gt_tree.Net_nodes[gt_node_i].num_descndnt_interior;interior_i++){
                        
                            //cout<<gt_tree.Net_nodes[gt_node_i].interior_nodes_below[interior_i]->num_descndnt_interior<<endl;
                            Node * dummy=gt_tree.Net_nodes[gt_node_i].interior_nodes_below[interior_i];
                            double num_int_des=double(dummy->num_descndnt_interior);
                            //cout<<"hea"<<endl;
                            //cout<<gt_tree.Net_nodes[gt_node_i].label<<"  "<<gt_tree.Net_nodes[gt_node_i].interior_nodes_below[interior_i]->label<<"  "  <<num_int_des<<endl;
                            num_int_des=num_int_des+1.0;
                            new_omega_prob=new_omega_prob/(1.0*(num_int_des));
                        }
                    }
                    
                    

                }
                //for (size_t interior_i=0;interior_i<gt_tree.Net_nodes[gt_node_i].interior_nodes_below.size();interior_i++){
                    //new_omega_prob=new_omega_prob/(1.0*(1+(gt_tree.Net_nodes[gt_node_i].interior_nodes_below[interior_i]->num_descndnt_interior)));
                    
                //}
                
            }
        }
    }
    return new_omega_prob;
}



