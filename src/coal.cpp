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

//#include"coal.hpp"
#include"net.hpp"

void Net::assign_bl_to_vec(){
	this->brchlens_vec = vector < double > ( this->NodeContainer.back().e_num(), 0 );
	this->max_num_brch_vec = vector < int > ( this->NodeContainer.back().e_num(),0);
		
	dout<<"branch lengths labellings"<<endl;
	for ( size_t node_i = 0; node_i < this->NodeContainer.size(); node_i++){
		Node * current_node = &this->NodeContainer[node_i];
		size_t index_enum = current_node->e_num();

		if ( index_enum != 0 ){			
			brchlens_vec[index_enum-1] = current_node->brchlen1();
			//cout<<current_node->label<<"  "<<current_node->brchlen1<<endl;
			max_num_brch_vec[index_enum-1] = this->samples_below[node_i].sum();
			//cout<<index_enum<<" "<<brchlens_vec[index_enum-1]<<" "<<max_num_brch_vec[index_enum-1]<<endl;
		}
		dout<<index_enum <<" "<<current_node->brchlen1()<<endl;
		
		if ( current_node->parent2 ){
			index_enum = current_node->e_num2();
			brchlens_vec[index_enum-1] = current_node->brchlen2();
			//cout<<current_node->label<<"  "<<current_node->brchlen1<<endl;
			max_num_brch_vec[index_enum-1] =  this->samples_below[node_i].sum();
			dout<<index_enum <<" "<<current_node->brchlen2()<<endl;
		}
	}
}


void Net::build_gijoe(){
	double v1,v2,v3,u1,u2;
	int u,v,h,k;
	for ( size_t b = 0; b < this->brchlens_vec.size(); b++){
		vector < vector < double > > gijoe_matrix_b;
		vector < double > empty_gijoe_matrix_b_1;
		gijoe_matrix_b.push_back(empty_gijoe_matrix_b_1);
		for ( u = 2; u <= this->max_num_brch_vec[b]; u++){
			
			vector < double > gijoe_matrix_b_u;
				for ( v = 1; v <= u; v++){
					double gi_joe_buv;
					//if (brchlens_vec[b]==0){
					if (b == brchlens_vec.size() - 1 ){
						gi_joe_buv = 1;
                    } else{
						gi_joe_buv=0;
						for( k = v; k <= u; k++) {
							v1 = 1.0;
							for( h = v; h <= v+k-2; h++ ){ v1 *= h; }
							v2 = 1.0; 
							for( h = v; h > 1; h-- ){ v2 *= h; }
							v3 = 1.0;
							for( h = k-v; h > 1; h-- ){ v3 *= h; }
							u1 = 1.0;
							for( h = u; h >= u - k + 1; h--){ u1 *= h; }
							u2 = 1.0;
							for( h = u; h <= u+k-1; h++){ u2 *= h; }
							gi_joe_buv += exp(.5*k*(1.0-k)*brchlens_vec[b])*(2.0*k-1.0)*pow(-1.0,(k-v)*1.0)*v1*u1/(v2*v3*u2);
						}
					}
					
					//cout<<"b:"<<b+1<<" u:"<<u<<" v:"<<v<<endl;
					//cout<<gi_joe_buv<<"  "<<gijoe(u,v,brchlens_vec[b])<<endl;
					//if (gi_joe_buv==gijoe(u,v,brchlens_vec[b])){cout<<"yes"<<endl;}
					//else{cout<<"no"<<endl;}
					//cout<<"g_["<<u<<","<<v<<"]("<<brchlens_vec[b]<<")="<<gi_joe_buv<<"  ";
					gijoe_matrix_b_u.push_back(gi_joe_buv);
					//cout<<"      "<<gijoe_matrix_b_u.size()<<endl;
				}
				//cout<<endl;
			gijoe_matrix_b.push_back(gijoe_matrix_b_u);
			//cout<<"  "<<gijoe_matrix_b.size()<<endl;
		}
		//cout<<endl<<endl;
		this->gijoemat.push_back(gijoe_matrix_b);
	}
    assert( print_gijoemat() );
}


bool Net::print_gijoemat(){
    for ( size_t branch_idx = 0; branch_idx < this->gijoemat.size(); branch_idx++ ){
        for ( size_t branch_in = 0; branch_in < this->gijoemat[branch_idx].size(); branch_in++){
            for (size_t branch_out = 0; branch_out < this->gijoemat[branch_idx][branch_in].size(); branch_out++){
                dout << "g_[" << branch_in << "," << branch_out << "](" << brchlens_vec[branch_idx] << ")=" << this->gijoemat[branch_idx][branch_in][branch_out] << "  ";
            }
            dout<<endl;
        }
        dout<<endl;
    }
    return true;
}

	
void Tree::building_M_matrix( Net & sp_net ) {
	int sp_max_enum = sp_net.NodeContainer.back().e_num();
	int gt_max_enum = this->NodeContainer.back().e_num();
	
	
	for ( int i_sp_enum = 0; i_sp_enum < sp_max_enum; i_sp_enum++ ){
		vector < int > M_matrix_row( ( gt_max_enum - 1 ), 1 );
		this->M_matrix.push_back( M_matrix_row );
	}
    
	for ( size_t i = 0; i < sp_net.NodeContainer.size() - 1; i++ ){
		if ( sp_net.NodeContainer[i].tip() ) continue;
        for (size_t j = 0; j < this->NodeContainer.size() - 1; j++){
            if ( this->NodeContainer[j].tip() ) continue;
            int des_diff_prod = 1;
            for ( size_t i_des = 0; i_des < this->tax_name.size(); i_des++ ){
                int des_diff = sp_net.descndnt[i][i_des] - this->descndnt[j][i_des]; /*! \todo this may not be right for multiple lineages per species*/
                if ( des_diff < 0 ){
                    des_diff_prod = 0;
                    M_matrix[sp_net.NodeContainer[i].e_num()-1][this->NodeContainer[j].e_num()-1] = 0;
                    break;
                }
                else{  /*! \todo check this else!!! dont think this is needed */
                    if (des_diff>0){
                        des_diff_prod=0;}
                }
            }
            
            if ( des_diff_prod == 1 ){
                M_matrix[sp_net.NodeContainer[i].e_num()-1][this->NodeContainer[j].e_num()-1] = 1;
            }
        }
	}
	
    dout<<"M matrix"<<endl;
    assert( print_matrix( M_matrix ) );
}


void Tree::building_R_matrix( ){
	size_t gt_max_enum = this->NodeContainer.back().e_num();	
	for ( size_t i_gt_enum=0; i_gt_enum < gt_max_enum; i_gt_enum++ ){
		//vector <int> R_matrix_row;
		//for (size_t j_gt_enum=0;j_gt_enum<i_gt_enum;j_gt_enum++){
			//R_matrix_row.push_back(1);
		//}
		vector <int> R_matrix_row( i_gt_enum, 1 );
		//vector <int> R_matrix_row(gt_max_enum,1);
		this->R_matrix.push_back(R_matrix_row);
	}
    
	for (size_t i = 0; i < this->NodeContainer.size()-1; i++ ){
		if ( this->NodeContainer[i].tip()) continue;
        for ( size_t j = 0; j < i; j++){
            if ( this->NodeContainer[j].tip() ) continue;
            for ( size_t i_des = 0; i_des < this->descndnt[i].size(); i_des++ ){
                int des_diff = this->descndnt[i][i_des] - this->descndnt[j][i_des];
                if ( des_diff < 0 ){
                    this->R_matrix[this->NodeContainer[i].e_num()-1][this->NodeContainer[j].e_num()-1]=0;
                    break;
                }
            }
        }
	}
	dout << "R matrix" << endl;
	assert(this->print_matrix(R_matrix));
}


void Net::building_S_matrix(){
	int sp_max_enum=this->NodeContainer.back().e_num();
	//vector < vector < int > > S_matrix;
	
	for ( int i_sp_enum = 0; i_sp_enum < sp_max_enum; i_sp_enum++ ){
		vector < int > S_matrix_row( sp_max_enum, 1 );
		S_matrix_row[i_sp_enum] = 0;
		this->S_matrix.push_back( S_matrix_row );
	}
    
	for ( size_t i = 0; i < this->NodeContainer.size()-1; i++ ){
		if ( this->NodeContainer[i].tip() ) continue;
            
        for ( size_t j = 0; j < this->NodeContainer.size()-1; j++){
            if ( this->NodeContainer[j].tip() ) continue;
                
            for ( size_t i_des = 0; i_des < this->descndnt[i].size(); i_des++ ){
                int des_diff = this->descndnt[i][i_des] - this->descndnt[j][i_des];
                if ( des_diff < 0 ){
                    S_matrix[this->NodeContainer[i].e_num()-1][this->NodeContainer[j].e_num()-1] = 0;
                    break;
                }
            }
        }
		
	}
	assert(print_matrix(S_matrix));
}


bool Tree::print_matrix( vector < vector < int > > & mat ){
	for ( size_t i = 0; i < mat.size(); i++ ){
		for ( size_t j = 0; j < mat[i].size(); j++){
			dout << mat[i][j];
		}
		dout << endl;
	}
	dout<<endl;
    return true;
}


void Tree::build_coal_hist_mat( Net& sp_net ){
	
	size_t sp_max_enum = sp_net.NodeContainer.back().e_num();
	size_t gt_max_enum = this->NodeContainer.back().e_num();
	
	if ( gt_max_enum == 1 ){
		vector <size_t> coal_hist_vec;
		coal_hist_vec.push_back(1);
		coal_hist_mat.push_back(coal_hist_vec);
	} else{
		int max_coal_hist_num = 1;
		for ( size_t i_gt_enum = 0; i_gt_enum < gt_max_enum-1; i_gt_enum++ ){
			vector <size_t> coal_hist_vec;
			for ( size_t i_sp_enum = 0; i_sp_enum < sp_max_enum; i_sp_enum++){
				if ( M_matrix[i_sp_enum][i_gt_enum] == 1 ){
					coal_hist_vec.push_back(i_sp_enum+1);
				}
			}
			max_coal_hist_num = max_coal_hist_num*(coal_hist_vec.size());
			this->coal_hist_mat.push_back(coal_hist_vec);
		}
	}
}

double calc_prob_of_hist_para(vector < vector <size_t> > coal_hist,
vector < vector <int> > all_w,
vector < vector <int> > all_d,
vector < vector <int> > num_enter,
vector < vector <int> > num_out,
Net & sp_net
){
	double probability=0.0;
	for (size_t i_coal_hist=0;i_coal_hist<coal_hist.size();i_coal_hist++){
		double current_prob_of_hist=1;
		for (size_t i=0;i<all_w[i_coal_hist].size();i++){
			size_t current_enum;
			for ( size_t node_i = 0; node_i < sp_net->NodeContainer.size(); node_i++ ){
				Node * current_node = sp_net->NodeContainer[node_i];
				if ( current_node->e_num() == (i+1) ){
					//current_branch_lengths=my_sp_net.NodeContainer[node_i].brchlen1;
					current_enum = current_node->e_num();
					break;
				}
			}	
			double current_gijoe = ( (i+1) == sp_net->NodeContainer.back()->e_num() ) ? 1 :
                                                                                           sp_net.gijoe_matrix->mat[current_enum-1][num_enter[i_coal_hist][i]-1][num_out[i_coal_hist][i]-1];
			//if ((i+1) == my_sp_net->NodeContainer.back()->e_num()){
				//current_gijoe=1;
				//}
				//else{
				//current_gijoe=gijoe_matrix->mat[current_enum-1][num_enter[i_coal_hist][i]-1][num_out[i_coal_hist][i]-1];
////!!! this works!!! //current_gijoe=gijoe(num_enter[i_coal_hist][i], num_out[i_coal_hist][i], current_branch_lengths); 				
			//}
			current_prob_of_hist *= all_w[i_coal_hist][i] / all_d[i_coal_hist][i] * current_gijoe;
            dout<<current_prob_of_hist<<endl;
		}
		probability += current_prob_of_hist;
		dout<<"gt prob = "<<probability<<endl;
	}
	return probability;
}



double calc_prob_of_hist(
vector < vector < int > > S_matrix,
vector < vector < int > > R_matrix,
vector < vector < int > > M_matrix,
vector < vector <size_t> > coal_hist_mat,
vector < vector <size_t> > coal_hist,
Net * my_sp_net,
gijoemat * gijoe_matrix
){
		vector < vector <int> > all_w;
		vector < vector <int> > all_d;
		vector < vector <int> > num_enter;
		vector < vector <int> > num_out;
    size_t sp_max_enum=my_sp_net->NodeContainer.back().e_num();

    for (size_t i_coal_hist=0;i_coal_hist<coal_hist.size();i_coal_hist++){	
        //for (size_t j_coal_hist=0;j_coal_hist<coal_hist[i_coal_hist].size();j_coal_hist++){
            //cout<<coal_hist[i_coal_hist][j_coal_hist];
        //}
        vector <int> num_enter_i_coal;
        vector <int> num_out_i_coal;
        vector <int> num_coal_in_branch_i_coal;
        vector <int> w_i_coal;
        vector <int> d_i_coal;
        for (size_t i=0;i<sp_max_enum;i++){
            for (size_t ith_node=0;ith_node<my_sp_net->NodeContainer.size();ith_node++){
                if (!my_sp_net->NodeContainer[ith_node]->tip() && my_sp_net->NodeContainer[ith_node].e_num()==i+1){
                    num_enter_i_coal.push_back(my_sp_net->NodeContainer[ith_node]->num_descndnt_tips());
                }
            }
        }	
        for (size_t i=0;i<sp_max_enum;i++){
            vector < int > clades_coal_in_branch;
            int num_coal_in_branch_i_coal_dummy=0;
            for (size_t coal_hist_position=0;coal_hist_position<coal_hist[i_coal_hist].size();coal_hist_position++){
                if (coal_hist[i_coal_hist][coal_hist_position]==i+1){
                    num_coal_in_branch_i_coal_dummy++;
                    clades_coal_in_branch.push_back(coal_hist_position);
                    for (size_t i_S_mat=i;i_S_mat<(sp_max_enum);i_S_mat++){
                        if (S_matrix[i_S_mat][i]==1){
                            num_enter_i_coal[i_S_mat]--;
                        }
                    }	
                }			
            }
            num_coal_in_branch_i_coal.push_back(num_coal_in_branch_i_coal_dummy);	
            num_out_i_coal.push_back(num_enter_i_coal[i]-num_coal_in_branch_i_coal[i]);
            int d_dummy=1;
            int coal_dummy=num_enter_i_coal[i];
            for (int y=0;y<num_coal_in_branch_i_coal_dummy;y++){
                //d_dummy=d_dummy*n_choose_k(coal_dummy,2);
                d_dummy=d_dummy*coal_dummy*(coal_dummy-1)/2;
                coal_dummy--;
            }
            d_i_coal.push_back(d_dummy);
            
            int w_dummy=factorial(num_coal_in_branch_i_coal_dummy);
            vector < vector < int > > updated_R_mat;
            if (clades_coal_in_branch.size()>1){
                for (size_t updated_R_i=0;updated_R_i<clades_coal_in_branch.size();updated_R_i++){
                    vector < int > updated_R_mat_row;
                    for (size_t updated_R_j=0;updated_R_j<updated_R_i;updated_R_j++){
                        updated_R_mat_row.push_back(R_matrix[clades_coal_in_branch[updated_R_i]][clades_coal_in_branch[updated_R_j]]);
                    }
                    updated_R_mat.push_back(updated_R_mat_row);
                }
                for (size_t updated_R_mat_i=0;updated_R_mat_i<updated_R_mat.size();updated_R_mat_i++){
                    int sum_r=1;
                    for (size_t updated_R_mat_j=0;updated_R_mat_j<updated_R_mat[updated_R_mat_i].size();updated_R_mat_j++){
                        sum_r=sum_r+updated_R_mat[updated_R_mat_i][updated_R_mat_j];
                    }
                    w_dummy=w_dummy/sum_r;
                }
            }	
            
            w_i_coal.push_back(w_dummy);
            //if (w_dummy/d_dummy!=1){
                //cout<<w_dummy<<"/"<<d_dummy;
            //}
            //if (i<sp_max_enum-1){
                //cout<<"p_{"<<num_enter_i_coal[i]<<num_out_i_coal[i]<<"}(lambda_{"<<i+1<<"})";	
            //}
        }
        //cout<<endl;
        num_enter.push_back(num_enter_i_coal);
        num_out.push_back(num_out_i_coal);
        //num_coal_in_branch.push_back(num_coal_in_branch_i_coal);
        all_w.push_back(w_i_coal);
        all_d.push_back(d_i_coal);	
        coal_hist[i_coal_hist].push_back(sp_max_enum);
    }	
	return calc_prob_of_hist_para(coal_hist,all_w,all_d,num_enter,num_out,my_sp_net,gijoe_matrix);
}			


vector < vector <size_t> >recur_coal_hist(
vector < vector <size_t> > coal_hist_mat,
vector < vector < int > > R_matrix,
//Net * my_sp_net,
Net * my_gt_tree
){
	vector < vector <size_t> > coal_hist;

	for (size_t first_coal_mat_i=0;first_coal_mat_i<coal_hist_mat[0].size();first_coal_mat_i++){
		vector <size_t> coal_hist_dummy;
		coal_hist_dummy.push_back(coal_hist_mat[0][first_coal_mat_i]);
		coal_hist.push_back(coal_hist_dummy);
	}
	int gt_max_enum=my_gt_tree->NodeContainer.back()->e_num();	
	if (gt_max_enum-1>1){
		coal_hist=build_coal_hist(coal_hist,1,coal_hist_mat,R_matrix);
	}

	//size_t sp_max_enum=my_sp_net->NodeContainer.back()->e_num();
	//for (size_t i_coal_hist=0;i_coal_hist<coal_hist.size();i_coal_hist++){
		//coal_hist[i_coal_hist].push_back(sp_max_enum);
	//}
	return coal_hist;
}


vector < vector <size_t> > build_coal_hist(
vector < vector <size_t> > coal_hist,
size_t node_i,
vector < vector <size_t> > coal_hist_mat, 
vector < vector < int > > R_matrix)
{
	//dout<<"Start of build_coal_hist"<<endl;
	
	vector < vector <size_t> > new_coal_hist;
	for (size_t coal_hist_i=0;coal_hist_i<coal_hist.size();coal_hist_i++){
		vector <size_t> coal_hist_dummy;
		coal_hist_dummy=coal_hist[coal_hist_i];
		size_t coal_hist_mat_node_i_i=0;
		for (size_t j_R_mat=0;j_R_mat<node_i;j_R_mat++){
			while (R_matrix[node_i][j_R_mat]==1 && coal_hist_mat[node_i][coal_hist_mat_node_i_i]< coal_hist_dummy[j_R_mat]){
				coal_hist_mat_node_i_i++;
			}
		}
		for (;coal_hist_mat_node_i_i<coal_hist_mat[node_i].size();coal_hist_mat_node_i_i++ ){
			vector <size_t> coal_hist_dummy_dummy;
			coal_hist_dummy_dummy=coal_hist_dummy;
			coal_hist_dummy_dummy.push_back(coal_hist_mat[node_i][coal_hist_mat_node_i_i]);
			new_coal_hist.push_back(coal_hist_dummy_dummy);
		}
	}
	
		//dout<<new_coal_hist.size()<<"new_coal_hist"<<endl;	
		//print_matrix(new_coal_hist);
		dout<<"Recurse to "<<new_coal_hist.size()<<" histories: "<<endl;	
		for (size_t i=0;i<new_coal_hist.size();i++){
			for (size_t j=0;j<new_coal_hist[i].size()-1;j++){
				dout<<new_coal_hist[i][j]<<",";
			}
			dout<<new_coal_hist[i].back()<<endl;
		}
	
	if (node_i<coal_hist_mat.size()-1){
		node_i++;
		new_coal_hist=build_coal_hist(new_coal_hist,node_i,coal_hist_mat,R_matrix);		
	}

		//dout<<"End of build_coal_hist"<<endl;

	return new_coal_hist;
}
