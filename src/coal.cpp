#include"coal.hpp"
#include"net.hpp"

gijoemat::gijoemat(Net* my_sp_net){
	size_t sp_max_enum=my_sp_net->Net_nodes.back()->e_num();
	vector < double > brchlens_vec(sp_max_enum,0);
	vector < int > max_num_brch_vec(sp_max_enum,0);
		
	dout<<"branch lengths labellings"<<endl;
	for (size_t node_i=0;node_i<my_sp_net->Net_nodes.size();node_i++){
		Node * current_node = my_sp_net->Net_nodes[node_i];
		size_t index_enum = current_node->e_num();

		if (index_enum!=0){
			
			brchlens_vec[index_enum-1]=current_node->brchlen1();
			//cout<<current_node->label<<"  "<<current_node->brchlen1<<endl;
			max_num_brch_vec[index_enum-1]=current_node->num_descndnt();
			//cout<<index_enum<<" "<<brchlens_vec[index_enum-1]<<" "<<max_num_brch_vec[index_enum-1]<<endl;
			
		}
		dout<<index_enum <<" "<<current_node->brchlen1()<<endl;
		
		if (current_node->parent2()){
			index_enum=current_node->e_num2();
			brchlens_vec[index_enum-1]=current_node->brchlen2();
			//cout<<current_node->label<<"  "<<current_node->brchlen1<<endl;
			max_num_brch_vec[index_enum-1]=current_node->num_descndnt();
			dout<<index_enum <<" "<<current_node->brchlen1()<<endl;
		}		

	}
	//if (coal_debug_bool){cout<<endl;}
	//vector < vector < vector < double > > > gijoe_matrix;	
	double v1,v2,v3,u1,u2;
	int u,v,h,k;
	for (size_t b=0;b<brchlens_vec.size();b++){
		vector < vector < double > > gijoe_matrix_b;
		vector < double > empty_gijoe_matrix_b_1;
		gijoe_matrix_b.push_back(empty_gijoe_matrix_b_1);
		for (u=2;u<=max_num_brch_vec[b];u++){
			
			vector < double > gijoe_matrix_b_u;
				for (v=1;v<=u;v++){
					double gi_joe_buv;
					//if (brchlens_vec[b]==0){
					if (b == brchlens_vec.size()-1){
						gi_joe_buv=1;
						}
					else{
						gi_joe_buv=0;
						for(k=v; k <= u; k++) {
							v1 = 1.0;
							for(h = v; h <= v+k-2; h++){
								v1 *= h;
							}
							v2 = 1.0; 
							for(h = v; h > 1; h--){         
								v2 *= h;
							}
							v3 = 1.0;
							for(h = k-v; h > 1; h--){
								v3 *= h;
							}
							u1 = 1.0;
							for(h = u; h >= u - k + 1; h--){
								u1 *= h;
							}
							u2 = 1.0;
							for(h = u; h <= u+k-1; h++){
								u2 *= h;
							}
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
		//gijoe_matrix.push_back(gijoe_matrix_b);
		mat.push_back(gijoe_matrix_b);
		//cout<<gijoe_matrix.size()<<endl;
	}
		for (size_t mat_b=0;mat_b<mat.size();mat_b++){
			for (size_t mat_i=0;mat_i<mat[mat_b].size();mat_i++){
				for (size_t mat_j=0;mat_j<mat[mat_b][mat_i].size();mat_j++){
					dout<<"g_["<<mat_i<<","<<mat_j<<"]("<<brchlens_vec[mat_b]<<")="<<mat[mat_b][mat_i][mat_j]<<"  ";
				
				}
			dout<<endl;
				
			}
			dout<<endl;
		}
	//return gijoe_matrix;
}

