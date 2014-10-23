size_t Tree::first_coal_rank(){
    size_t min_rank = nodes_.back()->rank();
    for ( auto it = nodes_.iterator(); it.good(); ++it){
        if ( (*it)->is_tip() ) continue;
        min_rank = ( (*it)->rank() < min_rank ) ?  (*it)->rank() : min_rank ;
    }
    return min_rank;
}


//size_t Tree::first_coal_index (){    
    //size_t min_rank = this->first_coal_rank();
    //size_t dummy_index = this->nodes_.size()-1;
    //double min_coal_time = this->nodes_[dummy_index].height();
    //for (size_t i = 0 ; i < nodes_.size(); i++){
        //if ( this->nodes_[i].rank() == min_rank &&  this->nodes_[i].height() < min_coal_time ){
            //dummy_index = i;
            //min_coal_time = this->nodes_[dummy_index].height();
        //}        
    //}
    //return dummy_index;
//}


bool Node::find_descndnt ( string &name, NAMETYPE type ){	
	if ( this->is_tip() ) {
        string tmp = ( type == TAXA ) ? this->name : this->label;
        return ( tmp == name ) ? true : false;
    }
	else {
        bool descndnt_found = false;
		for (size_t i = 0; i < this->child.size(); i++ ){
            descndnt_found = this->child[i]->find_descndnt( name, type );
            if ( descndnt_found ) break;
		}
        return descndnt_found;
	}	
}


///*! \brief label a node if its a tip node */
//void Node::find_tip(){
	//if ( this->child.size() == 0) this->tip_ = true;
	//else {
		//for ( size_t ith_child = 0; ith_child < this->child.size(); ith_child++ ){
			//(*this->child[ith_child]).find_tip();
            ////this->child[ith_child]->find_tip();
		//}
	//}
//}


//void Tree::check_isUltrametric(){
    //vector <int> remaining_node( nodes_.size(), 0 );
    //for ( size_t node_i = 0; node_i < nodes_.size(); node_i++ ){
        //remaining_node[node_i] = node_i;
    //}
    //size_t rank_i = 1;
    //size_t remaining_node_i=0;    
    //while ( remaining_node.size() > 0 ){
        //int node_i = remaining_node[remaining_node_i];
        //if ( nodes_[node_i].rank() == rank_i ){
            //if (rank_i == 1) nodes_[node_i].path_time.push_back(0.0);
            //else{
                //for (size_t child_i = 0; child_i < nodes_[node_i].child.size(); child_i++ ){

                    //double current_child_time = (nodes_[node_i].child[child_i]->parent1->label==nodes_[node_i].label)?                    
                                                //nodes_[node_i].child[child_i]->brchlen1():
                                                //nodes_[node_i].child[child_i]->brchlen2();
                    //for (size_t child_i_time_i=0;child_i_time_i<nodes_[node_i].child[child_i]->path_time.size();child_i_time_i++){
                        //nodes_[node_i].path_time.push_back(current_child_time+nodes_[node_i].child[child_i]->path_time[child_i_time_i]);
                    //}
                //}
            //}            
            //remaining_node.erase(remaining_node.begin()+remaining_node_i);
        //}
        //else{
            //remaining_node_i++;
        //}

        //if ( remaining_node_i == remaining_node.size()-1 ){
            //rank_i++;
            //remaining_node_i=0;
        //}
    //}

    //for (size_t node_i=0;node_i<nodes_.size();node_i++){
        //for (size_t path_time_i=0;path_time_i<nodes_[node_i].path_time.size();path_time_i++){
            //if (pow((nodes_[node_i].path_time[path_time_i]-nodes_[node_i].path_time[0]),2)>0.000001){
                //this->is_ultrametric = false;
                //break;
            //}
        //}
        //nodes_[node_i].set_height( nodes_[node_i].path_time[0] );
    //}
//}



///////////////////////////////////////////////////////////////////////////////////////////// consider removed



void Tree::init_descendant(){
    for ( auto it = nodes_.iterator(); it.good(); ++it ){
        valarray <int> descndnt_dummy(0,tax_name.size());
        descndnt.push_back(descndnt_dummy);
        valarray <int> samples_below_dummy(0,tip_name.size());
        samples_below.push_back(samples_below_dummy);
        for ( size_t tax_name_i = 0; tax_name_i < tax_name.size(); tax_name_i++ ) descndnt[i][tax_name_i] = it->find_descndnt( tax_name[tax_name_i], TAXA) ? 1:0;
        for ( size_t tip_name_i = 0; tip_name_i < tip_name.size(); tip_name_i++ ) samples_below[i][tip_name_i] = it->find_descndnt( tip_name[tip_name_i], TIP) ? 1:0;
        it->num_descndnt = descndnt[i].sum();
    }

    for ( size_t i = 0; i < nodes_.size(); i++){
        for (size_t j = 0; j < nodes_.size(); j++){
            if ( i == j ) continue;

            valarray <int> descndnt_diff=(descndnt[i]-descndnt[j]);
            if (descndnt_diff.min() >= 0 && nodes_[i].rank() > nodes_[j].rank() && nodes_[j].rank() >= 2){
                this->nodes_[i].num_descndnt_interior += 1 ;
                this->nodes_[i].interior_nodes_below.push_back( &this->nodes_[j] );
            }                
        }
    }
}



void Tree::init_node_clade(){
    for ( auto it = nodes_.iterator(); it.good(); ++it){
        //if ( this->descndnt[i].sum() == 0 ) break;

        it->clade.clear();
        for ( size_t tax_name_i = 0; tax_name_i < tax_name.size(); tax_name_i++ ){
            if ( descndnt[i][tax_name_i] != 1) continue;

            nodes_[i].clade = ( nodes_[i].clade.size() == 0 ) ? tax_name[tax_name_i]:
                                                                            nodes_[i].clade + tax_name[tax_name_i];                                                                                        
            nodes_[i].clade.push_back('&');
        }
        nodes_[i].clade.erase(nodes_[i].clade.size()-1,1);
    }
}


	check_and_remove("log_file");
	
	//debug_bool=false;
	
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
	
	maple_bool=false;
	string maple_F_name="maple_prob.mw";
	
	bool print_gene_topo_bool=false;
	string gtopo_F_name="GENE_topo";
	
	bool help=false;
	string net_str;
	vector <string> gt_tree_str_s;
	
	bool all_gt_tree_bool=true;
	bool symb_bool=false;
	bool print_tree=false;
		
	bool list_sub_network_bool=false;

	int plot_option;//=0;
	bool plot_label=false;
	bool plot_branch=false;
	bool plot_bool=false;
	string tex_fig_name="texfigure.tex";
	bool dot_bool=false;
	string dot_fig_name="figure.dot";
	bool out_bool=false;
	string out_file="out_coal";

		
	if (help || argc==1){
		print_help();
	}
	else{
		
		if (net_str.size()<0){
			exit(1);
		}
		
		//process net_str, check if gt have multiple lineages
		if (!all_gt_tree_bool){
			Net test_gt(gt_tree_str_s[0]);
			//Net test_sp(net_str);
			if (test_gt.tip_name.size()>test_gt.tax_name.size()){
				net_str=rewrite_net_str(net_str, test_gt.tax_name,test_gt.tip_name );
			}
			for (size_t i=0;i<gt_tree_str_s.size();i++){
				gt_tree_str_s[i]=rm_dash_sign(gt_tree_str_s[i]);
				//cout<<gt_tree_str_s[i]<<endl;
			}
			cout<<net_str<<endl;
		}
		
		
		
		
		Net net_dummy(net_str);
		if (print_tree){
			net_dummy.print_all_node();
			appending_log_file("Tree printed");
		}
		plot_option=set_plot_option(plot_label,plot_branch);

		if (plot_bool){
			plot_in_latex_file(tex_fig_name.c_str(), net_dummy,plot_option);	
		}
		
		if (dot_bool){
			plot_in_dot(dot_fig_name.c_str(), net_dummy,plot_option);			
		}
			
		if (all_gt_tree_bool){
			gt_tree_str_s=all_n_tax_gene_tree(net_dummy.tip_name);
			if (print_gene_topo_bool){
				print_all_gt_topo(gtopo_F_name.c_str(),gt_tree_str_s);
			}
		}

		if (list_sub_network_bool && gt_tree_str_s.size()==1){
			list_sub(net_str,gt_tree_str_s[0]);
		}

		if (list_sub_network_bool && all_gt_tree_bool){
			string 	empty_str="";
			list_sub(net_str,empty_str);
		}
				
		if (print_tree || plot_bool || dot_bool || print_gene_topo_bool || list_sub_network_bool && gt_tree_str_s.size()>1){
			//return 0;
			return my_exit();
		}		
					
				
		if (latex_bool){
			
			latex_header(latex_F_name.c_str());	
			plot_in_latex(latex_F_name.c_str(), net_str,1);		
		}
		
		check_and_remove(out_file.c_str());
		
		if (net_dummy.is_net){		
		
			double total_prob=0;
			for (size_t topo_i=0;topo_i<gt_tree_str_s.size();topo_i++){
				string current_gene_topo=gt_tree_str_s[topo_i];
				//check_gt_str(current_gene_topo);

				if (main_debug_bool){
					cout<<current_gene_topo<<endl;
				}
				
				vec_Net_wiz_prior_p net_str_s=simplify_Networks(net_str,current_gene_topo);
				//vec_Net_wiz_prior_p * net_str_s_pt;
				//net_str_s_pt= new vec_Net_wiz_prior_p;
				//vec_Net_wiz_prior_p net_str_s=simplify_Networks(net_str,current_gene_topo);
				//net_str_s_pt=&net_str_s;
				
				if (latex_bool){
					latex_nt_body1(latex_F_name.c_str(),topo_i+1,current_gene_topo);
				}
				double total_prob_topo_i=0;
				for (size_t i=0;i<net_str_s.Net_vec.size();i++){
					gt_coal_in_st gt_in_sp(current_gene_topo,net_str_s.Net_vec[i]);
					//gt_coal_in_st gt_in_sp(current_gene_topo,net_str_s_pt->Net_vec[i]);
					if (main_debug_bool){
						cout<<current_gene_topo<<endl;
						//cout<<net_str_s_pt->Net_vec[i].s_net_string<<endl;
						cout<<topo_i<<"  "<<i<<"  "<<net_str_s.Net_vec[i].omega <<"  "<<net_str_s.Net_vec[i].s_net_string<<endl;
					}
					total_prob_topo_i=total_prob_topo_i+net_str_s.Net_vec[i].omega*gt_in_sp.probability;
					//total_prob_topo_i=total_prob_topo_i+net_str_s_pt->Net_vec[i].omega*gt_in_sp.probability;
					if (latex_bool){
						latex_nt_body2(latex_F_name.c_str(),net_str_s, i);
						gt_in_sp.latex_print(latex_F_name.c_str());
					}
					gt_in_sp.clear();
					net_str_s.Net_vec[i].clear();
				}
				
				prob_out_body(out_file.c_str(),topo_i+1,gt_tree_str_s[topo_i],total_prob_topo_i,net_dummy.tax_name.size());
				
				total_prob=total_prob+total_prob_topo_i;
				
				net_str_s.clear();
			}
			prob_out_tail(out_file.c_str(),gt_tree_str_s[0].size(), total_prob);
			if (main_debug_bool){
				cout<<"End of gene tree probability given species network"<<endl;
			}
		}
		else{
			if (main_debug_bool){
				cout<<"Start of gene tree probability given species tree"<<endl;
			}
			net_dummy.clear();
			net_str=rm_one_child_root(net_str);
			
			Net new_net_dummy(net_str);
		
			net_dummy=new_net_dummy;
			double total_prob=0;
			vector < valarray <int> > prior_coal_clades_dummy;
			vector <int> prior_coal_hist;

			
			
			double total_prob_topo_i=0;

			if (prior_coal_clades_dummy.size()>0){
				Net_wiz_prior_p coaled_sp_wiz_prior;
	Net_wiz_prior_p initial_Networks;	
	Net old_net_enum=initialize_enum_net(net_str);
	coaled_sp_wiz_prior.s_net_string=net_str;
	coaled_sp_wiz_prior.omega=1; // fix this!!!!!!!!!!!!!!!
	coaled_sp_wiz_prior.lambda_sum=initialize_lambda_sum(old_net_enum);
	coaled_sp_wiz_prior.s_net_string_enum=construct_adding_new_Net_str(old_net_enum);
	coaled_sp_wiz_prior.root_enum=old_net_enum.Net_nodes.back().e_num;
	coaled_sp_wiz_prior.prior_clade_list=prior_coal_clades_dummy;//=initial_prior_clade_list;
	coaled_sp_wiz_prior.prior_coal_hist=prior_coal_hist;//=initial_prior_coal_hist;
				
				for (size_t topo_i=0;topo_i<gt_tree_str_s.size();topo_i++){
					//gt_coal_in_st gt_in_sp(current_gene_topo,net_str_s.Net_vec[i]);
					string current_gene_topo=gt_tree_str_s[topo_i];
					gt_coal_in_st gt_in_sp(current_gene_topo, coaled_sp_wiz_prior);
					//total_prob_topo_i=total_prob_topo_i+net_str_s.Net_vec[i].omega*gt_in_sp.probability;
					
					if (latex_bool){
						//latex_nt_body2(latex_F_name.c_str(),net_str_s, 0);
						gt_in_sp.latex_print(latex_F_name.c_str());
					}
					total_prob_topo_i=total_prob_topo_i+coaled_sp_wiz_prior.omega*gt_in_sp.probability;
					
					prob_out_body(out_file.c_str(),topo_i+1,gt_tree_str_s[topo_i],gt_in_sp.probability,net_dummy.tax_name.size());
					total_prob=total_prob+total_prob_topo_i;
					gt_in_sp.clear();
				}
			
				
			}
			else{
				for (size_t topo_i=0;topo_i<gt_tree_str_s.size();topo_i++){
					Net current_gt(gt_tree_str_s[topo_i]);
					gt_coal_in_st gt_in_sp(current_gt,net_dummy);
					if (latex_bool){
						latex_tre_body(latex_F_name.c_str(),gt_tree_str_s[topo_i],net_str);
						gt_in_sp.latex_print(latex_F_name.c_str());
					}
					
					prob_out_body(out_file.c_str(),topo_i+1,gt_tree_str_s[topo_i],gt_in_sp.probability,net_dummy.tax_name.size());
				
					total_prob=total_prob+gt_in_sp.probability;
					gt_in_sp.clear();
				}
				
			}
			prob_out_tail(out_file.c_str(),gt_tree_str_s[0].size(), total_prob);
			
			if (main_debug_bool){
				cout<<"End of gene tree probability given species tree"<<endl;
			}
		}

		if (latex_bool){
			latex_tail(latex_F_name.c_str());
		}
		if (maple_bool){
			produce_maple(maple_F_name.c_str(),net_str, symb_bool,  gt_tree_str_s);
		}
		
	}
	
	if (pdf_bool){
		string command="pdflatex "+latex_F_name;
		system(command.c_str());
		//system("pdflatex latex_prob.tex");
	}
		return my_exit();
	//return 0;
