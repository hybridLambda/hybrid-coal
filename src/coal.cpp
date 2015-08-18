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

#include "coal.hpp"

CoalST::CoalST ( string sp_str ) : GraphBuilder ( sp_str ){ 
    CoalSTdout << "Constrcut species tree: " << sp_str << endl; 
    this->which_taxa_is_below();
    this->which_sample_is_below();
    // TODO: consider, should these be done here??!!
    this->build_gijoe();
    this->building_S_matrix();
}


CoalST::CoalST ( const CoalST & spIn ) : GraphBuilder ( spIn ) {
    this->brchlens_vec = spIn.brchlens_vec;
    this->max_num_brch_vec = spIn.max_num_brch_vec;
    this->S_matrix = spIn.S_matrix;
    this->gijoemat = spIn.gijoemat;
}


void CoalST::assign_bl_to_vec(){
    CoalSTdout << " Start of assign_bl_to_vec()" << endl;
    assert( this->print_all_node_dout() );
    size_t numberOfInternalBranches = this->nodes_.back()->edge1.name();
    this->brchlens_vec = vector < double > ( numberOfInternalBranches, 0.0 );
    this->max_num_brch_vec = vector < int > ( numberOfInternalBranches,0);

    for ( auto it = nodes_.iterator(); it.good(); ++it ){
        //Node * (*it) = &this->NodeContainer[node_i];
        size_t tmpEdge = (*it)->edge1.name();

        if ( tmpEdge != 0 ){
            size_t edgeIndex = tmpEdge - 1;
            brchlens_vec[edgeIndex]     = (*it)->edge1.bl();
            max_num_brch_vec[edgeIndex] = (*it)->samples_below.sum();
        }
        CoalSTdout << " Node " << (*it) << "  " << tmpEdge << " " <<(*it)->edge1.bl()<<endl;
        
        if ( (*it)->parent2() ){
            tmpEdge = (*it)->edge2.name();
            size_t edgeIndex = tmpEdge - 1;
            brchlens_vec[edgeIndex]     = (*it)->edge2.bl();
            max_num_brch_vec[edgeIndex] = (*it)->samples_below.sum();
            dout << edgeIndex << " " << (*it)->edge2.bl() <<endl;
        }
    }
    CoalSTdout << "assign_bl_to_vec() finished" << endl;
}


void CoalST::build_gijoe(){
    this->assign_bl_to_vec();
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


bool CoalST::print_gijoemat(){
    CoalSTdout << "gijoe matrix"<<endl;
    for ( size_t branch_idx = 0; branch_idx < this->gijoemat.size(); branch_idx++ ){
        for ( size_t branch_in = 0; branch_in < this->gijoemat[branch_idx].size(); branch_in++){
            CoalSTdout << " ";
            for (size_t branch_out = 0; branch_out < this->gijoemat[branch_idx][branch_in].size(); branch_out++){
                dout << "g_[" << branch_in << "," << branch_out << "](" << brchlens_vec[branch_idx] << ")=" <<setw(10)<< this->gijoemat[branch_idx][branch_in][branch_out] << "  ";
            }
            dout<<endl;
        }
    }
    dout<<endl;
    return true;
}


CoalGT::CoalGT ( string gt_str ) : GraphBuilder ( gt_str ){
    CoalSTdout << "Constrcut gene tree: " << gt_str << endl;
    this->which_taxa_is_below();
    this->which_sample_is_below();
}


void CoalGT::building_M_matrix( CoalST & sp_net ) {
    size_t spNumberOfInternalBranches = sp_net.nodes_.back()->edge1.name();
    size_t gtNumberOfInteriorNode = this->nodes_.back()->edge1.name();

    for ( int i_sp_enum = 0; i_sp_enum < spNumberOfInternalBranches; i_sp_enum++ ){
        vector < int > M_matrix_row( ( gtNumberOfInteriorNode - 1 ), 1 );
        this->M_matrix.push_back( M_matrix_row );
    }

    for ( auto sp_it = sp_net.nodes_.iterator(); sp_it.good(); ++sp_it ){
        if ( (*sp_it)->isTip() || (*sp_it)->parent1() == NULL ) continue;
        for ( auto gt_it = this->nodes_.iterator(); gt_it.good(); ++gt_it){
            if ( (*gt_it)->isTip() || (*gt_it)->parent1() == NULL ) continue;
            int des_diff_prod = 1;
            for ( size_t i_des = 0; i_des < (*sp_it)->taxa_below.size(); i_des++ ){
                int des_diff = (*sp_it)->taxa_below[i_des] - (*gt_it)->taxa_below[i_des]; /*! \todo this may not be right for multiple lineages per species*/
                if ( des_diff < 0 ){
                    des_diff_prod = 0;
                    M_matrix[ (*sp_it)->edge1.name()-1][(*gt_it)->edge1.name()-1] = 0;
                    break;
                } else{  /*! \todo check this else!!! dont think this is needed */
                    if (des_diff>0){
                        des_diff_prod=0;}
                }
            }

            if ( des_diff_prod == 1 ){
                M_matrix[ (*sp_it)->edge1.name()-1][(*gt_it)->edge1.name()-1] = 1;
            }
        }
    }
    
    dout<<"M matrix"<<endl;
    assert( print_2D_matrix( M_matrix ) );
}


void CoalGT::building_R_matrix( ){
    size_t gtNumberOfInteriorNode = this->nodes_.back()->edge1.name();
    for ( size_t i_gt_enum=0; i_gt_enum < gtNumberOfInteriorNode; i_gt_enum++ ){
        vector <int> R_matrix_row( i_gt_enum, 1 );
        this->R_matrix.push_back(R_matrix_row);
    }
    
    for ( auto i = nodes_.iterator(); i.good(); ++i){
    //for ( size_t i = 0; i < this->nodes_.size()-1; i++ ){
        if ( (*i)->isTip()) continue;
        auto j = nodes_.iterator();
        for ( ; ; ++j){
            if ( (*j)->next() == (*i) ) break;
            if ( (*j)->isTip() ) continue;
            for ( size_t i_des = 0; i_des < (*i)->taxa_below.size(); i_des++ ){
                int des_diff = (*i)->taxa_below[i_des] - (*j)->taxa_below[i_des];
                if ( des_diff < 0 ){
                    this->R_matrix[(*i)->edge1.name()-1][(*j)->edge1.name()-1]=0;
                    break;
                }
            }
        }
    }
    dout << "R matrix" << endl;
    assert( print_2D_matrix(R_matrix) );
}


void CoalST::building_S_matrix(){
    int spNumberOfInternalBranches = this->nodes_.back()->edge1.name();
    
    for ( int i_sp_enum = 0; i_sp_enum < spNumberOfInternalBranches; i_sp_enum++ ){
        vector < int > S_matrix_row( spNumberOfInternalBranches, 1 );
        S_matrix_row[i_sp_enum] = 0;
        this->S_matrix.push_back( S_matrix_row );
    }
        //for ( auto it = nodes_.iterator(); it.good(); ++it){

    for ( auto i = nodes_.iterator(); i.good(); ++i ){
        if ( (*i)->isTip() || (*i)->parent1() == NULL ) continue;
            
        for ( auto j = nodes_.iterator(); j.good(); ++j){
            if ( (*j)->isTip() || (*j)->parent1() == NULL) continue;
                
            for ( size_t i_des = 0; i_des < (*i)->taxa_below.size(); i_des++ ){
                int des_diff = (*i)->taxa_below[i_des] - (*j)->taxa_below[i_des];
                if ( des_diff < 0 ){
                    S_matrix[(*i)->edge1.name()-1][(*j)->edge1.name()-1] = 0;
                    break;
                }
            }
        }
        
    }
    assert( print_2D_matrix(S_matrix) );
}


void CoalGT::initialize_possible_coal_hist( CoalST & sp_net ){
    size_t spNumberOfInternalBranches = sp_net.nodes_.back()->edge1.name();
    size_t gtNumberOfInteriorNode = this->nodes_.back()->edge1.name();
    
    if ( gtNumberOfInteriorNode == 1 ){
        vector <size_t> coal_hist_vec;
        coal_hist_vec.push_back(1);
        coal_hist_mat.push_back(coal_hist_vec);
    } else{
        int max_coal_hist_num = 1;
        for ( size_t i_gt_enum = 0; i_gt_enum < gtNumberOfInteriorNode-1; i_gt_enum++ ){
            vector <size_t> coal_hist_vec;
            for ( size_t i_sp_enum = 0; i_sp_enum < spNumberOfInternalBranches; i_sp_enum++){
                if ( M_matrix[i_sp_enum][i_gt_enum] == 1 ){
                    coal_hist_vec.push_back(i_sp_enum+1);
                }
            }
            max_coal_hist_num *= (coal_hist_vec.size());
            this->coal_hist_mat.push_back(coal_hist_vec);
        }
    }
    assert ( print_2D_matrix (coal_hist_mat));
}


void CoalGT::sum_coalescent_history_prob( CoalST & sp_net ){
    this->probability = 0.0;
    double current_prob_of_hist ;
    
    for ( size_t hist_i = 0; hist_i < this->valid_coal_hist.size(); hist_i++ ){
        dout << "The "<<hist_i<<"th history"<<endl;
        current_prob_of_hist = 1.0;
        for ( size_t i = 0; i < all_w[hist_i].size(); i++ ) {
            //cout << "all_w[hist_i].size()"<<all_w[hist_i].size()<<endl;
            size_t current_enum;
            for ( auto it = sp_net.nodes_.iterator(); it.good(); ++it){
                //Node * current_node = &sp_net.NodeContainer[node_i];
                if ( (*it)->edge1.name() == ( i + 1 ) ){
                    //current_branch_lengths=my_sp_net.NodeContainer[node_i].brchlen1;
                    current_enum = (*it)->edge1.name();
                    break;
                }
            }    
            double current_gijoe = ( (i+1) == sp_net.nodes_.back()->edge1.name() ) ? 1 :
                                                                              sp_net.gijoemat[current_enum-1][num_enter[hist_i][i]-1][num_out[hist_i][i]-1];
            dout << " current_gijoe ["<<num_enter[hist_i][i]<<"]["<<num_out[hist_i][i]<<"]" <<current_gijoe << " all_w[hist_i][i] " <<all_w[hist_i][i] <<" all_d[hist_i][i] "<<all_d[hist_i][i] << endl;
            current_prob_of_hist *= ( (double)all_w[hist_i][i] / (double)all_d[hist_i][i] * current_gijoe );            
        }
        dout<<current_prob_of_hist<<endl;
        this->probability += current_prob_of_hist;
        //dout<<"gt prob = "<<probability<<endl;
    }
}


void CoalGT::enumerate_coal_events( CoalST & sp_net ){
    size_t spNumberOfInternalBranches = sp_net.nodes_.back()->edge1.name();

    //size_t spNumberOfInternalBranches = sp_net.nodes_.back()->edge();

    for (size_t i_coal_hist = 0; i_coal_hist < valid_coal_hist.size(); i_coal_hist++){
        vector <int> num_enter_i_coal;
        vector <int> num_out_i_coal;
        vector <int> num_coal_in_branch_i_coal;
        vector <double> w_i_coal;
        vector <double> d_i_coal;

        for ( size_t i = 0; i < spNumberOfInternalBranches; i++ ){
            for ( auto it = sp_net.nodes_.iterator(); it.good(); ++it){
                if ( !(*it)->isTip() && (*it)->edge1.name() == i + 1 ){
                    num_enter_i_coal.push_back( (*it)->samples_below.sum() );
                }
            }
        }

        for ( size_t i = 0; i < spNumberOfInternalBranches; i++ ){
            vector < int > clades_coal_in_branch;
            int num_coal_in_branch_i_coal_dummy = 0;
            for ( size_t coal_hist_position=0; coal_hist_position < valid_coal_hist[i_coal_hist].size(); coal_hist_position++ ){
                if ( valid_coal_hist[i_coal_hist][coal_hist_position] != i+1 ) continue;
                num_coal_in_branch_i_coal_dummy++;
                clades_coal_in_branch.push_back( coal_hist_position );
                for ( size_t i_S_mat = i; i_S_mat < (spNumberOfInternalBranches); i_S_mat++ ){
                    if ( sp_net.S_matrix[i_S_mat][i] != 1 ) continue;
                    num_enter_i_coal[i_S_mat]--;
                }    
            }
            num_coal_in_branch_i_coal.push_back( num_coal_in_branch_i_coal_dummy );
            num_out_i_coal.push_back( num_enter_i_coal[i] - num_coal_in_branch_i_coal[i] );
            int d_dummy=1;
            int coal_dummy=num_enter_i_coal[i];
            for ( int y = 0; y < num_coal_in_branch_i_coal_dummy; y++ ){
                //d_dummy=d_dummy*n_choose_k(coal_dummy,2);
                d_dummy *= coal_dummy*(coal_dummy-1)/2;
                coal_dummy--;
            }
            d_i_coal.push_back(d_dummy);
            
            int w_dummy = factorial( num_coal_in_branch_i_coal_dummy );
            vector < vector < int > > updated_R_mat;
            if ( clades_coal_in_branch.size() > 1 ){
                for (size_t updated_R_i=0;updated_R_i<clades_coal_in_branch.size();updated_R_i++){
                    vector < int > updated_R_mat_row;
                    for (size_t updated_R_j=0;updated_R_j<updated_R_i;updated_R_j++){
                        updated_R_mat_row.push_back(R_matrix[clades_coal_in_branch[updated_R_i]][clades_coal_in_branch[updated_R_j]]);
                    }
                    updated_R_mat.push_back(updated_R_mat_row);
                }
                for (size_t updated_R_mat_i=0;updated_R_mat_i<updated_R_mat.size();updated_R_mat_i++){
                    int sum_r = 1;
                    for (size_t updated_R_mat_j=0;updated_R_mat_j<updated_R_mat[updated_R_mat_i].size();updated_R_mat_j++){
                        sum_r += updated_R_mat[updated_R_mat_i][updated_R_mat_j];
                    }
                    w_dummy /= sum_r; //w_dummy=w_dummy/sum_r;
                }
            }    
            
            w_i_coal.push_back( w_dummy );
            //if (w_dummy/d_dummy!=1){
                //cout<<w_dummy<<"/"<<d_dummy;
            //}
            //if (i<sp_max_enum-1){
                //cout<<"p_{"<<num_enter_i_coal[i]<<num_out_i_coal[i]<<"}(lambda_{"<<i+1<<"})";    
            //}
        }
        //cout<<endl;
        num_enter.push_back( num_enter_i_coal );
        num_out.push_back( num_out_i_coal );
        //num_coal_in_branch.push_back(num_coal_in_branch_i_coal);
        all_w.push_back( w_i_coal );
        all_d.push_back( d_i_coal );
        valid_coal_hist[i_coal_hist].push_back( spNumberOfInternalBranches );
    }    
    
    dout<<"valid_coal_hist"<<endl;    
    assert( print_2D_matrix(this->valid_coal_hist) );
}            


vector < vector < size_t > > CoalGT::recur_coal_hist( vector < vector <size_t > > coal_hist, size_t  node_i ) {
    vector < vector <size_t> > new_coal_hist;
    for ( size_t  coal_hist_i = 0; coal_hist_i < coal_hist.size(); coal_hist_i++ ){
        vector <size_t > coal_hist_dummy;
        coal_hist_dummy = coal_hist[coal_hist_i];
        int  coal_hist_mat_node_i_i=0;
        for ( size_t j_R_mat = 0; j_R_mat < node_i; j_R_mat++ ){
            while ( R_matrix[node_i][j_R_mat] == 1 && coal_hist_mat[node_i][coal_hist_mat_node_i_i] < coal_hist_dummy[j_R_mat] ){
                coal_hist_mat_node_i_i++;
            }
        }
        for ( ; coal_hist_mat_node_i_i < coal_hist_mat[node_i].size(); coal_hist_mat_node_i_i++ ){
            vector <size_t> coal_hist_dummy_dummy;
            coal_hist_dummy_dummy = coal_hist_dummy;
            coal_hist_dummy_dummy.push_back( coal_hist_mat[node_i][coal_hist_mat_node_i_i] );
            new_coal_hist.push_back( coal_hist_dummy_dummy );
        }
    }
        
    if ( node_i < coal_hist_mat.size() - 1 ){
        node_i++;
        new_coal_hist = recur_coal_hist( new_coal_hist, node_i );
    }

    return new_coal_hist;
}


void CoalGT::build_coal_hist ( ){
    for ( size_t first_coal_mat_i = 0; first_coal_mat_i < coal_hist_mat[0].size(); first_coal_mat_i++ ){
        vector <size_t> coal_hist_dummy;
        coal_hist_dummy.push_back( coal_hist_mat[0][first_coal_mat_i] );
        valid_coal_hist.push_back( coal_hist_dummy );
    }

    int gt_max_enum = this->nodes_.back()->edge1.name();
    if ( gt_max_enum-1 > 1){
        valid_coal_hist = recur_coal_hist( valid_coal_hist, (size_t)1 );
    }
}


void CoalGT::prob_given_sp_tree ( CoalST & sp_tree ){
    this->building_R_matrix();
    this->building_M_matrix( sp_tree );
    this->initialize_possible_coal_hist( sp_tree );
    this->build_coal_hist();
    this->enumerate_coal_events( sp_tree );
    this->sum_coalescent_history_prob( sp_tree );
}


TmpSN::TmpSN ( string tmpStr ) : GraphBuilder ( tmpStr ){
    this->chooseRemoveNode();
}


void TmpSN::chooseRemoveNode() {
    this->toBeRemovedNodeIndex_ = this->nodes()->size()-1;
    size_t nodeIndex = 0;
    for ( auto it = nodes_.iterator(); it.good(); ++it ){
        if ((((*it)->isHybrid() + (*it)->isBelowHybrid()) * ( 1 - (*it)->isTip())) >=1 && (*it)->rank() < nodes()->at(toBeRemovedNodeIndex_)->rank()){ /*! \todo Check if a tip node should be removed or not */
            toBeRemovedNodeIndex_ = nodeIndex;
        }
        nodeIndex++;
    }

    if ( toBeRemovedNodeIndex_ == this->nodes()->size() - 1 ){
        toBeRemovedNodeIndex_ = -1;
    }

    return;
}


CoalSN::CoalSN ( string sp_str ) : CoalST ( sp_str ){
    this->initializeNetStrWizPriorList( sp_str );
}


void CoalSN::initializeNetStrWizPriorList ( string spStr ) {
    this->maple_bool_local = false;
    NetStrWizPrior initNetStr ( spStr );
    this->NetStrWizPriorList.push_back( initNetStr );
    this->currentSubNetworkIndex = 0;
}


void CoalSN::simplifyNetworks( string gt_string, bool mapleSymbolic, bool latexSymbolic ){
    // parse in currentSubNetworkIndex
    while ( currentSubNetworkIndex < this->NetStrWizPriorList.size() ){
        TmpSN current_net ( this->NetStrWizPriorList[currentSubNetworkIndex].netStr );
        int rm_node_index = current_net.toBeRemovedNodeIndex();
        if ( rm_node_index >=0 ){
            if ( current_net.nodes_.at((size_t)rm_node_index)->isHybrid () ){ // Removing a hybrid node 
                this->removeHnode(current_net, rm_node_index, this->NetStrWizPriorList[currentSubNetworkIndex], true, true );  
            } else{ // Removing a internal node which is below a hybrid node
                dout << " Removie internal node " << current_net.nodes_.at(rm_node_index) << endl;
                assert ( current_net.nodes_.at(rm_node_index)->isBelowHybrid () );
                this->removeSnode( gt_string, current_net, rm_node_index, this->NetStrWizPriorList[currentSubNetworkIndex], true, true );
            }
            this->NetStrWizPriorList.erase(this->NetStrWizPriorList.begin() + currentSubNetworkIndex);
        } else{
            currentSubNetworkIndex++;
        }
    }
}


void CoalSN::removeHnode ( TmpSN &tmpSN, size_t rmNodeIndex, NetStrWizPrior netStrWizPrior, bool mapleSymbolic, bool latexSymbolic ){
    dout <<  " Start of rm_H_node from " << netStrWizPrior.netStr << endl;
    
    new_node_name = tmpSN.nodes_.at(rmNodeIndex)->nodeName;    
    
    //current_root_enum=new_Net_wiz_prior_p.root_enum;
    //current_lambda_sum=new_Net_wiz_prior_p.lambda_sum;
    
    // TODO size_t newly declared....
    hashSingIdx = hybrid_hash_index(new_node_name);

    assert( hashSingIdx < new_node_name.size() );
    
    // TODO string newly declared....
    string hybrid_parameter = extract_hybrid_para_str(new_node_name);
    // TODO string newly declared....
    
    //string left_hybrid_parameter = hybrid_parameter;
    //if ( sybolicMode == MAPLE ){
        //cout << "nothing here" <<endl;
        ////left_hybrid_parameter_str="gamma["+new_node_name.substr(1,hashSingIdx-1)+"]";
    //}
    //else if ( sybolicMode == SYMB_LATEX ){
        //cout << "nothing here" <<endl;

        ////left_hybrid_parameter_str="\\gamma"+new_node_name.substr(1,hashSingIdx-1);
    //}
    //else {
        //cout << "nothing here" <<endl;

        //assert ( sybolicMode == NONE);
    //}
    // TODO rework on the symbolic mode
    //right_hybrid_parameter="(1-"+hybrid_parameter+")";
    //right_hybrid_parameter_str="(1-"+left_hybrid_parameter_str+")";
    left_hybrid_parameter_num = extract_hybrid_para(new_node_name);
    right_hybrid_parameter_num = 1.0 - left_hybrid_parameter_num;
    
    if ( tmpSN.nodes_.at(rmNodeIndex)->child.size() > 1 ){
        //block_rm_H=nchild_gt_one(rmNodeIndex,maple_bool_local,tmpSN.Net_nodes[rmNodeIndex].num_child);
        this->removeHnodeManyChild( tmpSN, rmNodeIndex, netStrWizPrior, false, false ); // TODO, change the mapleBool and latexBool
    } else{
        assert ( tmpSN.nodes_.at(rmNodeIndex)->child.size() == 1 );
        this->removeHnodeOneChild( tmpSN, rmNodeIndex, netStrWizPrior, false, false ); // TODO, change the mapleBool and latexBool
    }

    dout<<"        End of rm_H_node "<< tmpSN.nodes_.at(rmNodeIndex)->nodeName << " from " << netStrWizPrior.netStr << endl;
}



string CoalSN::removeHnodeOneChildCore(TmpSN &tmpSN, size_t rmNodeIndex, size_t removingFromParentIndex){        
    GraphBuilder localTmpSN(tmpSN);
    Node * removingNode = localTmpSN.nodes_.at(rmNodeIndex);
    Node * removingFromParent[2];
    removingFromParent[0] = removingNode->parent1();
    removingFromParent[1] = removingNode->parent2();
    
    for ( size_t parentIndex = 0; parentIndex < 2; parentIndex++ ){
        int removingChildIndex = -1;
        dout << "From parent " << removingFromParent[parentIndex]->nodeName << "( " << removingFromParent[parentIndex]->child.size() << " ) child, ";
        for ( size_t childIndex = 0; removingFromParent[parentIndex]->child.size(); childIndex++){
            if ( removingFromParent[parentIndex]->child[childIndex] == removingNode ){
                dout << "removing " << removingChildIndex << "th child" << endl;
                removingChildIndex = childIndex;
                break;
            }
        }
        assert ( removingChildIndex != -1 );
        removingFromParent[parentIndex]->child.erase(removingFromParent[parentIndex]->child.begin() + (size_t)removingChildIndex);    
    }
    removingNode->child[0]->setIsBelowHybrid(false);
    removingNode->child[0]->set_parent1 (NULL);
    removingNode->child[0]->set_parent2 (NULL);
    removingNode->child[0]->edge1.setLength( removingNode->child[0]->edge1.bl() + removingNode->edge1.bl());
    removingFromParent[removingFromParentIndex]->add_child(removingNode->child[0]);
    removingNode->child.clear();

    //check hybrid node has zero kids
    assert ( removingNode->child.size() == 0);

    removingNode->set_parent1(NULL);
    removingNode->set_parent2(NULL);
    
    // Cleanning up the removing node, and remove one child internal node
    localTmpSN.nodes_.remove(removingNode);
    localTmpSN.rewrite_subTreeStr();    
    localTmpSN.removeOneChildInternalNode();
    return localTmpSN.reWritesubTreeStrAtRoot();
}


void CoalSN::removeHnodeOneChild( TmpSN &tmpSN, size_t rmNodeIndex, NetStrWizPrior netStrWizPrior, bool mapleSymbolic, bool latexSymbolic ){

    for (size_t removingFromParentIndex = 0; removingFromParentIndex < 2; removingFromParentIndex++){
        string adding_new_Net_string = removeHnodeOneChildCore(tmpSN, rmNodeIndex, removingFromParentIndex);
        
        // Comment out the following line for now...
        //vector < vector < int > > new_lambda_sum=rm_one_child_interior_lambda_sum(adding_new_Net_string_enum, current_lambda_sum);                 

        double uni_hybrid_paramter_num = removingFromParentIndex == 0 ? left_hybrid_parameter_num : 
                                                                        right_hybrid_parameter_num;

        int leftPower = removingFromParentIndex == 0 ? 0 : 1;

        double current_omega = ( netStrWizPrior.prior.omega() * uni_hybrid_paramter_num );

        NetStrWizPrior newNetStrWizPrior(netStrWizPrior);
        newNetStrWizPrior.netStr = adding_new_Net_string;
        newNetStrWizPrior.prior.setOmega ( current_omega * newNetStrWizPrior.prior.omega() );
        if ( mapleSymbolic ) newNetStrWizPrior.prior.setMapleOmegaAtH ( uni_hybrid_paramter_num, leftPower, uni_hybrid_paramter_num, 1 - leftPower );
        if ( latexSymbolic ) newNetStrWizPrior.prior.setLatexOmegaAtH ( uni_hybrid_paramter_num, leftPower, uni_hybrid_paramter_num, 1 - leftPower );
        this->NetStrWizPriorList.push_back(newNetStrWizPrior);
    }
}


string CoalSN::removeHnodeManyChildCore(TmpSN &tmpSN, size_t rmNodeIndex, NetStrWizPrior netStrWizPrior, int numberOfChildAtRemovingNode, vector < valarray <int> > &A_matrix_ith_row, size_t A_matrix_i ){
    GraphBuilder localTmpSN(tmpSN);
    Node * removingNode = localTmpSN.nodes_.at(rmNodeIndex);
    Node * removingFromParent[2];
    removingFromParent[0] = removingNode->parent1();
    removingFromParent[1] = removingNode->parent2();

    string nameToLeft  = new_node_name.substr(0,hashSingIdx) + "L";
    string nameToRight = new_node_name.substr(0,hashSingIdx) + "R";
    double blToLeft  = removingNode->edge1.bl();
    double blToRight = removingNode->edge2.bl();

    size_t max_of_taxa = removingNode->taxa_below.size();
    size_t max_of_sample = removingNode->samples_below.size();

    Node * newLeft = new Node(max_of_taxa, max_of_sample, nameToLeft, string(""), blToLeft, false );
    Node * newRight= new Node(max_of_taxa, max_of_sample, nameToRight, string(""), blToRight, false );
    localTmpSN.nodes_.add_before(newLeft, localTmpSN.nodes_.back());
    localTmpSN.nodes_.add_before(newRight, localTmpSN.nodes_.back());

    for ( size_t parentIndex = 0; parentIndex < 2; parentIndex++ ){
        int removingChildIndex = -1;
        for ( int childIndex = 0; removingFromParent[parentIndex]->child.size(); childIndex++ ){
            if (removingFromParent[parentIndex]->child[childIndex] == removingNode){                                
                removingChildIndex = childIndex;
                break;
            }
        }
        assert ( removingChildIndex != -1 );
        removingFromParent[parentIndex]->child.erase(removingFromParent[parentIndex]->child.begin()+(size_t)removingChildIndex);
        if ( A_matrix_i > 1 ){
            for ( int A_matrix_i_i_i = 0; A_matrix_i_i_i < numberOfChildAtRemovingNode; A_matrix_i_i_i++ ){
                removingNode->child[A_matrix_i_i_i]->set_parent1( NULL );
                if (A_matrix_ith_row[parentIndex][A_matrix_i_i_i] == 1){
                    if ( parentIndex == 0){
                        newLeft->add_child (removingNode->child[A_matrix_i_i_i]);
                    } else {
                        assert ( parentIndex == 1);
                        newRight->add_child (removingNode->child[A_matrix_i_i_i]);
                    }
                }
            }
            if ( parentIndex == 0){
                removingFromParent[parentIndex]->add_child( newLeft );
            } else {
                assert ( parentIndex == 1);
                removingFromParent[parentIndex]->add_child( newRight );
            }
        }
        else if (A_matrix_i == parentIndex){
            for ( int A_matrix_i_i_i = 0; A_matrix_i_i_i < numberOfChildAtRemovingNode; A_matrix_i_i_i++ ){
                removingNode->child[A_matrix_i_i_i]->set_parent1( NULL );
                if ( parentIndex == 0){
                    newLeft->add_child (removingNode->child[A_matrix_i_i_i]);
                } else {
                    assert ( parentIndex == 1);
                    newRight->add_child (removingNode->child[A_matrix_i_i_i]);
                }
            }
            if ( parentIndex == 0){
                removingFromParent[parentIndex]->add_child( newLeft );
            } else {
                assert ( parentIndex == 1);
                removingFromParent[parentIndex]->add_child( newRight );
            }
        }
        
    }

    removingNode->set_parent1(NULL);
    removingNode->set_parent2(NULL);

    // Cleanning up the removing node, and remove one child internal node
    localTmpSN.nodes_.remove(removingNode);      // 1. Remove the node which was supposed to be removed
    localTmpSN.nodes_.back()->CalculateRank();   // 2. Calculate the node rank including the new nodes
    //localTmpSN.rewrite_subTreeStr();             // 3. Rewrite the sub tree string at each node
    localTmpSN.removeOneChildInternalNode();     // 4. Remove internal nodes who has only one child
    localTmpSN.rewrite_subTreeStr();             // 5. Rewrite the sub tree string at each node
    return localTmpSN.reWritesubTreeStrAtRoot();
}


void CoalSN::removeHnodeManyChild( TmpSN &tmpSN, size_t rmNodeIndex, NetStrWizPrior netStrWizPrior, bool mapleSymbolic, bool latexSymbolic ){
    int numberOfChildAtRemovingNode = (int)tmpSN.nodes_.at(rmNodeIndex)->child.size();
    vector < vector < valarray <int> > > A_matrix = this->build_h_child(numberOfChildAtRemovingNode);
    dout << " numberOfChildAtRemovingNode = " << numberOfChildAtRemovingNode << endl;
    for ( size_t A_matrix_i = 0; A_matrix_i < A_matrix.size(); A_matrix_i++ ){
        string adding_new_Net_string = this->removeHnodeManyChildCore( tmpSN, rmNodeIndex, netStrWizPrior, numberOfChildAtRemovingNode, A_matrix[A_matrix_i], A_matrix_i);
        dout << adding_new_Net_string <<endl;
        //string adding_new_Net_string_enum=nchild_gt_one_core(current_removing_net_string_enum,new_node_name,rm_node_index,numberOfChildAtRemovingNode, A_matrix[A_matrix_i],A_matrix_i);
        vector <int> left_or_right = this->build_left_or_right_vec(numberOfChildAtRemovingNode, A_matrix[A_matrix_i], A_matrix_i);
        //vector < vector < int > > new_lambda_sum=rm_one_child_interior_lambda_sum(adding_new_Net_string_enum, current_lambda_sum);
                
        //bool enum_false=false;
        //string adding_new_Net_string_updated=rm_one_child_interior_node(adding_new_Net_string,enum_false);    
        
        //bool enum_true=true;
        //string adding_new_Net_string_updated_enum=rm_one_child_interior_node(adding_new_Net_string_enum,enum_true);                        

        double current_omega = ( netStrWizPrior.prior.omega() * 
                                 pow(1.0*left_hybrid_parameter_num, 1.0*left_or_right[0]) * 
                                 pow(1.0*right_hybrid_parameter_num,1.0*left_or_right[1]) );

        NetStrWizPrior newNetStrWizPrior(netStrWizPrior);
        newNetStrWizPrior.netStr = adding_new_Net_string;
        newNetStrWizPrior.prior.setOmega ( current_omega );
        if ( mapleSymbolic ) newNetStrWizPrior.prior.setMapleOmegaAtH (left_hybrid_parameter_num, left_or_right[0], right_hybrid_parameter_num, left_or_right[1]);
        if ( latexSymbolic ) newNetStrWizPrior.prior.setLatexOmegaAtH (left_hybrid_parameter_num, left_or_right[0], right_hybrid_parameter_num, left_or_right[1]);
        this->NetStrWizPriorList.push_back(newNetStrWizPrior);
    }
}


vector <int> CoalSN::build_left_or_right_vec( int n_child,vector < valarray <int> > A_matrix_ith_row, size_t A_matrix_i ){
    vector <int> left_or_right(2, 0);
    for (size_t parent_i = 0; parent_i < 2; parent_i++ ){
        if ( A_matrix_i > 1 ){
            for ( int A_matrix_i_i_i = 0; A_matrix_i_i_i < n_child; A_matrix_i_i_i++ ){
                if ( A_matrix_ith_row[parent_i][A_matrix_i_i_i] == 1 ){
                    left_or_right[parent_i]++;
                }
            }
        }
        else{
            if (A_matrix_i == parent_i){
                for (int A_matrix_i_i_i=0;A_matrix_i_i_i<n_child;A_matrix_i_i_i++){
                    left_or_right[parent_i]++;
                }
            }
        }
    }
    return left_or_right;
}


size_t hybrid_hash_index(string &in_str){
    return in_str.find('#');
}
string extract_hybrid_para_str(string &in_str){
    size_t hash_index=hybrid_hash_index(in_str);
    return in_str.substr(hash_index+1);//,in_str.size()-1);
}

double extract_hybrid_para(string in_str){
    return strtod(extract_hybrid_para_str(in_str).c_str(), NULL);
}

// This need to return:
// 1. net string
// 2. indicator of which child were coalesced
NetStrWizPrior CoalSN::removeSnodeCore( TmpSN &tmpSN, size_t rmNodeIndex, NetStrWizPrior netStrWizPrior, vector < valarray <int> > &A_matrix_ith_row, int numberOfChildAtRemovingNode ){
    NetStrWizPrior newNetStrWizPrior(netStrWizPrior);
    newNetStrWizPrior.tmpCladeList.clear();
    GraphBuilder localTmpSN(tmpSN);
    Node * removingNode = localTmpSN.nodes_.at(rmNodeIndex);
    Node * removingFromParent = removingNode->parent1();
    //vector < vector <int> > new_lambda_sum=current_lambda_sum;

    int removingChildIndex = -1;
    for ( int childIndex = 0; removingFromParent->child.size(); childIndex++ ){
        if (removingFromParent->child[childIndex] == removingNode){
            removingChildIndex = childIndex;
            break;
        }
    }
    assert ( removingChildIndex != -1 );
    removingFromParent->child.erase(removingFromParent->child.begin()+(size_t)removingChildIndex);

    vector < Node* > remainingNodes;

    //bool coal_ed=false;
    
    //int brch_lambda=sp_nodes_ptr_rm_enum[rm_node_index]->brchlen1;

    //vector <int> new_prior_coal_hist=current_prior_coal_hist;
    //vector < valarray < int > > new_current_prior_coal_clades=current_prior_coal_list;
    //vector < valarray < int > > brand_new_current_prior_coal_clades;
    for ( size_t A_matrix_i_i = 0; A_matrix_i_i < A_matrix_ith_row.size(); A_matrix_i_i++ ){
        Node* remainingChild = NULL;    
        int howmany_coaled = 0;
        bool A_matrix_i_i_coaled = false;

        for ( int A_matrix_i_i_i = 0; A_matrix_i_i_i < numberOfChildAtRemovingNode; A_matrix_i_i_i++ ){
            if ( A_matrix_ith_row[A_matrix_i_i][A_matrix_i_i_i] == 1 ){
                if (remainingChild){
                    remainingChild->nodeName += "&" + removingNode->child[A_matrix_i_i_i]->nodeName;
                    //coal_ed = true;
                    howmany_coaled++;
                    A_matrix_i_i_coaled = true;
                    localTmpSN.nodes_.remove (removingNode->child[A_matrix_i_i_i]);
                }
                else{
                    remainingChild = removingNode->child[A_matrix_i_i_i];
                }
            }
        }
        remainingNodes.push_back( remainingChild );
        if ( A_matrix_i_i_coaled ){
            valarray <size_t> newClade = convertNewCoalChildToClade ( remainingChild->nodeName );
            //newNetStrWizPrior.prior.priorCladeList.push_back(newClade);
            newNetStrWizPrior.tmpCladeList.push_back(newClade);
        }                    
    }

    double blOfRemovingNode = removingNode->edge1.bl();

    for (size_t nodeIdx = 0; nodeIdx < remainingNodes.size(); nodeIdx++ ){
        remainingNodes[nodeIdx]->edge1.setLength(remainingNodes[nodeIdx]->edge1.bl() + blOfRemovingNode );
        remainingNodes[nodeIdx]->set_parent1(NULL);
        remainingNodes[nodeIdx]->set_parent2(NULL);
        removingFromParent->add_child(remainingNodes[nodeIdx]);
    }

    //string adding_new_Net_string=construct_adding_new_Net_str(current_removing_net);
    localTmpSN.nodes_.remove(removingNode);      // 1. Remove the node which was supposed to be removed
    localTmpSN.nodes_.back()->CalculateRank();   // 2. Calculate the node rank including the new nodes
    localTmpSN.rewrite_subTreeStr();             // 3. Rewrite the sub tree string at each node
    localTmpSN.removeOneChildInternalNode();     // 4. Remove internal nodes who has only one child
    localTmpSN.rewrite_subTreeStr();             // 5. Rewrite the sub tree string at each node

    newNetStrWizPrior.netStr = localTmpSN.reWritesubTreeStrAtRoot();

    return newNetStrWizPrior;
}


bool CoalSN::checkSpCoalValid( GraphBuilder &tmpGt, vector < valarray <size_t> > new_coal_clade ){
    tmpGt.which_sample_is_below();
    bool spIsValid = true;
    for ( size_t coal_clade_i = 0; coal_clade_i < new_coal_clade.size(); coal_clade_i++ ){
        bool currentCladeIsValid = false;
        for ( auto it = tmpGt.nodes_.iterator(); it.good(); ++it){
            valarray <bool> comp = ( new_coal_clade[coal_clade_i] == (*it)->samples_below );
            if ( comp.min() == true ){
                currentCladeIsValid = true;
                break;
            }
        }

        if ( !currentCladeIsValid ){ 
            spIsValid = false;
            break;
        }
    }
    return spIsValid;
}


valarray <size_t> CoalSN::convertNewCoalChildToClade ( string nodeName ){
    valarray <size_t> A_matrix_i_i_valarray ( (size_t)0, this->tip_name.size() ); // We need to look at the actual tip of the network.
    for ( size_t i = 0; i < this->tip_name.size(); i++ ){
        size_t found = nodeName.find(this->tip_name[i]);
        if ( found != string::npos){
            A_matrix_i_i_valarray[i] = 1;
        }
    }
    return A_matrix_i_i_valarray;
}


double CoalSN::computeNumOfRepTopo ( size_t uBranchIn, size_t vBranchOut, vector < valarray < size_t > > tmpCladeList, GraphBuilder &tmpGt ){
// Refer to thesis page 75
    double w = 1.0 * factorial( uBranchIn - vBranchOut );
    for ( size_t cladeIndex = 0; cladeIndex < tmpCladeList.size(); cladeIndex++){
        if (tmpCladeList[cladeIndex].sum() == 1) 
            continue;

        for ( auto it = tmpGt.nodes_.iterator(); it.good(); ++it ){
            valarray <bool> comp = (tmpCladeList[cladeIndex] == (*it)->samples_below );
            if ( comp.min() == true ){
                // j = 1
                w /= ( 1.0 * ( 1 + (*it)->NumberOfInteriorNodesBelow() ) );

                if ( (*it)->NumberOfInteriorNodesBelow() > 1 ){
                    assert ( (uBranchIn - vBranchOut) > 2 );
                    for ( size_t interiorIndex = 0; interiorIndex < (*it)->interior_nodes_below.size(); interiorIndex++ ){
                        double aj = (*it)->interior_nodes_below[interiorIndex]->NumberOfInteriorNodesBelow();
                        w /= ( 1.0 + aj);
                    }
                }
            }
        }
    }
    return w;
}

double CoalSN::computeWaysToCoal( size_t uBranchIn, size_t vBranchOut ){
    // Refer to thesis, page 75, Eq.4.3
    double c = 1;
    for ( size_t i = uBranchIn; i > vBranchOut; i-- ){
        c *= (double)i * ((double)i-1.0) / 2.0;
    }
    return c;
}


void CoalSN::removeSnode( string gtStr, TmpSN &tmpSN, size_t rmNodeIndex, NetStrWizPrior netStrWizPrior, bool mapleSymbolic, bool latexSymbolic ){

    /*!
     * 1. check if the current tmpSN contains a coalesced Tip node
     *       If so, check the previous coalscent history of such node, whether it is consistant with the coalescet process of the gene tree.
     *          If so, merge the gene tree, then carry on.
     *          Otherwise, return
     *       Otherwise, carry on to step 2
     * 2. For each possibility of new topology, for the case of that coalsecent event has happened, check if gt coalsecnt could actually happen.
     *       If so, record the coalescent event, update the weight
     *       Othersiwe, return
     */

    //NetStrWizPrior new_Net_wiz_prior_p = this->NetStrWizPriorList[currentSubNetworkIndex];
    //vec_Net_wiz_prior_p my_rmd_networks;
    //Net current_net;
    
    //Net pre_current_net(new_Net_wiz_prior_p.s_net_string);

    
    //Net my_gt_tree(gt_str);
    //string updated_gt_str;
    
    //string net_str=new_Net_wiz_prior_p.s_net_string;
    //vector <string> original_tax_name=pre_current_net.tax_name;
    
    //bool multiple_sp_tip=false;
    //bool gt_tip_already_coaled=false;
    //if (gt_str.size()>0){
        //multiple_sp_tip=check_multiple_sp_tip(pre_current_net);
        //gt_tip_already_coaled=check_gt_tip_already_coaled(my_gt_tree);
    //}
    
    bool sp_valid_out = true;
    
    ////check if my_Net is valid    
    //if (multiple_sp_tip && gt_str.size()>0){
        //sp_valid_out=check_sp_valid(pre_current_net,my_gt_tree);    
    //}

    if (!sp_valid_out) return;

    //if (multiple_sp_tip && !gt_tip_already_coaled){
        
        //string dummy_string=construct_adding_new_Net_str(my_gt_tree);
        //Net load_gt_tree(dummy_string);
        ////vector <Node*> gt_nodes_ptr_rm;
        ////for (size_t i_link_ptr_rm=0;i_link_ptr_rm<load_gt_tree.Net_nodes.size();i_link_ptr_rm++){
            ////Node* new_node_ptr=NULL;
            ////gt_nodes_ptr_rm.push_back(new_node_ptr);
            ////gt_nodes_ptr_rm[i_link_ptr_rm]=&load_gt_tree.Net_nodes[i_link_ptr_rm];
        ////}
        //for (size_t i_gt_node=0;i_gt_node<load_gt_tree.descndnt.size();i_gt_node++){
            //for (size_t i_st_node=0;i_st_node<pre_current_net.descndnt.size();i_st_node++){
                //valarray <bool> comp = (load_gt_tree.descndnt[i_gt_node]==pre_current_net.descndnt[i_st_node]);
                //if (comp.min() && pre_current_net.Net_nodes[i_st_node].tip_bool ){
                    //gt_nodes_ptr_rm[i_gt_node]->num_child=0;
                    //gt_nodes_ptr_rm[i_gt_node]->child.clear();
                    ////gt_nodes_ptr_rm[i_gt_node]->rank=1;
                    //gt_nodes_ptr_rm[i_gt_node]->tip_bool=true;
                    ////cout<<gt_nodes_ptr_rm[i_gt_node]->label<<&(gt_nodes_ptr_rm[i_gt_node]->label)<<endl;
                    ////cout<<gt_nodes_ptr_rm[i_gt_node]->clade<<&(gt_nodes_ptr_rm[i_gt_node]->clade)<<endl;
                    //gt_nodes_ptr_rm[i_gt_node]->label=gt_nodes_ptr_rm[i_gt_node]->clade;
                    ////cout<<gt_nodes_ptr_rm[i_gt_node]->label<<&(gt_nodes_ptr_rm[i_gt_node]->label)<<endl;
                    //gt_nodes_ptr_rm[i_gt_node]->node_content.clear();
                //}
            //}        
        //}        
        ////load_gt_tree.print_all_node();
        //rewrite_node_content(gt_nodes_ptr_rm);
        //updated_gt_str=construct_adding_new_Net_str(load_gt_tree);            

        //string new_nt_str=construct_adding_new_Net_str(pre_current_net);
        //Net new_current_net(new_nt_str);
        //current_net=new_current_net;
    //}
    //else{
        //updated_gt_str=gt_str;
        //current_net=pre_current_net;
    //}
            
    //vector < valarray < int > > current_prior_coal_list=new_Net_wiz_prior_p.prior_clade_list;
    //vector <int> current_prior_coal_hist=new_Net_wiz_prior_p.prior_coal_hist;
    //vector < vector <int> > current_lambda_sum=new_Net_wiz_prior_p.lambda_sum;
    //double prior_prior_prob_num=new_Net_wiz_prior_p.omega;
    //string prior_prior_prob_string=new_Net_wiz_prior_p.omega_string;

    int numberOfChildAtRemovingNode = (int)tmpSN.nodes_.at(rmNodeIndex)->child.size();
    vector < vector < valarray <int> > > A_matrix = this->build_s_child( numberOfChildAtRemovingNode );
    
    //*** section of removing the descdent of the hybrid node.
    //string current_removing_net_string=construct_adding_new_Net_str(current_net);
    
    //string current_s_net_string_enum=new_Net_wiz_prior_p.s_net_string_enum;
    //int current_root_enum=new_Net_wiz_prior_p.root_enum;
    
        //gt_tree.print_all_node();
            //cout<<"gt_tree.print_all_node(); not changed?"<<endl;

    //cout << " A_matrix.size() = " << A_matrix.size()<<endl;
    
    
    for ( size_t A_matrix_i = 0; A_matrix_i < A_matrix.size(); A_matrix_i++ ){
        
        NetStrWizPrior newNetStrWizPrior = removeSnodeCore( tmpSN, rmNodeIndex, netStrWizPrior, A_matrix[A_matrix_i], numberOfChildAtRemovingNode );

        //if (gt_str.size()>0){
            ////sp_coal_valid=check_sp_coal_valid( coal_ed, gt_tree, new_current_prior_coal_clades);
            //sp_coal_valid=check_sp_coal_valid( coal_ed, my_gt_tree, new_current_prior_coal_clades);
        //}
        //else{
            //sp_coal_valid=true;
        //}

        GraphBuilder tmpGt(gtStr);
        bool spIsValid = checkSpCoalValid ( tmpGt, newNetStrWizPrior.tmpCladeList );
        
        // If the current sp is invalid, exit
        if ( !spIsValid ) continue;
        size_t numberBranchOut = A_matrix[A_matrix_i].size();
        double w = this->computeNumOfRepTopo ( numberOfChildAtRemovingNode, numberBranchOut, newNetStrWizPrior.tmpCladeList, tmpGt );
        double c = this->computeWaysToCoal ( numberOfChildAtRemovingNode, numberBranchOut );
        cout << " c = " << c <<endl;
        cout << " gijoe = " << gijoe (numberOfChildAtRemovingNode, numberBranchOut, tmpSN.nodes_.at(rmNodeIndex)->edge1.bl()) << endl; 
        double omega = w/c * gijoe (numberOfChildAtRemovingNode, numberBranchOut, tmpSN.nodes_.at(rmNodeIndex)->edge1.bl());
        cout << "omega is " << omega <<endl; 
        newNetStrWizPrior.prior.setOmega ( newNetStrWizPrior.prior.omega() * omega );
        // If it is valid, then compute omega
        

        this->NetStrWizPriorList.push_back(newNetStrWizPrior);
    }
    cout<<"        End of rm_S_node "<< tmpSN.nodes_.at(rmNodeIndex)->nodeName << " from " << netStrWizPrior.netStr << endl;
}

double CoalSN::gijoe(
size_t u, /*!< number of branch in */
size_t v, /*!< number of branch out */
double T) /*!< branch length*/
{
	//if (T==0){
		//return 1;
	//}
	//else{
	double sums=0;
	for (size_t k=v;k<=u;k++){
		double prods=exp(-k*(k-1)*T/2)*(2*k-1)*pow(-1.0,1.0*(k-v))/factorial(v*1.0)/factorial((k-v)*1.0)/(v+k-1);
		for (size_t y=0;y<k;y++){
			prods=prods*(v+y)*(u-y)/(u+y);
		}
		sums=sums+prods;
	}
	return sums;//}
}



vector <int> CoalSN::disjoint_list_s(int n, valarray <int> A_i,int i,vector <valarray <int> >A){
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

int CoalSN::disjoint_list_h( int n, int i, vector <valarray <int> >A ){
    valarray <int> A_i=A[i];
    int disjoint_list_h_return;
    for (size_t j=0;j<pow(2.0,1.0*n);j++){
        valarray <int> A_j = A[j];
        valarray <int> A_ij_sum = A_i + A_j;
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


vector < valarray <int> > CoalSN::all_possible_comb(int n){
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


// Todo, parse by reference, and use void type
vector < valarray <int> > CoalSN::rearrange_A( vector < valarray <int> > A, int n){
    vector <int> A_sum;
    int numLength = pow(1.0*2, 1.0*(n-1));
    for ( int i = 0; i < numLength; i++ ){
        A_sum.push_back(A[i].sum());
    }
     
    for ( int i = 0; i < numLength-1; i++){
        for ( int j = (i+1); j < numLength; j++){
            if (A_sum[i] < A_sum[j]){
                int temp = A_sum[i];                      // swap
                A_sum[i] = A_sum[j];
                A_sum[j] = temp;
                valarray<int> temp_valarray = A[i];       // swap
                A[i] = A[j];
                A[j] = temp_valarray;
            }
        }
    }
    return A;
}


vector < vector < valarray <int> > > CoalSN::build_s_child(int n){
    vector < valarray <int> > A_old = this->all_possible_comb(n);
    vector < valarray <int> > A = this->rearrange_A(A_old,n);
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
                all_compare=all_compare+compare; // all_compare += compare;
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



vector < vector < valarray < int > > > CoalSN::build_h_child( int n ){
    vector < vector < valarray < int > > > descdent_all_list;
    vector < valarray < int > > A = this->all_possible_comb( n );
    for (size_t i = 0; i < pow(1.0 * 2, 1.0 * n );i++){
        int i_comp = this->disjoint_list_h( n, i, A );
        vector < valarray < int > > current_descdent_list;
        valarray <int> current_descdent_list_i = A[i];
        valarray <int> current_descdent_list_i_comp = A[i_comp];
        current_descdent_list.push_back(current_descdent_list_i);
        current_descdent_list.push_back(current_descdent_list_i_comp);
        descdent_all_list.push_back(current_descdent_list);
    }
    return descdent_all_list;
}
