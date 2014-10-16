

//this is not used!!!
bool Net::is_ultrametric_func(){
	bool is_ultrametric_return=true;
	vector <int> remaining_node(Net_nodes.size(),0);
	for (unsigned int node_i=0;node_i<Net_nodes.size();node_i++){
		remaining_node[node_i]=node_i;
	
	}
	//cout<<remaining_node.size()<<endl;
	//for (int rank_i=1;rank_i<=max_rank;rank_i++){
		//for (unsigned int node_i=0;node_i<Net_nodes.size();node_i++){
	int rank_i=1;
	unsigned int remaining_node_i=0;	
	while (remaining_node.size()>0){
		//remaining_node_i=0;
		int node_i=remaining_node[remaining_node_i];
		//cout<<remaining_node[remaining_node_i]<<endl;
		//rank_i=1;
		if (Net_nodes[node_i]->rank()==rank_i){
			if (rank_i==1){
				Net_nodes[node_i]->path_time.push_back(0.0);
			}
			else{
				for (size_t child_i = 0; child_i<Net_nodes[node_i]->num_child(); child_i++){
					double current_child_time;
					//if ((Net_nodes[node_i]->child[child_i]->parent1()->nodeName_start() == Net_nodes[node_i]->nodeName_start()) && (Net_nodes[node_i]->child[child_i]->parent1()->nodeName_end() == Net_nodes[node_i]->nodeName_end())){
					if (Net_nodes[node_i]->child[child_i]->parent1() == Net_nodes[node_i]){	
						current_child_time=Net_nodes[node_i]->child[child_i]->brchlen1();
					}
					else{
						//if (Net_nodes[node_i].child[child_i]->parent2->label==Net_nodes[node_i].label){
							current_child_time=Net_nodes[node_i]->child[child_i]->brchlen2();
						//}
						//else{
							//cout<<"warning!!!!! check code again"<<endl;
						//}
							
					}
					for (unsigned int child_i_time_i=0;child_i_time_i<Net_nodes[node_i]->child[child_i]->path_time.size();child_i_time_i++){
						Net_nodes[node_i]->path_time.push_back(current_child_time+Net_nodes[node_i]->child[child_i]->path_time[child_i_time_i]);
					}
				}
			}
			
			remaining_node.erase(remaining_node.begin()+remaining_node_i);
			//cout<<rank_i<<" "<<Net_nodes[node_i].label<<" "<<Net_nodes[node_i].path_time.size()<<endl;
		}
		else{
			remaining_node_i++;
		}
		//
		if (remaining_node_i==remaining_node.size()-1){
			rank_i++;
			remaining_node_i=0;
		}
	}
	//}
	for (unsigned int node_i=0;node_i<Net_nodes.size();node_i++){
		//bool node_i_path_time_eq=true;
		//cout<<Net_nodes[node_i].label<<"  ";
		for (unsigned int path_time_i=0;path_time_i<Net_nodes[node_i]->path_time.size();path_time_i++){
			if (pow((Net_nodes[node_i]->path_time[path_time_i]-Net_nodes[node_i]->path_time[0]),2)>0.000001){
				is_ultrametric_return=false;
				//cout<<Net_nodes[node_i].label<<endl;
				break;
			}
			//cout<<Net_nodes[node_i].path_time[path_time_i]<<" ";
		}
		//cout<<Net_nodes[node_i].label<<endl;
		//if (!is_ultrametric){
			//break;
		//}
		Net_nodes[node_i]->set_height(Net_nodes[node_i]->path_time[0]);
	}
	return is_ultrametric_return;
}




/*! \brief Rank network node from the bottom.
 * 
 * Child node has lower rank than the parent node. Tip nodes have rank one, the root node has the highest rank
 */
int ranking(Node *current){
	//int current_rank;
	if (current->tip()){
		current->set_rank(1);}
	else
	{
		int child_max_rank=0;
		for (size_t i_num_child=0;i_num_child<current->num_child();i_num_child++){
			int child_rank=ranking(current->child[i_num_child]);
			child_max_rank=max(child_rank,child_max_rank);
		}
		current->set_rank(child_max_rank+1);
	}		
	//return current_rank=current->rank;
	return current->rank();
}



/*! \brief Add child node to parent node */
void add_node(
Node *parent_node /*! pointer to the parent node*/, 
Node *child_node /*! pointer to the child node*/){
    parent_node->child.push_back(child_node);
	if (child_node->parent1()){
		child_node->set_parent2(parent_node);
		child_node->set_hybrid(true);
	}
	else {
		child_node->set_parent1(parent_node);
		//cout<<child_node->label<<" parent1 is "<<child_node->parent1->label<<endl;
	}
	//Node *kidspintr=parent_node->child[parent_node->num_child];
	//parent_node->num_child++;
} 

/*! \brief Label a node if its a descendant of a hybrid node */
void find_hybrid_descndnt(Node *current/*! pointer to a node*/){
	//if (current->num_child==0){
		//current->tip_bool=false;}
	//else {
	if (!current->tip()){
		for (size_t i_num_child=0; i_num_child<current->num_child();i_num_child++){
			if (current->hybrid() || current->descndnt_of_hybrid()){
				current->child[i_num_child]->set_descndnt_of_hybrid(true);
			}
			find_hybrid_descndnt(current->child[i_num_child]);
		}
	}
}


void Node::clear(){
	//clade=" ";
	clade.clear();
	//~clade();
	label=" ";
	node_content=" ";
	num_child=0;
	num_descndnt=0;
	parent1=NULL;
	parent2=NULL;
	brchlen1=0.0;
	brchlen2=0.0;
	rank=0;
	e_num=0;
	//visit=0;
	hybrid=false;
	e_num2=0;
	descndnt_of_hybrid=false;
	tip_bool=false;
	//label.clear();
	//node_content.clear();
	//child.clear();
	//descndnt.clear();
}


/*! \brief Write a fixed parameter into a externed newick formatted network string*/
string write_para_into_tree(string in_str /*! Externed newick formatted network string*/, 
double para /*! Coalescent parameter or fixed population sizes */){
	if (in_str.size()==0){
		cout<<"Define input tree (network) first!"<<endl;
		exit(1);;
	}
	Net para_Net(in_str);
	vector <Node*> para_Net_node_ptr;
	if (utility_debug_bool){
		cout<<"write_para_into_tree flag1"<<endl;
	}
	for (unsigned int node_i=0;node_i<para_Net.Net_nodes.size();node_i++){
		Node* new_node_ptr=NULL;
        para_Net_node_ptr.push_back(new_node_ptr);
        para_Net_node_ptr[node_i]=&para_Net.Net_nodes[node_i];
		para_Net_node_ptr[node_i]->brchlen1=para;
		if (para_Net.Net_nodes[node_i].hybrid){
			para_Net_node_ptr[node_i]->brchlen2=para;
		}
	}
	if (utility_debug_bool){
		cout<<"write_para_into_tree flag2"<<para<<endl;
	}
	para_Net_node_ptr.back()->brchlen1=para;
	rewrite_node_content(para_Net_node_ptr);
	string para_string=construct_adding_new_Net_str(para_Net);
	return para_string;	
}

string construct_adding_new_Net_str(Net in_Net){
	string out_str;
	out_str=in_Net.Net_nodes.back().node_content;
	out_str=out_str+in_Net.Net_nodes.back().label;
	if (in_Net.Net_nodes.back().brchlen1!=0){
		ostringstream brchlen_str;
		brchlen_str<<in_Net.Net_nodes.back().brchlen1;
		out_str=out_str+":"+brchlen_str.str();
	}
	out_str.push_back(';');
	return out_str;
}

