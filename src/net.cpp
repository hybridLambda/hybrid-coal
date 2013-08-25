#include"net.hpp"

void Net::init(){
	this->set_is_net(false);
	this->set_max_rank(0);
}

Net::Net(){
	this->init();
}


/*! \brief Construct Net object from a (extended) Newick string */
Net::Net(string in_str /*! input (extended) newick form string */){
	this->init();
	//if (in_str.size()>0){
		checking_Parenthesis(in_str);
		net_str=checking_labeled(in_str);
		dout<<"Build Net from string: "<<net_str<<endl;
		dout<<"Extract nodes from net string"<<endl;
		this->create_node_list();	
		dout<<"Merge hybrid nodes:"<<endl;
		this->merge_repeated_nodes_as_one();
		dout<<"Linking nodes:"<<endl;
		this->linking_nodes();
		this->Net_nodes.back()->find_hybrid_descndnt();
		this->set_max_rank(Net_nodes.back()->ranking());
		this->is_net_func();
		assert(print_all_node());
	dout<<"Net of "<<net_str <<" is built"<<endl;

}

void Net::find_tips(){
	for (size_t i = 0; i < Net_nodes.size(); i++){
		if (Net_nodes[i]->num_child() == 0){Net_nodes[i]->set_tip(true);}
	}
}

void Net::create_node_list(){
	for (size_t i_str_len=1;i_str_len<net_str.size();){
		if (net_str[i_str_len]=='e' && (net_str[i_str_len+1]=='-' || net_str[i_str_len+1]=='+')){
			i_str_len++;
		}
		else{
			//if (isalpha(net_str[i_str_len])){
			if ( start_of_tax_name(net_str,i_str_len)){
			//if (isalpha(net_str[i_str_len]) || isdigit(net_str[i_str_len])){	
				size_t nodename_start_index = i_str_len;
				size_t nodename_end_index = end_of_label_or_bl(net_str, nodename_start_index);
				size_t node_content_start_idx;
				size_t node_content_end_idx;
				if (net_str[nodename_start_index-1]==')'){
					node_content_start_idx=Parenthesis_balance_index_backwards(net_str,nodename_start_index-1);
					node_content_end_idx=nodename_start_index-1;
					//size_t substr_len=str_start_index-rev_dummy_i;
					//node_content=net_str.substr(rev_dummy_i,substr_len);			
				}
				else {
					node_content_start_idx = nodename_start_index;
					node_content_end_idx = nodename_end_index;
				}
				i_str_len = nodename_end_index + 1;
				//i_str_len=label.size()+i_str_len;
				//cout<<node_content<<endl;
				//node_contents.push_back(node_content);
				//size_t found_bl=net_str.find(':');
				//if (found_bl!=string::npos){					
					//size_t found=min(min(net_str.find(",",i_str_len+1),net_str.find(")",i_str_len+1)),net_str.size());
					//brchlen=net_str.substr(i_str_len+1,found-i_str_len-1);
				
				//}
				//found_bl=net_str.find(":",found_bl+1);
				//brchlens.push_back(brchlen);
			
				Node * node_ptr ;
				if (net_str[i_str_len]==':'){
					i_str_len++;
					size_t bl_end_index = end_of_label_or_bl(net_str, i_str_len);
					istringstream bl_istrm(net_str.substr(i_str_len,bl_end_index-i_str_len-1));
					double bl;
					bl_istrm >> bl;
					node_ptr= new Node(nodename_start_index,nodename_end_index, node_content_start_idx, node_content_end_idx, bl);
				}
				else{
					node_ptr= new Node(nodename_start_index,nodename_end_index, node_content_start_idx, node_content_end_idx);
					}
				this->Net_nodes.push_back(node_ptr);
			}
			else {
				i_str_len++;
			}
		}
	}
	assert(this->print_all_node());
}

void Net::merge_repeated_nodes_as_one(){
		for (size_t i=1;i<Net_nodes.size()-1;i++){
			size_t j;
			//dout<<net_str.substr(Net_nodes[i]->nodeName_start(),Net_nodes[i]->nodeName_length())<<endl;
			for ( j=i+1;j<Net_nodes.size()-1;j++){
			//dout<<"    "<<net_str.substr(Net_nodes[j]->nodeName_start(),Net_nodes[j]->nodeName_length())<<endl;
			if (net_str.substr(Net_nodes[j]->nodeName_start(),Net_nodes[j]->nodeName_length()) == net_str.substr(Net_nodes[i]->nodeName_start(),Net_nodes[i]->nodeName_length())){
					//cout<<"reapeated"<<endl;
					dout<<net_str.substr(Net_nodes[j]->nodeName_start(),Net_nodes[j]->nodeName_length())<<endl;
					if (net_str[Net_nodes[j]->node_content_start()]=='('){
						Net_nodes[i]->set_node_content_start(Net_nodes[j]->node_content_start());
						Net_nodes[i]->set_node_content_end(Net_nodes[j]->node_content_end());
					}
	//				else{
					Net_nodes[i]->set_brchlen2(Net_nodes[j]->brchlen1());
					//Net_nodes[i].brchlen2=Net_nodes[j].brchlen1;
	//				}
					break;
				}
			}
			if (net_str.substr(Net_nodes[j]->nodeName_start(),Net_nodes[j]->nodeName_length()) == net_str.substr(Net_nodes[i]->nodeName_start(),Net_nodes[i]->nodeName_length())){
				delete Net_nodes[j];
				Net_nodes.erase(Net_nodes.begin()+j);
			//	Net_nodes_ptr.erase(Net_nodes_ptr.begin()+j);
			}
		}
}

/*! \brief free the meomory */
void Net::clear(){
	//tax_name.clear();
	Net_nodes.clear();
};


void Net::is_net_func(){
	//false stands for tree, true stands for net_work
	//bool is_net_return=false;
	for (unsigned int i=0;i<Net_nodes.size();i++){
		if (Net_nodes[i]->parent2()){
			//is_net_return=true;
			this->set_is_net(true);
			break;
		}
	}
	//return is_net_return;	/*!< if its net or not*/
}

bool Net::print_all_node(){
	//bool debug_switch=false;
	if ( this->is_net() ){
		dout<<"           label  hybrid hyb_des tpNode parent1  brchln1 parent2 brchln2 #child #dsndnt #id rank   e_num   Clade "<<endl;
		
		for (size_t i=0;i<Net_nodes.size();i++){
			//for (unsigned int j=0;j<descndnt[i].size();j++){
				//cout<<descndnt[i][j];
			//}
		assert(Net_nodes[i]->print_net_Node(net_str.c_str()));
			//Net_nodes[i].print_net_Node();
			//cout<<"  ";
			//for (unsigned int j=0;j<descndnt2[i].size();j++){
				//cout<<descndnt2[i][j];
			//}
			dout<<endl;
		}	
	}
	else{
		dout<<"label            tpNode   parent        brchln #child #dsndnt #id rank e_num   Clade "<<endl;
		for (size_t i=0;i<Net_nodes.size();i++){
			//for (unsigned int j=0;j<descndnt[i].size();j++){
				//cout<<descndnt[i][j];
			//}
			assert(Net_nodes[i]->print_tree_Node(net_str.c_str()));
			//Net_nodes[i].print_tree_Node();
						//cout<<"  ";
			//for (unsigned int j=0;j<descndnt2[i].size();j++){
				//cout<<descndnt2[i][j];
			//}
			dout<<endl;
		}	
	}
	return true;		
}


vector<string> extract_tax_names(vector<string> tip_name){
	
	
}


vector<string> Net::extract_tip_names(){
	vector <string> tip_name;
	for (size_t i=0;i<Net_nodes.size();i++){
		Node * current_node = Net_nodes[i];
		if(current_node->tip()){
			
			//if (Net_nodes[i].label.find("_")>0){
				////Net_nodes[i].name=Net_nodes[i].label.substr(0,Net_nodes[i].label.find("_"));
				//bool new_tax_bool=true;
				//for (size_t tax_i=0;tax_i<tax_name.size();tax_i++){
					//if (tax_name[tax_i]==Net_nodes[i].name){
						//new_tax_bool=false;
						//break;
					//}
				//}
				////if (new_tax_bool){
					////tax_name.push_back(Net_nodes[i].name);
				////}
			//}
			////else{
				////tax_name.push_back(Net_nodes[i].label);
			////}
			string tipname=net_str.substr(current_node->nodeName_start(), current_node->nodeName_length()+1);
			tip_name.push_back(tipname);
		}
	}
	
	//sort(tax_name.begin(), tax_name.end());
	sort(tip_name.begin(), tip_name.end());	
	return tip_name;
	
}

/*! \brief Checking Parenthesis of a (extended) Newick string */
void Net::checking_Parenthesis(string in_str){
	int num_b=0;
	for (size_t i=0;i<in_str.size();i++){
		if (in_str[i]=='('){
			num_b++;
		}
		if (in_str[i]==')'){
			num_b--;
		}
	}
	if (num_b!=0){
		cout<<"Error:"<<endl;
		cout<<in_str<<endl;
		cout<<"Parenthesis not balanced!"<<endl;
		exit (1);
	}
}

void Net::linking_nodes(){
	for (size_t i=0;i<Net_nodes.size();i++){
		if (net_str[Net_nodes[i]->node_content_start()] == '('){
			dout<<net_str.substr(Net_nodes[i]->nodeName_start(),Net_nodes[i]->nodeName_length())<<": ";			
			for (size_t ii = Net_nodes[i]->node_content_start()+1; ii < Net_nodes[i]->node_content_end(); ii++){
						//if (net_str[ii] =='(' || start_of_tax_name(net_str,ii)){
				if (net_str[ii] == '('){
					ii = Parenthesis_balance_index_forwards(net_str,ii)+1;
				}
				size_t ii_end = end_of_label_or_bl(net_str, ii);
				string child_node1_str = net_str.substr(ii, ii_end-ii+1);
				
				for (size_t j=0;j<Net_nodes.size();j++){
					if (child_node1_str == net_str.substr(Net_nodes[j]->nodeName_start(),Net_nodes[j]->nodeName_length())){
						dout<<child_node1_str<<" ";
						Net_nodes[j]->add_to_parent(Net_nodes[i]);
					}
				}
						//}			
			}
			dout<<endl;
		}
	}
	this->find_tips();
}





		//for (unsigned int i=0;i<Net_nodes.size();i++){
			//if(Net_nodes[i].label==Net_nodes[i].node_content){
				//if (Net_nodes[i].label.find("_")>0){
					//Net_nodes[i].name=Net_nodes[i].label.substr(0,Net_nodes[i].label.find("_"));
					////cout<<Net_nodes[i].name<<endl;
					//bool new_tax_bool=true;
					//for (unsigned int tax_i=0;tax_i<tax_name.size();tax_i++){
						//if (tax_name[tax_i]==Net_nodes[i].name){
							//new_tax_bool=false;
							//break;
						//}
					//}
					//if (new_tax_bool){
						//tax_name.push_back(Net_nodes[i].name);
					//}
					////cout<<tax_name.back()<<endl;
				//}
				//else{
					//tax_name.push_back(Net_nodes[i].label);
				//}
				//tip_name.push_back(Net_nodes[i].label);
			//}
		//}
		
		//sort(tax_name.begin(), tax_name.end());
		//sort(tip_name.begin(), tip_name.end());
		
			//if ( start_of_tax_name(net_str,i_str_len)){
			////if (isalpha(net_str[i_str_len]) || isdigit(net_str[i_str_len])){	
				//size_t nodename_start_index = i_str_len;
				//size_t nodename_end_index = end_of_label_or_bl(net_str, nodename_start_index);
			
			//}
			
			
			//char child_node1[Net_nodes[i]->node_content_length()];
			//size_t i_content_len;
			//size_t j_content_len;
			
			
			//for (i_content_len=1;i_content_len<Net_nodes[i]->node_content_length();){
				//if (Net_nodes[i]->node_content[i_content_len]=='(' ||  start_of_tax_name(Net_nodes[i].node_content,i_content_len) ){	
				////if (Net_nodes[i].node_content[i_content_len]=='(' ||  isalpha(Net_nodes[i].node_content[i_content_len]) || isdigit(Net_nodes[i].node_content[i_content_len]) ){	
					//if (Net_nodes[i].node_content[i_content_len]=='('){
						//j_content_len=Parenthesis_balance_index_forwards(Net_nodes[i].node_content,i_content_len)+1;
					//}
					//else{
						//j_content_len=i_content_len;
					//}
					//int child1_node_content_i=0;
					//for (;j_content_len<Net_nodes[i].node_content.length(); j_content_len++){
						//child_node1[child1_node_content_i]=Net_nodes[i].node_content[j_content_len];
						//char stop=Net_nodes[i].node_content[j_content_len+1];
						//if (stop==',' || stop==')' || stop==':'){
							//child_node1[child1_node_content_i+1]='\0';
							//break;}
						//child1_node_content_i++;
						//}
						//string child_node1_str=child_node1;		
						//i_content_len=j_content_len+2;
						
						
						
						//for (size_t j=0;j<Net_nodes.size();j++){
							//if (child_node1_str==Net_nodes[j].label){
								//add_node(Net_nodes_ptr[i],Net_nodes_ptr[j]);
							//}
						//}
				//}
				//else{i_content_len++;}
			//}	
	
		
		

		////cout<<descndnt.size()<<endl;
		
		//if (utility_debug_bool){
			//cout<<"Net::Net flag4"<<endl;
		//}
		//find_tip(Net_nodes_ptr.back());
		////cout<<"number of child " <<Net_nodes_ptr.back()->child.size()<<endl;
		////cout<<"node content "<<Net_nodes_ptr.back()->node_content<<endl;
		//find_hybrid_descndnt(Net_nodes_ptr.back());
		//max_rank=ranking(Net_nodes_ptr.back());
	
		//for (unsigned int i=0;i<Net_nodes.size();i++){
			//valarray <int> descndnt_dummy(0,tax_name.size());
			//descndnt.push_back(descndnt_dummy);
			//valarray <int> descndnt2_dummy(0,tip_name.size());
			//descndnt2.push_back(descndnt2_dummy);
			////cout<<Net_nodes_ptr[i]->name<<"  "<<Net_nodes_ptr[i]->label<<"  "<<Net_nodes_ptr[i]->tip_bool << " ";
			//for (unsigned int tax_name_i=0;tax_name_i<tax_name.size();tax_name_i++){
				////cout<<Net_nodes_ptr[i]->label<<"   "<<tax_name[tax_name_i]<<"  ";
				//if (find_descndnt(Net_nodes_ptr[i],tax_name[tax_name_i])){
					//descndnt[i][tax_name_i]=1;
				//}
				////cout<<descndnt[i][tax_name_i];
			//}
			////cout<<endl;
			//for (unsigned int tip_name_i=0;tip_name_i<tip_name.size();tip_name_i++){
				//if (find_descndnt2(Net_nodes_ptr[i],tip_name[tip_name_i])){
					//descndnt2[i][tip_name_i]=1;
				//}
			//}
			//Net_nodes[i].num_descndnt=descndnt[i].sum();
		//}
		
		//if (utility_debug_bool){
			//cout<<"Net::Net flag5"<<endl;
		//}
		
		////num_descndnt_interior(Net_nodes_ptr.back()); //replace by following
		//for (size_t i=0;i<Net_nodes_ptr.size();i++){
			//for (size_t j=0;j<Net_nodes_ptr.size();j++){
				//if (i!=j){
					//valarray <int> descndnt_diff=(descndnt[i]-descndnt[j]);
					////cout<<"hea"<<endl;
					//if (descndnt_diff.min() >= 0 && Net_nodes[i].rank > Net_nodes[j].rank && Net_nodes[j].rank>=2){
						//Net_nodes_ptr[i]->num_descndnt_interior=Net_nodes_ptr[i]->num_descndnt_interior+1;
						//Net_nodes_ptr[i]->descndnt_interior_node.push_back(Net_nodes_ptr[j]);
					//}
					
				//}
			
			//}
			////cout<<"checking #interior_des "<<Net_nodes_ptr[i]->label<<"  "<<Net_nodes_ptr[i]->num_descndnt_interior<<" "<<Net_nodes_ptr[i]->descndnt_interior_node.size()<<endl;
		//}
		
		////for (int rank_i=3;rank_i<=max_rank;rank_i++){
			
			
		////}
		
		
		
		
		//if (utility_debug_bool){
			//cout<<"Net::Net flag6"<<endl;
		//}	
		//int e_num_old=0;
		//enumerate_internal_branch(Net_nodes_ptr.back(),e_num_old);
		//if (utility_debug_bool){
			//cout<<"Net::Net flag6.5"<<endl;
		//}
		////cout<<descndnt.size()<<endl;
		////cout<<net_str<<endl;
		//for (unsigned int i=0;i<Net_nodes.size();i++){
			//for (unsigned int tax_name_i=0;tax_name_i<tax_name.size();tax_name_i++){
				////cout<<tax_name[tax_name_i]<<endl;
				//if (descndnt[i][tax_name_i] == 1){
					//if (Net_nodes[i].clade.size()==0){
						//Net_nodes[i].clade=tax_name[tax_name_i];
					//}
					//else{
						//Net_nodes[i].clade=Net_nodes[i].clade+tax_name[tax_name_i];
					//}
					//Net_nodes[i].clade.push_back('&');
				//}
			//}
			//Net_nodes[i].clade.erase(Net_nodes[i].clade.size()-1,1);
		//}
		
		//if (utility_debug_bool){
			//cout<<"Net::Net flag7"<<endl;
		//}
		////check for coaleased tips(& sign in the tips)
		//bool rewrite_descndnt=false;
		//for (unsigned int i=0;i<Net_nodes.size();i++){
			//if (Net_nodes[i].tip_bool ){
				//size_t found=Net_nodes[i].clade.find('&');
				//if (found!=string::npos){
					//rewrite_descndnt=true;
						//break;
				//}
			//}
			////if (rewrite_descndnt){
				////break;
			////}
		//}
		
		//if (rewrite_descndnt){
			//tax_name.clear();
			//size_t tax_name_start=0;
			//size_t tax_name_length=0;
			//for (unsigned int new_i_str=0;new_i_str<Net_nodes.back().clade.size();new_i_str++){
				//tax_name_length++;
				//if (Net_nodes.back().clade[new_i_str]=='&'){
					//tax_name_length--;
					//tax_name.push_back(Net_nodes.back().clade.substr(tax_name_start,tax_name_length));
					//tax_name_start=new_i_str+1;
					//tax_name_length=0;
				//}				
				//if (new_i_str==Net_nodes.back().clade.size()-1){
					//tax_name.push_back(Net_nodes.back().clade.substr(tax_name_start,tax_name_length));
				//}
			//}
			//sort(tax_name.begin(), tax_name.end());
			
			////cout<<descndnt.size()<<endl;
			//descndnt.clear();
		////	cout<<descndnt.size()<<endl;
			////~descndnt();
			//for (size_t i=0;i<Net_nodes.size();i++){
				//vector <string> contained_tips;
				//valarray <int> re_initial_descndnt(0,tax_name.size());
				//size_t tax_name_start=0;
				//size_t tax_name_length=0;
				//for (size_t new_i_str=0;new_i_str<Net_nodes[i].clade.size();new_i_str++){
					//tax_name_length++;
					//if (Net_nodes[i].clade[new_i_str] == '&'){
						//tax_name_length--;
						//contained_tips.push_back(Net_nodes[i].clade.substr(tax_name_start,tax_name_length));
						//tax_name_start=new_i_str+1;
						//tax_name_length=0;
					//}				
					//if (new_i_str==Net_nodes[i].clade.size()-1){
						//contained_tips.push_back(Net_nodes[i].clade.substr(tax_name_start,tax_name_length));
					//}
					
					////contained_tips.push_back(Net_nodes[i].clade.substr(tax_name_start,tax_name_length));
				//}
				
				//for (size_t tax_i=0;tax_i<tax_name.size();tax_i++){
					//for (size_t contained_tax_i=0;contained_tax_i<contained_tips.size();contained_tax_i++){
						//if (tax_name[tax_i]==contained_tips[contained_tax_i]){
							////descndnt[i][tax_i]=1;
							//re_initial_descndnt[tax_i]=1;
						//}
					//}
				//}
					
				//descndnt.push_back(re_initial_descndnt);
			//}
			
			//for (unsigned int i=0;i<Net_nodes.size();i++){
				////Net_nodes[i].clade="";
				//Net_nodes[i].clade.clear();
				//for (unsigned int tax_name_i=0;tax_name_i<tax_name.size();tax_name_i++){
					//if (descndnt[i][tax_name_i] == 1){
						////if (Net_nodes[i].clade == ""){
						//if (Net_nodes[i].clade.size()== 0){
							//Net_nodes[i].clade=tax_name[tax_name_i];
						//}
						//else{
							//Net_nodes[i].clade=Net_nodes[i].clade+tax_name[tax_name_i];
						//}
						//Net_nodes[i].clade.push_back('&');
					//}
				//}
				//Net_nodes[i].clade.erase(Net_nodes[i].clade.size()-1,1);				
			//}
		//}
		
		//is_net=is_net_func();
		
		
		////is_ultrametric=true;
		////is_ultrametric=is_ultrametric_func();

	////}
	////else{
		////descndnt.clear();
		////tax_name.clear();
		////Net_nodes.clear();
	////}
