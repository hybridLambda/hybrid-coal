#include"node.hpp"




Node::Node(){
	this->init();
}

Node::Node(size_t nodeName_start,	size_t nodeName_end, size_t node_content_start, size_t node_content_end){
	this->init(nodeName_start,nodeName_end, node_content_start, node_content_end); 
}


Node::Node(size_t nodeName_start,	size_t nodeName_end, size_t node_content_start, size_t node_content_end, double bl1){
	this->init(nodeName_start,nodeName_end, node_content_start, node_content_end, bl1); 
}
	
Node::Node(size_t nodeName_start,	size_t nodeName_end, size_t node_content_start, size_t node_content_end, double bl1, double bl2){
	this->init(nodeName_start,nodeName_end, node_content_start, node_content_end, bl1, bl2); 
}
		
void Node::init(size_t nodeName_start, size_t nodeName_end, size_t node_content_start, size_t node_content_end, double bl1, double bl2){
	this->set_nodeName_start(nodeName_start);
	this->set_nodeName_end(nodeName_end);
	this->set_node_content_start(node_content_start);
	this->set_node_content_end(node_content_end);
	
	
	//this->set_label(label);
	//this->set_nodeContent(nodeContent);
	//num_descndnt=0;
	//num_descndnt_interior=0;
	this->set_parent1(NULL);
	this->set_parent2(NULL);
	
	this->set_brchlen1(bl1);
	this->set_brchlen2(bl2);
	
	this->set_rank(0);
	this->set_e_num(0);
	this->set_e_num2(0);
	//visit=0;
	this->set_hybrid(false);
	
	this->set_descndnt_of_hybrid(false);
	this->set_tip(false);
	

	//this->set_height(0.0);
	//prob_to_hybrid_left=1.0;
	this->set_visited(false);
	this->set_num_descndnt(0);
	//vector<Node*> descndnt_interior_node; /*!< \brief list of pointers to its descndent interior nodes */
	//vector<Node*> child; 
	
}


bool Node::print_net_Node(const char* net_str) const{
	for (size_t i=0;i<=13-(nodeName_end_ - nodeName_start_ );i++){dout<<" ";}
	for (size_t i=nodeName_start();i<=nodeName_end();i++){dout<<net_str[i];}
	dout<<setw(6)<<hybrid();
	dout<<setw(8)<<descndnt_of_hybrid();
	dout<<setw(8)<<tip();
	if (this->parent1()){
		for (size_t i= 0;i<=9-this->parent1()->nodeName_length();i++){dout<<" ";}
		for (size_t i = this->parent1()->nodeName_start(); i <= this->parent1()->nodeName_end();i++){dout<<net_str[i];}
		}
	else dout<<"          ";
	//dout<<setw (8)<<height();
	dout<<setw (8)<<brchlen1();
	if (this->parent2()){
		for (size_t i= 0;i<=7-this->parent2()->nodeName_length();i++){dout<<" ";}
		for (size_t i = this->parent2()->nodeName_start(); i <= this->parent2()->nodeName_end();i++){dout<<net_str[i];}
		}
	else dout<<"        ";

	dout<<setw (8)<<brchlen2();
	dout<<setw (7)<<num_child();
	dout<<setw (8)<<num_descndnt();
	dout<<setw(4)<<num_descndnt_interior();
//	dout<<setw (7)<<current.e_num;
//	dout<<setw (3)<<current.e_num2;
	dout<<setw (6)<<rank()<<"   ";
	//for (unsigned int i=0;i<descndnt.size();i++){
		//dout<<setw (1)<<descndnt[i];
	//}
	dout<<setw(2)<<e_num();
	dout<<setw(3)<<e_num2();
	//dout<<"    "<<clade;
	//dout<<endl;
	return true;
}


bool Node::print_tree_Node(const char* tree_str) const{
	for (size_t i=0;i<=12-(nodeName_end_ - nodeName_start_ );i++){dout<<" ";}	
	for (size_t i=this->nodeName_start();i<=this->nodeName_end();i++){dout<<tree_str[i];}
	dout<<setw(5)<<tip();
	if (this->parent1()){
		for (size_t i= 0;i<=11-this->parent1()->nodeName_length();i++){dout<<" ";}
		for (size_t i = this->parent1()->nodeName_start(); i <= this->parent1()->nodeName_end();i++){dout<<tree_str[i];}
		}
	else dout<<"            ";
	//dout<<setw (11)<<height();
	dout<<setw (11)<<brchlen1();
	dout<<setw (7)<<num_child();
	dout<<setw (8)<<num_descndnt();
	dout<<setw(4)<<num_descndnt_interior();
	dout<<setw (6)<<rank()<<"   ";
	//for (unsigned int i=0;i<descndnt.size();i++){
		//dout<<setw (1)<<descndnt[i];
	//}	
	dout<<setw(3)<<e_num();
	//dout<<"    "<<clade;
	//dout<<setw(3)<<num_descndnt_interior;
	return true;
}

Node::~Node(){
	for (size_t i = 0; i < descndnt_interior_node.size(); i++){descndnt_interior_node[i] = NULL;}
	for (size_t i = 0; i < child.size(); i++){child[i] = NULL;}
	//for (size_t i = 0; i < descndnt.size(); i++){descndnt[i] = NULL;}
}



//int enumerate_internal_branch(Node *current,int e_num_old){
/*! \brief enumerate the internal branches */
size_t Node::enumerate_internal_branch(size_t e_num_old) /*!< enumerator which is about to be updated \todo change e_num_old to int* type */ 
{
//needs modification, this is not correct.	
	
	
	//dout<<"Net::enumerate_internal_branch start"<<endl;
	
	int e_num_new;
	if (this->tip()){
		e_num_new=e_num_old;}
	else {
		if (this->visited()){
			e_num_new=e_num_old+1;
			this->set_e_num2(e_num_new);
			}
		else{
			int e_num_dummy=e_num_old;
			for (size_t i_child = 0; i_child < this->num_child(); i_child++){
				e_num_dummy=this->child[i_child]->enumerate_internal_branch(e_num_dummy);
			}
			this->set_visited(true);
			e_num_new=e_num_dummy+1;
			////if (current->e_num>0 && ){
				////current->e_num2=e_num_new;
			////}
			////else{
			//current->e_num=e_num_new;
			this->set_e_num(e_num_new);
			////}
			
		}
	}

	//dout<<"Net::enumerate_internal_branch end"<<endl;
	return e_num_new;
}


void Node::print_all_child(/*! pointer to the parent node*/){
    cout << this << " has " << this->num_child() << " kids" << endl;
    for (size_t child_i =0; child_i < this->num_child(); child_i++){
        cout<<this->child[child_i]<<endl;
    }
}


void Node::print_parent( /*! pointer to the child node*/){
	cout<< this <<" has parents "<<endl;
	cout<<this->parent1()<<" and "<<this->parent2() <<endl;
}


void Node::add_to_parent(Node* parent){
	parent->child.push_back(this);
	if (this->parent1()){
		this->set_parent2(parent);
		this->set_hybrid(true);
	}
	else {
		this->set_parent1(parent);
		//cout<<child_node->label<<" parent1 is "<<child_node->parent1->label<<endl;
	}
	
}


/*! \brief Rank network node from the bottom.
 * 
 * Child node has lower rank than the parent node. Tip nodes have rank one, the root node has the highest rank
 */
int Node::ranking() {
	if (this->rank() == 0){
		//int current_rank;
		if (this->tip()){
			this->set_rank(1);
			}
		else
		{
			int child_max_rank=0;
			for (size_t child_i = 0; child_i < this->num_child(); child_i++){
				int child_rank = this->child[child_i]->ranking();
				child_max_rank = max(child_rank,child_max_rank);
			}
			this->set_rank(child_max_rank+1);
			//current->rank=child_max_rank+1;
		}		
		//return current_rank=current->rank;
	}

	return this->rank();
}





/*! \brief Label a node if its a descendant of a hybrid node */
void Node::find_hybrid_descndnt(){
	//if (current->num_child==0){
		//current->tip_bool=false;}
	//else {
	//if (!current->tip()){
		for (size_t i_num_child=0; i_num_child<this->num_child();i_num_child++){
			if (this->hybrid() || this->descndnt_of_hybrid()){
				this->child[i_num_child]->set_descndnt_of_hybrid(true);
			}
			this->child[i_num_child]->find_hybrid_descndnt();
		}
	//}
}

bool find_descndnt(Node* current, string taxname){
	bool descndnt_found=false;
	if (current->tip()){
		//if (current->label==taxname){
			//cout<<current->name<<endl;
			//cout<<taxname<<endl;
		if (current->name==taxname){
			descndnt_found=true;
		}
		//else{
			//descndnt_found=false;
		//}
	}
	else {//int i;
		for (int i=0;i<current->num_child;i++){
			if (find_descndnt(current->child[i],taxname)){
				descndnt_found=true;
				break;
			}
			//else descndnt_found=false;
		}
	}	
	return descndnt_found;
}
