#include"node.hpp"



void Node::init(){
	this->set_label(label);
	this->set_nodeContent(nodeContent);
	//num_descndnt=0;
	//num_descndnt_interior=0;
	this->set_parent1(NULL);
	this->set_parent2(NULL);
	
	this->set_brchlen1(0.0);
	this->set_brchlen2(0.0);
	
	this->set_rank(0);
	e_num=0;
	e_num2=0;
	//visit=0;
	hybrid=false;
	e_num2=0;
	descndnt_of_hybrid=false;
	tip_bool=false;
	//clade=" ";
	absolute_time=0.0;
	//prob_to_hybrid_left=1.0;
	visited = false;
	//double pAC=0.0;
	vector<Node*> descndnt_interior_node; /*!< \brief list of pointers to its descndent interior nodes */
	vector<Node*> child; 
	
}

Node::Node(){
	label=" ";
	node_content=" ";
	num_child=0;
	num_descndnt=0;
	num_descndnt_interior=0;
	parent1=NULL;
	parent2=NULL;
	brchlen1=0.0;
	brchlen2=0.0;
	rank=0;
	e_num=0;
	e_num2=0;
	//visit=0;
	hybrid=false;
	e_num2=0;
	descndnt_of_hybrid=false;
	tip_bool=false;
	//clade=" ";
	absolute_time=0.0;
	//prob_to_hybrid_left=1.0;
	visited = false;
	//double pAC=0.0;
	vector<Node*> descndnt_interior_node; /*!< \brief list of pointers to its descndent interior nodes */
	vector<Node*> child; 
}




void Node::print_net_Node(){
	cout<<setw(12)<<label;
	cout<<setw(6)<<hybrid;
	cout<<setw(8)<<descndnt_of_hybrid;
	cout<<setw(5)<<tip_bool;
	if (parent1){cout<<setw (11)<<(parent1->label);}
	else cout<<"           ";
	cout<<setw (8)<<absolute_time;
	cout<<setw (8)<<brchlen1;
	if (parent2){cout<<setw (9)<<parent2->label;}
	else cout<<"         ";
	cout<<setw (8)<<brchlen2;
	cout<<setw (7)<<num_child;
	cout<<setw (8)<<num_descndnt;
	cout<<setw(4)<<num_descndnt_interior;
//	cout<<setw (7)<<current.e_num;
//	cout<<setw (3)<<current.e_num2;
	cout<<setw (6)<<rank<<"   ";
	//for (unsigned int i=0;i<descndnt.size();i++){
		//cout<<setw (1)<<descndnt[i];
	//}
	cout<<setw(2)<<e_num;
	cout<<setw(3)<<e_num2;
	cout<<"    "<<clade;
	//cout<<endl;
}


void Node::print_tree_Node(){
	cout<<setw (12)<<label;
	cout<<setw(5)<<tip_bool;
	if (parent1){cout<<setw (11)<<(parent1->label);}
	else cout<<"           ";
	cout<<setw (11)<<absolute_time;
	cout<<setw (11)<<brchlen1;
	cout<<setw (7)<<num_child;
	cout<<setw (8)<<num_descndnt;
	cout<<setw(4)<<num_descndnt_interior;
	cout<<setw (6)<<rank<<"   ";
	//for (unsigned int i=0;i<descndnt.size();i++){
		//cout<<setw (1)<<descndnt[i];
	//}	
	cout<<setw(3)<<e_num;
	cout<<"    "<<clade;
	//cout<<setw(3)<<num_descndnt_interior;
}

	//~Node(){
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
	label.clear();
	node_content.clear();
	child.clear();
	//descndnt.clear();
}
