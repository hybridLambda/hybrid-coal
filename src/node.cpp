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

#include "node.hpp"

Node::Node ( size_t max_of_descendant, // number of tip
             string label,
             string content,
             double bl,
             bool tip ){
	this->init();
    this->label = label;
	this->node_content = content;
	this->brchlen1_ = bl;    
    this->is_tip_ = tip;
	//clade=" ";
    
    this->descendant = valarray < size_t > ( (size_t)0, max_of_descendant );
}

void Node::init(){
	this->parent1_ = NULL;
	this->parent2_ = NULL;
    this->previous_ = NULL;
    this->next_     = NULL;
	this->brchlen2_ = 0.0;
	this->rank_     = 0;
	this->e_num_    = 0;
	this->e_num2_   = 0;
	this->height_ = 1.0/0.0;
	//prob_to_hybrid_left=1.0;
	this->visited_ = false;    
    this->is_below_hybrid_ = false;
	num_descndnt=0;
	num_descndnt_interior=0;
}


//Node::Node(){
	//label="";
	//node_content="";
	////num_child=0;
	//num_descndnt=0;
	//num_descndnt_interior=0;
	//this->parent1_ = NULL;
	//this->parent2_ = NULL;
    //this->previous_ = NULL;
    //this->next_     = NULL;
	//this->brchlen1_ = 0.0;
	//this->brchlen2_ = 0.0;
	//this->rank_     = 0;
	//this->e_num_    = 0;
	//this->e_num2_   = 0;
	////visit=0;
	//descndnt_of_hybrid=false;
	//this->is_tip()=false;
	////clade=" ";
	//this->height_ = 1.0/0.0;
	////prob_to_hybrid_left=1.0;
	//this->visited_ = false;
//}


void Node::print( bool is_Net ){
    cout << setw(7) << this->label;
    cout << setw(12) << (this);
	if ( is_Net ) cout << setw(6) << this->hybrid();
    if ( is_Net ) cout << setw(8) << this->is_below_hybrid();
	cout << setw(5) << this->is_tip();
    if ( this->parent1() != NULL ) cout << setw (11) << ( this->parent1() );//if (this->parent1) cout << setw (11) << (parent1->label);
    else cout << "           ";
	cout << setw (12) << this->height();
	cout << setw (12) << this->brchlen1();
    if (is_Net){
        if ( this->parent2() != NULL) cout << setw (11) << ( this->parent2() ); //if (this->parent2) cout << setw (11) << (parent2->label);
        else cout << "           ";
        cout<<setw (12) << this->brchlen2();
    }
	cout << setw (7) << this->child.size();
	cout << setw (8) << num_descndnt;
	cout << setw(4) << num_descndnt_interior;
	cout << setw(6) << this->rank() << "   ";
	
	cout << setw(2)<<this->e_num();
	if ( is_Net ) cout << setw(3) << this->e_num2();
	cout << "    " << this->clade;
    for ( size_t i = 0; i < this->descendant.size(); i++ ){
		cout<<this->descendant[i];
	}
	//cout<<endl;
}

/*! \brief Add child node to parent node */
void Node::add_child( Node *child_node /*! pointer to the child node*/){
    this->child.push_back(child_node);
	if ( child_node->parent1() != NULL ){
		child_node->set_parent2 ( this );
	}
	else child_node->set_parent1 ( this );
} 


/*! \brief Rank network node from the bottom.
 * 
 * Child node has lower rank than the parent node. Tip nodes have rank one, the root node has the highest rank
 */
void Node::CalculateRank(){
	if ( this->is_tip() ) this->rank_ = 1;
	else {
		size_t child_max_rank = 0;
		for ( size_t ith_child = 0; ith_child < this->child.size(); ith_child++ ){
            this->child[ith_child]->CalculateRank();			
			child_max_rank = max( child_max_rank, this->child[ith_child]->rank() );
		}
		this->rank_ = child_max_rank + 1;
	}
}

///*! \brief Label a node if its a descendant of a hybrid node */
//void Node::find_hybrid_descndnt(){
	//if ( this->is_tip() ) return;
    //for ( size_t ith_child = 0; ith_child < this->child.size(); ith_child++){
        //if ( this->hybrid() || this->descndnt_of_hybrid ) this->child[ith_child]->descndnt_of_hybrid = true;
        //this->child[ith_child]->find_hybrid_descndnt();
    //}
//}


Node::~Node(){}
