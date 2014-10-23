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

#include <iomanip>  // std::setw
#include <iostream> //std::cout
#include <vector>
#include <string>
#include <cassert>
#include <stdexcept>
#include <valarray>

using namespace std;

//Unless compiled with options NDEBUG, we will produce a debug output using 
//'dout' instead of cout and execute (expensive) assert statements.
#ifndef NDEBUG
#define dout std::cout
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && std::cout
#endif

#ifndef NODE
#define NODE
/*! \brief Node of a tree or network, it also represent the branch between this node and its parent node
 */

enum NAMETYPE { TAXA, TIP };

// use ChildContainer as a new class for the child nodes? we will need to include pointers to indicate the previous child and the next child.
//class NodeContainer;

class Node {
    friend class NodeIterator;
    friend class NodeContainer;
    friend class GraphBuilder;
    friend class Net;
    friend class CoalGT;
    friend class CoalST;
    friend class simTree;
    friend class HybridLambda;
    friend class Figure;
	public:
        ~Node();
        // Getters and Setters
        double height() const { return this->height_;} // This has no use for hybrid-coal
        void set_height ( double h ){ this->height_ = h; }// This has no use for hybrid-coal
        
        double brchlen1() const { return this->brchlen1_;}
        void set_brchlen1 ( double bl ){ this->brchlen1_ = bl; }
        
        double brchlen2() const { return this->brchlen2_;}
        void set_brchlen2 ( double bl ){ this->brchlen2_ = bl; }

        string label; /*!< \brief String label of a node, each node has unique label */
        //size_t node_index; /*!< \brief node index in the array, \todo use this more often!!!*/
        string node_content; /*!< \brief node content, the subtree string at this node */
        bool hybrid() const { return ( this->parent2() != NULL ) ;} /*!< \brief Hybrid node only, indicator of a hybrid node */
        
    private:    
        // Members
        size_t rank_;     /*!< \brief rank of the node, tip node has rank one, the root has the highest rank */        
        bool visited_;
        size_t e_num_;    /*!< \brief numbering the branch */
        size_t e_num2_;   /*!< \brief Hybrid node only, numbering the branch between the node and its second parent */
        double brchlen1_; /*!< \brief Branch length */
        double brchlen2_;/*!< \brief Hybrid node only, Branch length to the second parent*/

        valarray < size_t > taxa_below;
        valarray < size_t > samples_below; // tips_below
        vector < Node* > interior_nodes_below; /*!< \brief list of pointers to its descndent interior nodes */
        vector < Node* > child; /*!< \brief list of pointers to its child nodes */	
        Node* parent1_; /*!< \brief pointer to its parent node. */
        Node* previous_;
        Node* next_;

        string clade; /*!< \brief clade at this node, \todo this should be modified to a vector <string> */
        
        int num_descndnt; /*!< \brief number of the tip nodes, that are descendant from this node */
        int num_descndnt_interior; /*!< \brief number of the interior nodes, that are descendant from this node \todo to be replaced by interior_nodes_below.size()? */
        size_t NumberOfInteriorNodesBelow() const { return this->interior_nodes_below.size(); }
        vector <double> path_time; 
        double height_; /*!< \brief distance to the bottom of the tree */  // This has no use for hybrid-coal
            
        
        bool is_tip_; /*!< \brief Indicator of tip nodes. It's true, if it is a tip node, otherwise it is false. */
        bool is_tip() const { return this->is_tip_ ;}
        bool is_below_hybrid_; //bool descndnt_of_hybrid; /*!< \brief Indicator of descendant of hybrid nodes. It's true, if it is a descendant of hybrid nodes; false, otherwise. */
        bool is_below_hybrid() const { return this->is_below_hybrid_; }        
        
        Node* parent2_; /*!< \brief Hybrid node only, pointer to its second parent node. */
        //double prob_to_hybrid_left; /*!< \brief Hybrid node only, the probability that a lineage goes to the left */
        
        string name; /*!< \brief Name of a node, this is not unique for nodes. e.g. if its label is A_1, name is A */
        
        //vector <size_t> Net_node_contains_gt_node1; /*!< Used while simulation, check if a Network node contains a gene tree node */
        //vector <size_t> Net_node_contains_gt_node2; /*!< Used while simulation, check if a Network node contains a gene tree node */
                
        size_t e_num() const {return this->e_num_;}
        void set_enum( size_t num ) { this->e_num_ = num; }
    
        size_t e_num2() const {return this->e_num2_;}
        void set_enum2( size_t num ) { this->e_num2_ = num; }
    
        bool visited() const { return this->visited_; }
        void set_visited ( bool TorF ){ this->visited_ = TorF; }

        size_t rank() const { return this->rank_; }
        
        //bool tip() const { return this->tip_bool; }
        
        Node* parent1() const { return this->parent1_ ; }
        void set_parent1 ( Node * node ) { this->parent1_ = node; }

        Node* parent2() const { return this->parent2_ ; }
        void set_parent2 ( Node * node ) { this->parent2_ = node; }

        // NodeIterator
        Node* previous() const { return this->previous_ ; }
        void set_previous ( Node * node ) { this->previous_ = node; }
        
        Node* next() const { return this->next_; }
        void set_next ( Node * node ) { this->next_ = node; }
        
        bool is_first() const { return( this->previous_ == NULL ); }
        bool is_last()  const { return( this->next_ == NULL ); }
        
        // Methods    
        //Node (); /*!< \brief Initialize Node class*/
        Node ( size_t max_of_taxa,
               size_t max_of_sample, // number of tip
               string label,
               string content,
               double bl,
               bool tip );
                          
        void init();
        void add_child( Node *child_node /*! pointer to the child node*/);
        void CalculateRank();
        void print( bool is_Net = false );
        bool print_dout( bool is_Net = false );
        //void find_tip();
        void find_hybrid_descndnt();
        //bool find_descndnt ( string &name, NAMETYPE type );
        
        double extract_hybrid_para(){
            size_t hash_index = this->label.find('#');
            return strtod( this->label.substr(hash_index+1).c_str(), NULL) ;
        }   
};


#endif    
