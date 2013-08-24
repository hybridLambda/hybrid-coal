#include<vector>
//#include<string>
#include<iostream>
#include<iomanip>
#include<cassert>

#ifndef NDEBUG
#define dout std::cout
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && std::cout
#endif


// do not use string here!!!

using namespace std;

#ifndef NODE
#define NODE

/*! \brief Node of a tree or network, it also represent the branch between this node and its parent node
 */
class Node {
	public:
		
		//vector<int> descndnt;
		vector<int> des_num_descndnt_interior;
		vector<Node*> descndnt_interior_node; /*!< \brief list of pointers to its descndent interior nodes */
		vector<Node*> child; /*!< \brief list of pointers to its child nodes */		
		//vector<Node*> descndnt;
		vector <double> path_time;
//valarray <int>	descndnt;
		
		size_t num_descndnt() const {return this->num_descndnt_;}; /*!< \brief number of the tip nodes, that are descendant from this node */
		void set_num_descndnt(size_t num){this->num_descndnt_ = num;};
		
		
		size_t num_descndnt_interior() const {return this->descndnt_interior_node.size();}; /*!< \brief number of the interior nodes, that are descendant from this node \todo to be replaced by descndnt_interior_node.size()? */
		size_t num_child() const {return this->child.size();}; /*!< \brief number of child \todo this can be replaced by child.size */


		
		//class Net * SubNetwork; /*!< \brief \todo pointer to the subtree of this node */

// Getter and setter
		Node* parent1() const {return this->parent1_;} /*!< \brief pointer to its parent node. */
		void set_parent1(Node* parent1){this->parent1_ = parent1;}
		
		size_t e_num() const {return this->e_num_;}; /*!< \brief numbering the branch */
		void set_e_num(size_t e_num){this->e_num_ = e_num;}
		
		int ranking();
		int rank() const {return this->rank_;}; /*!< \brief rank of the node, tip node has rank one, the root has the highest rank */
		void set_rank(int rank){this->rank_ = rank;}
		 
		double height() const {return this->height_;}
		void set_height(double height){this->height_ = height;}
		
		double brchlen1() const {return this->brchlen1_;}; /*!< \brief Branch length */
		void set_brchlen1(double bl){this->brchlen1_ = bl;}
		
		bool visited() const {return this->visited_;};
		void set_visited (bool yesorno){this->visited_ = yesorno;}
		
		bool descndnt_of_hybrid() const {return this->descndnt_of_hybrid_;}; /*!< \brief Indicator of descendant of hybrid nodes. It's true, if it is a descendant of hybrid nodes; false, otherwise. */
		void set_descndnt_of_hybrid(bool yesorno){this->descndnt_of_hybrid_ = yesorno;}
		
		bool tip() const {return this->tip_;}; /*!< \brief Indicator of tip nodes. It's true, if it is a tip node, otherwise it is false. */
		void set_tip(bool yesorno){this->tip_ = yesorno;}
		
		/* These members apply to only hybrid nodes */
		bool hybrid() const {return this->hybrid_;}; /*!< \brief Hybrid node only, indicator of a hybrid node */
		void set_hybrid(bool yesorno){this->hybrid_ = yesorno;};
		
		Node* parent2() const {return this->parent2_;} /*!< \brief pointer to its parent node. */
		void set_parent2(Node* parent2){this->parent2_ = parent2;}

		size_t e_num2() const { return this->e_num2_;}; /*!< \brief Hybrid node only, numbering the branch between the node and its second parent */
		void set_e_num2(size_t e_num2){this->e_num2_ = e_num2;}
		
		
		double brchlen2() const {return this->brchlen2_;};/*!< \brief Hybrid node only, Branch length to the second parent*/
		void set_brchlen2(double bl){this->brchlen2_ = bl;};
		//double prob_to_hybrid_left; /*!< \brief Hybrid node only, the probability that a lineage goes to the left */
		
		size_t nodeName_start() const {return this->nodeName_start_;};
		void set_nodeName_start(size_t index){this->nodeName_start_ = index;};
		
		size_t nodeName_end() const {return this->nodeName_end_;};
		void set_nodeName_end(size_t index){this->nodeName_end_ = index;};
		
		size_t node_content_start() const {return this->node_content_start_;};
		void set_node_content_start(size_t index){this->node_content_start_ = index;};
		
		size_t node_content_end() const {return this->node_content_end_;};
		void set_node_content_end(size_t index){this->node_content_end_ = index;};
		
		void add_to_parent(Node* parent);
		
		
		void init(
			size_t nodeName_start = 0,
			size_t nodeName_end = 0, 
			size_t node_content_start = 0, 
			size_t node_content_end = 0, 
			double bl1=0.0,
			double bl2=0.0);
			
		Node(); /*!< \brief Initialize Node class*/
		Node(size_t nodeName_start,	size_t nodeName_end, size_t node_content_start, size_t node_content_end);
		Node(size_t nodeName_start,	size_t nodeName_end, size_t node_content_start, size_t node_content_end, double bl1);
		Node(size_t nodeName_start,	size_t nodeName_end, size_t node_content_start, size_t node_content_end, double bl1, double bl2);
		
		~Node();
		
		
		size_t enumerate_internal_branch(size_t e_num_old);  /*!< enumerator which is about to be updated \todo change e_num_old to int* type */ 
		void print_all_child();
		void print_parent();
		
		bool print_net_Node(const char* net_str) const; 
		bool print_tree_Node(const char* tree_str ) const;
		void clear();
	
	private:
		size_t num_descndnt_;
	
		Node* parent1_; /*!< \brief pointer to its parent node. */
		Node* parent2_; /*!< \brief Hybrid node only, pointer to its second parent node. */
		
		size_t e_num_; /*!< \brief numbering the branch */
		size_t e_num2_; /*!< \brief Hybrid node only, numbering the branch between the node and its second parent */

		int rank_; /*!< \brief rank of the node, tip node has rank one, the root has the highest rank */
		
		size_t node_index_; /*!< \brief node index in the array, \todo use this more often!!!*/

		double height_; /*!< \brief distance to the bottom of the tree */

		double brchlen1_; /*!< \brief Branch length */
		double brchlen2_;/*!< \brief Hybrid node only, Branch length to the second parent*/

		bool visited_;
		bool descndnt_of_hybrid_; /*!< \brief Indicator of descendant of hybrid nodes. It's true, if it is a descendant of hybrid nodes; false, otherwise. */
		bool tip_; /*!< \brief Indicator of tip nodes. It's true, if it is a tip node, otherwise it is false. */
		bool hybrid_; /*!< \brief Hybrid node only, indicator of a hybrid node */

		size_t nodeName_start_;
		size_t nodeName_end_;
		
		size_t node_content_start_;
		size_t node_content_end_;
		
		//size_t name_;
		
		//string clade; /*!< \brief clade at this node, \todo this should be modified to a vector <string> */		
		//string label; /*!< \brief String label of a node, each node has unique label */
		//string node_content; /*!< \brief node content, the subtree string at this node */
		//string name; /*!< \brief Name of a node, this is not unique for nodes. e.g. if its label is A_1, name is A */

		
};



void find_hybrid_descndnt(Node *current);
void find_tip(Node *current);


#endif

