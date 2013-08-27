//coal.hpp
#include<vector>
#include"net.hpp"
class gijoemat{
	public:
	gijoemat();
	gijoemat(Net* sp_net);
	~gijoemat(){};
	vector < vector < vector < double > > > mat;
};

void print_matrix(vector < vector < int > > mat);
vector < vector < int > > building_S_matrix(Net* my_sp_net);
vector < vector < int > > building_M_matrix(Net* my_gt_tree, Net* my_sp_net);
vector < vector < int > > building_R_matrix(Net* gt_tree);
vector < vector <unsigned int> > build_coal_hist_mat(
vector < vector < int > > M_matrix,
Net * my_gt_tree, Net * my_sp_net
);

vector < vector <unsigned int> >recur_coal_hist(
vector < vector <unsigned int> > coal_hist_mat,
vector < vector < int > > R_matrix,
Net * my_sp_net,
Net * my_gt_tree
);
