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
vector < vector <size_t> > build_coal_hist_mat(
vector < vector < int > > M_matrix,
Net * my_gt_tree, Net * my_sp_net
);

vector < vector <size_t> >recur_coal_hist(
vector < vector <size_t> > coal_hist_mat,
vector < vector < int > > R_matrix,
//Net * my_sp_net,
Net * my_gt_tree
);


double calc_prob_of_hist_para(vector < vector <size_t> > coal_hist,
vector < vector <int> > all_w,
vector < vector <int> > all_d,
vector < vector <int> > num_enter,
vector < vector <int> > num_out,
Net * my_sp_net,
gijoemat * gijoe_matrix
//vector < vector < vector < double > > > gijoe_matrix
);

double calc_prob_of_hist(
vector < vector < int > > S_matrix,
vector < vector < int > > R_matrix,
vector < vector < int > > M_matrix,
vector < vector <size_t> > coal_hist_mat,
vector < vector <size_t> > coal_hist,
Net * my_sp_net,
gijoemat * gijoe_matrix
);
