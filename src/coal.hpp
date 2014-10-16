/* 
 * hybrid_coal is used to compute gene tree probabilities given species network under coalescent process.
 * 
 * Copyright (C) 2010 -- 2014 Sha (Joe) Zhu
 * 
 * This file is part of hybrid_coal
 * 
 * hybrid_coal is free software: you can redistribute it and/or modify
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

    /*! \brief Compute factorial of a \return double a! */
    template < class T > T factorial ( T a ){
        if (a > 1) return (a * factorial (a-1));
        else       return (1);
    }
    
    /*! \brief Compute a permutations of n \return double */
    template < class T > T n_permu_a ( T n, T a ){
        if   ( a > 1 ) return (n*n_permu_a(n-1,a-1));
        else if (a==1) return (n);
        else           return (1);
    }
    
    /*! \brief Compute n choose k \return double */
    template < class T > T n_choose_k ( T n, T k ){
        if ( k < ( n/2 ) ) return (n_choose_k(n,n-k));
        else               return (n_permu_a(n,k)/factorial(k));
    }
