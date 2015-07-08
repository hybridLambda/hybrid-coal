/* 
 * hybrid-coal is used to compute gene tree probabilities given species network under coalescent process.
 * 
 * Copyright (C) 2010 -- 2015 Sha (Joe) Zhu
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

#include "graph.hpp"

#ifndef NDEBUG
#define CoalSTdout (std::cout << "      CoalST ")
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define CoalSTdout 0 && (std::cout << "      CoalST ")
#endif

#ifndef COAL
#define COAL

enum SYMBOLIC_MODE { NONE, MAPLE, SYMB_LATEX };

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


template < class T > bool print_2D_matrix( vector < vector < T > > & mat ){
    for ( size_t i = 0; i < mat.size(); i++ ){
        CoalSTdout << " " ;
        for ( size_t j = 0; j < mat[i].size(); j++){
            dout << mat[i][j];
        }
        dout << endl;
    }
    dout<<endl;
    return true;
}

class CoalST;
class CoalGT: public GraphBuilder {
    friend class HybridCoal;
    friend class TestCoal;
    // Members
    vector < vector < int > >    R_matrix;
    vector < vector < int > >    M_matrix;
    vector < vector < size_t > > coal_hist_mat;
    vector < vector < size_t > > valid_coal_hist;
    vector < vector < double > > all_w;
    vector < vector < double > > all_d;
    vector < vector < int > >    num_enter;
    vector < vector < int > >    num_out;    
    double probability;
    
    // Methods
    void prob_given_sp_tree ( CoalST & sp_tree );
    void initialize_possible_coal_hist( CoalST & sp_tree );
    void building_R_matrix();
    void building_M_matrix( CoalST & sp_tree ) ;
    void sum_coalescent_history_prob( CoalST & sp_tree );
    void enumerate_coal_events( CoalST & sp_tree );
    vector < vector < size_t > > recur_coal_hist ( vector < vector <size_t > > coal_hist, size_t node_i);
    void build_coal_hist( );
    
  public:
    CoalGT ( string gt_str );
    ~CoalGT(){};
};


class CoalST: public GraphBuilder {
  friend class HybridCoal;
  friend class CoalGT;
  friend class TestCoal;
  friend class CoalSN;
    // Members
    vector < double > brchlens_vec;
    vector < int > max_num_brch_vec;
    vector < vector < int > > S_matrix;
    vector < vector < vector < double > > > gijoemat;
    
    // Methods
    void assign_bl_to_vec();
    void build_gijoe();
    bool print_gijoemat();
    void building_S_matrix();

    CoalST ( string sp_str );
    CoalST ( const CoalST & spIn );
    ~CoalST(){}
};


class NetStrWizPrior {
  friend class CoalSN;
  friend class TestRm;
  friend class HybridCoal;
    string netStr;
    //string netStrLabelled;
    vector < valarray < int > > priorCladeList;
    vector < int > priorCoalList;
    //vector < vector < int > > lambdaSum;
    string omegaStr(){
        return to_string(this->omega_);
    }
    void setOmega ( double omega ){
        this->omega_ = omega;
    }
    double omega() const {
        return this->omega_;
    }
    double omega_;
  public:
    NetStrWizPrior( string firstStr ){
        this->netStr = firstStr;
        this->setOmega(1.0);
    }

    ~NetStrWizPrior(){}
};


class TmpSN : public GraphBuilder {
  friend class CoalSN;
  friend class TestRm;
    int toBeRemovedNodeIndex_;
    int toBeRemovedNodeIndex() const {return this->toBeRemovedNodeIndex_;}

    TmpSN ( string tmpStr );
    ~TmpSN(){}

    void chooseRemoveNode();

};


class CoalSN : public CoalST {
  friend class CoalST;
  friend class TestRm;
  friend class HybridCoal;
    bool maple_bool_local; // TODO, need to rework on this
    vector < NetStrWizPrior > NetStrWizPriorList;
    size_t currentSubNetworkIndex;
    
    void initializeNetStrWizPriorList ( string spStr );
    void simplifyNetworks( string gt_string );
    void removeSnode( size_t rmNodeIndex, string gtStr );
    
    double left_hybrid_parameter_num;
    double right_hybrid_parameter_num;
    string new_node_name;
    size_t hashSingIdx;
    void removeHnode ( TmpSN &tmpSN, size_t rmNodeIndex, NetStrWizPrior netStrWizPrior, SYMBOLIC_MODE sybolicMode = NONE );
       void removeHnodeOneChild( TmpSN &tmpSN, size_t rmNodeIndex, NetStrWizPrior netStrWizPrior, SYMBOLIC_MODE sybolicMode );
         string removeHnodeOneChildCore( TmpSN &tmpSN, size_t rmNodeIndex, size_t removingFromParentIndex);
       void removeHnodeManyChild( TmpSN &tmpSN, size_t rmNodeIndex, NetStrWizPrior netStrWizPrior, SYMBOLIC_MODE sybolicMode );
         string removeHnodeManyChildCore(TmpSN &tmpSN, size_t rmNodeIndex, NetStrWizPrior netStrWizPrior, int numberOfChildAtRemovingNode, vector < valarray <int> > &A_matrix_ith_row, size_t A_matrix_i );
         vector <int> build_left_or_right_vec(int n_child,vector < valarray <int> > A_matrix_ith_row,size_t A_matrix_i);

    CoalSN ( string sp_str );
    ~CoalSN(){}
};

vector < valarray <int> > all_possible_comb(int n);
vector < vector < valarray < int > > > build_h_child(int n);
vector < vector < valarray < int > > > build_s_child(int n);
vector < valarray <int> > rearrange_A(vector < valarray <int> > A,int n);
size_t hybrid_hash_index(string &in_str);
string extract_hybrid_para_str(string &in_str);
double extract_hybrid_para(string in_str);

#endif
