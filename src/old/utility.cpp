/* 
 * hybrid_coal is used to compute gene tree probabilities given species network under coalescent process.
 * 
 * hybrid_sim is used to simulate gene trees given species network under 
 * coalescent process.
 * 
 * Copyright (C) 2011 Sha (Joe) Zhu
 * 
 * This file is part of hybrid_sim and hybrid_coal
 * 
 * hybrid_sim is free software: you can redistribute it and/or modify
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

/*! \file utility.cpp
 *  \brief Core function of converting a Newick (extended Newick) format string into a species tree (network), and simple string manipulation for tree strings
 */



#include"utility.hpp"
//extern bool utility_debug_bool;
//extern bool log_file_bool;
bool utility_debug_bool=false;







/*! \brief Identify if its the start of the taxon name in a newick string, should be replaced by using (isalpha() || isdigit())  */
bool start_of_tax_name(string in_str,size_t i){
	bool start_bool=false;
	if ( (in_str[i]!='(' && in_str[i-1]=='(') || (in_str[i-1]==',' && in_str[i]!='(') || ( (in_str[i-1]==')') && ( in_str[i]!=')' || in_str[i]!=':' || in_str[i]!=',' || in_str[i]!=';' ) ) ) {
		start_bool=true;	
	}
	
	return 	start_bool;
}


size_t Parenthesis_balance_index_backwards(string in_str,size_t i){
	size_t j=i;
	int num_b=0;
	for (;j>0;j--){
		if (in_str[j]=='('){
			num_b--;
		}
		if (in_str[j]==')'){
			num_b++;
		}
		if (num_b==0){
			break;
		}
	}
	return j;
}

size_t Parenthesis_balance_index_forwards(string in_str,size_t i){
	size_t j=i;
	int num_b=0;
	for (;j<in_str.size();j++){
		if (in_str[j]=='('){
			num_b++;
		}
		if (in_str[j]==')'){
			num_b--;
		}
		if (num_b==0){
			break;
		}
	}
	return j;
}



string checking_labeled(string in_str){
	bool labeled_bool=true;
	string out_str;
	for (size_t i=0;i<in_str.size();i++){
		if (in_str[i]==')' && i==end_of_label_or_bl(in_str, i) ){
			labeled_bool=false;
			//break;
		}
	}
	
	if (labeled_bool){
		out_str=in_str;
	}
	else{
		out_str=label_interior_node(in_str);
	}	
	return out_str;
}







//bool find_descndnt2(Node* current, string taxname){
	//bool descndnt_found=false;
	//if (current->tip_bool){
		////if (current->label==taxname){
		//if (current->label==taxname){
			//descndnt_found=true;
		//}
		//else{
			//descndnt_found=false;
		//}
	//}
	//else {//int i;
		//for (int i=0;i<current->num_child;i++){
			//if (find_descndnt2(current->child[i],taxname)){
				//descndnt_found=true;
				//break;
			//}
			//else descndnt_found=false;
		//}
	//}	
	//return descndnt_found;
//}





/*! \brief Label interior node if the interior nodes of the tree string are not labeled */
string label_interior_node(string in_str /*!< input newick form string */){
	vector <string> in_str_partition;
	int interior_node_counter=0;
	int sub_str_start_index=0;			
	size_t i=in_str.find(')');
	while ( i<in_str.size() ){
		interior_node_counter++;
		string current_string;
		size_t found_next_bracket=min(in_str.find(")",sub_str_start_index),in_str.size());
		current_string=in_str.substr(sub_str_start_index,found_next_bracket - sub_str_start_index +1);
		if (in_str[i+1]==';' || i==(in_str.size()-1)){
			current_string=current_string+"root;";
			in_str_partition.push_back(current_string);
		}
		else{
			ostringstream interior_node_counter_str;
			interior_node_counter_str<<interior_node_counter;
			current_string=current_string+"Int";
			current_string=current_string+interior_node_counter_str.str();
			in_str_partition.push_back(current_string);
			sub_str_start_index=i+1;
		}
		i=in_str.find(")",i+1);
	}
	string out_str;
	for (size_t i=0;i<in_str_partition.size();i++){
		//cout<<in_str_partition[i]<<endl;
		out_str=out_str+in_str_partition[i];
	}
	return out_str;
}

void appending_debug_file(string debug_file_input){
	ofstream debug_file;
	debug_file.open ("debug_file", ios::out | ios::app | ios::binary); 
	debug_file << debug_file_input << "\n";
	debug_file.close();
}

/*! \brief Add more information to log_file */
void appending_log_file(string log_file_input /*! Information added*/){
	ofstream log_file;
	log_file.open ("log_file", ios::out | ios::app | ios::binary); 
	log_file << log_file_input << "\n";
	log_file.close();
}


///*! \brief rewrite node content of nodes */
//void rewrite_node_content(vector <Node*> Net_ptr /*! vector of pointers pointing to nodes */){
	//int highest_i=0;
	//for (unsigned int i =0; i<Net_ptr.size();i++){
		//if (Net_ptr[i]->num_descndnt > Net_ptr[highest_i]->num_descndnt ){highest_i=i;}
	//}
	
	////for (unsigned int node_i=0;node_i<Net_ptr.size();node_i++){
		////Net_ptr[node_i]->print_net_Node();
		////cout<<endl;
	////}
	
	//ranking(Net_ptr[highest_i]);
	////cout<<Net_ptr[highest_i]->node_content<<endl;
	////for (unsigned int node_i=0;node_i<Net_ptr.size();node_i++){
		////Net_ptr[node_i]->print_net_Node();
		////cout<<endl;
	////}
	//for (int rank_i=1;rank_i<=Net_ptr.back()->rank;rank_i++){
		//for (unsigned int node_i=0;node_i<Net_ptr.size();node_i++){
			//if (Net_ptr[node_i]->rank==1){
				//Net_ptr[node_i]->node_content=Net_ptr[node_i]->label;
			//}
			//else{
			//if (Net_ptr[node_i]->rank==rank_i){
				//string new_node_content="(";
				//for (int child_i=0;child_i<Net_ptr[node_i]->num_child;child_i++){
					//ostringstream brchlen_str;
					//ostringstream brchlen_str2;
					//brchlen_str<<Net_ptr[node_i]->child[child_i]->brchlen1;
					//if (Net_ptr[node_i]->child[child_i]->node_content==Net_ptr[node_i]->child[child_i]->label){
						//new_node_content=new_node_content+Net_ptr[node_i]->child[child_i]->label+":"+brchlen_str.str();}
					//else{
						//bool new_hybrid_node=false;
						//for (unsigned int node_ii=0;node_ii<node_i;node_ii++){
							//for (int node_ii_child_i=0;node_ii_child_i<Net_ptr[node_ii]->num_child;node_ii_child_i++){
								//if (Net_ptr[node_ii]->child[node_ii_child_i]->node_content==Net_ptr[node_i]->child[child_i]->node_content){
									//new_hybrid_node=true;
									//brchlen_str2<<Net_ptr[node_i]->child[child_i]->brchlen2;
								//break;}
							//}
							////if (new_hybrid_node==1){break;}
						//}
						//if (new_hybrid_node){
							//new_node_content=new_node_content+Net_ptr[node_i]->child[child_i]->label+":"+brchlen_str2.str();
						//}
						//else{
							//new_node_content=new_node_content+Net_ptr[node_i]->child[child_i]->node_content+Net_ptr[node_i]->child[child_i]->label+":"+brchlen_str.str();
						//}
					//}
					////new_node_content=new_node_content+sp_nodes_ptr_rm1[node_i]->child[child_i]->node_content+sp_nodes_ptr_rm1[node_i]->child[child_i]->label+":"+brchlen_str.str();
					//if (child_i<Net_ptr[node_i]->num_child-1){
						//new_node_content=new_node_content+",";
					//}
				//}
				//new_node_content=new_node_content+")";
				//Net_ptr[node_i]->node_content=new_node_content;
				//}
			//}
		//}	
	//}
	////for (unsigned int i =0; i<Net_ptr.size();i++){
		////cout<<Net_ptr[i]->label<<" "<<Net_ptr[i]->node_content<<endl;//<<"!!!!";
	////}
//}

/*! \brief Remove interior nodes label of a string */
string remove_interior_label(string in_str/*!< input newick form string */){
	string out_str;
	out_str=in_str;
	size_t found_bracket=out_str.find(')');
	while ( found_bracket<out_str.size() ){
		if (isalpha(out_str[found_bracket+1]) || isdigit(out_str[found_bracket+1])){
			size_t char_j=end_of_label_or_bl(out_str, found_bracket+1);
			out_str.erase(out_str.begin()+found_bracket+1,out_str.begin()+char_j+1);
		}
		found_bracket=out_str.find(")",found_bracket+1);
		
	}
	return out_str;
}


/*! \brief Remove branch length in a tree string, gives the tree topology */
string remove_brchlen(string in_str /*!< input newick form string */){
	string out_str = in_str;
	size_t found_col=out_str.find(':');
	while ( found_col<out_str.size() ){
		size_t char_j=end_of_label_or_bl(out_str,found_col)+1;
		//out_str.erase(out_str.begin()+found_col,out_str.begin()+char_j);
		out_str.erase(found_col,char_j-found_col);
		found_col=out_str.find(":",found_col+1);
	}
	return out_str;
}

/*! \brief Remove branch length in a tree string, gives the tree topology \todo combine with remove_brchlen, only keep one of them*/
string tree_topo(string in_str /*!< input newick form string */){
	return remove_brchlen(remove_interior_label(in_str));
}






/*! \brief Compute factorial of a \return double a! */
double factorial(double factnum)
{
    double a=0;
    for(double i = 0;i<=factnum;i++)
    {
        if(i==0)
        {
            a = 1;
        }
        else
        {
            a=(a*i);
        }
    }

    return (a);

}

//double factorial (double a){
	//if (a > 1){
		//return (a * factorial (a-1));}
	//else{
		//return (1);}
//}





/*! \brief Compute a permutations of n \return double */
double n_permu_a (double n, double a){
	if (a>1){
		return (n*n_permu_a(n-1,a-1));
	}
	else{
		if (a==1){
			return (n);
		}
		else{
			return (1);
		}
	}
}

double n_choose_k(double n, double k){
		if (n>=k){
		double n_minus_k;//=n-k;
		if (n/2>=k){
			n_minus_k=n-k;
		}
		else{
			n_minus_k=k;
			k=n-k;
		}
		double above=1;
		double below=1;
		for (double i=n;i>n_minus_k;i--){
			above=above*i;
		}
		for(double i = 1;i<=k;i++){
	        below=(below*i);
	    }	
		return (above/below);
	}
	else{ return 0;}
}


/*! \brief Compute n choose k \return double */
double n_choose_k_old(double n, double k){
	if (k<(n/2)){
		return (n_choose_k_old(n,n-k));}
	else{
		return (n_permu_a(n,k)/factorial(k));}
}



/*! \brief Compute factorial of a \return int a! */
int factorial_int (int a){
	if (a > 1){
		return (a * factorial_int (a-1));}
	else{
		return (1);}
}

/*! \brief Compute a permutations of n \return int */
int n_permu_a_int (int n, int a){
	if (a>1){
		return (n*n_permu_a_int(n-1,a-1));}
	else{
		if (a==1){
			return (n);}
		else{
			return (1);}
	}
}

/*! \brief Compute n choose k \return int */
int n_choose_k_int(int n, int k){
	if (k<(n/2)){
		return (n_choose_k_int(n,n-k));}
	else{
		return (n_permu_a_int(n,k)/factorial_int(k));}
}

/*! \brief Remove the '&' and '#' signs from a string \return string */
string rm_and_hash_sign(string in_str){
	return rm_hash_sign(rm_and_sign(in_str));
}


/*! \brief Remove '#' signs and the gamma parameter from a string \return string */
string rm_hash_sign(string in_str){
	while (int(in_str.find('#'))>0 && in_str.find('#')!=string::npos) {
		size_t i=end_of_label_or_bl(in_str, in_str.find('#'));
		in_str.erase(in_str.find('#'),i-in_str.find('#')+1);
		//in_str.erase(in_str.find('#'),1);
	}
	return in_str;
}

string rm_and_sign(string in_str){
	while (int(in_str.find('&'))>0 && in_str.find('&')!=string::npos) {
		in_str.erase(in_str.find('&'),1);
	}
	return in_str;
}

string rm_dash_sign(string in_str){
	while (int(in_str.find('_'))>0 && in_str.find('_')!=string::npos) {
		in_str.erase(in_str.find('_'),1);
	}
	return in_str;
}

/*! \brief Terminate the program and print out the log file*/
int my_exit(){ // change my_exit() to return 0, and terminate program
	int sys=system("cat log_file");		
	//exit(1);
	return 0;
}




/*! Check and remove files*/
void check_and_remove(const char* file_name){
	ifstream my_file(file_name);
	if (my_file.good())
	{
	  remove(file_name);
	}
}

size_t end_of_label_or_bl(string in_str, size_t i){
	size_t j ;
	for (j=i;j<in_str.size();j++){
		char stop=in_str[j+1];
		if (stop==',' || stop==')' || stop==':' || stop==';'){
			break;
		}
	}
	return j;
}

string extract_label(string in_str, size_t i){
	size_t j=end_of_label_or_bl(in_str, i);
	//cout<<"i="<<i<<", j="<<j<<endl;
	return in_str.substr(i,j+1-i);
}


size_t hybrid_hash_index(string in_str){
	return in_str.find('#');
}

string extract_hybrid_label(string in_str){
	size_t hash_index=hybrid_hash_index(in_str);
	return in_str.substr(0,hash_index);
}

string extract_hybrid_para_str(string in_str){
	size_t hash_index=hybrid_hash_index(in_str);
	return in_str.substr(hash_index+1);//,in_str.size()-1);
}

double extract_hybrid_para(string in_str){
	double para;
	istringstream para_istr(extract_hybrid_para_str(in_str));
	para_istr>>para;
	return para;
}


string read_input_line(char inchar[]){
	ifstream in_file;
	string out_str;
	in_file.open(inchar);
	if (in_file.good()){
		getline (in_file,out_str);}
	else{
		string dummy_str(inchar);
		if (dummy_str.find('(')!=string::npos && dummy_str.find(')')!=string::npos){
		out_str=dummy_str;
		}else{
			cout<<"Error: check input '"<<inchar<<"'"<<endl;
			my_exit();
		}
	}
	in_file.close();
	
	//check_sp_str_dash_sign(out_str);
			
	return 	out_str;
}

vector <string> read_input_lines(char inchar[]){
	vector <string> out_vec;
	ifstream in_file;
	in_file.open(inchar);
	string out_str;
	if (in_file.good()){
		getline (in_file,out_str);
		while (out_str.size()>0){   
			check_gt_str_and_sign(out_str);
			out_vec.push_back(out_str);
			getline (in_file,out_str);
		}
	}	
	else{
		string dummy_str(inchar);
		if (dummy_str.find('(')!=string::npos && dummy_str.find(')')!=string::npos){
			out_str=dummy_str;
			check_gt_str_and_sign(out_str);
			out_vec.push_back(out_str);
		}else{
			cout<<"Error: check input '"<<inchar<<"'"<<endl;
			my_exit();
		}
	}
	in_file.close();
	return out_vec;	
}

bool is_num(char inchar[]){
	bool is_num_return=true;
	string in_str(inchar);
	for (size_t i=0;i<in_str.size();i++){
		if (isalpha(in_str[i]) && in_str[i]!='e'){
			is_num_return=false;
			break;
		}
	}
	return is_num_return;
}

//string read_input_para(char inchar[],string in_str){
	
	//string out_str;
	//if (is_num(inchar)){
		//istringstream para_istrm(inchar);
		//double para;
		//para_istrm>>para;
		//out_str=write_para_into_tree(in_str, para);
	//}
	//else{
		//out_str=read_input_line(inchar);
	//}
	//return out_str;
//}
		
		

		
		
string rewrite_net_str(string net_str, vector <string> tax_name, vector <string> tip_name){
	//vector < vector <string> > new_names;
	net_str=checking_labeled(net_str);
	for (size_t i=0;i<tax_name.size();i++){
		vector <string> new_names;
		for (size_t j=0;j<tip_name.size();){
			size_t found_tip=tip_name[j].find(tax_name[i]);
			if (found_tip!=string::npos){
			//if (){
				new_names.push_back(tip_name[j]);
				tip_name.erase(tip_name.begin()+j);
			}
			else{
				j++;
			}
		}
		if (new_names.size()>1){
			size_t label_start=net_str.find(tax_name[i]);
			size_t label_end=end_of_label_or_bl(net_str, label_start);
			string new_sub="(";
			for (size_t k=0;k<new_names.size();k++){
				new_sub=new_sub+new_names[k]+":0,";
				//if (k<new_names.size()-1){
					
				//}
				
			}
			string int_label=")"+tax_name[i];
			new_sub.replace(new_sub.size()-1,1,int_label);
			net_str.replace(label_start,label_end-label_start+1,new_sub);
		}
		
		//new_names.clear();
	}
	
	
	
	return rm_dash_sign(net_str);
}



void check_gt_str_and_sign(string gt_str){
	size_t found=gt_str.find('&');
	if (found!=string::npos){
		string gt_log="Error: gene tree "+gt_str+" can not be a gene tree string";
		appending_log_file(gt_log);
		my_exit();
		exit(1);
	}
	
}

//void check_sp_str_dash_sign(string sp_str){
	//Net sp_net(sp_str);
	//for (size_t i=0;i<sp_net.Net_nodes.size();i++){
		//if (sp_net.Net_nodes[i].tip_bool){
			//size_t found=sp_net.Net_nodes[i].label.find('_');
			//if (found!=string::npos){
				//string sp_log="Error:  "+sp_str+" can not be a species string";
				//appending_log_file(sp_log);
				//my_exit();
				//exit(1);				
			//}
		//}
	//}

	
//}
