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
