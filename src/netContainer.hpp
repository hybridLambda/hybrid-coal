#include"net.hpp"

/*! \brief Collection of Networks in a vector container*/
class Network_s{
	public:
	vector <class Net *> Net_vec;
	
	~Network_s(){for(size_t i = 0; i < Net_vec.size();i++){delete Net_vec[i];}}
	//void clear(){Net_vec.clear();};
	
	//Network_s(){
		//vector <class Net *> Net_vec;
	//}
	
};
