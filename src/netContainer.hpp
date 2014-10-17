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
