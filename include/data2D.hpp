#ifndef DATA2D_HPP_
#define DATA2D_HPP_

#include "boost/multi_array.hpp"

class arr2D {

	private:
		boost::multi_array<float,2> data;

	public:
		arr2D() : data() {};

		void resize ( const unsigned int nR, const unsigned int nC ) {

			//boost::multi_array<int, 3> extents; 
			//data.resize(extents[nR][nC]);

		}	
};

#endif
