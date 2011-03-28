#ifndef DATA2D_HPP_
#define DATA2D_HPP_

#include "boost/multi_array.hpp"

class arr2D {

	private:
		typedef boost::multi_array<float,2> boost2D;
		boost2D data_;

	public:
		arr2D() : data_() {};

		void resize ( const unsigned int nR, const unsigned int nC ) {

			data_.resize(boost::extents[nR][nC]);

		}	

		float operator() (unsigned int row, unsigned int col) {
			return data_[row][col];
		}
};

#endif
