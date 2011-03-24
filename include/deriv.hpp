#ifndef DERIV_HPP_
#define DERIV_HPP_

#include <vector>

template <typename T>
std::vector<T> deriv ( const std::vector<T> &data, const T dx ) {

	int n = data.size ();

	std::vector<T> result (n);

	for (int i=1;i<n-1;i++)
	{
    	result[i] = 1.0 / ( 2.0 * dx ) *
			( data[i+1] - data[i-1] );

	}

	result[0] = 1.0 / ( 2.0 * dx ) * 
		( -3.0 * data[0] + 4.0 * data[1] - data[2] );

	result[n-1] = 1.0 / ( 2.0 * dx ) * 
		( 3.0 * data[n-1] - 4.0 * data[n-2] + data[n-3] );


	return result;
}

#endif

