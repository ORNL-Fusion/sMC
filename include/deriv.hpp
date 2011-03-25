#ifndef DERIV_HPP_
#define DERIV_HPP_

#include <vector>

template <typename T>
std::vector<T> deriv ( const std::vector<T> &data, const T dx, const int order ) {

	int n = data.size ();

	std::vector<T> result (n);

	if(order==2) {

		for (int i=1;i<n-1;i++)
		{
    		result[i] = 1.0 / ( 2.0 * dx ) *
				( data[i+1] - data[i-1] );

		}

		result[0] = 1.0 / ( 2.0 * dx ) * 
			( -3.0 * data[0] + 4.0 * data[1] - data[2] );

		result[n-1] = 1.0 / ( 2.0 * dx ) * 
			( 3.0 * data[n-1] - 4.0 * data[n-2] + data[n-3] );

	}

	if(order==4) {

		for(int i=2;i<n-2;i++)
		{
    		result[i] = 1.0 / dx *
				( data[i-2]*1.0/12.0 
				  - data[i-1]*2.0/3.0 
				  + data[i+1]*2.0/3.0
				  - data[i+2]*1.0/12.0 );

		}

		for(int i=0;i<2;i++)
			result[i] = 1.0 / dx * 
				( -data[i]*25.0/12.0
				  +data[i+1]*4.0
				  -data[i+2]*3.0
				  +data[i+3]*4.0/3.0
				  -data[i+4]*1.0/4.0 );

		for(int i=n-2;i<n;i++)
			result[i] = -1.0 / dx * 
				( -data[i]*25.0/12.0
				  +data[i-1]*4.0
				  -data[i-2]*3.0
				  +data[i-3]*4.0/3.0
				  -data[i-4]*1.0/4.0 );

	}


	return result;
}

#endif

