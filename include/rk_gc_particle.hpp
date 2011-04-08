#ifndef RK_PARTICLE_HPP_
#define RK_PARTICLE_HPP_

#include "constants.hpp"

class Crk {

	private:

	public:
		REAL r, p, z, vPar;
		
		// Default constructor
		__host__ __device__
			Crk () {r=0.0;p=0.0;z=0.0;vPar=0.0;}; 

		// Copy constructor
		__host__ __device__
			Crk ( const Crk &K ) {*this = K;}	

		__host__ __device__
    		Crk& operator+=(const Crk &K);
		__host__ __device__
	   		Crk& operator-=(const Crk &K);
		__host__ __device__
			Crk& operator=(const Crk &K);

 		__host__ __device__
   			friend Crk operator+(const Crk &K1, const Crk &K2);
    	__host__ __device__
			friend Crk operator-(const Crk &K1, const Crk &K2);
    	__host__ __device__
			friend Crk operator*(const REAL &f, const Crk &K);
    	__host__ __device__
			friend Crk operator*(const Crk &K, const REAL &f);
    	__host__ __device__
			friend Crk operator/(const REAL &f, const Crk &K);
    	__host__ __device__
			friend Crk operator/(const Crk &K, const REAL &f);

	   	__host__ __device__
			friend Crk Kabs(const Crk &K);
	   	__host__ __device__
			friend REAL Kmax(const Crk &K);
		
		//void print();
};

#endif
