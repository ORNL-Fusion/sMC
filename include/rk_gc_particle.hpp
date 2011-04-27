#ifndef RK_PARTICLE_HPP_
#define RK_PARTICLE_HPP_

#include "constants.hpp"
#include <cmath>

namespace {

class Crk {

	private:

	public:
		REAL r, p, z, vPar;
		
		// Default constructor
#ifdef __CUDACC__
		__host__ __device__
#endif
			Crk () {r=0.0;p=0.0;z=0.0;vPar=0.0;}; 

#ifdef __CUDACC__
		__host__ __device__
#endif
			Crk (const REAL r_, const REAL p_, const REAL z_, const REAL vPar_) {
					r=r_;
					p=p_;
					z=z_;
					vPar=vPar_;
			}; 

		// Copy constructor
#ifdef __CUDACC__
		__host__ __device__
#endif
			Crk ( const Crk &K ) {*this = K;}	

#ifdef __CUDACC__
		__host__ __device__
#endif
    		Crk& operator+=(const Crk &K);
#ifdef __CUDACC__
		__host__ __device__
#endif
	   		Crk& operator-=(const Crk &K);
#ifdef __CUDACC__
		__host__ __device__
#endif
			Crk& operator=(const Crk &K);

#ifdef __CUDACC__
 		__host__ __device__
#endif
   			friend Crk operator+(const Crk &K1, const Crk &K2);
#ifdef __CUDACC__
    	__host__ __device__
#endif
			friend Crk operator-(const Crk &K1, const Crk &K2);
#ifdef __CUDACC__
    	__host__ __device__
#endif
			friend Crk operator*(const REAL &f, const Crk &K);
#ifdef __CUDACC__
    	__host__ __device__
#endif
			friend Crk operator*(const Crk &K, const REAL &f);
#ifdef __CUDACC__
    	__host__ __device__
#endif
			friend Crk operator/(const REAL &f, const Crk &K);
#ifdef __CUDACC__
    	__host__ __device__
#endif
			friend Crk operator/(const Crk &K, const REAL &f);

#ifdef __CUDACC__
	   	__host__ __device__
#endif
			friend Crk Kabs(const Crk &K);
#ifdef __CUDACC__
	   	__host__ __device__
#endif
			friend REAL Kmax(const Crk &K);
#ifdef __CUDACC__
	   	__host__ __device__
#endif
			bool isValid();
};


#ifdef __CUDACC__
__host__ __device__ 
#endif
Crk& Crk::operator=(const Crk &K) {
	r = K.r;
	p = K.p;
	z = K.z;
	vPar = K.vPar;
	return *this;
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
Crk& Crk::operator+=(const Crk &K) {
    r += K.r;
    p += K.p;
    z += K.z;
	vPar += K.vPar;
    return *this;
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
Crk& Crk::operator-=(const Crk &K) {
    r -= K.r;
    p -= K.p;
    z -= K.z;
	vPar -= K.vPar;
    return *this;
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
Crk operator+(const Crk &K1, const Crk &K2) {
    // Calls copy constructor.
	Crk K = K1;
	return K+=K2;
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
Crk operator-(const Crk &K1, const Crk &K2) {
    // Calls copy constructor.
	Crk K = K1;
	return K-=K2;
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
Crk operator*(const REAL &f, const Crk &K) {
	Crk Kout;
	Kout.r = f * K.r;
	Kout.p = f * K.p;
	Kout.z = f * K.z;
	Kout.vPar = f * K.vPar;
	return Kout;
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
Crk operator*(const Crk &K, const REAL &f) {
	return (f * K);
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
Crk operator/(const REAL &f, const Crk &K) {
	Crk Kout;
	Kout.r = f / K.r;
	Kout.p = f / K.p;
	Kout.z = f / K.z;
	Kout.vPar = f / K.vPar;
	return Kout;
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
Crk operator/(const Crk &K, const REAL &f) {
	Crk Kout;
	Kout.r = K.r / f;
	Kout.p = K.p / f;
	Kout.z = K.z / f;
	Kout.vPar = K.vPar / f;
	return Kout;
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
Crk Kabs (const Crk &K) {
	Crk Kout;
	Kout.r = std::abs(K.r);
	Kout.p = std::abs(K.p);
	Kout.z = std::abs(K.z);
	return Kout;
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
REAL Kmax (const Crk &K) {
	REAL maxVal = K.r;
	if(K.p>maxVal) maxVal = K.p;
	if(K.z>maxVal) maxVal = K.z;
	return maxVal;
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
bool Crk::isValid () {
	bool good = true;
	if(r!=r) good = false;
	if(p!=p) good = false;
	if(z!=z) good = false;
	if(vPar!=vPar) good = false;
	return good;
}

} // end namespace

#endif
