#include "constants.hpp"
#include "rk_gc_particle.hpp"
#include <cmath>

__host__ __device__ Crk& Crk::operator=(const Crk &K) {
	r = K.r;
	p = K.p;
	z = K.z;
	vPar = K.vPar;
	return *this;
}

__host__ __device__ Crk& Crk::operator+=(const Crk &K) {
    r += K.r;
    p += K.p;
    z += K.z;
	vPar += K.vPar;
    return *this;
}

__host__ __device__ Crk& Crk::operator-=(const Crk &K) {
    r -= K.r;
    p -= K.p;
    z -= K.z;
	vPar -= K.vPar;
    return *this;
}

__host__ __device__ Crk operator+(const Crk &K1, const Crk &K2) {
    // Calls copy constructor.
	Crk K = K1;
	return K+=K2;
}

__host__ __device__ Crk operator-(const Crk &K1, const Crk &K2) {
    // Calls copy constructor.
	Crk K = K1;
	return K-=K2;
}

__host__ __device__ Crk operator*(const REAL &f, const Crk &K) {
	Crk Kout;
	Kout.r = f * K.r;
	Kout.p = f * K.p;
	Kout.z = f * K.z;
	Kout.vPar = f * K.vPar;
	return Kout;
}

__host__ __device__ Crk operator*(const Crk &K, const REAL &f) {
	return (f * K);
}

__host__ __device__ Crk operator/(const REAL &f, const Crk &K) {
	Crk Kout;
	Kout.r = f / K.r;
	Kout.p = f / K.p;
	Kout.z = f / K.z;
	Kout.vPar = f / K.vPar;
	return Kout;
}

__host__ __device__ Crk operator/(const Crk &K, const REAL &f) {
	Crk Kout;
	Kout.r = K.r / f;
	Kout.p = K.p / f;
	Kout.z = K.z / f;
	Kout.vPar = K.vPar / f;
	return Kout;
}

__host__ __device__ Crk Kabs (const Crk &K) {
	Crk Kout;
	Kout.r = std::abs(K.r);
	Kout.p = std::abs(K.p);
	Kout.z = std::abs(K.z);
	return Kout;
}

__host__ __device__ REAL Kmax (const Crk &K) {
	REAL maxVal = K.r;
	if(K.p>maxVal) maxVal = K.p;
	if(K.z>maxVal) maxVal = K.z;
	return maxVal;
}

/*
void Crk::print () {
	std::cout <<"\t"<<r<<"\t"<<p<<"\t"<<z<<"\t"<<vPar<<std::endl;
}
*/
