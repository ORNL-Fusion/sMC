#include "constants.hpp"
#include "particle.hpp"

Crk& Crk::operator+=(const Crk &K) {
    this->r = this->r + K.r;
    this->p = this->p + K.r;
    this->z = this->z + K.r;
    return *this;
}

Crk operator+(const Crk &K1, const Crk &K2){
    // shallow copy
	Crk K = K1;
	return K+=K2;
}

Crk operator*(const REAL &f, const Crk &K){
	Crk Kout;
	Kout.r = f * K.r;
	Kout.p = f * K.p;
	Kout.z = f * K.z;
	return Kout;
}

Crk operator*(const Crk &K, const REAL &f){
	return (f * K);
}

Crk operator/(const REAL &f, const Crk &K){
	Crk Kout;
	Kout.r = f / K.r;
	Kout.p = f / K.p;
	Kout.z = f / K.z;
	return Kout;
}

Crk operator/(const Crk &K, const REAL &f){
	Crk Kout;
	Kout.r = K.r / f;
	Kout.p = K.p / f;
	Kout.z = K.z / f;
	return Kout;
}
