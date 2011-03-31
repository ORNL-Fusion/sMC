#ifndef PARTICLE_HPP_
#define PARTICLE_HPP_

#include <iostream>
#include "constants.hpp"
#include <iomanip>

class C_GCparticle {

	public:
		REAL r, p, z, vPer, vPar, weight, mu;
		int status;
};

class C_rkGCparticle {

	public:
		REAL r, p, z;
		REAL v_r, v_p, v_z; 
		REAL vPer, vPar, dvPar_dt;
		REAL mu;
		REAL bDotGradB, bmag,
			b_r, b_p, b_z,
			bCurv_r, bCurv_p, bCurv_z, 
			bGrad_r, bGrad_p, bGrad_z, 
		   	unitb_r, unitb_p, unitb_z;

		void print () {

			std::cout << std::endl;
			std::cout << "\t   r: " << r   <<"\t\t    p: "<< p    <<"\t    z: " << z << std::endl;
			std::cout << "\t v_r: " << v_r <<"\t\t  v_p: "<< v_p  <<"\t  v_z: " << v_z << std::endl;
			std::cout << "\t  mu: " << mu  <<"\t vPer: "<< vPer <<"\t vPar: " << vPar << std::endl;

		}
};

class Crk {

	public:
		REAL r, p, z;

    Crk& operator+=(const Crk &K);

    friend Crk operator+(const Crk &K1, const Crk &K2);
    friend Crk operator*(const REAL &f, const Crk &K);
    friend Crk operator*(const Crk &K, const REAL &f);
    friend Crk operator/(const REAL &f, const Crk &K);
    friend Crk operator/(const Crk &K, const REAL &f);
};

#endif
