#ifndef PARTICLE_HPP_
#define PARTICLE_HPP_

class C_GCparticle {

	public:
		float r, p, z, vPer, vPar, weight, mu;
		int status;
};

class C_rkGCparticle {

	public:
		float r, p, z;
		float v_r, v_p, v_z; 
		float vPer, vPar, dvPar_dt;
		float mu;
		float bDotGradB, bmag,
			b_r, b_p, b_z,
			bCurv_r, bCurv_p, bCurv_z, 
			bGrad_r, bGrad_p, bGrad_z, 
		   	unitb_r, unitb_p, unitb_z;
};

#endif
