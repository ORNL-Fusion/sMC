#ifndef CUDA_CLASSES_
#define CUDA_CLASSES_

#include "constants.hpp"

class Ccu_ptr_pitch {

    public:
        REAL *ptr;
        size_t pitch;
    
        Ccu_ptr_pitch(): ptr(NULL), pitch(0) {};
};

class Ccu_ptrs {

    public:
        REAL *r, *z;
        Ccu_ptr_pitch 
            bmag, b_r, b_p, b_z, 
            bCurv_r, bCurv_p, bCurv_z,
            bGrad_r, bGrad_p, bGrad_z,
            bDotGradB;

        Ccu_ptrs(): r(NULL), z(NULL) {};

//#ifdef __CUDACC__
//        __host__ __device__
//#endif
//        bool isValid() const {
//       
//            bool err = true; 
//
//            if(!bmag.ptr || bmag.pitch<=1) 
//                err = false;
//
//            if(!b_r.ptr || b_r.pitch<=1) 
//                err = false;
//            if(!b_p.ptr || b_p.pitch<=1) 
//                err = false;
//            if(!b_z.ptr || b_z.pitch<=1) 
//                err = false;
//
//            if(!bCurv_r.ptr || bCurv_r.pitch<=1) 
//                err = false;
//            if(!bCurv_p.ptr || bCurv_p.pitch<=1) 
//                err = false;
//            if(!bCurv_z.ptr || bCurv_z.pitch<=1) 
//                err = false;
//
//            if(!bGrad_r.ptr || bGrad_r.pitch<=1) 
//                err = false;
//            if(!bGrad_p.ptr || bGrad_p.pitch<=1) 
//                err = false;
//            if(!bGrad_z.ptr || bGrad_z.pitch<=1) 
//                err = false;
//
//            if(!bDotGradB.ptr || bDotGradB.pitch<=1) 
//                err = false;
//
//            return err;
//
//        }
};

#endif
