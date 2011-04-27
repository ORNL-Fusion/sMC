//Mostly stolen form Ben Bales' Data.hpp ;)

#ifndef ARRAY2D_HPP_
#define ARRAY2D_HPP_

#include "util.hpp"
#include <iostream>
#include <cstring>

#define BCHECK true

template <class T, bool boundscheck = false>
class array2D {

    private:

    public:
 
        T *ptr;
        unsigned int M, N;
        size_t pitchBytes;

        // Constructors
#ifdef __CUDACC__
        __host__ __device__
#endif
        array2D(): M(0), N(0), ptr(NULL), pitchBytes(0) {};
        array2D(const unsigned int _M, const unsigned int _N) : ptr(NULL), pitchBytes(0) {
           resize(_M,_N); 
        }
        // copy constructor
        array2D(const array2D &source) {
			pitchBytes = 0;
			ptr = NULL;
            *this = source;
        }

        // conversion constructor 
        template<class T2>
        array2D(const array2D<T2,boundscheck>&source) {
            pitchBytes = 0;
            ptr = NULL;
            resize(source.M,source.N);
            for(size_t n=0;n<N;++n) {
                for(size_t m=0;m<M;++m) {
                    ptr[ n + N * m ] = (T)source(m,n);
                }
            }
        }

        // assignment 
        array2D &operator=(const array2D &source) {
            // this resize
            resize ( source.M, source.N );
            // copy data pointer to
            memcpy ( ptr, source.ptr, source.M * source.N * sizeof(T) );
            return *this;
        }

        // cuda assignment 
        array2D &operator<<(const array2D &source) {
            this->ptr = source.ptr;
            this->pitchBytes = source.pitchBytes;
            return *this;
        }


        // indexing
#ifdef __CUDACC__
        __host__ __device__
#endif
        T &operator()( const unsigned int m, const unsigned int n ) const {
#ifndef __CUDA_ARCH__
            if(boundscheck)
                ASSERT_MSG((n < N) && (m < M), "array2D bounds check failed, !((%u < %u) && (%u < %u))", n, N, m, M);

            return ptr[n + N * m];
#else
			T *row = (T*)( (char*)ptr + pitchBytes * m);
            return row[n];
#endif
        } 

        void resize(int _M, int _N) {

            M = _M;
            N = _N;
           
            if(ptr) 
                delete ptr;
            
            if(M != 0 && N != 0) {
              ptr = new T[M * N];
              memset ( ptr, 0, M * N * sizeof(T) );
            }
        }

}; 

#endif
