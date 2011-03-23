#ifndef DATA_HPP_
#define DATA_HPP_

#include "util.hpp"
#include <iostream>

using namespace std;

template <class T, bool boundscheck = false>
class Data1D {
private:
  T *ptr;
  unsigned int N;

public:
  Data1D() : N(0), ptr(NULL) {};

  Data1D(const unsigned int _N) : ptr(NULL) {
    resize(_N);
  }

  Data1D(const Data1D &other) {
    *this = other;
  }

  ~Data1D() {
    if(ptr)
      delete ptr;
  }

  Data1D &operator=(const Data1D &other) {
    resize(other.N);
    memcpy(ptr, other.ptr, other.N * sizeof(T));

    return *this;
  }

  T &operator()(const unsigned int n) {
    if(boundscheck)
      ASSERT_MSG(n < N, "Data1D bounds check failed, !(%u < %u)", n, N);

    return ptr[n];
  }

  void resize(int _N) {
    N = _N;
    
    if(ptr)
      delete ptr;

    if(N != 0) {
      ptr = new T[N];
      memset(ptr, 0, N * sizeof(T));
    }
  }

  void print () {
	for(int i=0;i<N;i++)
		cout << ptr[i] << endl;
  }

};

template <class T, bool boundscheck = false>
class Data2D {
private:
  T *ptr;
  unsigned int N, M;

public:
  Data2D() : N(0), M(0), ptr(NULL) {};

  Data2D(const unsigned int _N, const unsigned int _M) : ptr(NULL) {
    resize(_N, _M);
  }

  Data2D(const Data2D &other) {
    *this = other;
  }

  ~Data2D() {
    if(ptr)
      delete ptr;
  }

  Data2D &operator=(const Data2D &other) {
    resize(other.N, other.M);
    memcpy(ptr, other.ptr, other.N * other.M * sizeof(T));

    return *this;
  }

  T &operator()(const unsigned int n, const unsigned int m) {
    if(boundscheck)
      ASSERT_MSG((n < N) && (m < M), "Data2D bounds check failed, !((%u < %u) && (%u < %u))", n, N, m, M);

    return ptr[n + N * m];
  }

  void resize(int _N, int _M) {
    N = _N;
    M = _M;
    
    if(ptr)
      delete ptr;

    if(N != 0 && M != 0) {
      ptr = new T[N * M];
      memset(ptr, 0, N * M * sizeof(T));
    }
  }
};

#endif
