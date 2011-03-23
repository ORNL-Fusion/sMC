#ifndef UTIL_HPP_
#define UTIL_HPP_

#include <execinfo.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <csignal>
#include <algorithm>

#define BACKTRACE(N)                                            \
  {                                                             \
  void* tracePtrs[N];                                         \
  int count = backtrace(tracePtrs, N);                        \
  char** funcNames = backtrace_symbols(tracePtrs, count);     \
  printf("Backtrace:\n");                                     \
  if(funcNames != NULL)                                       \
    for(int i = 0; i < count; i++)                            \
      if(funcNames[i][0] != '\0')                             \
	printf("%s\n", funcNames[i]);                         \
  free(funcNames);                                            \
  }

#define ASSERT_MSG(s, msg...)                                       \
  {								      \
    if(!(s)) {                                                          \
      printf("Failed to assert '%s' in %s on line %d\n\n", #s, __FILE__, __LINE__); \
      printf(msg);                                                      \
      printf("\n");                                                     \
      BACKTRACE(20);                                                    \
      raise(SIGINT);                                                    \
    }                                                                   \
  }

#define ASSERT(s)                                                       \
  {                                                                     \
    if(!(s)) {                                                          \
      printf("Failed to assert '%s' in %s on line %d\n\n", #s, __FILE__, __LINE__); \
      BACKTRACE(20);                                                    \
      raise(SIGINT);                                                    \
    }                                                                   \
  }

#endif

