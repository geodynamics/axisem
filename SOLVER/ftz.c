#include <xmmintrin.h>

// A Fortran callable function to activate flush to zero for denormal float handling
// http://software.intel.com/en-us/articles/how-to-avoid-performance-penalties-for-gradual-underflow-behavior
// http://stackoverflow.com/questions/9314534/why-does-changing-0-1f-to-0-slow-down-performance-by-10x

#if defined(__GNUC__) || defined(linux)
  #define set_ftz set_ftz_
#endif

void set_ftz(){
  _MM_SET_FLUSH_ZERO_MODE (_MM_FLUSH_ZERO_ON);
}

