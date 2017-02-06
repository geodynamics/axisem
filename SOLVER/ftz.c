//
//    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
//                    Simon Stähler, Kasra Hosseini, Stefanie Hempel
//
//    This file is part of AxiSEM.
//    It is distributed from the webpage <http://www.axisem.info>
//
//    AxiSEM is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    AxiSEM is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with AxiSEM.  If not, see <http://www.gnu.org/licenses/>.
//
#include <stdint.h>


#if defined(__i386__) || defined(__x86_64__)
#ifndef _CRAYC
#include <xmmintrin.h>
#endif
#endif

#if defined(__PPC__) || defined(__PPC64__)
#include <altivec.h>    // -maltivec still required ?
#endif


// A Fortran callable function to activate flush to zero for denormal float handling
// http://software.intel.com/en-us/articles/how-to-avoid-performance-penalties-for-gradual-underflow-behavior
// http://stackoverflow.com/questions/9314534/why-does-changing-0-1f-to-0-slow-down-performance-by-10x

#if defined(__GNUC__) || defined(linux)
  #define set_ftz set_ftz_
#endif


void set_ftz(){

#if defined(__i386__) || defined(__x86_64__)
#ifndef _CRAYC
  _MM_SET_FLUSH_ZERO_MODE (_MM_FLUSH_ZERO_ON);
#endif

#elif defined(__PPC__) || defined(__PPC64__)

//    Altivec non-IEEE mode for subnormal (denormalized) values.
//  m*vscr requires vector types even for writing to registers (disturbing)
//  so the high order bits are index'd.
  
  vector unsigned short vscr = vec_mfvscr();
  vscr[1] |= 1;   // (1<<16) in reg
  vec_mtvscr(vscr);

#endif
}

