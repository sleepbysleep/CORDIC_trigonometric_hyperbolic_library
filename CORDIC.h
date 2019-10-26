/*
 Copyright 2011 Hoyoung Lee.

 This program is free software; you can redistribute it and/or modify it
 under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, please visit www.gnu.org.
*/
#ifndef __CORDIC_H__
#define __CORDIC_H__

#include <floatpt/realtype.h>
//#include "realtype.h"

#ifdef __cplusplus
extern "C" {
#endif

  // In linear mode
  real_t CORDIC_multiply(real_t a, real_t b, int iteration);
  real_t CORDIC_divide(real_t a, real_t b, int iteration);

  // In circular mode
  void CORDIC_compute_cos_and_sin(real_t *cos_t, real_t *sin_t, real_t t, int iteration);
  real_t CORDIC_atan2(real_t y, real_t x, int iteration);
  real_t CORDIC_enorm2(real_t x, real_t y, int iteration);

  // In hyperbolic mode
  void CORDIC_compute_cosh_and_sinh(real_t *cosh_t, real_t *sinh_t, real_t t, int iteration);
  real_t CORDIC_exp(real_t t, int iteration);
  real_t CORDIC_atanh(real_t x, int iteration);
  real_t CORDIC_atanh2(real_t y, real_t x, int iteration);
  real_t CORDIC_ln(real_t v, int iteration);
  real_t CORDIC_sqrt(real_t v, int iteration);

#ifdef __cplusplus
}
#endif

#endif /* __CORDIC_H__ */
