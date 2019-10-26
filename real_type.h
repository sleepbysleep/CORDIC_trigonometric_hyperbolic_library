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
#ifndef __REAL_TYPE_H__
#define __REAL_TYPE_H__

#include <stdbool.h>
#include <float.h>

#include <common.h>

#if defined(USE_FLOAT_AS_REAL)
typedef float real_t;
# define REAL_EPSILON FLT_EPSILON

///////////// Basic operations
INLINE real_t real_abs(real_t x) { return fabsf(x); }
INLINE real_t real_mod(real_t x, real_t y) { return fmodf(x); }
INLINE real_t real_remainder(real_t x, real_t y) { return remainderf(x, y); }
INLINE real_t real_remquo(real_t x, real_t y, int *quo) { return remquof(x, y, quo); }
INLINE real_t real_multiply_add(real_t x, real_t y, real_t z) { return fmaf(x, y, z); }
INLINE real_t real_max(real_t x, real_t y) { return fmaxf(x, y); }
INLINE real_t real_min(real_t x, real_t y) { return fminf(x, y); }
// the positive difference: real_max(x-y, 0)
INLINE real_t real_dim(real_t x, real_t y) { return fdimf(x, y); }
INLINE real_t real_NaN(const char *arg) { return nanf(arg); }
///////////// Exponential functions
INLINE real_t real_exp(real_t arg) { return expf(arg); }
INLINE real_t real_exp2(real_t n) { return exp2f(n); }
// exp(arg) - 1
INLINE real_t real_expm1(real_t arg) { return expm1f(arg); }
INLINE real_t real_log(real_t arg) { return logf(arg); }
INLINE real_t real_log10(real_t arg) { return log10f(arg); }
INLINE real_t real_log2(real_t arg) { return log2f(arg); }
// log(1+arg)
INLINE real_t real_log1p(real_t arg) { return log1pf(arg); }
///////////// Power functions
INLINE real_t real_pow(real_t base, real_t exponent) { return powf(base, exponent); }
INLINE real_t real_sqrt(real_t arg) { return sqrtf(arg); }
INLINE real_t real_cbrt(real_t arg) { return cbrtf(arg); }
INLINE real_t real_hypot(real_t x, real_t y) { return hypotf(x, y); }
///////////// Trigonometric functions
INLINE real_t real_sin(real_t x) { return sinf(x); }
INLINE real_t real_cos(real_t x) { return cosf(x); }
INLINE real_t real_tan(real_t x) { return tanf(x); }
INLINE real_t real_asin(real_t x) { return asinf(x); }
INLINE real_t real_acos(real_t x) { return acosf(x); }
INLINE real_t real_atan(real_t x) { return atanf(x); }
INLINE real_t real_atan2(real_t y, real_t x) { return atan2f(y, x); }
// TODO:
INLINE real_t real_tan2(real_t y, real_t x) { return tan2f(y, x); }
///////////// Hyperbolic functions
INLINE real_t real_sinh(real_t arg) { return sinhf(arg); }
INLINE real_t real_cosh(real_t arg) { return coshf(arg); }
INLINE real_t real_tanh(real_t arg) { return tanhf(arg); }
INLINE real_t real_asinh(real_t arg) { return asinhf(arg); }
INLINE real_t real_acosh(real_t arg) { return acoshf(arg); }
INLINE real_t real_atanh(real_t arg) { return atanhf(arg); }
///////////// Gamma functions
INLINE real_t real_tgamma(real_t arg) { return tgammaf(arg); }
INLINE real_t real_lgamma(real_t arg) { return lgammaf(arg); }
///////////// Nearest integer floating-point operations
INLINE real_t real_ceil(real_t arg) { return ceilf(arg); }
INLINE real_t real_floor(real_t arg) { return floorf(arg); }
INLINE real_t real_trunc(real_t arg) { return truncf(arg); }
INLINE real_t real_round(real_t x) { return roundf(x); }
INLINE long real_lround(real_t x) { return lroundf(x); }
INLINE long long real_llround(real_t x) { return llroundf(x); }
INLINE real_t real_nearbyint(real_t arg) { return nearbyintf(arg); }
INLINE real_t real_rint(real_t arg) { return rintf(arg); }
INLINE long real_lrint(real_t arg) { return lrintf(arg); }
INLINE long long real_llrint(real_t arg) { return llrintf(arg); }



INLINE bool real_is_zero(real_t x) { return (fabsf(x) < FLT_EPSILON); }
INLINE bool real_is_equal(real_t x, real_t y) { return (fabsf(x-y) < FLT_EPSILON); }

#elif defined(USE_DOUBLE_AS_REAL)
typedef double real_t;
# define REAL_EPSILON DBL_EPSILON
INLINE real_t real_cos(real_t x) { return cos(x); }
INLINE real_t real_acos(real_t x) { return acos(x); }
INLINE real_t real_sin(real_t x) { return sin(x); }
INLINE real_t real_asin(real_t x) { return asin(x); }
INLINE real_t real_tan(real_t x) { return tan(x); }
INLINE real_t real_atan(real_t x) { return atan(x); }
INLINE real_t real_tan2(real_t y, real_t x) { return tan2(y, x); }
INLINE real_t real_atan2(real_t y, real_t x) { return atan2(y, x); }
INLINE real_t real_floor(real_t x) { return floor(x); }
INLINE real_t real_ceil(real_t x) { return ceil(x); }
INLINE real_t real_round(real_t x) { return round(x); }
INLINE long real_lround(real_t x) { return lround(x); }
INLINE long long real_llround(real_t x) { return llround(x); }
INLINE real_t real_abs(real_t x) { return fabs(x); }
INLINE bool real_is_zero(real_t x) { return (fabs(x) < DBL_EPSILON); }
INLINE bool real_is_equal(real_t x, real_t y) { return (fabs(x-y) < DBL_EPSILON); }
#elif defined(USE_LONG_DOUBLE_AS_REAL)
typedef long double real_t;
# define REAL_EPSILON LDBL_EPSILON
INLINE real_t real_cos(real_t x) { return cosl(x); }
INLINE real_t real_acos(real_t x) { return acosl(x); }
INLINE real_t real_sin(real_t x) { return sinl(x); }
INLINE real_t real_asin(real_t x) { return asinl(x); }
INLINE real_t real_tan(real_t x) { return tanl(x); }
INLINE real_t real_atan(real_t x) { return atanl(x); }
INLINE real_t real_tan2(real_t y, real_t x) { return tan2l(y, x); }
INLINE real_t real_atan2(real_t y, real_t x) { return atan2l(y, x); }
INLINE real_t real_floor(real_t x) { return floorl(x); }
INLINE real_t real_ceil(real_t x) { return ceill(x); }
INLINE real_t real_round(real_t x) { return roundl(x); }
INLINE long real_lround(real_t x) { return lroundl(x); }
INLINE long long real_llround(real_t x) { return llroundl(x); }
INLINE real_t real_abs(real_t x) { return fabsl(x); }
INLINE bool real_is_zero(real_t x) { return (fabsl(x) < LDBL_EPSILON); }
INLINE bool real_is_equal(real_t x, real_t y) { return (fabsl(x-y) < LDBL_EPSILON); }
#else
# error "Undefined real_t!: data-type"
#endif


#endif /* __REAL_TYPE_H__ */
