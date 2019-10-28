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
// CORDIC(COordinate Rotation DIgital Computer)
// - efficient algorithm to calculate trigoometric functions.
// Reference http://en.wikibooks.org/wiki/Digital_Circuits/CORDIC
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <assert.h>
#include <stdbool.h>

#include "cordic.h"

#ifdef USE_FLOAT_AS_REAL

#define NR_CIRCULAR_DELTA_RADIANS 13

const real_t circular_delta_radians[NR_CIRCULAR_DELTA_RADIANS] = {
	7.853981852531433105e-01,
	4.636476039886474609e-01,
	2.449786663055419922e-01,
	1.243549957871437073e-01,
	6.241881102323532104e-02,
	3.123983368277549744e-02,
	1.562372874468564987e-02,
	7.812341209501028061e-03,
	3.906230209395289421e-03,
	1.953122555278241634e-03,
	9.765622089616954327e-04,
	4.882812208961695433e-04,
	2.441406250000000000e-04
};

#define NR_CIRCULAR_K_PRODUCTS 11

const real_t circular_K_products[NR_CIRCULAR_K_PRODUCTS] = {
	7.071067690849304199e-01,
	6.324555277824401855e-01,
	6.135720014572143555e-01,
	6.088339090347290039e-01,
	6.076482534408569336e-01,
	6.073517799377441406e-01,
	6.072776317596435547e-01,
	6.072590947151184082e-01,
	6.072544455528259277e-01,
	6.072533130645751953e-01,
	6.072530150413513184e-01
};

#define NR_HYPERBOLIC_DELTA_RADIANS 12

const real_t hyperbolic_delta_radians[NR_HYPERBOLIC_DELTA_RADIANS] = {
	5.493061542510986328e-01,
	2.554128170013427734e-01,
	1.256572157144546509e-01,
	6.258156895637512207e-02,
	3.126017749309539795e-02,
	1.562627218663692474e-02,
	7.812659256160259247e-03,
	3.906270023435354233e-03,
	1.953127561137080193e-03,
	9.765628492459654808e-04,
	4.882813082076609135e-04,
	2.441406250000000000e-04
};

#define NR_HYPERBOLIC_K_PRODUCTS 11

const real_t hyperbolic_K_products[NR_HYPERBOLIC_K_PRODUCTS] = {  
	1.154700517654418945e+00,
	1.192569613456726074e+00,
	1.201997160911560059e+00,
	1.204351663589477539e+00,
	1.204940199851989746e+00,
	1.205087304115295410e+00,
	1.205124139785766602e+00,
	1.205133318901062012e+00,
	1.205135583877563477e+00,
	1.205136179924011230e+00,
	1.205136299133300781e+00
};

#else //if USE_DOUBLE_AS_REAL
#define NR_CIRCULAR_DELTA_RADIANS 28

const real_t circular_delta_radians[NR_CIRCULAR_DELTA_RADIANS] = {
	7.853981633974482790e-01,
	4.636476090008060935e-01,
	2.449786631268641435e-01,
	1.243549945467614382e-01,
	6.241880999595735002e-02,
	3.123983343026827744e-02,
	1.562372862047683129e-02,
	7.812341060101111144e-03,
	3.906230131966971757e-03,
	1.953122516478818758e-03,
	9.765621895593194594e-04,
	4.882812111948982899e-04,
	2.441406201493617712e-04,
	1.220703118936702079e-04,
	6.103515617420877259e-05,
	3.051757811552609573e-05,
	1.525878906131576154e-05,
	7.629394531101969981e-06,
	3.814697265606496142e-06,
	1.907348632810186965e-06,
	9.536743164059608441e-07,
	4.768371582030888423e-07,
	2.384185791015579737e-07,
	1.192092895507806809e-07,
	5.960464477539055221e-08,
	2.980232238769530257e-08,
	1.490116119384765460e-08,
	7.450580596923828125e-09
};

#define NR_CIRCULAR_K_PRODUCTS 26

const real_t circular_K_products[NR_CIRCULAR_K_PRODUCTS] = {
	7.071067811865474617e-01,
	6.324555320336757713e-01,
	6.135719910778962838e-01,
	6.088339125177524291e-01,
	6.076482562561682510e-01,
	6.073517701412960434e-01,
	6.072776440935261366e-01,
	6.072591122988928447e-01,
	6.072544793325624912e-01,
	6.072533210898752865e-01,
	6.072530315291344571e-01,
	6.072529591389449477e-01,
	6.072529410413972650e-01,
	6.072529365170102889e-01,
	6.072529353859135171e-01,
	6.072529351031393796e-01,
	6.072529350324458175e-01,
	6.072529350147723992e-01,
	6.072529350103540446e-01,
	6.072529350092494838e-01,
	6.072529350089733713e-01,
	6.072529350089043154e-01,
	6.072529350088871070e-01,
	6.072529350088827771e-01,
	6.072529350088816669e-01,
	6.072529350088814448e-01
};

#define NR_HYPERBOLIC_DELTA_RADIANS 26

const real_t hyperbolic_delta_radians[NR_HYPERBOLIC_DELTA_RADIANS] = {
	5.493061443340547800e-01,
	2.554128118829953609e-01,
	1.256572141404530274e-01,
	6.258157147700300904e-02,
	3.126017849066699272e-02,
	1.562627175205220931e-02,
	7.812658951540421212e-03,
	3.906269868396826210e-03,
	1.953127483532549752e-03,
	9.765628104410359420e-04,
	4.882812888051127689e-04,
	2.441406298506386082e-04,
	1.220703131063298193e-04,
	6.103515632579122063e-05,
	3.051757813447390088e-05,
	1.525878906368423676e-05,
	7.629394531398029172e-06,
	3.814697265643503435e-06,
	1.907348632814812823e-06,
	9.536743164065390500e-07,
	4.768371582031611048e-07,
	2.384185791015669999e-07,
	1.192092895507818059e-07,
	5.960464477539069117e-08,
	2.980232238769531912e-08,
	1.490116119384765625e-08
};

#define NR_HYPERBOLIC_K_PRODUCTS 26

const real_t hyperbolic_K_products[NR_HYPERBOLIC_K_PRODUCTS] = {
	1.154700538379251684e+00,
	1.192569587999887881e+00,
	1.201997162280556752e+00,
	1.204351713336804908e+00,
	1.204940206757381338e+00,
	1.205087321122956912e+00,
	1.205124099152999895e+00,
	1.205133293625433977e+00,
	1.205135592241350251e+00,
	1.205136166895192318e+00,
	1.205136310558644341e+00,
	1.205136346474506848e+00,
	1.205136355453472419e+00,
	1.205136357698213700e+00,
	1.205136358259399021e+00,
	1.205136358399695462e+00,
	1.205136358434769628e+00,
	1.205136358443538169e+00,
	1.205136358445730194e+00,
	1.205136358446278200e+00,
	1.205136358446415201e+00,
	1.205136358446449396e+00,
	1.205136358446458056e+00,
	1.205136358446460276e+00,
	1.205136358446460720e+00,
	1.205136358446460942e+00 
};

//#else
//# error "Undefined real_t!: data-type"
#endif

#define CORDIC_HYPERBOLIC_THETA_MIN (-2.0*hyperbolic_delta_radians[0])
#define CORDIC_HYPERBOLIC_THETA_MAX (2.0*hyperbolic_delta_radians[0])

// In linear mode
// As the iteration goes, *x to *x, *y to *y + *x * *z, and *z to zero.
#define CORDIC_linear_rotating(x, y, z, n) CORDIC_linear_process(x, y, z, n, true)
#define CORDIC_linear_vectoring(x, y, z, n) CORDIC_linear_process(x, y, z, n, false)

static void CORDIC_linear_process(real_t *x, real_t *y, real_t *z, int iteration, bool rotate_or_vector)
{
  int n, dir;
  real_t delta, power_of_2;

  power_of_2 = 1;
  for (n = 0; n < iteration; n++) {
    if (rotate_or_vector) {
      if (*z < 0) dir = -1;
      else dir = +1;
    } else {
      if (*y > 0) dir = -1;
      else dir = +1;
    }

    delta = dir * power_of_2;

    *y += *x * delta;
    *z -= delta;

    // No more process due to precise limitation
    if (power_of_2 < REAL_EPSILON) break;

    power_of_2 /= 2;
  }
}

real_t CORDIC_multiply(real_t a, real_t b, int iteration)
{
  real_t c = 0.0;

#if 0
  assert(b >= -2.0 && b <= 2.0);

  CORDIC_linear_rotating(&a, &c, &b, iteration);

  return c;
#else
  int shift = 0;

  while (!(b >= -2.0 && b <= 2.0)) {
    b /= 2.0;
    shift++;
  }
  assert(shift <= 32);

  CORDIC_linear_rotating(&a, &c, &b, iteration);

  return c * (real_t)(1UL << shift);
#endif
}

real_t CORDIC_divide(real_t a, real_t b, int iteration)
{
  real_t c = 0.0;

#if 0
  assert(a >= -2.0*b && a <= 2.0*b);

  CORDIC_linear_vectoring(&b, &a, &c, iteration);

  return c;
#else
  int sign = +1, mul = 0;

  if ((a > 0 && b < 0) || (a < 0 && b > 0)) sign = -1;

  a = ABS(a), b = ABS(b);
  while (!(a >= -2.0*b && a <= 2.0*b)) {
    a  -= b;
    mul++;
  }
  CORDIC_linear_vectoring(&b, &a, &c, iteration);
  return (c + mul) * sign;
#endif
}

// In circular mode
/**
   The operation in circular mode
   x' = cos(z) * x - sin(z) * y
   y' = sin(z) * x + cos(z) * y
   where cos(z) = 1 / sqrt(1 + tan(z)^2), and sin(z) = tan(z) / sqrt(1 + tan(z)^2)
   x' = 1/sqrt(1 + tan(z)^2) * (x - tan(z) * y)
   y' = 1/sqrt(1 + tan(z)^2) * (y + tan(z) * x)
   Let tan(z) be a radix 2 number.
*/
static void CORDIC_generate_circular_delta_radians(int max_iteration)
{
  int i;
  real_t value, *delta_angles;

  delta_angles = (real_t *)malloc(max_iteration * sizeof(*delta_angles));

  for (value = 1.0, i = 0; i < max_iteration; i++) {
    delta_angles[i] = atan(value);
    value /= 2.0;
  }

  for (i = 0; i < max_iteration-1; i++) {
    printf("%.18e", delta_angles[i]);
    if (fabs((delta_angles[i]/delta_angles[i+1]) - 2.0) < REAL_EPSILON) {
      printf("\nFurther things are meaningless due to having a value divided by 2 of the last.\n");
      break;
    }
    printf(",\n");
  }
  free(delta_angles);
}

static void CORDIC_generate_circular_K_cumulative_products(int max_iteration)
{
  int i;
  real_t *k_products, value;

  k_products = (real_t *)malloc(max_iteration * sizeof(*k_products));

  for (value = 1.0, i = 0; i < max_iteration; i++) {
    value *= 1.0 / sqrt(1.0 + pow(2.0, -2.0*i));
    k_products[i] = value;
  }

  for (i = 0; i < max_iteration-1; i++) {
    printf("%.18e", k_products[i]);
    if (fabs(k_products[i] - k_products[i+1]) < REAL_EPSILON) {
      printf("\nFurther things are meaningless due to having the same value with the last.\n");
      break;
    }
    printf(",\n");
  }
  free(k_products);
}

/* CORDIC circular rotating operation
   x[i+1] <= (x[i] - dir * 2^-i * y[i]) / sqrt(1 + 2^(-2*i))
   y[i+1] <= (y[i] + dir * 2^-i * x[i]) / sqrt(1 + 2^(-2*i))
   z[i+1] <= z[i] - dir * atan(2^-i)
   where dir = sgn(z[i])
   As z[n] goes to zero, these become
   x[i] = x[0] * cos(z[0]) - y[0] * sin(z[0]), and
   y[i] = y[0] * cos(z[0]) + x[0] * sin(z[0]), respectively.
   If x[0] = 1, and y[0] = 0, and z[0] = t, then x[n] = cos(t), y[n] = sin(t).

   CORDIC circular vectoring operation
   x[i+1] <= (x[i] - dir * 2^-i * y[i]) / sqrt(1 + 2^(-2*i))
   y[i+1] <= (y[i] + dir * 2^-i * x[i]) / sqrt(1 + 2^(-2*i))
   z[i+1] <= z[i] - dir * atan(2^-i)
   where dir = -sgn(y[i])
   As y[n] goes to zero, these become
   x[n] = sqrt(x[0] * x[0] + y[0] * y[0]) / Kn, and
   z[n] = atan(y[0] / x[0]), respectively.
   Here Kn is the cumulative products of scale factors belonging to atan(2^-i).
   Kn = 1/(sqrt(1+1) * sqrt(1+2^-1) * sqrt(1+2^-2) ... sqrt(1+2^-n)).
*/
// 'z' has to be placed in (-M_PI_2, M_PI_2).
#define CORDIC_circular_rotating(x, y, z, n) CORDIC_circular_process(x, y, z, n, true)
#define CORDIC_circular_vectoring(x, y, z, n) CORDIC_circular_process(x, y, z, n, false)

static void CORDIC_circular_process(real_t *x, real_t *y, real_t *z, int iteration, bool rotate_or_vector)
{
  int n, dir;
  real_t x_nth, y_nth;
  real_t dt, tan_dt, power_of_2;
  real_t K_n = circular_K_products[MIN(iteration, NR_CIRCULAR_K_PRODUCTS-1)];

  dt = circular_delta_radians[0]; // M_PI_4
  power_of_2 = 1; // 2^0 = tan(M_PI_4)

  for (n = 0; n < iteration; n++) {
    if (rotate_or_vector) { // rotating
      if (*z < 0) dir = -1;
      else dir = +1;
    } else { // vectoring
      if (*y > 0) dir = -1;
      else dir = +1;
    }

    // Applying the rotation matrix with angle of dt
    tan_dt = dir * power_of_2;
    x_nth = *x - tan_dt * *y;
    y_nth = *y + tan_dt * *x;
    *x = x_nth, *y = y_nth;

    // Update the current angle changed by dt.
    *z -= dir * dt;

    // Update the next delta-angle for rotating.
    if ((n+1) < NR_CIRCULAR_DELTA_RADIANS) dt = circular_delta_radians[n+1];
    else dt /= 2; // atan(dt) is very closer to the half of the previous if dt is very small.

    power_of_2 /= 2;

    // No more process due to precise limitation.
    if (dt < REAL_EPSILON || power_of_2 < REAL_EPSILON) break;
  }

  // Putting pre-calculated nth K value and results together.
  *x *= K_n;
  *y *= K_n;
}

void CORDIC_compute_cos_and_sin(real_t *cos_t, real_t *sin_t, real_t t, int iteration)
{
  real_t sign_change;

  // 'theta' manages to be ranged between -M_PI_2 and M_PI_2.
  // Corresponding sign changes of x and y has to be calculated.
  sign_change = +1;
  while (!(t >= -M_PI_2 && t <= M_PI_2)) {
    if (t < 0) t += M_PI;
    else t -= M_PI;
    sign_change *= -1;
  }

  *cos_t = 1.0, *sin_t = 0.0;
  CORDIC_circular_rotating(cos_t, sin_t, &t, iteration);

  *cos_t *= sign_change, *sin_t *= sign_change;
}

real_t CORDIC_atan2(real_t y, real_t x, int iteration)
{
  real_t theta, t;
#if 0
  assert(y >= -2.0 && y <= 2.0);
#else
  int shift = 0;

  while (!(y >= -2.0 && y <= 2.0)) {
    y /= 2.0;
    x /= 2.0;
    shift++;
  }
  assert(shift <= 32);
#endif

  // Reflecting in the origin to make the position toward positive X-axis.
  theta = 0;
  if (x < 0) {
    if (y > 0) theta += M_PI;
    else theta -= M_PI;
    x *= -1, y *= -1;
  }

  t = 0.0;
  CORDIC_circular_vectoring(&x, &y, &t, iteration);
  // which makes that x becomes the distance from origin, y to zero, and theta to the atan2(y, x)

  return theta + t;
}

real_t CORDIC_enorm2(real_t x, real_t y, int iteration)
{
  real_t theta;
  int shift = 0;
#if 0
  assert(y >= -2.0 && y <= 2.0);
#else
  while (!(y >= -2.0 && y <= 2.0)) {
    y /= 4.0;
    x /= 4.0;
    shift += 2;
  }
#endif

  // Reflecting in the origin to make the position toward positive X-axis.
  if (x < 0) {
    x *= -1, y *= -1;
  }
  theta = 0;
  CORDIC_circular_vectoring(&x, &y, &theta, iteration);
  // which makes that x becomes the distance from origin, y to zero, and theta to the atan2(y, x).

  return x * (real_t)(1UL<<shift);
} 

// In hyperbolic mode
/* The operation in hyperbolic mode
   x' <= cosh(z) * x + sinh(z) * y
   y' <= sinh(z) * x + cosh(z) * y
   where cosh(z) = 1 / sqrt(1 - tanh(z)^2), and sinh(z) = tanh(z) / sqrt(1 - tanh(z)^2).
   x' <= 1/sqrt(1 - tanh(z)^2) * (x + tanh(z) * y)
   y' <= 1/sqrt(1 - tanh(z)^2) * (y + tanh(z) * z)
   Let tanh(z) be a radix 2 number.
*/
static void CORDIC_generate_hyperbolic_delta_radians(int max_iteration)
{
  int i;
  real_t value, *delta_angles;

  delta_angles = (real_t *)malloc(max_iteration * sizeof(*delta_angles));

  // since tanh(1) = 0, 1/tanh(1) = INF, it is need to be excluded by starting from i = 1.
  for (value = 0.5, i = 1; i <= max_iteration; i++) {
    delta_angles[i-1] = atanh(value);
    value /= 2.0;
  }

  for (i = 0; i < max_iteration-1; i++) {
    printf("%.18e", delta_angles[i]);

    if (ABS((delta_angles[i]/delta_angles[i+1]) - 2.0) < REAL_EPSILON) {
      printf("\nFurther things are meaningless due to having the value divided by 2 of the previous.\n");
      break;
    }
    printf(",\n");
  }
  printf("\n");
  free(delta_angles);
}

static void CORDIC_generate_hyperbolic_K_cumulative_products(int max_iteration)
{
  int i;
  real_t *k_products, value;

  k_products = (real_t *)malloc(max_iteration * sizeof(*k_products));
  
  // Since tanh(1) = 0, 1/tanh(1) = INF, it is needed to be excluded as starting from i =1.
  for (value = 1.0, i = 1; i < max_iteration; i++) {
    value *= 1.0 / sqrt(1.0 - pow(2.0, -2.0 * i));
    k_products[i-1] = value;
  }

  for (i = 0; i < max_iteration-1; i++) {
    printf("%.18e", k_products[i]);
    if (ABS(k_products[i] - k_products[i+1]) < REAL_EPSILON) {
      printf("\nFurther things are meaningless due to having the same value with the last.\n");
      break;
    }
    printf(",\n");
  }
  printf("\n");
  free(k_products);
}

/* CORDIC hyperbolic rotating operation
   x[i+1] <= (x[i] + dir * 2^-i * y[i]) / sqrt(1 - 2^(-2*i))
   y[i+1] <= (y[i] + dir * 2^-i * x[i]) / sqrt(1 - 2^(-2*i))
   where dir = sgn(z[i])
   As z[n] = 0, these become
   x[n] = x[1] * cosh(z[1]) + y[1] * sinh(z[1]), and
   y[n] = y[1] * cosh(z[1]) + x[1] * sinh(z[1]), respectively.
   If x[1] =1, y[1] = 0, and z[1] = t, then x[n] = cosh(t), y[n] = sinh(t).

   CORDIC hyperbolic vectoring operation
   x[i+1] <= (x[i] + dir * 2^-i * y[i]) / sqrt(1 - 2^(-2*i))
   y[i+1] <= (y[i] + dir * 2^-i * x[i]) / sqrt(1 - 2^(-2*i))
   z[i+1] <= z[i] - dir * atanh(2^-i)
   where dir = -sgn(y[i])
   As y[n] = 0 under y[1] < x[1], these become
   x[n] = sqrt(x[1] * x[1] - y[1] * y[1]) / Kn, and
   z[n] = atanh(y[1]/x[1]), respectively.
   Here Kn is the cumulative products of scale factors belonging to atan(2^-i).
   Kn = 1 / (sqrt(1 - 2^-1) * sqrt(1 - 2^-2) ... sqrt(1 - 2^-n)).
*/
/* 'z' has to be placed in (-2*delta_theta[0], 2*delta_theta[0]),
   which is (-1.0986122886681096, 1.0986122886681096)
   actually covering to (-62.95 degree, 62.95 degree).
   And, iterations 4, 13, 40, 121, ..., j, 3j+1, ... must be repeated twice.
*/
#define CORDIC_hyperbolic_rotating(x, y, z, n) CORDIC_hyperbolic_process(x, y, z, n, true)
// Given y < x,
#define CORDIC_hyperbolic_vectoring(x, y, z, n) CORDIC_hyperbolic_process(x, y, z, n, false)

static void CORDIC_hyperbolic_process(real_t *x, real_t *y, real_t *z, int iteration, bool rotate_or_vector)
{
  int n, dir, k;
  real_t x_nth, y_nth;
  real_t dt, tanh_dt, power_of_2;
  real_t K_n = hyperbolic_K_products[MIN(iteration, NR_HYPERBOLIC_K_PRODUCTS-1)];

  dt = hyperbolic_delta_radians[0]; // atanh(2^-1)
  power_of_2 = 0.5; // 2^-1

  k = 4;
  for (n = 0; n < iteration; n++) {
    if (rotate_or_vector) { // rotating
      if (*z < 0) dir = -1;
      else dir = +1;
    } else { // vectoring
      if (*y > 0) dir = -1;
      else dir = +1;
    }

    // Rotating with the angle of dt along with 1st-hyperbolic curve.
    tanh_dt = dir * power_of_2;
    x_nth = *x + tanh_dt * *y;
    y_nth = *y + tanh_dt * *x;
    *x = x_nth, *y = y_nth;

    // Update the current angle of rotation.
    *z -= dir *dt;

    // For more precise result, iteration 4, 13, 40, 121, ..., j, 3j+1, ... must be repeated twice.
    if ((n+1) == k) {
      if (rotate_or_vector) { // rotating
	if (*z < 0) dir = -1;
	else dir = +1;
      } else { // vectoring
	if (*y > 0) dir = -1;
	else dir = +1;
      }

      // Rotating with the angle of dt along with 1st-hyperbolic curve.
      tanh_dt = dir * power_of_2;
      x_nth = *x + tanh_dt * *y;
      y_nth = *y + tanh_dt * *x;
      *x = x_nth, *y = y_nth;

      // Update the current rotating angle
      *z -= dir * dt;

      k = 3*k + 1;
    }

    // Update the next rotating delta-angle
    if ((n+1) < NR_HYPERBOLIC_DELTA_RADIANS) dt = hyperbolic_delta_radians[n+1];
    else dt /= 2; // atan(dt) is very closer to the half of the previous if dt is vary small.

    power_of_2 /= 2;

    // No more process due to precise limitation.
    if (dt < REAL_EPSILON || power_of_2 < REAL_EPSILON) break;
  }

  // Putting pre-calculated nth K value and results together.
  *x *= K_n;
  *y *= K_n;
}

// http://www.nahee.com/spanky/pub/fractals/docs/cordic.math
/*
  tan(z) = sin(z)/cos(z)
  tanh(z) = sinh(z) / cosh(z)
  exp(z) = sinh(z) + cosh(z)
  ln(w) = 2 * atanh(y/x) where x = w+1, and y = w-1.
  sqrt(w) = sqrt(x*x - y*y) where x = w + 1/4, and y = w - 1/4.
*/
/* For expansion
   sinh(Q*Ln(2)+D) = 2^(Q-1) * (cosh(D) + sinh(D) - (cosh(D) - sinh(D)) / 2^(2*Q))
   cosh(Q*Ln(2)+D) = 2^(Q-1) * (cosh(D) + sinh(D) + (cosh(D) - sinh(D)) / 2^(2*Q))
   tanh(Q*Ln(2)+D) = sinh(Q*Ln(2)+D) / cosh(Q*Ln(2)+D)
   exp(Q*Ln(2)+D) = 2^Q * (cosh(D) + sinh(D)) or 2^-Q * (cosh(D)  sinh(D))

   atanh(1-M/2^E) = atanh(T) + Ln(2) * E/2, where T = (2 - M - M/2^E) / (2 + M - M/2^E).
   Ln(M*2^E) = Ln(M) + E*Ln(2), where 0.5 <= M < 1.0.
   sqrt(M*2^E) = 2^(E/2) * sqrt(M) if E mod 2 = 0, 0.5 <= M < 1.0
   2^((E+1)/2) * sqrt(M/2) if E mod 2 =1, 0.25 <= M < 0.4.
*/

void CORDIC_compute_cosh_and_sinh(real_t *cosh_t, real_t *sinh_t, real_t t, int iteration)
{
  /*
    assert(t >= CORDIC_HYPERBOLIC_THETA_MIN);
    assert(t <= CORDIC_HYPERBOLIC_THETA_MAX);

    CORDIC_hyperbolic_rotating(cosh_t, sinh_t, &t, iteration);
  */

  // sinh(Q*Ln(2)+D) = 2^(Q-1) * (cosh(D) + sinh(D) - (cosh(D) - sinh(D)) / 2^(2*Q))
  // cosh(Q*Ln(2)+D) = 2^(Q-1) * (cosh(D) + sinh(D) + (cosh(D) - sinh(D)) / 2^(2*Q))
  // sinh(Q*Ln(2)+D) = ((cosh(D) + sinh(D)) * 2^Q - (cosh(D) - sinh(D))/2^Q) / 2
  // cosh(Q*Ln(2)+D) = ((cosh(D) + sinh(D)) * 2^Q + (cosh(D) - sinh(D))/2^Q) / 2

  // Ln(2) is M_LN2 about 0.69 in C
  int Q, sign;
  real_t D, cosh_D, sinh_D;

  assert(t >= CORDIC_HYPERBOLIC_THETA_MIN - 16 * M_LN2);
  assert(t <= CORDIC_HYPERBOLIC_THETA_MAX + 16 * M_LN2);

  sign = SGN(t);
  Q = 0, D = ABS(t);
  while (!(D >= CORDIC_HYPERBOLIC_THETA_MIN && D <= CORDIC_HYPERBOLIC_THETA_MAX)) {
    D -= M_LN2;
    Q++;
  }
  assert(Q <= 16);

  if (Q == 0) {
    *cosh_t = 1.0, *sinh_t = 0.0;
    CORDIC_hyperbolic_rotating(cosh_t, sinh_t, &t, iteration);
  } else {
    cosh_D = 1.0, sinh_D = 0.0;
    CORDIC_hyperbolic_rotating(&cosh_D, &sinh_D, &D, iteration);
    *cosh_t = ((cosh_D + sinh_D)*(1UL<<Q) - (cosh_D - sinh_D)/(1UL<<Q))/(1UL<<1);
    *sinh_t = ((cosh_D + sinh_D)*(1UL<<Q) + (cosh_D - sinh_D)/(1UL<<Q))/(1UL<<1);
    //*cosh_t = (cosh_D + sinh_D - (cosh_D - sinh_D) / (1UL<<(2*Q))) * (1UL<<(Q-1));
    //*sinh_t = (cosh_D + sinh_D + (cosh_D - sinh_D) / (1UL<<(2*Q))) * (1UL<<(Q-1));
    *sinh_t *= sign;
  }
}

real_t CORDIC_exp(real_t t, int iteration)
{
  /*
    real_t cosh_t, sinh_t;

    assert(t >= CORDIC_HYPERBOLIC_THETA_MIN);
    assert(t <= CORDIC_HYPERBOLIC_THETA_MAX);

    CORDIC_hyperbolic_rotating(&cosh_t, &sinh_t, &t, iteration);

    return cosh_t + sinh_t;
  */
  // exp(Q*Ln(2)+D) = 2^Q * (cosh(D) + sinh(D))
  int Q, sign;
  real_t D, cosh_D = 1.0, sinh_D = 0.0;

  sign = SGN(t);
  Q = 0, D = ABS(t);
  while (!(D >= CORDIC_HYPERBOLIC_THETA_MIN && D <= CORDIC_HYPERBOLIC_THETA_MAX)) {
    D -= M_LN2;
    Q++;
  }

  CORDIC_hyperbolic_rotating(&cosh_D, &sinh_D, &D, iteration);

  if (sign > 0) return (cosh_D + sinh_D) * (1UL<<Q);

  return (cosh_D  - sinh_D) / (1UL<<Q);
}

real_t CORDIC_atanh(real_t x, int iteration)
{
  int sign, E;
  real_t x1, y1, z = 0, M;

  /* atanh(1 - M/2^E) = atanh(T) + Ln(2) * E^2,
     where T = (2 - M - M/2^E) / (2 + M - M/2^E), 0.5 <= M < 1, and E >= 1 integer.
     Thus, if 3/4 <= y/x = 1-M/2^E < 1,
     then we obtain M be repeatedly multiplying 1-y/x by 2 until we get
     0.5 <= 2^E * (1 - y/x) = M < 1.
     To compute atanh(T) we use the given x and y to compute new values.
     x1 <= 1 + M + y/x, and y1 <= 1 - M + y/x. Or
     x1 <= x + y + M*x, and y1 <= x + y - M*x. */

  // ABS(x) < 1
  assert(ABS(x) < 1);

  sign = SGN(x);
  x = ABS(x);

  E = 0, M = 1.0 - x;
  while (!(M >= 0.5 && M < 1)) {
    M *= 2.0;
    E++;
  }

  x1 = 1 + M + x, y1 = 1 - M + x;

  z = 0.0;
  CORDIC_hyperbolic_vectoring(&x1, &y1, &z, iteration);

  return (z + M_LN2 * ((real_t)E/2)) * sign;
}

real_t CORDIC_atanh2(real_t y, real_t x, int iteration)
{
  int sign, E;
  real_t x1, y1, z = 0, M;

  /* atanh(1 - M/2^E) = atanh(T) + Ln(2) * E^2,
     where T = (2 - M - M/2^E) / (2 + M - M/2^E), 0.5 <= M < 1, and E >= 1 integer.
     Thus, if 3/4 <= y/x = 1-M/2^E < 1,
     then we obtain M be repeatedly multiplying 1-y/x by 2
     until we get 0.5 <= 2^E * (1 - y/x) = M < 1.
     To compute atanh(T) we use the given x and y to compute new values.
     x1 <= 1 + M + y/x, and y1 <= 1 - M + y/x. Or
     x1 <= x + y + M*x, and y1 <= x + y - M*x. */
  // ABS(x) < 1
  assert(ABS(y) < ABS(x));

  if ((x > 0 && y < 0) || (x < 0 && y > 0)) sign = -1;
  else sign = +1;

  y = ABS(y), x = ABS(x);
  E = 0, M = 1.0 - y/x;
  while (!(M >= 0.5 && M < 1)) {
    M *= 2.0;
    E++;
  }

  // x1 = 1 + M + y/x, y1 = 1 - M + y/x, z = 0;
  x1 = x + y + M*x, y1 = x + y - M*x, z = 0.0;
  CORDIC_hyperbolic_vectoring(&x1, &y1, &z, iteration);

  return (z + M_LN2 * ((real_t)E/2)) * sign;
}  

real_t CORDIC_ln(real_t v, int iteration)
{
  int E;
  real_t M, x, y, z;

  assert(v >= 0);

  E = 0, M = v;
  while (!(M >= 0.5 && M < 1.0)) {
    if (M >= 1.0) {
      M /= 2.0;
      E++;
    } else {
      M *= 2.0;
      E--;
    }
  }

  x = M + 1, y = M - 1, z = 0.0;
  CORDIC_hyperbolic_vectoring(&x, &y, &z, iteration);

  return 2.0 * z + M_LN2 * E;
}

real_t CORDIC_sqrt(real_t v, int iteration)
{
  int E;
  real_t M, x, y, z = 0;

  assert(v >= 0);

  E = 0, M = v;
  while (!(M >= 0.5 && M < 1.0)) {
    if (M >= 1.0) {
      M /= 2.0;
      E++;
    } else {
      M *= 2.0;
      E--;
    }
  }

  if ((E%2) == 0) {
    x = M + 0.25, y = M - 0.25, z = 0.0;
    CORDIC_hyperbolic_vectoring(&x, &y, &z, iteration);
    if (E > 0) return x * (real_t)(1UL << (E/2));
    return x / (real_t)(1UL << (E/2));
  }

  x = M/2 + 0.25, y = M/2 -0.25, z = 0.0;
  CORDIC_hyperbolic_vectoring(&x, &y, &z, iteration);

  if (E > 0) return x * (real_t)(1UL << ((E+1)/2));
  return x / (real_t)(1UL << ((E+1)/2));
}

#if 1
#include <time.h>
#include <stdlib.h>
#include <inttypes.h>

  // generate a random floating point number from min to max
  static real_t randfrom(real_t vmin, real_t vmax)
  {
    real_t range = (vmax - vmin);
    real_t div = RAND_MAX / range; 

    return vmin + (rand() / div);
  }

  static void CORDIC_linear_test(int n)
  {
    real_t a, b, c;

    srand((unsigned)time(NULL));

    printf("=== Test of CORDIC in linear mode ===\n");

    a = randfrom(INT32_MIN, INT32_MAX);
    b = randfrom(INT32_MIN, INT32_MAX);
    c = CORDIC_multiply(a, b, n);
    printf("CORDIC_multiply in iteration (%d):\n", n);
    printf(" %.8e * %.8e = %.8e goes to %.8e\n", a, b, a*b, c);

    a = randfrom(INT32_MIN, INT32_MAX);
    b = randfrom(INT32_MIN, INT32_MAX);
    c = CORDIC_divide(a, b, n);
    printf("CORDIC_divide in iteration (%d):\n", n);
    printf(" %.8e / %.8e = %.8e goes to %.8e\n", a, b, a/b, c);

    printf("=== End of the test of CORDIC in linear mode ===\n");
  }

  static void CORDIC_circular_test(int n)
  {
    real_t x, y, z;

    srand((unsigned)time(NULL));

    printf("=== Test of CORDIC in circular mode ===\n");

    z = randfrom(-2.0*M_PI, 2.0*M_PI);
    CORDIC_compute_cos_and_sin(&x, &y, z, n);
    printf("CORDIC_compute_cos_and_sin in iteration(%d).\n", n);
    printf(" (cos(%.8e), sin(%.8e) = (%.8e, %.8e) goes to (%.8e, %.8e).\n", z, z, cos(z), sin(z), x, y);

    x = randfrom(INT32_MIN, INT32_MAX);
    y = randfrom(INT32_MIN, INT32_MAX);
    z = CORDIC_atan2(y, x, n);
    printf("CORDIC_atan2 in iteration(%d).\n", n);
    printf(" atan2(%.8e, %.8e) = %.8e goes to %.8e.\n", y, x, atan2(y, x), z);

    x = randfrom(INT32_MIN, INT32_MAX);
    y = randfrom(INT32_MIN, INT32_MAX);
    z = CORDIC_enorm2(x, y, n);
    printf("CORDIC_enorm2 in iteration(%d).\n", n);
    printf(" enorm2(%.8e, %.8e) = %.8e goes to %.8e.\n", x, y, sqrt(x*x + y*y), z);

    printf("=== End of the test of CORDIC in circular mode ===\n");
  }

  static void CORDIC_hyperbolic_test(int n)
  {
    real_t x, y, z;

    srand((unsigned)time(NULL));

    printf("=== Test of CORDIC in hyperbolic mode ===\n");

    z = randfrom(CORDIC_HYPERBOLIC_THETA_MIN - 16 * M_LN2, CORDIC_HYPERBOLIC_THETA_MAX + 16 * M_LN2);
    CORDIC_compute_cosh_and_sinh(&x, &y, z, n);
    printf("CORDIC_compute_cosh_and_sinh in iteration(%d).\n", n);
    printf(" (cosh(%.8e), sinh(%.8e)) = (%.8e, %.8e) goes to (%.8e, %.8e).\n", z, z, cosh(z), sinh(z), x, y);

    x = randfrom(CORDIC_HYPERBOLIC_THETA_MIN - 31 * M_LN2, CORDIC_HYPERBOLIC_THETA_MAX + 31 * M_LN2);
    y = CORDIC_exp(x, n);
    printf("CORDIC_exp in iteration(%d).\n", n);
    printf(" exp(%.8e) = %.8e goes to %.8e.\n", x, exp(x), y);

    x = randfrom(-0.9999, 0.9999);
    y = CORDIC_atanh(x, n);
    printf("CORDIC_atanh in iteration(%d).\n", n);
    printf(" atanh(%.8e) = %.8e goes to %.8e.\n", x, atanh(x), y);

    x = randfrom(INT32_MIN, INT32_MAX);
    y = randfrom(0, ABS(x));
    z = CORDIC_atanh2(y, x, n);
    printf("CORDIC_atanh2 in iteration(%d).\n", n);
    printf(" atanh2(%.8e, %.8e) = %.8e goes to %.8e.\n", y, x, atanh(y/x), z);

    x = randfrom(0, INT32_MAX);
    y = CORDIC_ln(x, n);
    printf("CORDIC_ln in iteration(%d).\n", n);
    printf(" ln(%.8e) = %.8e goes to %.8e.\n", x, log(x), y);

    x = randfrom(0, INT32_MAX);
    y = CORDIC_sqrt(x, n);
    printf("CORDIC_sqrt in iteration(%d).\n", n);
    printf(" sqrt(%.8e) = %.8e goes to %.8e.\n", x, sqrt(x), y);

    printf("=== End of the test of CORDIC in hyperbolic mode ===\n");
  }

  int main(void)
  {
    CORDIC_linear_test(20);
	
    CORDIC_generate_circular_delta_radians(100);
    CORDIC_generate_circular_K_cumulative_products(100);
    CORDIC_circular_test(20);

    CORDIC_generate_hyperbolic_delta_radians(100);
    CORDIC_generate_hyperbolic_K_cumulative_products(100);
    CORDIC_hyperbolic_test(20);
  }
#endif
