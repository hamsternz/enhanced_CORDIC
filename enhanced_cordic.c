///////////////////////////////////////////////////////////////////////////
// enhanced_cordic.c : An enhanced CORDIC test program
//
// Author: Mike Field <hamster@snap.net.nz>
//
// This is a test for an CORDIC SIN()/COS() implementation that is
// optimized for implementation in FPGA hardware. The differences are:
//
// - That the 'z value (outstanding angle error) is doubled during each
//   iteration of the CORDIC algorithm.
//
// - A block RAM lookup table is used to implement some of the first
//   set of CORDIC iterations. If sufficent block RAM is used, the 
//   values in angles[] can be replaced with a constant, as it quickly
//   approaches a constant value
//
// The input is a phase (not degrees). The outputs are signed SIN 
// and COS values
//
// INDEX_BITS     How many bits are resolved using a lookup table
// CORDIC_BITS    How many bits are resolved using CORDIC
// INPUT_BITS     The total size of the input paramter, in bits
// 
// CORDIC_REPS    How many CORDIC iterations are to be performed
// OUTPUT_SCALE   The positive range of the CORDIC output 
//
// OUTPUT_EXTRA_BITS Scaling factor for the results in progress
// Z_EXTRA_BITS      Scaling vactor for the 'z' (angle yet to be resolved
//
// MAX_ERROR         The limit where the working will be printed out, for
//                   debugging
//
// The benefits of this optimizations are lower latency, lower resource 
// usage, and maybe allow higher Fmax performance 
// 
// Say thanks: This is released under MIT license and can be openly used,
// but feel free to reward me for my efforts with a PayPal donation at the
// above email address if you find this useful.
///////////////////////////////////////////////////////////////////////////
//MIT License
//
//Copyright (c) 2018 Mike Field
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.
///////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include <stdint.h>

/* How the parameter is made broken up */
#define INDEX_BITS     (11)
#define CORDIC_BITS    (19)
#define INPUT_BITS     (2+INDEX_BITS+CORDIC_BITS)

/* Details for the output */
#define CORDIC_REPS       (24)
#define OUTPUT_SCALE      ((int64_t)1<<31)

/* Scaling factor or the results in progress */
#define OUTPUT_EXTRA_BITS (4)

/* Scaling vactor for the z */
#define Z_EXTRA_BITS   (2)

/* Limit where we print out errors */
#define MAX_ERROR  (3.0)

#define PI                (3.14159265358979323846)
#define FULL_CIRCLE       ((int64_t)1<<INPUT_BITS)


/* Masks for extracting parts */
#define QUADRANT_MASK (3 <<(INDEX_BITS+CORDIC_BITS))
#define CORDIC_MASK  ((1<<CORDIC_BITS)-1)
#define INDEX_MASK   (((1<<INDEX_BITS)-1) << CORDIC_BITS)
#define TABLE_SIZE   (1<<INDEX_BITS)
#define TARGET       (1<<(CORDIC_BITS+Z_EXTRA_BITS-1))


int32_t angles[CORDIC_REPS];
int32_t shifts[CORDIC_REPS];
int64_t initial[TABLE_SIZE];
/****************************************************************
 * Calculate the values required for CORDIC sin()/cos() function
 ***************************************************************/
void setup(void) {
   int i, start_shifts;
   double scale = pow(0.5,0.5);
   double table_angle, half_table_angle;
   double cordic_start;
   double table_magnitude;

   table_angle      = PI / 2.0 / TABLE_SIZE;
   half_table_angle = table_angle / 2.0;

   cordic_start     = log(atan(half_table_angle))/log(2.0);
   start_shifts     = ceil(cordic_start);
   printf("Starting CORDIC at lest %13.11f => %i shifts\n", cordic_start, start_shifts);

   scale = 1.0;
   for(i = 0; i < CORDIC_REPS; i++ ) {
     double angle = atan(1.0/pow(2,i-start_shifts));
     angles[i]  = FULL_CIRCLE * angle / (2*PI) * ((int64_t)1<<(Z_EXTRA_BITS+i))+1;
     shifts[i]  = INDEX_BITS+i;
     scale     *= cos(angle);
     printf("angle[%i] = %i\n",i, angles[i]);
   }
   table_magnitude = (OUTPUT_SCALE * scale)*pow(2,OUTPUT_EXTRA_BITS);

   for(i = 0; i < TABLE_SIZE; i++) {
     initial[i] = (int64_t)(table_magnitude * sin(table_angle * i + half_table_angle)-pow(2,OUTPUT_EXTRA_BITS-1));
   }
   if(angles[0] == angles[CORDIC_REPS-1]) {
      printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      printf("!! NOTE = All entries in 'angles' are the same, so a constant can be used     !!!\n");
      printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
   }
}

/***************************************************************
 * Cordic routine to calculate Sine and Cosine for angles
 * with 2^INPUT_BITS representing the full circle
 **************************************************************/
void cordic_sine_cosine(int64_t z, int64_t *s, int64_t *c, int show) {
   int8_t i, flip_sin_sign, flip_cos_sign, quadrant_bit0, quadrant_bit1;
   int32_t index;
   int64_t x, y; 

   /* Split into sections */
   quadrant_bit1 = (z >> (CORDIC_BITS+INDEX_BITS+1)) & 1;
   quadrant_bit0 = (z >> (CORDIC_BITS+INDEX_BITS  )) & 1;
   index         = (z &  INDEX_MASK) >> CORDIC_BITS;
   z             = (z & CORDIC_MASK) << Z_EXTRA_BITS;

   /* Sort out hot to respond to the quadrant we are in */
   flip_sin_sign = quadrant_bit1;
   flip_cos_sign = quadrant_bit1 ^ quadrant_bit0;

   if(quadrant_bit0) 
      z = (1<<(CORDIC_BITS+Z_EXTRA_BITS)) -z; 

   z -= TARGET;

   /* Subtract half the sector angle from Z */
   /* Use Dual Port memory for this */
   if(quadrant_bit0) {
     x = initial[index];
     y = initial[TABLE_SIZE-1-index];
   } else {
     x = initial[TABLE_SIZE-1-index];
     y = initial[index];
   }

   if(show) {
     printf("      SIN        COS        Z\n");
     printf("%10li, %10li, %10li\n", y, x, z);
   }

   for(i = 0; i < CORDIC_REPS; i++ ) {
     int64_t tx = x >> shifts[i];
     int64_t ty = y >> shifts[i];

     x -= (z < 0) ?       -ty :        ty;
     y += (z < 0) ?       -tx :        tx;
     z += (z < 0) ? angles[i] : -angles[i];
     z <<= 1;

     if(show) {
       printf("%10li, %10li, %10li\n", y, x, z);
     }
   }
   *c = (flip_cos_sign  ? -x : x)>>OUTPUT_EXTRA_BITS;
   *s = (flip_sin_sign ? -y : y)>>OUTPUT_EXTRA_BITS;
}

/**************************************************************/
int main(int argc, char *argv[]) {
  int64_t a = 0.8;
  double max = 0.0;
  double total_e = 0.0;
  int64_t count = 0;
  int64_t out_of_range = 0;
  setup();

  if(FULL_CIRCLE > 20000000) {
    printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("!! INPUT_BITS is very large, so this may take a long time to prove all test cases\n");
    printf("!! Please wait........................\n");
    printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  }
  for(a = 0; a < FULL_CIRCLE; a++) {
    int64_t s, c;
    double es,ec;

    cordic_sine_cosine(a, &s, &c, 0);
    es = s-(int64_t)(sin(a*(2*PI/FULL_CIRCLE))*(OUTPUT_SCALE)-0.5);
    ec = c-(int64_t)(cos(a*(2*PI/FULL_CIRCLE))*(OUTPUT_SCALE)-0.5);

    if(es >= MAX_ERROR || es <= -MAX_ERROR || ec >= MAX_ERROR || ec <= -MAX_ERROR) {
      out_of_range++;
      cordic_sine_cosine(a, &s, &c, 1);
      printf("%10li  => %10li, %10li  (error %10f, %10f)\n\n", a, s, c, es, ec);
    }

    if(es > 0) total_e += es;
    else       total_e -= es;
    
    if(ec > 0) total_e += ec;
    else       total_e -= ec;

    if(max < es)  max =  es;
    if(max < -es) max = -es;
    if(max < ec)  max =  ec;
    if(max < -ec) max = -ec;
    count++;
  }
  printf("Error is %13.11f per calcuation out of +/-%li\n",total_e/count, OUTPUT_SCALE);
  printf("Max error is %13.11f, occured %li times\n",max, out_of_range);
  return 0;
}
/**************************************************************/

