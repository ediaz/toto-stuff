/* This file is automatically generated. DO NOT EDIT! */

#ifndef _fft1_h
#define _fft1_h


void fft1_init (int n, float d, float o, bool opt, bool symin);
/*< >*/


void fft1 (float *field, float *Field, sf_file in, bool inv, bool sym, bool opt);
/*< Fast Fourier Transform along the first axis >*/


void fft1_2D (float *field, float *Field, int n2, bool inv);
/*< Fast Fourier Transform along the first axis >*/


void fft1_2D_fwd (float *input, float *output, int n2);
/*< Fast Fourier Transform along the first axis forward>*/


void fft1_2D_inv (float *input, float *output, int n2);
/*< Fast Fourier Transform along the first axis >*/

#endif