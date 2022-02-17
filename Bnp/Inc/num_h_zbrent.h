/*--------------------------------------------------------------
	FILE: num_h_zbrent.h
	PURPOSE: Brent method one-dimensional root finding from NR (Chapter 9.3)
	AUTHOR: Dimitri Mayevski
	DATE: 11/07/2002
  --------------------------------------------------------------*/

#ifndef __NUM_H_ZBRENT_H__
#define __NUM_H_ZBRENT_H__

Err num_f_zbrac(Err (*func)(double, double *, void *), double *x1, double *x2,
				double *f1, double *f2, void *static_data);

Err num_f_zbrent(Err (*func)(double, double *, void *), double x1, double x2,
				 double f1, double f2, double tol, void *static_data, double *answer);

#endif  /* #ifndef __NUM_H_ZBRENT_H__ */