/* ======================================================
   FILENAME:  num_h_golden.h

   PURPOSE:   Golden section for root finding
   ====================================================== */

#ifndef UTGOLDEN_H
#define UTGOLDEN_H

double golden_section(
    double ax, double bx, double cx, double (*function)(double), double tol, double* xmin);

double golden_section_va_list(
    double ax,
    double bx,
    double cx,
    double (*function)(double, va_list),
    double  tol,
    double* xmin,
    ...);

#endif
