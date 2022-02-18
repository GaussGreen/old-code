/*  More NR utils for use with SDP. */

//#include "math.h"
#include "utallhdr.h"

/*  Balances a matrix. */
void balanc(double** a, int n);

/*  Transform into Hessenberg. */
void elmhes(double** a, int n);

/*  QR for Hessenberg matrixes. */
Err hqr(double** a, int n, double wr[], double wi[]);