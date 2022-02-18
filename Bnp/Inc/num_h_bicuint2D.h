/*************************************************
 * Declaration for bicubic interpolation functions
 * Copied from NR book
 * C. Godart
 * Date: 24/01/2000
 **************************************************/

void bcucof(double y[], double y1[], double y2[], double y12[], double d1, double d2, double** c);

Err bcuint(
    double  y[],
    double  y1[],
    double  y2[],
    double  y12[],
    double  x1l,
    double  x1u,
    double  x2l,
    double  x2u,
    double  x1,
    double  x2,
    double* ansy,
    double* ansy1,
    double* ansy2);

Err grid_fd_12(
    double   x1a[],
    double   x2a[],
    double** ya,
    int      m,
    int      n,
    double** y1a,
    double** y2a,
    double** y12a);
