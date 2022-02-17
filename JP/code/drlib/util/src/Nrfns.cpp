//----------------------------------------------------------------------------
//
//      File           : Nrfns.cpp
//
//      Description    : Some Numerical Recipes functions.
//                                              Use arrays with ptr-1.
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#define NRANSI
//#include <math.h>
#include "edginc/Nrfns.hpp"
#include "edginc/Nrutil.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Addin.hpp"
//#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE

// utility method for determining if the exception or
// cause of an exception is an NRException
bool causeIsNR(exception* e)
{
    //is the exception itself an NRException
    NRException* nre = dynamic_cast<NRException*>(e);

    if (nre) {
        return true;
    }
    else {
        //is it even a ModelException
        ModelException* me = dynamic_cast<ModelException*>(e);

        if (me) {
            //is there a sub-cause ?
            ModelException* mec = me->getCause();

            if (mec) {
                //recurse to check the cause
                return causeIsNR(mec);
            }
            else {
                //no cause, from a non-NR exception
                return false;
            }
        }
        else {
            //not even a ModelException
            return false;
        }
    }
}

#define nrerror(text) throw NRException(text)
#include "edginc/Nrutil.hpp"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0 || !iy) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ;
            *idum=IA*(*idum-k*IQ)-IR*k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    if (*idum < 0) *idum += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = *idum;
    if ((temp=AM*iy) > RNMX) return (float) RNMX;
    else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */

float gasdev(long *idum)
{
    static int iset=0;
    static float gset;
    float fac,rsq,v1,v2;

	if ( *idum < 0 ) iset = 0; /* Re initialize */

    if  (iset == 0) {
        do {
            v1=2.0*ran1(idum)-1.0;
            v2=2.0*ran1(idum)-1.0;
            rsq=v1*v1+v2*v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac=sqrt(-2.0*log(rsq)/rsq);
        gset=v1*fac;
        iset=1;
        return (v2*fac);
    } else {
        iset=0;
        return gset;
    }
}
/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */

float gamdev(int ia, long *idum)
/** Returns a deviate distributed as a gamma distribution of integer order ia, i.e., a waiting time
    to the iath event in a Poisson process of unit mean, using ran1(idum) as the source of
    uniform deviates. */
{
    int j;
    float am,e,s,v1,v2,x,y;
    if (ia < 1) nrerror("Error in routine gamdev");
    if (ia < 6) { // Use direct method, adding waiting times. 
        x=1.0;
        for (j=1;j<=ia;j++) x *= ran1(idum);
        x = -log(x);
    } else { // Use rejection method.
        do {
            do {
                do { // These four lines generate the tangent
                     // of a random angle, i.e., they
                     // are equivalent to
                    v1=ran1(idum);
                    v2=2.0*ran1(idum)-1.0;
                } while (v1*v1+v2*v2 > 1.0); // y = tan(ð * ran1(idum)).
                y=v2/v1;
                am=ia-1;
                s=sqrt(2.0*am+1.0);
                x=s*y+am; // We decide whether to reject x:
            } while (x <= 0.0); // Reject in region of zero probability.
            e=(1.0+y*y)*exp(am*log(x/am)-s*y); // Ratio of prob. fn. to comparison fn.
        } while (ran1(idum) > e); // Reject on basis of a second uniform
                                  // deviate. 
    }
    return x;
}

// all floats below are converted to doubles
#define float double
#define fvector dvector
#define matrix dmatrix
#define free_vector free_dvector
#define free_matrix free_dmatrix

float gammln(float xx)
/** Returns the value ln[Ã(xx)] for xx > 0. */
{
    // Internal arithmetic will be done in double precision, a nicety that you can omit if .ve-.gure
    // accuracy is good enough.
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
                          24.01409824083091,-1.231739572450155,
                          0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return (-tmp+log(2.5066282746310005*ser/x));
}

void polint(float xa[], float ya[], int n, float x, float *y, float *dy)
{
    int i,m,ns=1;
    float den,dif,dift,ho,hp,w;
    float *c,*d;

    dif=fabs(x-xa[1]);
    c=fvector(1,n);
    d=fvector(1,n);
    for (i=1;i<=n;i++) {
        if ( (dift=fabs(x-xa[i])) < dif) {
            ns=i;
            dif=dift;
        }
        c[i]=ya[i];
        d[i]=ya[i];
    }
    *y=ya[ns--];
    for (m=1;m<n;m++) {
        for (i=1;i<=n-m;i++) {
            ho=xa[i]-x;
            hp=xa[i+m]-x;
            w=c[i+1]-d[i];
            if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
            den=w/den;
            d[i]=hp*den;
            c[i]=ho*den;
        }
        *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
    }
    free_vector(d,1,n);
    free_vector(c,1,n);
}
/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */

void polcof(float xa[], float ya[], int n, float cof[])
/** Given arrays xa[0..n] and ya[0..n] containing a tabulated function ya i = f(xai), this
routine returns an array of coe.cients cof[0..n] such that ya i =
P
j cofjxa
j
i . */
{
    int k,j,i;
    float xmin,dy,*x,*y;

    x=fvector(0,n);
    y=fvector(0,n);
    for (j=0;j<=n;j++) {
        x[j]=xa[j];
        y[j]=ya[j];
    }
    
    for (j=0;j<=n;j++) {
        polint(x-1,y-1,n+1-j,0.0,&cof[j],&dy);
        /* Subtract 1 from the pointers to x and y because polint uses dimensions [1..n]. We
           extrapolate to x = 0. */
        xmin=1.0e38;
        k = -1;
        for (i=0;i<=n-j;i++) { // Find the remaining xi of smallest absolute value, 
            if (fabs(x[i]) < xmin) {
                xmin=fabs(x[i]);
                k=i;
            }
            if (x[i]) y[i]=(y[i]-cof[j])/x[i]; // (meanwhile reducing all the terms)
        }
        for (i=k+1;i<=n-j;i++) { // and eliminate it.
            y[i-1]=y[i];
            x[i-1]=x[i];
        }
    }
    free_vector(y,0,n);
    free_vector(x,0,n);
}
/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */


#if defined(__GNUC__) && ( __GNUC__ >= 3)
void spline(vector<double>::iterator x,
            vector<double>::iterator y, 
            int                      n, 
            double                   yp1,
            double                   ypn, 
            vector<double>::iterator y2){
    spline(&*x, &*y, n, yp1, ypn, &*y2);
}
#endif
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[])
{
    int i,k;
    float p,qn,sig,un,*u;

    u=fvector(1,n-1);
    if (yp1 > 0.99e30)
        y2[1]=u[1]=0.0;
    else {
        y2[1] = -0.5;
        u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
    }
    for (i=2;i<=n-1;i++) {
        sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
        p=sig*y2[i-1]+2.0;
        y2[i]=(sig-1.0)/p;
        u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
        u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    if (ypn > 0.99e30)
        qn=un=0.0;
    else {
        qn=0.5;
        un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
    }
    y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
    for (k=n-1;k>=1;k--)
        y2[k]=y2[k]*y2[k+1]+u[k];
    free_vector(u,1,n-1);
}
/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */


#if defined(__GNUC__) && ( __GNUC__ >= 3)
void splint(vector<double>::iterator xa, 
            vector<double>::iterator ya, 
            vector<double>::iterator y2a,
            int n, 
            double x,
            double *y){
    splint(&*xa, &*ya, &*y2a, n, x, y);
}
#endif

void splint(float xa[], float ya[], float y2a[], int n, float x, float *y)
{
    int klo,khi,k;
    float h,b,a;

    klo=1;
    khi=n;
    while (khi-klo > 1) {
        k=(khi+klo) >> 1;
        if (xa[k] > x) khi=k;
        else klo=k;
    }
    h=xa[khi]-xa[klo];
    if (h == 0.0) nrerror("Bad xa input to routine splint");
    a=(xa[khi]-x)/h;
    b=(x-xa[klo])/h;
    *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}
/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */

void hpsort(unsigned long n, float ra[])
{
    unsigned long i,ir,j,l;
    float rra;

    if (n < 2) return;
    l=(n >> 1)+1;
    ir=n;
    for (;;) {
        if (l > 1) {
            rra=ra[--l];
        } else {
            rra=ra[ir];
            ra[ir]=ra[1];
            if (--ir == 1) {
                ra[1]=rra;
                break;
            }
        }
        i=l;
        j=l+l;
        while (j <= ir) {
            if (j < ir && ra[j] < ra[j+1]) j++;
            if (rra < ra[j]) {
                ra[i]=ra[j];
                i=j;
                j <<= 1;
            } else j=ir+1;
        }
        ra[i]=rra;
    }
}
/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */

/* Computes forward-di.erence approximation to Jacobian. On input, x[1..n] is the point at
which the Jacobian is to be evaluated, fvec[1..n] is the vector of function values at the
point, and vecfunc(n,x,f) is a user-supplied routine that returns the vector of functions at
x. On output, df[1..n][1..n] is the Jacobian array. */
#define EPS 1.0e-4 // Approximate square root of the machine precision.
void fdjac(int n, float x[], float fvec[], float **df,
           void (*vecfunc)(int, float [], float []))
{
    int i,j;
    float h,temp,*f;
    f=fvector(1,n);
    for (j=1;j<=n;j++) {
        temp=x[j];
        h=EPS*fabs(temp);
        if (h == 0.0) h=EPS;
        x[j]=temp+h; // Trick to reduce finite precision error.
        h=x[j]-temp;
        (*vecfunc)(n,x,f);
        x[j]=temp;
        for (i=1;i<=n;i++) df[i][j]=(f[i]-fvec[i])/h; // Forward difference formula.
    }
    free_vector(f,1,n);
}
#undef EPS
/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */


#define MAXIT 100
float rtsafe(void (*funcd)(float, float *, float *), float x1, float x2, float xacc)
{
    int j;
    float df,dx,dxold,f,fh,fl;
    float temp,xh,xl,rts;

    (*funcd)(x1,&fl,&df);
    (*funcd)(x2,&fh,&df);
    if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
        nrerror("Root must be bracketed in rtsafe");
    if (fl == 0.0) return x1;
    if (fh == 0.0) return x2;
    if (fl < 0.0) {
        xl=x1;
        xh=x2;
    } else {
        xh=x1;
        xl=x2;
    }
    rts=0.5*(x1+x2);
    dxold=fabs(x2-x1);
    dx=dxold;
    (*funcd)(rts,&f,&df);
    for (j=1;j<=MAXIT;j++) {
        if ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)
            || (fabs(2.0*f) > fabs(dxold*df))) {
            dxold=dx;
            dx=0.5*(xh-xl);
            rts=xl+dx;
            if (xl == rts) return rts;
        } else {
            dxold=dx;
            dx=f/df;
            temp=rts;
            rts -= dx;
            if (temp == rts) return rts;
        }
        if (fabs(dx) < xacc) return rts;
        (*funcd)(rts,&f,&df);
        if (f < 0.0)
            xl=rts;
        else
            xh=rts;
    }
    nrerror("Maximum number of iterations exceeded in rtsafe");
    return 0.0;
}
#undef MAXIT
/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */


//////////////////////////////////////////////////////////////////////
// ZBRENT removed from this file
// If required, use the implementation in util/edginc/RootFinder.hpp
//////////////////////////////////////////////////////////////////////
//
// #define ITMAX 100
// #define EPS 3.0e-8
//
// float zbrent(float (*func)(float), float x1, float x2, float tol)
// {
//     int iter;
//     float a=x1,b=x2,c=x2,min1,min2;
//     float d=0;
//     float e=0;
//     float fa=(*func)(a),fb=(*func)(b),fc,p,q,r,s,tol1,xm;

//     if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
//         nrerror("Root must be bracketed in zbrent");
//     fc=fb;
//     for (iter=1;iter<=ITMAX;iter++) {
//         if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
//             c=a;
//             fc=fa;
//             e=d=b-a;
//         }
//         if (fabs(fc) < fabs(fb)) {
//             a=b;
//             b=c;
//             c=a;
//             fa=fb;
//             fb=fc;
//             fc=fa;
//         }
//         tol1=2.0*EPS*fabs(b)+0.5*tol;
//         xm=0.5*(c-b);
//         if (fabs(xm) <= tol1 || fb == 0.0) 
//             return b;
//         if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
//             s=fb/fa;
//             if (a == c) {
//                 p=2.0*xm*s;
//                 q=1.0-s;
//             } else {
//                 q=fa/fc;
//                 r=fb/fc;
//                 p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
//                 q=(q-1.0)*(r-1.0)*(s-1.0);
//             }
//             if (p > 0.0) q = -q;
//             p=fabs(p);
//             min1=3.0*xm*q-fabs(tol1*q);
//             min2=fabs(e*q);
//             if (2.0*p < (min1 < min2 ? min1 : min2)) {
//                 e=d;
//                 d=p/q;
//             } else {
//                 d=xm;
//                 e=d;
//             }
//         } else {
//             d=xm;
//             e=d;
//         }
//         a=b;
//         fa=fb;
//         if (fabs(d) > tol1)
//             b += d;
//         else
//             b += SIGN(tol1,xm);
//         fb=(*func)(b);
//     }
//     nrerror("Maximum number of iterations exceeded in zbrent");
//     return 0.0;
// }
// #undef ITMAX
// #undef EPS

#define ITMAX 100
#define EPS 3.0e-16
/* zbrent is not very useful as it doesn't let you pass any data. Even more useful
   would be if you could specify a target. */
double zbrentUseful(double (*func)(double , void *), void *params, double x1, double x2, double tol)
{
    int iter;
    double a=x1,b=x2,c=x2,min1,min2;
    double d=0;
    double e=0;
    double fa=(*func)(a, params);
    double fb=(*func)(b, params);
    double fc,p,q,r,s,tol1,xm;

    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
        nrerror("Root must be bracketed in zbrent");
    fc=fb;
    for (iter=1;iter<=ITMAX;iter++) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c=a;
            fc=fa;
            e=d=b-a;
        }
        if (fabs(fc) < fabs(fb)) {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        tol1=2.0*EPS*fabs(b)+0.5*tol;
        xm=0.5*(c-b);
        if (fabs(xm) <= tol1 || fb == 0.0) 
            return b;
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            s=fb/fa;
            if (a == c) {
                p=2.0*xm*s;
                q=1.0-s;
            } else {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            if (p > 0.0) q = -q;
            p=fabs(p);
            min1=3.0*xm*q-fabs(tol1*q);
            min2=fabs(e*q);
            if (2.0*p < (min1 < min2 ? min1 : min2)) {
                e=d;
                d=p/q;
            } else {
                d=xm;
                e=d;
            }
        } else {
            d=xm;
            e=d;
        }
        a=b;
        fa=fb;
        if (fabs(d) > tol1)
            b += d;
        else
            b += SIGN(tol1,xm);
        fb=(*func)(b, params);
    }
    nrerror("Maximum number of iterations exceeded in zbrent");
    return 0.0;
}
#undef ITMAX
#undef EPS
/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */


#define ITMAX 100
#define EPS 3.0e-16
/* Copy of zbrent. Identical except that it accepts two extra parameters, 
 * with the value of the evaluating function in the two range points x1 
 * and x2 (which may have been obtained by other means, prior to calling 
 * this function) */
double zbrentUsefulBoundary(double (*func)(double , void *), 
                            void *params, 
                            double x1, 
                            double x2, 
                            double tol,
                            double fx1,
                            double fx2)
{
    int iter;
    double a=x1,b=x2,c=x2,min1,min2;
    double d=0;
    double e=0;
    double fa = fx1;
    double fb = fx2;
    double fc,p,q,r,s,tol1,xm;

    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
        nrerror("Root must be bracketed in zbrent");
    fc=fb;
    for (iter=1;iter<=ITMAX;iter++) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c=a;
            fc=fa;
            e=d=b-a;
        }
        if (fabs(fc) < fabs(fb)) {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        tol1=2.0*EPS*fabs(b)+0.5*tol;
        xm=0.5*(c-b);
        if (fabs(xm) <= tol1 || fb == 0.0) 
            return b;
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            s=fb/fa;
            if (a == c) {
                p=2.0*xm*s;
                q=1.0-s;
            } else {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            if (p > 0.0) q = -q;
            p=fabs(p);
            min1=3.0*xm*q-fabs(tol1*q);
            min2=fabs(e*q);
            if (2.0*p < (min1 < min2 ? min1 : min2)) {
                e=d;
                d=p/q;
            } else {
                d=xm;
                e=d;
            }
        } else {
            d=xm;
            e=d;
        }
        a=b;
        fa=fb;
        if (fabs(d) > tol1)
            b += d;
        else
            b += SIGN(tol1,xm);
        fb=(*func)(b, params);
    }
    nrerror("Maximum number of iterations exceeded in zbrent");
    return 0.0;
}
#undef ITMAX
#undef EPS
/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */


//////////////////////////////////////////////////////////////////////
// ZBRAC removed from this file
// If required, use the implementation in util/edginc/RootFinder.hpp
//////////////////////////////////////////////////////////////////////
//
// #define FACTOR 1.6
// #define NTRY 50
//
// int zbrac(float (*func)(float), float *x1, float *x2)
// {
//     int j;
//     float f1,f2;
//
//     if (*x1 == *x2) nrerror("Bad initial range in zbrac");
//     f1=(*func)(*x1);
//     f2=(*func)(*x2);
//     for (j=1;j<=NTRY;j++) {
//         if (f1*f2 < 0.0) return 1;
//         if (fabs(f1) < fabs(f2))
//             f1=(*func)(*x1 += FACTOR*(*x1-*x2));
//         else
//             f2=(*func)(*x2 += FACTOR*(*x2-*x1));
//     }
//     return 0;
// }
// #undef FACTOR
// #undef NTRY


// JLH - Work in progress. Will migrate to the implementation in 
// util/edginc/RootFinder.hpp and remove this
#define NTRY 40
#define NTRYFUNC 15
#define FACTOR 2
/** Similar to zbrac with a few differences:
 *  - Just like in zbrentUseful, it is convenient to be able to pass extra data
 *    to the evaluating function. 
 *  - Also, change floats for doubles. 
 *  - Accept two extra doubles in order to return by reference the value in the 
 *    two final bracket points. 
 *  - Instead of defining the range being calculated so that f(x) is zero for an
 *    x in (x1, x2), allow x to be in [x1, x2] - ie, if fx1 or fx2 are zero
 *    return immediately. 
 *  - Try to return the smallest possible interval bracketting the zero: remember
 *    the smallest range where fx changes sign.
 *  - Detect NRExceptions thrown during the evaluation of func - if any are 
 *    thrown, do not fail but rather reduce the range under evaluation and           
 *    halve the factor.
 *  - Return values modified: 
 *    + ZBRAC_SUCCESS: success
 *    + ZBRAC_BRAC_ERROR: couldn't bracket the root
 *    + ZBRAC_FUNC_ERROR: error evaluating the function in initial range
 */
ZbracReturn zbracUseful(double (*func)(double, void*),
                        void *params,
                        double *x_1,
                        double *x_2,
                        double *fx_1,
                        double *fx_2)
{
    int j;
    int i;
    double factor = FACTOR;
    double x1  = *x_1;  // For convenience
    double x2  = *x_2;  // For convenience
    double fx1 = *fx_1; // For convenience
    double fx2 = *fx_2; // For convenience

    double lastX1 = x1;
    double lastX2 = x2;
    double lastFx1;
    double lastFx2;

    if (x1 == x2) {
        nrerror("Bad initial range in zbracUseful");
    }
    try {
        lastFx1 = fx1 = (*func)(x1, params);
        lastFx2 = fx2 = (*func)(x2, params);
    }
    catch (exception&) {
        // There was an error evaluating the function in the initial range
        return ZBRAC_FUNC_ERROR;
    }

    for (j=1; j<=NTRY; j++) {
        if (fx1 * fx2 > 0.0) {
            // Both have the same sign, so expand the range

            if (fabs(fx1) < fabs(fx2)) {
                // Shifth x1 - And, as we know the root is NOT between x1 and x2,
                // move lastX2 to x1 too (to narrow the range as much as possible)
                lastX2  = x1;
                lastFx2 = lastFx1;

                for (i=0; i<NTRYFUNC; i++) {
                    // If an NRException is raised when evaluating func, halve
                    // the factor and try again
                    x1 = lastX2 + factor * (lastX2 - x2);

                    try {
                        fx1 = (*func)(x1, params);
                        // No exceptions were raised so got it -> get out of inner loop
                        lastX1  = x1;
                        lastFx1 = fx1;
                        break;
                    }
                    catch (exception& e) {
                        // Determine whether the exception is actually an NRException
                        // (derives from ModelException, which in turn derives from exception)
                        // in which case this tells us that we need to amend the trial & continue
                        // any other type of exception should just cause a failure
                        // Of course this requires that catch blocks between the throw and here
                        // do not change the type 
                        if (causeIsNR(&e))
                        {
                            // Something failed in the evaluation of "func", probably
                            // because x1 was out of boundaries -> Reduce its value
                            factor /= 1.6;
                        }
                        else
                        {
                            //just throw the original exception
                            throw e;
                        }
                    }
                }
            }
            else {
                // Shifth x2 - And, as we know the root is NOT between x1 and x2,
                // move lastXa to x2 too (to narrow the range as much as possible)
                lastX1  = x2;
                lastFx1 = fx2;

                for (i=0; i<NTRYFUNC; i++) {
                    // If an NRException is raised when evaluating func, halve
                    // the factor and try again
                    x2 = lastX1 + factor * (lastX1 - x1);

                    try {
                        fx2 = (*func)(x2, params);
                        // No exceptions were raised so got it -> get out of inner loop
                        lastX2  = x2;
                        lastFx2 = fx2;
                        break;
                    }
                    catch (exception& e) {
                        // Determine whether the exception is actually an NRException
                        // (derives from ModelException, which in turn derives from exception)
                        // in which case this tells us that we need to amend the trial & continue
                        // any other type of exception should just cause a failure
                        // Of course this requires that catch blocks between the throw and here
                        // do not change the type 
                        if (causeIsNR(&e))
                        {
                            // Something failed in the evaluation of "func", probably
                            // because x1 was out of boundaries -> Reduce its value
                            factor /= 1.6;
                        }
                        else
                        {
                            //just throw the original exception
                            throw e;
                        }
                    }
                }
            }
            if (i == NTRYFUNC) {
                // Failed to evaluate the function successfully and
                // therefore to bracket the root
                break; // Get out of the for loop (will return BRAC_ERROR)
            }
        }
        else {
            // Either fx1 and fx2 have different signs or one of them is 0
            // so we are done. Output the narrowest possible range
            *x_1  = lastX1;
            *x_2  = lastX2;
            *fx_1 = lastFx1;
            *fx_2 = lastFx2;
            return ZBRAC_SUCCESS;
        }
    }
    *x_1  = x1;
    *x_2  = x2;
    *fx_1 = fx1;
    *fx_2 = fx2;
    return ZBRAC_BRAC_ERROR;
}
#undef FACTOR
#undef NTRYFUNC
#undef NTRY

/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

float nrselect(unsigned long k, unsigned long n, float arr[])
{
    unsigned long i,ir,j,l,mid;
    float a,temp;

    l=1;
    ir=n;
    for (;;) {
        if (ir <= l+1) {
            if (ir == l+1 && arr[ir] < arr[l]) {
                SWAP(arr[l],arr[ir])
                    }
            return arr[k];
        } else {
            mid=(l+ir) >> 1;
            SWAP(arr[mid],arr[l+1])
                if (arr[l+1] > arr[ir]) {
                    SWAP(arr[l+1],arr[ir])
                        }
            if (arr[l] > arr[ir]) {
                SWAP(arr[l],arr[ir])
                    }
            if (arr[l+1] > arr[l]) {
                SWAP(arr[l+1],arr[l])
                    }
            i=l+1;
            j=ir;
            a=arr[l];
            for (;;) {
                do i++; while (arr[i] < a);
                do j--; while (arr[j] > a);
                if (j < i) break;
                SWAP(arr[i],arr[j])
                    }
            arr[l]=arr[j];
            arr[j]=a;
            if (j >= k) ir=j-1;
            if (j <= k) l=i;
        }
    }
}
#undef SWAP
/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */

#define ALF 1.0e-4 // Ensures su.cient decrease in function value.
#define TOLX_default 1.0e-7 // Convergence criterion on .x.
void lnsrch(int n, float xold[], float fold, float g[], float p[], float x[],
            float *f, float stpmax, int *check, float (*func)(float []),
            float TOLX = TOLX_default)
{
    int i;
    float a,alam,alamin,b,disc,rhs1,rhs2,slope,sum,temp,
        test,tmplam;
    float alam2=0.0;
    float f2=0.0;
    float fold2=0.0;

    *check=0;
    for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
    sum=sqrt(sum);
    if (sum > stpmax)
        for (i=1;i<=n;i++) p[i] *= stpmax/sum;
    for (slope=0.0,i=1;i<=n;i++)
        slope += g[i]*p[i];
    test=0.0;
    for (i=1;i<=n;i++) {
        temp=fabs(p[i])/Maths::max(fabs(xold[i]),1.0);
        if (temp > test) test=temp;
    }
    alamin=TOLX/test;
    alam=1.0;
    for (;;) {
        for (i=1;i<=n;i++) 
            x[i]=xold[i]+alam*p[i];
        *f=(*func)(x);
        if (alam < alamin) {
            for (i=1;i<=n;i++) x[i]=xold[i];
            *check=1;
            return;
        } else if (*f <= fold+ALF*alam*slope) return;
        else {
            if (alam == 1.0)
                tmplam = -slope/(2.0*(*f-fold-slope));
            else {
                rhs1 = *f-fold-alam*slope;
                rhs2=f2-fold2-alam2*slope;
                a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
                b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
                if (a == 0.0) tmplam = -slope/(2.0*b);
                else {
                    disc=b*b-3.0*a*slope;
                    if (disc<0.0) nrerror("Roundoff problem in lnsrch.");
                    else tmplam=(-b+sqrt(disc))/(3.0*a);
                }
                if (tmplam>0.5*alam)
                    tmplam=0.5*alam;
            }
        }
        alam2=alam;
        f2 = *f;
        fold2=fold;
        alam=Maths::max(tmplam,0.1*alam);
    }
}
#undef ALF
#undef TOLX_default
/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */

#define ITMAX 200
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0

#define FREEALL free_vector(xi,1,n);free_vector(pnew,1,n); \
free_matrix(hessin,1,n,1,n);free_vector(hdg,1,n);free_vector(g,1,n); \
free_vector(dg,1,n);
void dfpmin(float p[], int n, float gtol, int *iter, float *fret,
            float(*func)(float []), void (*dfunc)(float [], float []))
{
    int check,i,its,j;
    float den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
    float *dg,*g,*hdg,**hessin,*pnew,*xi;

    dg=fvector(1,n);
    g=fvector(1,n);
    hdg=fvector(1,n);
    hessin=matrix(1,n,1,n);
    pnew=fvector(1,n);
    xi=fvector(1,n);
    fp=(*func)(p);
    (*dfunc)(p,g);
    for (i=1;i<=n;i++) {
        for (j=1;j<=n;j++) hessin[i][j]=0.0;
        hessin[i][i]=1.0;
        xi[i] = -g[i];
        sum += p[i]*p[i];
    }
    stpmax=STPMX*Maths::max(sqrt(sum),(float)n);
    for (its=1;its<=ITMAX;its++) {
        *iter=its;
        lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,&check,func);
        fp = *fret;
        for (i=1;i<=n;i++) {
            xi[i]=pnew[i]-p[i];
            p[i]=pnew[i];
        }
        test=0.0;
        for (i=1;i<=n;i++) {
            temp=fabs(xi[i])/Maths::max(fabs(p[i]),1.0);
            if (temp > test) test=temp;
        }
        if (test < TOLX) {
            FREEALL
                return;
        }
        for (i=1;i<=n;i++) dg[i]=g[i];
        (*dfunc)(p,g);
        test=0.0;
        den=Maths::max(*fret,1.0);
        for (i=1;i<=n;i++) {
            temp=fabs(g[i])*Maths::max(fabs(p[i]),1.0)/den;
            if (temp > test) test=temp;
        }
        if (test < gtol) {
            FREEALL
                return;
        }
        for (i=1;i<=n;i++) dg[i]=g[i]-dg[i];
        for (i=1;i<=n;i++) {
            hdg[i]=0.0;
            for (j=1;j<=n;j++) hdg[i] += hessin[i][j]*dg[j];
        }
        fac=fae=sumdg=sumxi=0.0;
        for (i=1;i<=n;i++) {
            fac += dg[i]*xi[i];
            fae += dg[i]*hdg[i];
            sumdg += Maths::square(dg[i]);
            sumxi += Maths::square(xi[i]);
        }
        if (fac*fac > EPS*sumdg*sumxi) {
            fac=1.0/fac;
            fad=1.0/fae;
            for (i=1;i<=n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
            for (i=1;i<=n;i++) {
                for (j=1;j<=n;j++) {
                    hessin[i][j] += fac*xi[i]*xi[j]
                        -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
                }
            }
        }
        for (i=1;i<=n;i++) {
            xi[i]=0.0;
            for (j=1;j<=n;j++) xi[i] -= hessin[i][j]*g[j];
        }
    }
    nrerror("too many iterations in dfpmin");
    FREEALL
        }
#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX
#undef FREEALL
/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */

#define TOL 2.0e-4

int ncom;
float *pcom,*xicom,(*nrfunc)(float []);

float f1dim(float x) // local function used by linmin
{
    int j;
    float f,*xt;

    xt=fvector(1,ncom);
    for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
    f=(*nrfunc)(xt);
    free_vector(xt,1,ncom);
    return f;
}

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc,
            float (*func)(float))
{
    float ulim,u,r,q,fu,dum;

    *fa=(*func)(*ax);
    *fb=(*func)(*bx);
    if (*fb > *fa) {
        SHFT(dum,*ax,*bx,dum)
            SHFT(dum,*fb,*fa,dum)
            }
    *cx=(*bx)+GOLD*(*bx-*ax);
    *fc=(*func)(*cx);
    while (*fb > *fc) {
        r=(*bx-*ax)*(*fb-*fc);
        q=(*bx-*cx)*(*fb-*fa);
        u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
            (2.0*SIGN(Maths::max(fabs(q-r),TINY),q-r));
        ulim=(*bx)+GLIMIT*(*cx-*bx);
        if ((*bx-u)*(u-*cx) > 0.0) {
            fu=(*func)(u);
            if (fu < *fc) {
                *ax=(*bx);
                *bx=u;
                *fa=(*fb);
                *fb=fu;
                return;
            } else if (fu > *fb) {
                *cx=u;
                *fc=fu;
                return;
            }
            u=(*cx)+GOLD*(*cx-*bx);
            fu=(*func)(u);
        } else if ((*cx-u)*(u-ulim) > 0.0) {
            fu=(*func)(u);
            if (fu < *fc) {
                SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
                    SHFT(*fb,*fc,fu,(*func)(u))
                    }
        } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
            u=ulim;
            fu=(*func)(u);
        } else {
            u=(*cx)+GOLD*(*cx-*bx);
            fu=(*func)(u);
        }
        SHFT(*ax,*bx,*cx,u)
            SHFT(*fa,*fb,*fc,fu)
            }
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

float brent(float ax, float bx, float cx, float (*f)(float), float tol,
            float *xmin)
{
    int iter;
    float a,b,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
    float d=0.0;
    float e=0.0;

    a=(ax < cx ? ax : cx);
    b=(ax > cx ? ax : cx);
    x=w=v=bx;
    fw=fv=fx=(*f)(x);
    for (iter=1;iter<=ITMAX;iter++) {
        xm=0.5*(a+b);
        tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
        if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
            *xmin=x;
            return fx;
        }
        if (fabs(e) > tol1) {
            r=(x-w)*(fx-fv);
            q=(x-v)*(fx-fw);
            p=(x-v)*q-(x-w)*r;
            q=2.0*(q-r);
            if (q > 0.0) p = -p;
            q=fabs(q);
            etemp=e;
            e=d;
            if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
                d=CGOLD*(e=(x >= xm ? a-x : b-x));
            else {
                d=p/q;
                u=x+d;
                if (u-a < tol2 || b-u < tol2)
                    d=SIGN(tol1,xm-x);
            }
        } else {
            d=CGOLD*(e=(x >= xm ? a-x : b-x));
        }
        u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
        fu=(*f)(u);
        if (fu <= fx) {
            if (u >= x) a=x; else b=x;
            SHFT(v,w,x,u)
                SHFT(fv,fw,fx,fu)
                } else {
                    if (u < x) a=u; else b=u;
                    if (fu <= fw || w == x) {
                        v=w;
                        w=u;
                        fv=fw;
                        fw=fu;
                    } else if (fu <= fv || v == x || v == w) {
                        v=u;
                        fv=fu;
                    }
                }
    }
    nrerror("Too many iterations in brent");
    *xmin=x;
    return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */

void linmin(float p[], float xi[], int n, float *fret, float (*func)(float []))
{
    int j;
    float xx,xmin,fx,fb,fa,bx,ax;

    ncom=n;
    pcom=fvector(1,n);
    xicom=fvector(1,n);
    nrfunc=func;
    for (j=1;j<=n;j++) {
        pcom[j]=p[j];
        xicom[j]=xi[j];
    }
    ax=0.0;
    xx=1.0;
    mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
    *fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
    for (j=1;j<=n;j++) {
        xi[j] *= xmin;
        p[j] += xi[j];
    }
    free_vector(xicom,1,n);
    free_vector(pcom,1,n);
}
#undef TOL
/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */

#define ITMAX 200

void powell(float p[], float **xi, int n, float ftol, int *iter, float *fret,
            float (*func)(float []))
{
    int i,ibig,j;
    float del,fp,fptt,t,*pt,*ptt,*xit;

    pt=fvector(1,n);
    ptt=fvector(1,n);
    xit=fvector(1,n);
    *fret=(*func)(p);
    for (j=1;j<=n;j++) pt[j]=p[j];
    for (*iter=1;;++(*iter)) {
        fp=(*fret);
        ibig=0;
        del=0.0;
        for (i=1;i<=n;i++) {
            for (j=1;j<=n;j++) xit[j]=xi[j][i];
            fptt=(*fret);
            linmin(p,xit,n,fret,func);
            if (fabs(fptt-(*fret)) > del) {
                del=fabs(fptt-(*fret));
                ibig=i;
            }
        }
        if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
            free_vector(xit,1,n);
            free_vector(ptt,1,n);
            free_vector(pt,1,n);
            return;
        }
        if (*iter == ITMAX) nrerror("powell exceeding maximum iterations.");
        for (j=1;j<=n;j++) {
            ptt[j]=2.0*p[j]-pt[j];
            xit[j]=p[j]-pt[j];
            pt[j]=p[j];
        }
        fptt=(*func)(ptt);
        if (fptt < fp) {
            t=2.0*(fp-2.0*(*fret)+fptt)*
                Maths::square(fp-(*fret)-del)-
                del*Maths::square(fp-fptt);
            if (t < 0.0) {
                linmin(p,xit,n,fret,func);
                for (j=1;j<=n;j++) {
                    xi[j][ibig]=xi[j][n];
                    xi[j][n]=xit[j];
                }
            }
        }
    }
}
#undef ITMAX

/** Given an array xx[1..n], and given a value x, returns a value j such that x is between xx[j]
and xx[j+1]. xx must be monotonic, either increasing or decreasing. j=0 or j=n is returned
to indicate that x is out of range. */
void locate(double xx[], unsigned long n, double x, unsigned long *j) {
    unsigned long ju,jm,jl;
    int ascnd;
    jl=0; // Initialize lower and upper limits.
    ju=n+1; 
    ascnd=(xx[n] >= xx[1]);
    while (ju-jl > 1) { // If we are not yet done, compute a midpoint,
        jm=(ju+jl) >> 1; 
        if (x >= xx[jm] == ascnd)
            jl=jm; // and replace either the lower limit
        else
            ju=jm; // or the upper limit, as appropriate.
    } // Repeat until the test condition is satised.
    if (x == xx[1]) *j=1; // Then set the output
    else if(x == xx[n]) *j=n-1;
    else *j=jl;
}

/** Given an array xx[1..n], and given a value x, returns a value jlo such that x is between
xx[jlo] and xx[jlo+1]. xx[1..n] must be monotonic, either increasing or decreasing.
jlo=0 or jlo=n is returned to indicate that x is out of range. jlo on input is taken as the
initial guess for jlo on output. */
void hunt(double xx[], unsigned long n, double x, unsigned long *jlo) {
    unsigned long jm,jhi,inc;
    int ascnd;
    ascnd=(xx[n] >= xx[1]); // True if ascending order of table, false otherwise.
    if (*jlo <= 0 || *jlo > n) { // Input guess not useful. Go immediately to bisection
        *jlo=0;
        jhi=n+1;
    } 
    else {
        inc=1; // Set the hunting increment.
        if (x >= xx[*jlo] == ascnd) { // Hunt up:
            if (*jlo == n) return;
            jhi=(*jlo)+1;
            while (x >= xx[jhi] == ascnd) { // Not done hunting,
                *jlo=jhi;
                inc += inc; // so double the increment
                jhi=(*jlo)+inc;
                if (jhi > n) { // Done hunting, since off end of table.
                    jhi=n+1;
                    break;
                } // Try again.
            } //Done hunting, value bracketed.
        }
        else { // Hunt down:
            if (*jlo == 1) {
                *jlo=0;
                return;
            }
            jhi=(*jlo)--;
            while (x < xx[*jlo] == ascnd) { // Not done hunting,
                jhi=(*jlo);
                inc <<= 1; // so double the increment
                if (inc >= jhi) { // Done hunting, since off end of table.
                    *jlo=0;
                    break;
                }
                else *jlo=jhi-inc;
            } // and try again.
        } // Done hunting, value bracketed.
    } // Hunt is done, so begin the final bisection phase:
    while (jhi-(*jlo) != 1) {
        jm=(jhi+(*jlo)) >> 1;
        if (x >= xx[jm] == ascnd)
            *jlo=jm;
        else
            jhi=jm;
    }
    if (x == xx[n]) *jlo=n-1;
    if (x == xx[1]) *jlo=1;
}

//----------------------------------------------------------------------------
#ifndef ROTATE
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
a[k][l]=h+s*(g-h*tau);
#endif

/* Computes all eigenvalues and eigenvectors of a real symmetric fmatrix a[1..n][1..n]. On
output, elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a.
v[1..n][1..n] is a fmatrix whose columns contain, on output, the normalized eigenvectors of
a. nrot returns the number of Jacobi rotations that were required. */
void jacobi(double **a, int n, double d[], double **v, int *nrot)
{
    int j,iq,ip,i;
    double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

    b = fvector(1,n);
    z = fvector(1,n);

    for (ip=1;ip<=n;ip++) { // Initialize to the identity fmatrix. 
        for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
        v[ip][ip]=1.0;
    }
    for (ip=1;ip<=n;ip++) { // Initialize b and d to the diagonal of a. 
        b[ip]=d[ip]=a[ip][ip];
        z[ip]=0.0;
    }
    *nrot=0;
    for (i=1;i<=50;i++) {
        sm=0.0;
        for (ip=1;ip<=n-1;ip++) {   // Sum off-diagonal elements. 
            for (iq=ip+1;iq<=n;iq++)
                sm += fabs(a[ip][iq]);
        }
        if (sm == 0.0) {        // The normal return, which relies on quadratic convergence 
                                // to machine under ow.
            free_vector(z,1,n);
            free_vector(b,1,n);
            // delete [] z;
            // delete [] b;
            return;
        }
        if (i < 4)
            tresh=0.2*sm/(n*n); // ...on the first three sweeps. 
        else
            tresh=0.0;          // ...thereafter. 
        for (ip=1;ip<=n-1;ip++) {
            for (iq=ip+1;iq<=n;iq++) {
                g=100.0*fabs(a[ip][iq]);
                // After four sweeps, skip the rotation if the off-diagonal element is small. 
                if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
                    && (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
                    a[ip][iq]=0.0;
                else if (fabs(a[ip][iq]) > tresh) {
                    h=d[iq]-d[ip];
                    if ((double)(fabs(h)+g) == (double)fabs(h)) {
                        t= (a[ip][iq])/h; 
                    }else {
                        theta=0.5*h/(a[ip][iq]);
                        t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                        if (theta < 0.0) t = -t;
                    }
                    c=1.0/sqrt(1+t*t);
                    s=t*c;
                    tau=s/(1.0+c);
                    h=t*a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq]=0.0;
                    for (j=1;j<=ip-1;j++) {     // Case of rotations 1 . j < p. 
                        ROTATE(a,j,ip,j,iq)
                    }
                    for (j=ip+1;j<=iq-1;j++) {  // Case of rotations p < j < q. 
                        ROTATE(a,ip,j,j,iq)
                    }
                    for (j=iq+1;j<=n;j++) {     // Case of rotations q < j . n. 
                        ROTATE(a,ip,j,iq,j)
                    }
                    for (j=1;j<=n;j++) {
                        ROTATE(v,j,ip,j,iq)
                    }
                    ++(*nrot);
                }
            }
        }
        for (ip=1;ip<=n;ip++) {
            b[ip] += z[ip];
            d[ip]=b[ip];    // Update d with the sum of tapq, 
            z[ip]=0.0;      // and reinitialize z. 
        }
    }
    nrerror("Too many iterations in routine jacobi");
}

/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */

//----------------------------------------------------------------------------
/* Given the eigenvalues d[1..n] and eigenvectors v[1..n][1..n] as output from jacobi
( x 11.1) or tqli ( x 11.3), this routine sorts the eigenvalues into descending order, and rearranges
the columns of v correspondingly. The method is straight insertion. */
void eigsrt(double d[], double **v, int n)
{
    int k,j,i;
    double p;
    for (i=1;i<n;i++) {
        p=d[k=i];
        for (j=i+1;j<=n;j++)
            if (d[j] >= p) p=d[k=j];
        if (k != i) {
            d[k]=d[i];
            d[i]=p;
            for (j=1;j<=n;j++) {
                p=v[j][i];
                v[j][i]=v[j][k];
                v[j][k]=p;
            }
        }
    }
}

//----------------------------------------------------------------------------
#ifndef MIN     /* minimum */
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

/* Computes (a 2 + b 2 ) 1=2 without destructive under ow or over ow. */
double pythag(double a, double b)
{
    double absa,absb;
    absa=fabs(a);
    absb=fabs(b);
    if (absa > absb) return absa*sqrt(1.0+Maths::square(absb/absa));
    else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+Maths::square(absa/absb)));
}

//----------------------------------------------------------------------------
/* Given a matrix a[1..m][1..n], this routine computes its singular value decomposition, A =
U . W . V T . The matrix U replaces a on output. The diagonal matrix of singular values W is out-
put as a vector w[1..n]. The matrix V (not the transpose V T ) is output as v[1..n][1..n]. */
void svdcmp(double **a, int m, int n, double w[], double **v)
{
    double pythag(double a, double b);
    int flag,i,its,j,jj,k,l,nm;
    double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

    rv1=dvector(1,n);
    g=scale=anorm=0.0; // Householder reduction to bidiagonal form.
    for (i=1;i<=n;i++) {
        l=i+1;
        rv1[i]=scale*g;
        g=s=scale=0.0;
        if (i <= m) {
            for (k=i;k<=m;k++) scale += fabs(a[k][i]);
            if (scale) {
                for (k=i;k<=m;k++) {
                    a[k][i] /= scale;
                    s += a[k][i]*a[k][i];
                }
                f=a[i][i];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                a[i][i]=f-g;
                for (j=l;j<=n;j++) {
                    for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
                    f=s/h;
                    for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
                }
                for (k=i;k<=m;k++) a[k][i] *= scale;
            }
        }
        w[i]=scale *g;
        g=s=scale=0.0;
        if (i <= m && i != n) {
            for (k=l;k<=n;k++) scale += fabs(a[i][k]);
            if (scale) {
                for (k=l;k<=n;k++) {
                    a[i][k] /= scale;
                    s += a[i][k]*a[i][k];
                }
                f=a[i][l];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                a[i][l]=f-g;
                for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
                for (j=l;j<=m;j++) {
                    for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
                    for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
                }
                for (k=l;k<=n;k++) a[i][k] *= scale;
            }
        }
        anorm=Maths::max(anorm,(fabs(w[i])+fabs(rv1[i])));
    }
    for (i=n;i>=1;i--) { // Accumulation of right-hand transformations.
        if (i < n) {
            if (g) {
                for (j=l;j<=n;j++) // Double division to avoid possible under ow.
                    v[j][i]=(a[i][j]/a[i][l])/g;
                for (j=l;j<=n;j++) {
                    for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
                    for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
                }
            }
            for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
        }
        v[i][i]=1.0;
        g=rv1[i];
        l=i;
    }
    for (i=MIN(m,n);i>=1;i--) { // Accumulation of left-hand transformations.
        l=i+1;
        g=w[i];
        for (j=l;j<=n;j++) a[i][j]=0.0;
        if (g) {
            g=1.0/g;
            for (j=l;j<=n;j++) {
                for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
                f=(s/a[i][i])*g;
                for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
            }
            for (j=i;j<=m;j++) a[j][i] *= g;
        } else for (j=i;j<=m;j++) a[j][i]=0.0;
        ++a[i][i];
    }
    for (k=n;k>=1;k--) {    // Diagonalization of the bidiagonal form: Loop over
                            //singular values, and over allowed iterations. 
        for (its=1;its<=30;its++) {
            flag=1;
            for (l=k;l>=1;l--) { // Test for splitting.
                nm=l-1; // Note that rv1[1] is always zero.
                if ((double)(fabs(rv1[l])+anorm) == anorm) {
                    flag=0;
                    break;
                }
                if ((double)(fabs(w[nm])+anorm) == anorm) break;
            }
            if (flag) {
                c=0.0; // Cancellation of rv1[l], if l > 1.
                s=1.0;
                for (i=l;i<=k;i++) {
                    f=s*rv1[i];
                    rv1[i]=c*rv1[i];
                    if ((double)(fabs(f)+anorm) == anorm) break;
                    g=w[i];
                    h=pythag(f,g);
                    w[i]=h;
                    h=1.0/h;
                    c=g*h;
                    s = -f*h;
                    for (j=1;j<=m;j++) {
                        y=a[j][nm];
                        z=a[j][i];
                        a[j][nm]=y*c+z*s;
                        a[j][i]=z*c-y*s;
                    }
                }
            }
            z=w[k];
            if (l == k) { // Convergence.
                if (z < 0.0) { // Singular value is made nonnegative.
                    w[k] = -z;
                    for (j=1;j<=n;j++) v[j][k] = -v[j][k];
                }
                break;
            }
            if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
            x=w[l]; // Shift from bottom 2-by-2 minor.
            nm=k-1;
            y=w[nm];
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
            c=s=1.0; // Next QR transformation:
            for (j=l;j<=nm;j++) {
                i=j+1;
                g=rv1[i];
                y=w[i];
                h=s*g;
                g=c*g;
                z=pythag(f,h);
                rv1[j]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g = g*c-x*s;
                h=y*s;
                y *= c;
                for (jj=1;jj<=n;jj++) {
                    x=v[jj][j];
                    z=v[jj][i];
                    v[jj][j]=x*c+z*s;
                    v[jj][i]=z*c-x*s;
                }
                z=pythag(f,h);
                w[j]=z; // Rotation can be arbitrary if z = 0.
                if (z) {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;
                for (jj=1;jj<=m;jj++) {
                    y=a[jj][j];
                    z=a[jj][i];
                    a[jj][j]=y*c+z*s;
                    a[jj][i]=z*c-y*s;
                }
            }
            rv1[l]=0.0;
            rv1[k]=f;
            w[k]=x;
        }
    }
    free_vector(rv1,1,n);
}

/* Given a positive-definite symmetric matrix a[1..n][1..n], this routine constructs its Cholesky
decomposition, A = L . L T . On input, only the upper triangle of a need be given; it is not
modified. The Cholesky factor L is returned in the lower triangle of a, except for its diagonal
elements which are returned in p[1..n]. */
void choldc(double **a, int n, double p[])
{
    int i,j,k;
    double sum;
    for (i=1;i<=n;i++) {
        for (j=i;j<=n;j++) {
            for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
            if (i == j) {
                if (sum <= 0.00000000001) // a, with rounding errors, is not positive definite.
                    // XXX fix me! exactly as was in MonteCarlo
                {
                    // Try to be helpful. NB i==j, and >1 here!
                    throw ModelException("Nrfns::choldc failed at element [" + 
                                         Format::toString(i) + "," + 
                                         Format::toString(i) + 
                                         "].\nTry decreasing the size of elements in the lower triangle [" +
                                         Format::toString(i) + "..1," +
                                         Format::toString(i-1) + "..1] : between " +
                                         Format::toString(a[i-1][i]) + ", " +
                                         Format::toString(a[1][i]) + " and " +
                                         Format::toString(a[1][2]));
                }
                p[i]=sqrt(sum);
            } else a[j][i]=sum/p[i];
        }
    }
}

/* Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise
permutation of itself. a and n are input. a is output, arranged as n equation (2.3.14) above;
indx[1..n] is an output vector that records the row permutation effected by the partial
pivoting; d is output as ± 1 depending on whether the number of row interchanges was even
or odd, respectively. This routine is used in combination with lubksb to solve linear equations
or invert a matrix.*/
#define TINY 1.0e-20

void ludcmp(double **a, int n, int *indx, double *d)
{
    int i,imax,j,k;
    double big,dum,sum,temp;
    double *vv; // vv stores the implicit scaling of each row.
    vv=fvector(1,n);
    *d=1.0; // No row interchanges yet.
    for (i=1;i<=n;i++) { // Loop over rows to get the implicit scaling information. 
        big=0.0;
        for (j=1;j<=n;j++)
            if ((temp=fabs(a[i][j])) > big) big=temp;
        if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
        // No nonzero largest element.
        vv[i]=1.0/big; // Save the scaling.
    }
    for (j=1;j<=n;j++) { // This is the loop over columns of Crout’s method.
        for (i=1;i<j;i++) { // This is equation (2.3.12) except for i = j.
            sum=a[i][j];
            for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
        }
        big=0.0; // Initialize for the search for largest pivot element.
        for (i=j;i<=n;i++) {// This is i = j of equation (2.3.12) and i = j +1... N
                            // of equation (2.3.13). 
            sum=a[i][j];
            for (k=1;k<j;k++) 
                sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
            if ( (dum=vv[i]*fabs(sum)) >= big) {
            // Is the figure of merit for the pivot better than the best so far?
                big=dum;
                imax=i;
            }
        }
        if (j != imax) { // Do we need to interchange rows?
            for (k=1;k<=n;k++) { // Yes, do so...
                dum=a[imax][k];
                a[imax][k]=a[j][k];
                a[j][k]=dum;
            }
            *d = -(*d); // ...and change the parity of d.
            vv[imax]=vv[j]; // Also interchange the scale factor.
        }
        indx[j]=imax;
        if (a[j][j] == 0.0) a[j][j]=TINY;
        // If the pivot element is zero the matrix is singular (at least to the precision of the
        // algorithm). For some applications on singular matrices, it s desirable to substitute
        // TINY for zero.
        if (j != n) { // Now, finally, divide by the pivot element.
            dum=1.0/(a[j][j]);
            for (i=j+1;i<=n;i++) a[i][j] *= dum;
        }
    } // Go back for the next column in the reduction.
    free_vector(vv,1,n);
}

#undef TINY

/* Solves the set of n linear equations A · X = B.Here a[1..n][1..n] is input, not as the matrix
A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input
as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector
B, and returns with the solution vector X. a, n,and indx are not modified by this routine
and can be left in place for successive calls with different right-hand sides b. This routine takes
into account the possibility that b will begin with many zero elements, so it is efficient for use
in matrix inversion.*/
void lubksb(double **a, int n, int *indx, double b[])
{
    int i,ii=0,ip,j;
    double sum;
    for (i=1;i<=n;i++) {// When ii is set to a positive value, it will become the
                        // index of the first nonvanishing element of b.We now
                        // do the forward substitution, equation (2.3.6). The
                        // only new wrinkle is to unscramble the permutation
                        // as we go.
        ip=indx[i];
        sum=b[ip];
        b[ip]=b[i];
        if (ii)
            for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
        else if (sum) ii=i; // A nonzero element was encountered, so from now on we
                            // will have to do the sums in the loop above. 
        b[i]=sum;
    }
    for (i=n;i>=1;i--) { // Now we do the backsubstitution, equation (2.3.7).
        sum=b[i];
        for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
        b[i]=sum/a[i][i]; // Store a component of the solut on vector X.
    } // All done!
}
/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */


#define FREERETURN {free_matrix(fjac,1,n,1,n);free_vector(fvec,1,n);\
free_vector(p,1,n);free_ivector(indx,1,n);return;}
/** Given an initial guess x[1..n] for a root in n dimensions, take ntrial Newton-Raphson steps
    to improve the root. Stop if the root converges in either summed absolute variable increments
    tolx or summed absolute function values tolf */
void mnewt(int ntrial, float x[], int n, float tolx, float tolf,
           int* k, float* errx, float* errf,
           void (*usrfun)(float *x,int n,float *fvec,float **fjac)){
    int i,*indx;
    float d,*fvec,**fjac,*p;

    indx=ivector(1,n);
    p=fvector(1,n);
    fvec=fvector(1,n);
    fjac=matrix(1,n,1,n);

    for (*k=1;*k<=ntrial;(*k)++) {
        // User function supplies function values at x in
        // fvec and Jacobian matrix in fjac.
        usrfun(x,n,fvec,fjac);    
        // Check function convergence.
        *errf=0.0;
        for (i=1;i<=n;i++) *errf += fabs(fvec[i]);
        if (*errf <= tolf) FREERETURN
        for (i=1;i<=n;i++) p[i] = -fvec[i]; // Right-hand side of linear equations.
        ludcmp(fjac,n,indx,&d); // Solve linear equations using LU decomposition.
        lubksb(fjac,n,indx,p);
        // Check root convergence.
        *errx=0.0; 
        for (i=1;i<=n;i++) { // Update solution.
            *errx += fabs(p[i]);
            x[i] += p[i];
        }
        if (*errx <= tolx) FREERETURN
    }
    FREERETURN
}
#undef FREERETURN
/* (C) Copr. 1986-92 Numerical Recipes Software oK.V'. */

int nn; // Global variables to communicate with fmin.
float *fvec;
void (*nrfuncv)(int n, float v[], float f[]);
#define FREERETURN {free_vector(fvec,1,n);free_vector(xold,1,n);\
free_vector(p,1,n);free_vector(g,1,n);free_matrix(fjac,1,n,1,n);\
free_ivector(indx,1,n);return;}
#define FMAX(a,b) Maths::max((a),(b))
#define SQR(a) Maths::square((a))

float fmin(float x[])
/* Returns f = 1/2 F · F at x. The global pointer *nrfuncv points to a routine that returns the
fvector of functions at x. It is set to point to a user-supplied routine in the calling program.
Global variables also communicate the function values back to the calling program. */
{
    int i;
    float sum;
    (*nrfuncv)(nn,x,fvec);
    for (sum=0.0,i=1;i<=nn;i++) sum += SQR(fvec[i]);
    return 0.5*sum;
}

void newt(float x[], int n,
          int *its, double *errf, double *errx, int *check,
          void (*vecfunc)(int, float [], float []),
          void (*jac)(int, float [], float [], float **,
                      void (*)(int, float [], float [])),
          int MAXITS, float TOLF, float TOLMIN, float TOLX, float STPMX)
/* Given an initial guess x[1..n] for a root in n dimensions, .nd the root by a globally convergent
Newton’s method. The fvector of functions to be zeroed, called fvec[1..n] in the routine
below, is returned by the user-supplied routine vecfunc(n,x,fvec). The output quantity
check is false (0) on a normal return and true (1) if the routine has converged to a local
minimum of the function fmin de.ned below. In this case try restarting from a di.erent initial
guess. */
{
    int i,j,*indx;
    float d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold;
    indx=ivector(1,n);
    fjac=matrix(1,n,1,n);
    g=fvector(1,n);
    p=fvector(1,n);
    xold=fvector(1,n);
    fvec=fvector(1,n); // De.ne global variables.
    nn=n;
    nrfuncv=vecfunc;
    f=fmin(x); // fvec is also computed by this call.
    (*errf)=0.0; // Test for initial guess being a root. Use
    // more stringent test than simply TOLF. 
    for (i=1;i<=n;i++)
        if (fabs(fvec[i]) > (*errf)) (*errf)=fabs(fvec[i]);
    if ((*errf) < 0.01*TOLF) {
        *check=0;
        FREERETURN
    }
    for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]); // Calculate stpmax for line searches.
    stpmax=STPMX*FMAX(sqrt(sum),(float)n);
    for ((*its)=1;(*its)<=MAXITS;(*its)++) { // Start of iteration loop.
        // If analytic Jacobian is available, you can replace the routine jac below with your
        // own routine.
        jac(n,x,fvec,fjac,vecfunc);
        for (i=1;i<=n;i++) { // Compute .f for the line search.
            for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*fvec[j];
            g[i]=sum;
        }
        for (i=1;i<=n;i++) xold[i]=x[i]; // Store x,
        fold=f; // and f.
        for (i=1;i<=n;i++) p[i] = -fvec[i]; // Right-hand side for linear equations.
        ludcmp(fjac,n,indx,&d); // Solve linear equations by LU decomposition.
        lubksb(fjac,n,indx,p);
        lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fmin, /* additional parameter here */ TOLX);
        // lnsrch returns new x and f. It also calculates fvec at the new x when it calls fmin.
        (*errf)=0.0; // Test for convergence on function values.
        for (i=1;i<=n;i++)
            if (fabs(fvec[i]) > (*errf)) (*errf)=fabs(fvec[i]);
        if ((*errf) < TOLF) {
            *check=0;
            FREERETURN
        }
        if (*check) { // Check for gradient of f zero, i.e., spurious
            // convergence. 
            test=0.0;
            den=FMAX(f,0.5*n);
            for (i=1;i<=n;i++) {
                temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
                if (temp > test) test=temp;
            }
            *check=(test < TOLMIN ? 1 : 0);
            FREERETURN
        }
        (*errx)=0.0; // Test for convergence on äx.
        for (i=1;i<=n;i++) {
            temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
            if (temp > (*errx)) (*errx)=temp;
        }
        if ((*errx) < TOLX) FREERETURN
    }
    nrerror("MAXITS exceeded in newt");
}
#undef FREERETURN
#undef FMAX

#define BIG 1.0e30

void mprove(float **a, float **alud, int n, int indx[], float b[], float x[]) {
/** Improves a solution vector x[1..n] of the linear set of equations A · X = B. The matrix
    a[1..n][1..n], and the vectors b[1..n] and x[1..n] are input, as is the dimension n.
    Also input is alud[1..n][1..n], the LU decomposition of a as returned by ludcmp, and
    the vector indx[1..n] also returned by that routine. On output, only x[1..n] is modified,
    to an improved set of values. */
    int j,i;
    double sdp;
    float *r;
    r=fvector(1,n);
    for (i=1;i<=n;i++) {    // Calculate the right-hand side, accumulating
        sdp = -b[i];        //the residual in double precision. 
        for (j=1;j<=n;j++) sdp += a[i][j]*x[j];
        r[i]=sdp;
    }
    lubksb(alud,n,indx,r);  // Solve for the error term,
    for (i=1;i<=n;i++) x[i] -= r[i]; // and subtract it from the old solution.
    free_vector(r,1,n);
}

/** Pade approximation */
void pade(double cof[], int n, double *resid) {
/* Given cof[0..2*n], the leading terms in the power series expansion of a function, solve the
    linear Pad´e equations to return the coe.cients of a diagonal rational function approximation to
    the same function, namely (cof[0] + cof[1]x+· · ·+ cof[n]xN)/(1 + cof[n+1]x +· · ·+
    cof[2*n]xN). The value resid is the norm of the residual vector; a small value indicates a
    well-converged solution. Note that cof is double precision for consistency with ratval. */

    int j,k,*indx;
    float d,rr,rrold,sum,**q,**qlu,*x,*y,*z;
    
    indx=ivector(1,n);
    q=matrix(1,n,1,n);
    qlu=matrix(1,n,1,n);
    x=fvector(1,n);
    y=fvector(1,n);
    z=fvector(1,n);
    for (j=1;j<=n;j++) { //Set up matrix for solving.
        y[j]=x[j]=cof[n+j];
        for (k=1;k<=n;k++) {
            q[j][k]=cof[j-k+n];
            qlu[j][k]=q[j][k];
        }
    }
    ludcmp(qlu,n,indx,&d); //Solve by LU decomposition and backsubstitution.
    lubksb(qlu,n,indx,x);
    rr=BIG;
    do {            //Important to use iterative improvement, since
        rrold=rr;   //the Pad´e equations tend to be ill-conditioned. 
        for (j=1;j<=n;j++) z[j]=x[j];
        mprove(q,qlu,n,indx,y,x);
        for (rr=0.0,j=1;j<=n;j++) //Calculate residual.
            rr += SQR(z[j]-x[j]);
    } while (rr < rrold);   // If it is no longer improving, call it quits.
    *resid=sqrt(rrold);
    for (k=1;k<=n;k++) { // Calculate the remaining coe.cients.
        for (sum=cof[k],j=1;j<=k;j++) sum -= z[j]*cof[k-j];
        y[k]=sum;
    } // Copy answers to output.
    for (j=1;j<=n;j++) {
        cof[j]=y[j];
        cof[j+n] = -z[j];
    }
    free_vector(z,1,n);
    free_vector(y,1,n);
    free_vector(x,1,n);
    free_matrix(qlu,1,n,1,n);
    free_matrix(q,1,n,1,n);
    free_ivector(indx,1,n);
}

double ratval(double x, double cof[], int mm, int kk) {
/** Given mm, kk, and cof[0..mm+kk], evaluate and return the rational function (cof[0] +
    cof[1]x + · · · + cof[mm]xmm)/(1 + cof[mm+1]x + · · · + cof[mm+kk]xkk). */
    int j;
    double sumd,sumn;   // Note precision! Change to float if desired.
    for (sumn=cof[mm],j=mm-1;j>=0;j--) sumn=sumn*x+cof[j];
    for (sumd=0.0,j=mm+kk;j>=mm+1;j--) sumd=(sumd+cof[j])*x;
    return (sumn/(1.0+sumd));
}

#undef SQR

/** Computes Pade coefficients from an array of Taylor coefficients.
    The number of Taylor coefficients must be odd and the first element
    must be the 0 order expansion i.e. f(0) */
class PadeCoeffsAddin: public CObject {
public:
    static CClassConstSP const TYPE;
    DoubleArray coeffs;

    static IObjectSP padeCoefficients(PadeCoeffsAddin* params){
        static const string routine = "PadeCoeffsAddin::padeCoefficients";
        try {
            DoubleArray& coeffs = params->coeffs;
            
            int size = coeffs.size();
            if(size % 2 == 0) {
                throw ModelException(routine, "Number of coefficients must be odd.");
            }
            int nParams = (size - 1) / 2;
            float residual = 0.0;
            
            pade(&coeffs[0], nParams, &residual);            
            
            return IObjectSP(coeffs.clone());

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }


    /** for reflection */
    PadeCoeffsAddin(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PadeCoeffsAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultPadeCoeffsAddin);
        FIELD(coeffs, "Taylor coefficients");

        Addin::registerClassObjectMethod("PADE_COEFFICIENTS",
                                         Addin::UTILITIES,
                                         "Compute Pade approximation coefficients",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)padeCoefficients);

    }

    static IObject* defaultPadeCoeffsAddin(){
        return new PadeCoeffsAddin();
    }
};

CClassConstSP const PadeCoeffsAddin::TYPE = CClass::registerClassLoadMethod(
"PadeCoeffsAddin", typeid(PadeCoeffsAddin), PadeCoeffsAddin::load);

bool PadeCoeffsAddinLinkIn() { 
    return (PadeCoeffsAddin::TYPE != 0);
}


/** Evaluates the rational function based on the provided
    Pade coefficients at an array of values. */
class PadeAddin: public CObject {
public:
    static CClassConstSP const TYPE;
    DoubleArray padeCoeffs;
    DoubleArray xValues;

    static IObjectSP padeValue(PadeAddin* params){
        static const string routine = "PadeAddin::padeValue";
        try {
            DoubleArray& coeffs = params->padeCoeffs;
            if(coeffs.size() % 2 == 0) {
                throw ModelException(routine, "Number of coefficients must be odd.");
            }
            
            int nParams = (coeffs.size() - 1) / 2;
            // int nParamsDen = (coeffs.size() - 1) / 2;
            // int nParamsNum = coeffs.size() - nParamsDen;
            
            DoubleArray& xValues = params->xValues;
            int n = xValues.size();
            DoubleArray yValues(n);
            for(int i = 0; i < n; i++) {
                yValues[i] = ratval(xValues[i], &coeffs[0], nParams, nParams);
            }
            return IObjectSP(yValues.clone());

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }


    /** for reflection */
    PadeAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PadeAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultPadeAddin);
        FIELD(padeCoeffs, "Pade coefficients");
        FIELD(xValues, "Arrays of values at which to evaluate function");

        Addin::registerClassObjectMethod("PADE_VALUE",
                                         Addin::UTILITIES,
                                         "Compute Pade approximation coefficients",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)padeValue);

    }

    static IObject* defaultPadeAddin(){
        return new PadeAddin();
    }
};

CClassConstSP const PadeAddin::TYPE = CClass::registerClassLoadMethod(
"PadeAddin", typeid(PadeAddin), PadeAddin::load);

bool PadeAddinLinkIn() { 
    return (PadeAddin::TYPE != 0);
}

#undef nrerror
DRLIB_END_NAMESPACE
