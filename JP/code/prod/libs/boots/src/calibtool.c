/************************************************************************
 * Module:    BOOTS 
 * Name  :  Optimizer
 * Function:Minimizes the objective function
 * Author:    CA from Numrec
 * Revision:    
 ************************************************************************/
#include "drlstd.h"        /* platform compatibility */
#include "drlerr.h"        /* DrlErrMsg */
#include "drllinsm.h"

#include <math.h>
#include <float.h>
#include "boots.h"         /* Objective function definition*/


#define NRANSI
#define TOL 2.0e-4
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double	maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))


#define NR_END 1

static	double *dvector(long nl, long nh)
{
	double *v;
	v= (double *)MALLOC((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) DrlErrMsg("allocation failure in dvector()\n");
	return v-nl+NR_END;
}
static	void free_dvector(double *v, long nl, long nh)
{
	FREE((void*) (v+nl-NR_END));
}



int ncom;
double *pcom,*xicom;



/*f-------------------------------------------------------------
 * Value of the function along the line.
 * 
 */
int
fdim1(
    BootsData *that,  /*(I) BootsData parameters*/
    Target* target,   /*(I) Target correlations input*/
    double x,         /*(I) Point at which to evaluate*/
    double *result    /*(O) Function value*/
)
{
    int j;
    double *xt;
    int errCode = 1;

    xt=dvector(1,ncom);
    for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
    if(BootsObjFunc(that, target, xt, result)!=0){
    goto done;
    }

    errCode = 0;
done:
    free_dvector(xt,1,ncom);
    if (errCode != 0) {
        DrlErrMsg("BootsBuildTimeLineandInterpolate: failed\n");
    }
    return(errCode);
}


/*f-------------------------------------------------------------
 * Given ax and bx starting points, searches in the downhill
 * direction and finds three new points ax, bx, cx
 * 
 */
int
mnbrak1(
    BootsData *that,    /* (I) BootsData parameters*/
    Target* target,     /* (I) Target correlations input*/
    double *ax,         /* (I)(O) Initial and Ouput points*/
    double *bx,         /* (I)(O) Iniital and Ouput points*/
    double *cx,         /* (O) New Output point*/
    double *fa,         /* (O) Value of the function at ouput point A*/
    double *fb,         /* (O) Value of the function at ouput point B*/
    double *fc          /* (O) Value of the function at ouput point C*/
)
{
    double ulim,u,r,q,fu,dum,shift;
    int errCode = 1;

    if((fdim1)(that, target, *ax, fa)!=0){
        goto done;
    }
    if((fdim1)(that, target, *bx, fb)!=0){
        goto done;
    }
    if (*fb > *fa) {
        SHFT(dum,*ax,*bx,dum)
        SHFT(dum,*fb,*fa,dum)
    }
    *cx=(*bx)+GOLD*(*bx-*ax);
    if((fdim1)(that, target, *cx, fc)!=0){
    goto done;
    }
    while (*fb > *fc) {
        r=(*bx-*ax)*(*fb-*fc);
        q=(*bx-*cx)*(*fb-*fa);
        u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
            (2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
        ulim=(*bx)+GLIMIT*(*cx-*bx);
        if ((*bx-u)*(u-*cx) > 0.0) {
            if((fdim1)(that, target, u, &fu)!=0){
            goto done;
            }
            if (fu < *fc) {
                *ax=(*bx);
                *bx=u;
                *fa=(*fb);
                *fb=fu;
                errCode = 0;
                goto done;
            } else if (fu > *fb) {
                *cx=u;
                *fc=fu;
                errCode =0;
                goto done;
            }
            u=(*cx)+GOLD*(*cx-*bx);
            if((fdim1)(that, target, u, &fu)!=0){
            goto done;
            }
        } else if ((*cx-u)*(u-ulim) > 0.0) {
            if((fdim1)(that, target, u, &fu)!=0){
            goto done;
            }
            if (fu < *fc) {
                SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
                if((fdim1)(that, target, u, &shift)!=0){
                goto done;
                }
                SHFT(*fb,*fc,fu,shift)
            }
        } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
            u=ulim;
            if((fdim1)(that, target, u, &fu)!=0){
            goto done;
            }
        } else {
            u=(*cx)+GOLD*(*cx-*bx);
            if((fdim1)(that, target, u, &fu)!=0){
            goto done;
            }
        }
        SHFT(*ax,*bx,*cx,u)
        SHFT(*fa,*fb,*fc,fu)
    }

    errCode = 0;
done:
    if (errCode != 0) {
        DrlErrMsg("mnbrak1: failed\n");
    }
    return(errCode);
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT

#define ITMAX 10000
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

/*f-------------------------------------------------------------
 * Isolates locally a minimum, given three initial points.
 * bx should stand between ax and cx.
 * x result position as xmin.
 */
int
brent1(
       BootsData *that, /*(I) BootsData parameters*/
       Target* target,  /*(I) Target correlations input*/
       double ax,       /*(I) Initial point*/
       double bx,       /*(I) Initial point*/
       double cx,       /*(I) Initial point*/
       double tol,      /*(I) Precision*/
       double *xmin,    /*(O) Position of the minimum*/
       double *result   /*(O) Value of the minimum*/
)
{
    int iter;
    int errCode = 1;
    double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
    double e=0.0;

    a=(ax < cx ? ax : cx);
    b=(ax > cx ? ax : cx);
    x=w=v=bx;
    if((fdim1)(that, target, x, &fw)!=0){
        goto done;
    }
    fx = fv = fw;
    for (iter=1;iter<=ITMAX;iter++) {
        xm=0.5*(a+b);
        tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
        if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
            *xmin=x;
            errCode = 0;
            (*result) = fx;
            goto done;
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
        if((fdim1)(that, target, u, &fu)!=0){
        goto done;
        }
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
    DrlErrMsg("Too many iterations in brent1");
    *xmin=x;
    (*result) = fx;
done:
    if (errCode != 0) {
        DrlErrMsg("brent1: failed\n");
    }
    return(errCode);
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT

/*f-------------------------------------------------------------
 * Gives the minimum from point P in direction xi (ndimensional).
 * 
 * 
 */
int
linmin1(BootsData* that,  /*(I) BootsData parameters*/
        Target* target,   /*(I) Target correlations input*/
        double* p,        /*(I) Initial point*/
        double* xi,       /*(I) Initial direction*/
        int n,            /*(I) Dimension*/
        double *fret      /*(O) Minimum obtained*/
)
{
    int brent1(BootsData *that, 
                Target* target,
                double ax,
                double bx,
                double cx,
                double tol,
                double *xmin,
                double *result
                );

    int mnbrak1(BootsData *that,
                Target* target,
                double *ax,
                double *bx,
                double *cx,
                double *fa,
                double *fb,
                double *fc
                );
    int j;
    double xx,xmin,fx,fb,fa,bx,ax;
    int errCode =1;
    
    ncom=n;
    pcom=dvector(1,n);
    xicom=dvector(1,n);

    for (j=1;j<=n;j++) {
        pcom[j]=p[j];
        xicom[j]=xi[j];
    }
    ax=0.0;
    xx=1.0;
    if(mnbrak1(that, target, &ax, &xx, &bx, &fa, &fx, &fb)!=0){
        goto done;
    }
    if(brent1(that, target, ax, xx, bx, TOL, &xmin, fret)!=0){
        goto done;
    }
    for (j=1;j<=n;j++) {
        xi[j] *= xmin;
        p[j] += xi[j];
    }
    errCode = 0;
done:
    free_dvector(xicom,1,n);
    free_dvector(pcom,1,n);
    if (errCode != 0) {
        DrlErrMsg("linmin1: failed\n");
    }
    return(errCode);
}
#undef TOL
#undef NRANSI

#define ITMAX 12000
#define EPS 1.0e-10
#define FREEALL free_dvector(xi,1,n);free_dvector(h,1,n);free_dvector(g,1,n);


/*f-------------------------------------------------------------
 * Fletcher-Reeves-Polak-Ribiere minimization.
 * Starting point P.
 * Convergence tolerance ftol.
 * Location of the minimum, p as returned quantity.
 * Fret minimum value of the function.
 */
int
frprmn1(BootsData *that,    /*(I) BootsData parameters*/
        Target* target,     /*(I) Target correlations input*/
        double* p,          /*(I) Starting point*/
        int n,              /*(I) Dimension*/
        double  ftol,       /*(I) Tolerance on convergence*/
        int *iter,          /*(I) Maximum number of iterations*/
        double  *fret       /*(O) Minimum value of the function*/
)
{
    int linmin1(BootsData* that,
                Target* target,
                double*  p,
                double*  xi,
                int n,
                double* fret
                );
    int j,its;
    int errCode = 1;
    double  gg,gam,fp,dgg;
    double  *g,*h,*xi;

    g=dvector(1,n);
    h=dvector(1,n);
    xi=dvector(1,n);


    if(BootsObjFunc(that, target, p, &fp)!=0){
        goto done;
    }
    if(BootsdObjFunc(that, target, p, xi)!=0){
        goto done;
    }
    for (j=1;j<=n;j++) {
        g[j] = -xi[j];
        xi[j]=h[j]=g[j];
    }
    for (its=1;its<=ITMAX;its++) {
        *iter=its;
        if(linmin1(that, target, p, xi, n, fret)!=0){
            goto done;    
        }
        if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)){
            errCode = 0;
            goto done;
        }
        if(BootsObjFunc(that, target, p, &fp)!=0){
            goto done;
        }
        if(BootsdObjFunc(that, target, p,xi)!=0){
            goto done;
        }
        dgg=gg=0.0;
        for (j=1;j<=n;j++) {
            gg += g[j]*g[j];
            dgg += (xi[j]+g[j])*xi[j];
        }
        if (gg == 0.0) {
            errCode = 0;
            goto done;
        }
        gam=dgg/gg;
        for (j=1;j<=n;j++) {
            g[j] = -xi[j];
            xi[j]=h[j]=g[j]+gam*h[j];
        }
    }
    DrlErrMsg("Too many iterations in frprmn \n");

done:

    free_dvector(xi,1,n);
    free_dvector(h,1,n);
    free_dvector(g,1,n);
    if (errCode != 0) {
        DrlErrMsg("frprmn1: failed \n");
    }
    return(errCode);
}
#undef ITMAX
#undef EPS
#undef FREEALL
#undef NRANSI









