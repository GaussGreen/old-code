/*
***************************************************************************
** FILENAME: mdmin.c
**
** Multi-dimensional minimization routine - based on Powell from NumRec
***************************************************************************
*/

#include "mdmin.h"

#include <math.h>
#include <stdio.h>
#include <time.h>

#include "cxmacros.h"

static int linmin
(int              n,       /* (I) */
 double          *p,       /* (I/O) */
 double          *xi,      /* (I/O) */
 CxTMultiDimObjFunc func,    /* (I) */
 void            *data,    /* (I) */
 double          *value);  /* (O) */

/* Context infomation used to pass data to NumRecF1Dim */
typedef struct _CONTEXT
{
    CxTMultiDimObjFunc nrfunc;
    int              ncom;
    double          *pcom;
    double          *xicom;
    double          *xt;
    void            *data;
} CONTEXT;

static int powell
(FILE              *log,
 double            *p,
 double           **xi,
 int                n,
 double             ftol,
 int                maxiter,
 double             maxtime,
 CxTMultiDimObjFunc func,
 void              *data,
 int               *iter,
 double            *fret,
 double            *vdiff);

static int mnbrak
(double      *ax,
 double      *bx,
 double      *cx,
 double      *fa,
 double      *fb,
 double      *fc,
 CONTEXT *context);

static int brent
(double       ax,
 double       bx,
 double       cx,
 CONTEXT *context,
 double       tol,
 double      *xmin,
 double      *fret);

static int f1dim
(double       x,
 CONTEXT *context,
 double      *result);

int CxMultiDimMinStateValidate (CxTMultiDimMinState *state)
{
    static char routine[] = "CxMultiDimMinStateValidate";
    int         status    = FAILURE;

    int i;

    REQUIRE (state != NULL);
    if (state->direction == NULL)
    {
        state->direction = GtoMatrixNewEmpty(state->n, state->n);
        if (state->direction == NULL) goto done; /* failure */
        for (i = 0; i < state->n; ++i)
            state->direction->data[i][i] = 1.0;
    }
    else
    {
        if (state->direction->numDim1 != state->n ||
            state->direction->numDim2 != state->n)
        {
            GtoErrMsg ("%s: direction matrix has wrong dimensions [%d,%d] "
                       "instead of [%d,%d]\n", routine,
                       state->direction->numDim1, state->direction->numDim2, 
                       state->n, state->n);
            goto done; /* failure */
        }
    }

    status = SUCCESS;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}

static void logging (FILE *fp, time_t startTime, int n, double *x, double val,
                     int iter)
{
    if (fp != NULL)
    {
        time_t currentTime = time(NULL);
        int    i;
        char  *sep = "[ ";

        fprintf (fp, "\nIteration = %d. Elapsed time = %f\n",
                 iter, (double)(currentTime - startTime));
        fprintf (fp, "Value = %f\n", val);
        fprintf (fp, "x = ");
        for (i = 0; i < n; ++i)
        {
            fprintf (fp, "%s%f", sep, x[i]);
            sep = ", ";
        }
        fprintf (fp, " ]\n");
        fflush (fp);
    }
}
 
CxTMultiDimMinState* CxMultiDimMinimization
(CxTMultiDimMinState *initState,
 CxTMultiDimObjFunc   func,
 CxTMultiDimObjFunc   dfunc,
 void                *data,
 double               ftol,
 int                  maxiter,
 double               maxtime,
 char                *logfilename)
{
    static char routine[] = "CxMultiDimMinimization";
    int         status    = FAILURE;

    CxTMultiDimMinState *state = NULL;

    FILE       *log = NULL;
    
    if (dfunc != NULL)
    {
        GtoErrMsg ("%s: Methods using derivative function not yet supported\n",
                   routine);
        goto done; /* failure */
    }

    if (logfilename != NULL && *logfilename != '\0')
    {
        log = fopen (logfilename, "w");
        if (log == NULL)
        {
            GtoErrMsg ("%s: Could not open logfilename %s\n", routine, 
                       logfilename);
            goto done; /* failure */
        }
    }

    state = CxMultiDimMinStateCopy (initState);
    if (state == NULL)
        goto done; /* failure */

    /* using the Powell direction set method */
    if (powell (log,
                state->x,
                state->direction->data,
                state->n,
                ftol,
                state->iter + maxiter,
                maxtime,
                func,
                data,
                &state->iter,
                &state->value,
                &state->vdiff) != SUCCESS)
        goto done; /* failure */

    status = SUCCESS;

 done:

    if (log)
    {
        fclose(log);
        log = NULL;
    }

    if (status != SUCCESS)
    {
        CxMultiDimMinStateFree (state);
        state = NULL;
        GtoErrMsgFailure (routine);
    }

    return state;
}


/* Cconstants for the numerical recipes based private functions */
#define SHFT(a,b,c,d)             (a)=(b);(b)=(c);(c)=(d);
#define GOLD                      1.618034
#define CGOLD                     0.3819660
#define SIGN(a,b)                 ((b) > 0.0 ? fabs(a) : -fabs(a))
#define TINY                      1.0e-20
#define GLIMIT                    100.0
#define ZEPS                      1.0e-10

static int powell
(FILE              *log,
 double            *p,
 double           **xi,
 int                n,
 double             ftol,
 int                maxiter,
 double             maxtime,
 CxTMultiDimObjFunc func,
 void              *data,
 int               *iter,
 double            *fret,
 double            *vdiff)
{
    static char routine[] = "powell";
    int         status    = FAILURE;

    int         i;
    int         ibig;
    int         j;
    double      t;
    double      fptt;
    double      fp;
    double      del;
    double     *pt = NULL;
    double     *ptt = NULL;
    double     *xit = NULL;

    time_t      startTime;

    if (maxtime <= 0.0)
        maxtime = 100.0;

    pt=NEW_ARRAY(double,n);
    ptt=NEW_ARRAY(double,n);
    xit=NEW_ARRAY(double,n);

    if (pt == NULL || ptt == NULL || xit == NULL)
        goto done;

    startTime = time(NULL);
    if (log)
    {
        fprintf (log, "Powell algorithm started.\n");
        fflush (log);
    }
    if (func (n, p, data, fret) != SUCCESS)
        goto done;
    logging (log, startTime, n, p, *fret, *iter);

    /* save the initial point */
    COPY_ARRAY (pt, p, double, n);

    while (*iter < maxiter)
    {
        ++(*iter);

        fp = (*fret);
        ibig=0;
        del=0.0;  /* will be the biggest function decrease */
        for (i=0;i<n;++i)
        { /* loop over all directions in the set */
            for (j=0;j<n;j++)
            { /* copy the direction */
                xit[j]=xi[j][i];
            }
            fptt=(*fret);
            if (linmin(n, p, xit, func, data, fret) != SUCCESS)
                goto done; /* failure */
            logging (log, startTime, n, p, *fret, *iter);

            if (fptt-(*fret) > del)
            { /* this is the largest decrease so far */
                del=fptt-(*fret);
                ibig=i;
            }
        }
        *vdiff = fp-(*fret);
        if (2.0*(*vdiff) <= ftol*(fabs(fp)+fabs(*fret))+TINY)
            break;

        for (j=0;j<n;j++)
        {
            /* Construct the extrapolated point and the average direction
               moved. Save the old starting point. */
            ptt[j]=2.0*p[j] - pt[j];
            xit[j]=p[j] - pt[j];
            pt[j]=p[j];
        }

        /* fptt is function value at exptrapolated point */
        if (func (n, ptt, data, &fptt) != SUCCESS)
            goto done;

        if (fptt < fp)
        {
            double sq1 = fp-(*fret)-del;
            double sq2 = fp-fptt;
            t=2.0 * (fp - 2.0*(*fret) + fptt) * sq1 * sq1 - del*sq2*sq2;
            if (t < 0.0)
            {
                /* move the the minimum of the new direction and save the new
                   direction */
                if (linmin(n,p,xit,func,data,fret) != SUCCESS)
                    goto done;
                logging (log, startTime, n, p, *fret, *iter);

                for (j=0;j<n;j++)
                {
                    xi[j][ibig]=xi[j][n-1];
                    xi[j][n-1]=xit[j];
                    /* xi[j][ibig]=xit[j]; */
                }
            }
        }
        if (maxtime > 0.0 && (double)(time(NULL) - startTime) >= maxtime)
            break;
    }

    status = SUCCESS;

 done:

    FREE (xit);
    FREE (ptt);
    FREE (pt);

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}


/*
***************************************************************************
** Given an n dimensional point p[1..n] and an n dimensional direction 
** xi[1..n], moves and resets p to where the function func takes on a 
** minimum along the direction xi from p, and replaces xi by the actual
** vector displacement that p was moved.
** Also returns the value of func at the location p.
** (See Numerical Recipes page 316).
***************************************************************************
*/
static int linmin
(int              n,       /* (I) */
 double          *p,       /* (I/O) */
 double          *xi,      /* (I/O) */
 CxTMultiDimObjFunc func,    /* (I) */
 void            *data,    /* (I) */
 double          *value)   /* (O) */
{
    static char routine[] = "linmin";
    int         status = FAILURE;
    int         j;
    double      xx;
    double      xmin;
    double      fx;
    double      fb;
    double      fa;
    double      bx;
    double      ax;
    double     *xt = NULL;
    CONTEXT context;

    xt = NEW_ARRAY(double,n);
    if (xt == NULL)
        goto done;

    context.nrfunc=func;
    context.ncom=n;
    context.pcom=p;
    context.xicom=xi;
    context.xt=xt;
    context.data=data;

    /*----- Initial guesses -----*/
    ax=0.0;
    xx=1.0; 
    bx=2.0;
    if (mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,&context) != SUCCESS)
        goto done; /* failure */
    
    /*----- Use Brent or Golden search -----*/
    if (brent(ax,xx,bx,&context,2e-4,&xmin,value) != SUCCESS)
        goto done; /* failure */
    
    for (j=0;j<n;++j)
    {
        xi[j]*= xmin;
        p[j] += xi[j];
    }
    status = SUCCESS;

 done:

    FREE(xt);
    if (status != SUCCESS)
        GtoErrMsgFailure(routine);

    return(status);
}


/*
***************************************************************************
** Given a function func, and given distinct initial points ax and bx, this
** routine searches in the downhill direction (defined by the function as
** evaluated at the initial points) and returns new points, ax, bx cx, which
** bracket a minimum of the function. Also returned are the function values
** at the three points, fa, fb and fc. (See Numerical Recipes page 297).
***************************************************************************
*/
static int mnbrak
(double      *ax,
 double      *bx,
 double      *cx,
 double      *fa,
 double      *fb,
 double      *fc,
 CONTEXT *context)
{
    static char routine[] = "mnbrak";
    int         status = FAILURE;
    double      ulim;
    double      u;
    double      r;
    double      q;
    double      fu;
    double      dum;

    if (f1dim(*ax, context, fa) != SUCCESS ||
        f1dim(*bx, context, fb) != SUCCESS)
    {
        goto done;
    }

    if (*fb > *fa)
    {
        SHFT(dum,*ax,*bx,dum)
        SHFT(dum,*fb,*fa,dum)
    }
    *cx=(*bx)+GOLD*(*bx-*ax);
    if (f1dim(*cx, context, fc) != SUCCESS)
    {
        goto done;
    }
    while (*fb > *fc)
    {
        r=(*bx-*ax)*(*fb-*fc);
        q=(*bx-*cx)*(*fb-*fa);
        u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
        ulim=(*bx)+GLIMIT*(*cx-*bx);
        if ((*bx-u)*(u-*cx) > 0.0)
        {
            if (f1dim(u, context, &fu) != SUCCESS)
            {
                goto done;
            }
            if (fu < *fc)
            {
                *ax=(*bx);
                *bx=u;
                *fa=(*fb);
                *fb=fu;
                status = SUCCESS;
                goto done;
            }
            else
            {
                if (fu > *fb)
                {
                    *cx=u;
                    *fc=fu;
                    status = SUCCESS;
                    goto done;
                }
            }
            u=(*cx)+GOLD*(*cx-*bx);
            if (f1dim(u, context, &fu) != SUCCESS)
            {
                goto done;
            }
        }
        else
        {
            if ((*cx-u)*(u-ulim) > 0.0)
            {
                if (f1dim(u, context, &fu) != SUCCESS)
                {
                    goto done;
                }
                if (fu < *fc)
                {
                    double d;
                    SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
                    if (f1dim(u, context, &d) != SUCCESS)
                    {
                        goto done;
                    }
                    SHFT(*fb,*fc,fu,d)
                }
            }
            else
            {
                if ((u-ulim)*(ulim-*cx) >= 0.0)
                {
                    u=ulim;
                    if (f1dim(u, context, &fu) != SUCCESS)
                    {
                        goto done;
                    }
                }
                else
                {
                    u=(*cx)+GOLD*(*cx-*bx);
                    if (f1dim(u, context, &fu) != SUCCESS)
                    {
                        goto done;
                    }
                }
            }
        }
        SHFT(*ax,*bx,*cx,u)
        SHFT(*fa,*fb,*fc,fu)
    }
    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsgFailure(routine);
    }
    return(status);
}

/*
***************************************************************************
** Used by linmin. (See Numerical Recipes page 317).
***************************************************************************
*/
static int f1dim
(double       x,
 CONTEXT *context,
 double      *result)
{
    int     j;
    int     status = FAILURE;

    int     ncom = context->ncom;
    double *pcom = context->pcom;
    double *xicom = context->xicom;
    double *xt = context->xt;

    for (j=0;j<ncom;j++)
    {
        xt[j]=pcom[j]+x*xicom[j];
    }

    status = (context->nrfunc)(ncom, xt, context->data, result);
    
    if (status != SUCCESS)
    {
        GtoErrMsgFailure("f1dim");
    }
    return(status);
}


/*
***************************************************************************
** Given a function, and given a bracketing triplet of abscissas ax, bx, cx
** (such that bx is between ax and cx, and f(bx) is less than both f(ax) and
** f(cx)), this routine isolates the minimum to a fractional precision of about
** tol using Brent's method. The abscissa of the minimum returned as xmin.
** The function returns the value of the minimum funtion value.
***************************************************************************
*/
static int brent
(double       ax,
 double       bx,
 double       cx,
 CONTEXT *context,
 double       tol,
 double      *xmin,
 double      *fret)
{
    static char  routine[] = "brent";
    int          status = FAILURE;
    int          iter;
    double       a;
    double       b;
    double       etemp;
    double       fu;
    double       fv;
    double       fw;
    double       fx;
    double       p;
    double       q;
    double       r;
    double       tol1;
    double       tol2;
    double       u;
    double       v;
    double       w;
    double       x;
    double       xm;
    double       d=0.0;
    double       e=0.0;
    int          brentMaxIter = 50;

    a=((ax < cx) ? ax : cx);
    b=((ax > cx) ? ax : cx);
    x=w=v=bx;
    if (f1dim(x, context, &fx) != SUCCESS)
    {
        goto done;
    }
    fw=fv=fx;
    for (iter=1;iter<=brentMaxIter;iter++)
    {
        xm=0.5*(a+b);
        tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
        if (fabs(x-xm) <= (tol2-0.5*(b-a)))
        {
            *xmin=x;
            *fret = fx;
            status = SUCCESS;
            goto done;
        }
        if (fabs(e) > tol1)
        {
            r=(x-w)*(fx-fv);
            q=(x-v)*(fx-fw);
            p=(x-v)*q-(x-w)*r;
            q=2.0*(q-r);
            if (q > 0.0)
            {
                p = -p;
            }
            q=fabs(q);
            etemp=e;
            e=d;
            if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
            {
                d=CGOLD*(e=(x >= xm ? a-x : b-x));
            }
            else
            {
                d=p/q;
                u=x+d;
                if (u-a < tol2 || b-u < tol2)
                {
                    d=SIGN(tol1,xm-x);
                }
            }
        }
        else
        {
            d=CGOLD*(e=(x >= xm ? a-x : b-x));
        }
        u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
        if (f1dim(u, context, &fu) != SUCCESS)
        {
            goto done;
        }
        if (fu <= fx)
        {
            if (u >= x)
            {
                a=x;
            }
            else
            {
                b=x;
            }
            SHFT(v,w,x,u)
            SHFT(fv,fw,fx,fu)
        }
        else
        {
            if (u < x)
            {
                a=u;
            }
            else
            {
                b=u;
            }
            if (fu <= fw || IS_EQUAL(w,x))
            {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            }
            else
            {
                if (fu <= fv || IS_EQUAL(v,x) || IS_EQUAL(v,w))
                {
                    v=u;
                    fv=fu;
                }
            }
        }
    }

done:
    if (status != SUCCESS)
    {
        GtoErrMsgFailure(routine);
    }
    return(status);
}
