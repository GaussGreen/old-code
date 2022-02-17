/************************************************************************
 * Module:	DRL
 * Submodule:	SORT
 * File:	
 * Function:	Root Finding
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdarg.h>

#include "drlsort.h"		/* Prototype consistency */

#define	__DEBUG__
#undef	__DEBUG__

/*----------------------------------------------------------------------
 *
 */

#define	USEB

#if defined(USEA)

#define MAXSTACK (sizeof(size_t) * CHAR_BIT)

static void exchange(void *a, void *b, size_t size) {
    size_t i;

    /******************
     *  exchange a,b  *
     ******************/

    for (i = sizeof(int); i <= size; i += sizeof(int)) {
        int t = *((int *)a);
        *(((int *)a)++) = *((int *)b);
        *(((int *)b)++) = t;
    }
    for (i = i - sizeof(int) + 1; i <= size; i++) {
        char t = *((char *)a);
        *(((char *)a)++) = *((char *)b);
        *(((char *)b)++) = t;
    }
}

DLL_EXPORT(void)
DrlQSort(
	void *base,
	size_t nmemb,
	size_t size,
        int (*compar)(const void *, const void *))
{
    void *lbStack[MAXSTACK], *ubStack[MAXSTACK];
    int sp;
    unsigned int offset;

    /********************
     *  ANSI-C qsort()  *
     ********************/

    lbStack[0] = (char *)base;
    ubStack[0] = (char *)base + (nmemb-1)*size;
    for (sp = 0; sp >= 0; sp--) {
        char *lb, *ub, *m;
        char *P, *i, *j;

        lb = lbStack[sp];
        ub = ubStack[sp];

        while (lb < ub) {

            /* select pivot and exchange with 1st element */
            offset = (ub - lb) >> 1;
            P = lb + offset - offset % size;
            exchange (lb, P, size);

            /* partition into two segments */
            i = lb + size;
            j = ub;
            while (1) {
                while (i < j && compar(lb, i) > 0) i += size;
                while (j >= i && compar(j, lb) > 0) j -= size;
                if (i >= j) break;
                exchange (i, j, size);
                j -= size;
                i += size;
            }

            /* pivot belongs in A[j] */
            exchange (lb, j, size);
            m = j;

            /* keep processing smallest segment, and stack largest */
            if (m - lb <= ub - m) {
                if (m + size < ub) {
                    lbStack[sp] = m + size;
                    ubStack[sp++] = ub;
                }
                ub = m - size;
            } else {
                if (m - size > lb) {
                    lbStack[sp] = lb; 
                    ubStack[sp++] = m - size;
                }
                lb = m + size;
            }
        }
    }
}

#elif defined (USEB)




/* qsort -- qsort interface implemented by faster quicksort.
   J. L. Bentley and M. D. McIlroy, SPE 23 (1993) 1249-1265.
   Copyright 1993, John Wiley.
*/

    /*assume sizeof(long) is a power of 2 */
#define SWAPINIT(a, es) swaptype =         \
    (a-(char*)0 | es) % sizeof(long) ? 2 : es > sizeof(long);
#define swapcode(TYPE, parmi, parmj, n) {  \
    register TYPE *pi = (TYPE *) (parmi);  \
    register TYPE *pj = (TYPE *) (parmj);  \
    do {                                   \
        register TYPE t = *pi;             \
        *pi++ = *pj;                       \
        *pj++ = t;                         \
    } while ((n -= sizeof(TYPE)) > 0);     \
}
#include <stddef.h>
static void swapfunc(char *a, char *b, size_t n, int swaptype)
{   if (swaptype <= 1) swapcode(long, a, b, n)
    else swapcode(char, a, b, n)
}
#define swap(a, b)                         \
    if (swaptype == 0) {                   \
        t = *(long*)(a);                   \
        *(long*)(a) = *(long*)(b);         \
        *(long*)(b) = t;                   \
    } else                                 \
        swapfunc(a, b, es, swaptype)

#define PVINIT(pv, pm)                                \
    if (swaptype != 0) { pv = a; swap(pv, pm); }      \
    else { pv = (char*)&v; *(long*)pv = *(long*)pm; }

#define vecswap(a, b, n) if (n > 0) swapfunc(a, b, n, swaptype)

#define min(x, y) ((x)<=(y) ? (x) : (y))

static char *med3(char *a, char *b, char *c, int (*cmp)())
{	return cmp(a, b) < 0 ?
		  (cmp(b, c) < 0 ? b : cmp(a, c) < 0 ? c : a)
		: (cmp(b, c) > 0 ? b : cmp(a, c) > 0 ? c : a);
}

DLL_EXPORT(void)
DrlQSort(
	void *base,
	size_t n,
	size_t es,
        int (*cmp)(const void *, const void *))
{
	char	*a = (char*) base;
	char *pa, *pb, *pc, *pd, *pl, *pm, *pn, *pv;
	int r, swaptype;
	long t, v;
	size_t s;

	SWAPINIT(a, es);
	if (n < 7) {	 /* Insertion sort on smallest arrays */
		for (pm = a + es; pm < a + n*es; pm += es)
			for (pl = pm; pl > a && cmp(pl-es, pl) > 0; pl -= es)
				swap(pl, pl-es);
		return;
	}
	pm = a + (n/2)*es;    /* Small arrays, middle element */
	if (n > 7) {
		pl = a;
		pn = a + (n-1)*es;
		if (n > 40) {    /* Big arrays, pseudomedian of 9 */
			s = (n/8)*es;
			pl = med3(pl, pl+s, pl+2*s, cmp);
			pm = med3(pm-s, pm, pm+s, cmp);
			pn = med3(pn-2*s, pn-s, pn, cmp);
		}
		pm = med3(pl, pm, pn, cmp); /* Mid-size, med of 3 */
	}
	PVINIT(pv, pm);       /* pv points to partition value */
	pa = pb = a;
	pc = pd = a + (n-1)*es;
	for (;;) {
		while (pb <= pc && (r = cmp(pb, pv)) <= 0) {
			if (r == 0) { swap(pa, pb); pa += es; }
			pb += es;
		}
		while (pb <= pc && (r = cmp(pc, pv)) >= 0) {
			if (r == 0) { swap(pc, pd); pd -= es; }
			pc -= es;
		}
		if (pb > pc) break;
		swap(pb, pc);
		pb += es;
		pc -= es;
	}
	pn = a + n*es;
	s = min(pa-a,  pb-pa   ); vecswap(a,  pb-s, s);
	s = min(pd-pc, pn-pd-es); vecswap(pb, pn-s, s);
	if ((s = pb-pa) > es) DrlQSort(a,    s/es, es, cmp);
	if ((s = pd-pc) > es) DrlQSort(pn-s, s/es, es, cmp);
}




#endif







/*----------------------------------------------------------------------
 * Sorting : revers order in vector.
 * 
 * <br><br>
 * Inverts the order of elements in a vector x_0,...,x_{n-1}
 * to x_{n-1},...,x_0. To be used with <i> DoubleVectSort</i>
 * for descending order.
 */

DLL_EXPORT(int)
DrlDoubleVectRevert(double *x, int n)
{
	register int i, k;
	double	y;

	k = n/2;
	for (i=0; i<=k-1; i++) {
		y = x[i];
		x[i] = x[n-i-1];
		x[n-i-1] = y;
	}
	return(SUCCESS);
}


/*f---------------------------------------------------------------------
 * Sorting : double vector sorting using Quicksort.
 * 
 * <br><br>
 * Sorts a vector x_0,...,x_{n-1} of type <i> double</i>
 * in ascending order.
 */

DLL_EXPORT(int)
DrlDoubleVectSort(double *x, int n)
{
static	char	routine[] = "DrlDoubleVectSort";

#ifdef	_OLDVERSION
#define M 7
#define NSTACK 50
#define FM 7875
#define FA 211
#define FC 1663
	/*
	 * pour que ca commence en 0
	 */
	double	*arr = x - 1 ;

	int	l=1, jstack=0, j, ir, iq, i ;
	int	istack[NSTACK+1] ;
	long	int fx=0L ;
	double	a ;

	ir=n;
	for (;;) {
	    if (ir-l < M) {
		for (j=l+1;j<=ir;j++) {
			a=arr[j];
			for (i=j-1;arr[i]>a && i>0;i--) arr[i+1]=arr[i];
			arr[i+1]=a;
		}
		if (jstack == 0) return(SUCCESS);
		ir=istack[jstack--];
		l=istack[jstack--];
	    } else {
		i=l;
		j=ir;
		fx=(fx*FA+FC) % FM;
		iq= l+((ir-l+1)*fx)/FM;
		a=arr[iq];
		arr[iq]=arr[l];
		for (;;) {
			while (j > 0 && a < arr[j]) j--;
			if (j <= i) {
				arr[i]=a;
				break;
			}
			arr[i++]=arr[j];
			while (a > arr[i] && i <= n) i++;
			if (j <= i) {
				arr[(i=j)]=a;
				break;
			}
			arr[j--]=arr[i];
		}
		if (ir-i >= i-l) {
			istack[++jstack]=i+1;
			istack[++jstack]=ir;
			ir=i-1;
		} else {
			istack[++jstack]=l;
			istack[++jstack]=i-1;
			l=i+1;
		}
		/*NSTACK too small in QCKSRT*/
		if (jstack > NSTACK) {
			GtoErrMsg("%s: NSTACK too small.\n",
				routine);
			return(FAILURE);
		}
	    }
	}

#undef M
#undef NSTACK
#undef FM
#undef FA
#undef FC

#else /*_OLDVERSION*/

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

	unsigned	long i,ir=n,j,k,l=1;
	int		jstack=0,
			istack[NSTACK+3];
	double		*arr, a, temp;

	/* NRC convention */
	arr = x - 1;

	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				for (i=j-1;i>=1;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
				}
				arr[i+1]=a;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(arr[k],arr[l+1])
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
				SWAP(arr[i],arr[j]);
			}
			arr[l]=arr[j];
			arr[j]=a;
			jstack += 2;
			if (jstack > NSTACK) {
				GtoErrMsg("%s: NSTACK too small.\n",
					routine);
				return(FAILURE);
			}
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}

	return(SUCCESS);

#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI

#endif	/*_OLDVERSION*/
}



/*----------------------------------------------------------------------
 * static routine used in next routine.
 */

#define	NMAX	32
static	double	*arrays[NMAX],
		tol;
static	int	nArrays;

static	int
_DoubleVectLexSortCompare(const void *x, const void *y)
{
	double		*xd = (double*)x,
			*yd = (double*)y;
	register int	j, status = 0;


	for (j=0; j<=nArrays-1; j++) {
	    if ((xd[j] + tol)  <= yd[j]) {
		status = -1;
		break;
	    } else if ((xd[j] - tol)  >= yd[j]) {
		status =  1;
		break;
	    }
	}

#if defined(__DEBUG__)
	printf("_DoubleVectLexSortCompare: %lf\n", tol);
	for (j=0; j<=nArrays-1; j++) printf("\t%lf", xd[j]); printf("\n");
	for (j=0; j<=nArrays-1; j++) printf("\t%lf", yd[j]); printf("\n");
	printf("\tstatus=%d\n", status);
#endif

	return(status);
}


/*f---------------------------------------------------------------------
 * Sorting : lexicographic sorting of mutiple arrays.
 * 
 * <br><br>
 * Sorts <i> numArrays</i> of vectors of same length <i> numItems</i>
 * in lexigographic order. If the flags <i> remove\-Double\-Items</i>
 * is TRUE, the removes all multiple occurence of items 
 * whose distance (in max norm) is larger than <i> tolerance</i>.
 * <i> numItems</i> is then changed on exit.
 * Returns SUCCESS/FAILURE.\\
 * <b> Example</b>
 * <tt>
 * double  tolerance = 1e-4; 
 *         *vect1, *vect2, *vect3;
 * ...
I* status = DrlDoubleVectLexSort(
 *              tolerance,
 *              TRUE,
 *              vect1, vect2, vect3, NULL);
 * ...
 * </tt>
 */

DLL_EXPORT(int)
DrlDoubleVectLexSort(
	double tolerance,	/* (I) tolerance */
	int removeDoubleItems,	/* (I) TRUE=remove double items within tol */
	int *numItems,		/* (I) # of elements in each array */
	/* double *array # 0,
	 * ...
	 * double *array # numArrays-1,
	 * NULL (must terminate by NULL!)
	 */
	...)
{
static	char	routine[] = "DrlDoubleVectLexSort";
	int	status = FAILURE;

	va_list		ap;
	double		*x, *y, *base = NULL;
	register int	i, i1, imax = *numItems, j;


	/* if nothing to do */
	imax = *numItems;
	if (imax <= 0) return(SUCCESS);

	tol = (tolerance > DBL_EPSILON ? tolerance : DBL_EPSILON*1e1);

	/* get variable # of arguments (last MUST be NULL!) */
	nArrays = 0;
	va_start(ap, numItems);
	while ((arrays[nArrays] = (double*) va_arg(ap, double*)) != NULL) {
	    if (nArrays >= NMAX) {
		GtoErrMsg("%s: too many arrays.\n", routine);
		goto done;
	    }
	    nArrays++;
	}
	va_end(ap);

	/* copy into sort vector */
	if ((base = NEW_ARRAY(double, nArrays*imax)) == NULL)
		goto done;
	for (i=0; i<=imax-1; i++)
	for (j=0; j<=nArrays-1; j++) {
		base[nArrays*i+j] = arrays[j][i];
	}

#if defined(__DEBUG__)
	printf("%s: INPUT (numItems=%d):\n", routine, *numItems);
	for (i=0; i<=*numItems-1; i++) {
	    double	*x = base + (i    *nArrays);
	    printf("# %2d: ", i);
	    for (j=0; j<=nArrays-1; j++) printf("\t%lf", x[j]);
	    printf("\n");
	    fflush(stdout);
	}
#endif

	/* sort */
	qsort((void*) base, (size_t) imax,
		(size_t) nArrays*sizeof(double), 
		_DoubleVectLexSortCompare);

	/* remove identical elements */
	if ((removeDoubleItems) && (!IS_ALMOST_ZERO(tolerance))) {

	    for(i=0; i<=(*numItems)-2; ) {
		x = base + ( i   *nArrays);
		y = base + ((i+1)*nArrays);

		if (_DoubleVectLexSortCompare((void*)x, (void*)y) == 0) {
		    for (i1=i; i1<=*numItems-2; i1++)
		    for (j=0; j<=nArrays-1; j++) {
			base[nArrays*i1+j] = base[nArrays*(i1+1)+j];
		    }
		    --(*numItems);
		} else {
		    i++;
		}
	    }

#ifdef	_SKIP
	    for (i=0; i<=imax-2; i++) {
		x = base + ( i   *nArrays);
		y = base + ((i+1)*nArrays);
		if (_DoubleVectLexSortCompare((void*)x, (void*)y) == 0) {
		    for (i1=i; i1<=imax-2; i1++)
		    for (j=0; j<=nArrays-1; j++) {
			base[nArrays*i1+j] = base[nArrays*(i1+1)+j];
		    }
		    imax = --(*numItems);
		}
	    }
#endif
	}


#if defined(__DEBUG__)
	for (i=0; i<=imax-2; i++) {
	    double	*x = base + (i    *nArrays), 
	    		*y = base + ((i+1)*nArrays);
	    printf("---> # (%2d, %2d):\n", i, i+1); fflush(stdout);
	    for (j=0; j<=nArrays-1; j++) printf("\t%lf", x[j]); printf("\n");
	    for (j=0; j<=nArrays-1; j++) printf("\t%lf", y[j]); printf("\n");
	    printf("O: %d\n",
		_DoubleVectLexSortCompare((void*)x, (void*)y));
	    fflush(stdout);
	}
#endif


	/* copy back to arrays */
	for (i=0; i<=imax-1; i++)
	for (j=0; j<=nArrays-1; j++) {
		arrays[j][i] = base[nArrays*i+j];
	}

#if defined(__DEBUG__)
	printf("%s: OUTPUT (numItems=%d):\n", routine, *numItems);
	for (i=0; i<=*numItems-1; i++) {
	    printf("# %2d: ", i);
	    for (j=0; j<=nArrays-1; j++) printf("\t%lf", arrays[j][i]);
	    printf("\n");
	    fflush(stdout);
	}
#endif



	/* made it through */
	status = SUCCESS;
done:
	FREE(base);
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}

