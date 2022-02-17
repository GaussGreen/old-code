/* ==========================================================
   FILENAME :   num_f_golden.c

   PURPOSE:		Search min or max of a function via golden
                section
   ========================================================== */

#include        "utallhdr.h"
#include		<math.h"
#include		<stdarg.h"    

/* =========================================================================
	FUNC: golden_section & golden_section_va_list
	DESC: Searching min or max of the function. 
		See "Numerical recipies in C" p.401
	MODIFIES:
	DECLARATION:
   ========================================================================= */

#define R 0.61803399
#define C (1.0-R)
#define SHIFT2(a,b,c) (a)=(b);(b)=(c);
#define SHIFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double golden_section(  double ax, 
						double bx, 
						double cx, 
						double (*function)(double ),
						double tol, 
						double *xmin)
{
	double f1,f2,x0,x1,x2,x3;


	x0=ax; 
	x3=cx;
	if(fabs(cx-bx)>fabs(bx-ax))
	{
		x1=bx; 
		x2=bx+C*(cx-bx);
    }
    else
    {
		x2=bx; 
		x1=bx-C*(bx-ax);
	}
	f1=function(x1);
	f2=function(x2);
	
	while ( fabs(x3-x0) > tol*(fabs(x1)+fabs(x2)) )
    {
		if(f2<f1) 
		{ 
			SHIFT3(x0,x1,x2,R*x1+C*x3); 
	  		SHIFT2(f1,f2,function(x2));
		}
		else
		{ 
			SHIFT3(x3,x2,x1,R*x2+C*x0); 
			SHIFT2(f2,f1,function(x1));
		}
    }
    if(f1<f2)
    {
		*xmin=x1; 
		return f1;
    }
    else
    {
		*xmin=x2; return f2;
    }
}

/* ========================================================================= */

/*<startcom>*/
double golden_section_va_list(  double ax, 
								double bx, 
								double cx, 
								double (*function)(double, va_list),
								double tol, 
								double *xmin,
								...)
/*<endcom>*/
{
	va_list argptr;
	double f1,f2,x0,x1,x2,x3;

	va_start(argptr,xmin);

	x0=ax; 
	x3=cx;
	if(fabs(cx-bx)>fabs(bx-ax))
	{
		x1=bx; 
		x2=bx+C*(cx-bx);
    }
    else
    {
		x2=bx; 
		x1=bx-C*(bx-ax);
	}
	f1=function(x1,argptr);
	f2=function(x2,argptr);
	
	while ( fabs(x3-x0) > tol*(fabs(x1)+fabs(x2)) )
    {
		if(f2<f1) 
		{ 
			SHIFT3(x0,x1,x2,R*x1+C*x3); 
	  		SHIFT2(f1,f2,function(x2,argptr));
		}
		else
		{ 
			SHIFT3(x3,x2,x1,R*x2+C*x0); 
			SHIFT2(f2,f1,function(x1,argptr));
		}
    }
    if(f1<f2)
    {
		*xmin=x1; 
		return f1;
    }
    else
    {
		*xmin=x2; return f2;
    }
}

#undef R 
#undef C 
#undef SHIFT2
#undef SHIFT3
