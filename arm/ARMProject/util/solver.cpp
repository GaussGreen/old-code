/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: solver.cpp,v $
 * Revision 6 & 7 (continuus) 2004/06/17 09:50:00  jmprie
 * Ajout d'un throw si pente trop faible dans solver NewtonRaphson
 *
 * Revision 1.23  2004/05/24 15:28:41  jmprie
 * dans la recherche initiale oubli du ++iter. Ajout du throw si
 * iter depasse alors iter_max
 *
 * Revision 1.22  2004/05/14 09:43:13  emezzine
 * Modif and added random researsh.
 *
 * Revision 1.21  2004/03/02 14:33:42  emezzine
 * Improvment and devoloppment
 *
 * Revision 1.20  2004/01/27 12:46:34  emezzine
 * improvment and developpment
 *
 * Revision 1.19  2004/01/22 17:28:22  emezzine
 * improvment
 *
 * Revision 1.18  2004/01/21 12:05:05  emezzine
 * added bounadries into solver
 *
 * Revision 1.13  2004/01/07 07:41:23  jmprie
 * ajout du zero Brent solver
 *
 * Revision 1.12  2003/12/16 14:21:45  emezzine
 *  change Solve() if d_f = 0
 *
 * Revision 1.11  2003/12/09 17:27:01  emezzine
 * improvement
 *
 * Revision 1.10  2003/12/09 17:21:00  emezzine
 * added an exeption
 *
 * Revision 1.9  2003/11/18 16:06:58  emezzine
 * improvment
 *
 * Revision 1.8  2003/11/18 13:35:20  emezzine
 * Cut absolute bound in ModifiedNRSolver:: Solve()
 *
 * Revision 1.5  2003/10/17 14:03:02  jmprie
 * chgt de fonctionnement du NewtonRaphsonSolver pour pouvoir modifier
 * le contexte de la f(x) sans toucher au solveur
 *
 * Revision 1.4  2003/10/13 07:56:01  jmprie
 * ajout pour Newton, du SetInitialGuess()
 *
 * Revision 1.3  2003/08/18 11:19:34  ebenhamou
 * explicit message
 *
 * Revision 1.2  2003/08/18 11:16:16  ebenhamou
 * exception is sent in the proper format
 *
 */

           
#include "rand-gen.h"
#include "solver.h"

DichotomySolver::DichotomySolver( const RealValuedFunction& f_, 
	double minorant_,
	double majorant_,
	double b_,
	double precision_, 
	long max_iters_ )
:	f(f_),minorant(minorant_),majorant(majorant_),b(b_),precision(precision_),max_iters(max_iters_)
{}


double DichotomySolver::Solve()
{
	double dx,fmin,fmid,xmid;
	long i;
	
	fmin=f(minorant)-b;
	fmid=f(majorant)-b;	
	if (fmin*fmid>0.0)
	{
		char msg[100];
		sprintf( msg,"Could not bracket any solution between %f and %f in DichotomySolver::Solve",minorant,majorant);
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,msg );
	}	

	double result = fmin < 0.0 ? (dx=majorant-minorant,minorant) : (dx=minorant-majorant,majorant);
	for (i=1;i<=max_iters;i++) 
	{
		fmid=f(xmid=result+(dx *= 0.5))-b;
		if (fmid <= 0.0) result=xmid;
		if (fabs(dx) < precision || fmid == 0.0) return result;
	}

	char msg[100];
	sprintf( msg, "Not enough iterations to find root in <DichotomySolver::Solve>, max iteration: %ld", max_iters);
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		msg );
}




NewtonRaphsonSolver::NewtonRaphsonSolver(const FunctionWithKnownDerivative* f_,
	                                    double initial_guess_,
	                                    double b_,
	                                    double precision_, 
	                                    long max_iters_)
                                        :	f(f_),fprime(f_->Derivative()), 
                                            initial_guess(initial_guess_),
                                            b(b_),
                                            precision(precision_),
	max_iters(max_iters_)
{}

NewtonRaphsonSolver::~NewtonRaphsonSolver()
{}

void NewtonRaphsonSolver::setInitialGuess(double guess)
{
    initial_guess=guess;
}

double NewtonRaphsonSolver::Solve(){
  double df_x,dx,f_x,result;
  
  result=initial_guess;
  for (long i=1;i<=max_iters;i++)
  {
    f_x=(*f)(result)-b;
    df_x=(*fprime)(result);

    if(fabs(df_x)<K_DOUBLE_TOL)
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "Slope to low to find root in <NewtonRaphsonSolver>");

    dx=f_x/df_x;
    result -= dx;
    if (fabs(dx) < precision) 
    {
        f_x=(*f)(result);
        return result;
    }
  }

  throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "Not enough iterations (max_iters) to find root in <NewtonRaphsonSolver>");


  return 0.0;
}


ModifiedNRSolver::ModifiedNRSolver(const FunctionWithKnownDerivative* f_,
	                                    double initial_guess_,
	                                    double b_,
	                                    double precision_, 
	                                    long max_iters_)
                                        :	itsFunction(f_),itsFunctionDerivative(f_->Derivative()), 
                                            initial_guess(initial_guess_),
                                            b(b_),
                                            precision(precision_),
	max_iters(max_iters_)
{}

ModifiedNRSolver::~ModifiedNRSolver()
{}
void ModifiedNRSolver::SetInitialParams(double guess,
                                        double lowerbound,
                                        double upperbound)
{
    initial_guess=guess;
    lower_bound = lowerbound;
    upper_bound = upperbound;

}

double ModifiedNRSolver::Solve()
{
  double df_x,dx,f_x;  
  double result=initial_guess;
  long  FirstUse = 1;
  double unif = MRGK5Random(FirstUse);  
  double x0=result;

  int iter = 0.0;
  for (long i=1;i<=max_iters;i++)
  {
      f_x=(*itsFunction)(result)-b;
      df_x=(*itsFunctionDerivative)(result);

      double rtn0 = result;
      if(fabs(df_x) < K_FRM_TOL)
      {
          iter++;
          result = lower_bound + unif* upper_bound;
          FirstUse = 0;
          unif = MRGK5Random(FirstUse);
          if (iter == 5) 
              return x0;
      }
      else
      {
          dx=f_x/df_x;
          result -= dx;
          iter = 0;
          x0 = result;

          if (fabs(dx) < precision) 
          {
              f_x=(*itsFunction)(result);
              return result;
          }
      } 
      
      if (result > 1.25*rtn0)
          result = 1.25*rtn0;     

      else if (result < 0.75*rtn0)
          result = 0.75*rtn0;  
  }
  
  return result;
  
}

  



BrentSolver::BrentSolver( const RealValuedFunction& f_, 
		                    double minGuess_,
		                    double maxGuess_,
	                        double target_,
	                        double precision_, 
	                        long max_iters_)
                            :	f(f_),
                                minGuess(minGuess_),
                                maxGuess(maxGuess_),
                                target(target_),
                                precision(precision_),
                                max_iters(max_iters_)
{}

BrentSolver::~BrentSolver()
{}

double BrentSolver::Solve()
{
	int iter;
	double a=minGuess,b=maxGuess,c,d,e,min1,min2;
	double fa=f(a)-target,fb=f(b)-target,fc,p,q,r,s,tol1,xm;

    iter=0;
    double df;
    while( ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) && iter<max_iters)
    {
        /// Locate the solution
        df=(fb-fa)/(b-a);
        if( (df>K_NEW_DOUBLE_TOL && fa<0) || (df<-K_NEW_DOUBLE_TOL && fa>0) )
        {
            /// Shift interval towards positive values
            c=b;
            fc=fb;
            b += (b-a);
            fb=f(b)-target;
            a=c;
            fa=fc;

        }
        else if( (df>K_NEW_DOUBLE_TOL && fa>0) || (df<-K_NEW_DOUBLE_TOL && fa<0) )
        {
            /// Shift interval towards negative values
            c=a;
            fc=fa;
            a -= (b-a);
            fa=f(a)-target;
            b=c;
            fb=fc;
        }
        else
        {        
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Slope to low to do the bracketing in <BrentSolver>");
        }
        ++iter;
    }

    if(iter>=max_iters)
    {        
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "Can'y locate a valid initial solution in <BrentSolver>");
    }

    c=b;
    fc=fb;
	for (iter=1;iter<=max_iters;iter++)
    {
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
		tol1=2.0*3.0e-8*fabs(b)+0.5*precision;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
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
            b += (xm>0 ? tol1 : -tol1);
		fb=f(b)-target;
	}

    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "Not enough iterations (max_iters) to find root in <BrentSolver>");


  return 0.0;
}

