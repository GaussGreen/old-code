
/*
*$Log: armquasinewton.cpp,v $
*Revision 1.17  2004/05/12 14:03:06  emezzine
*Added a new parameter to control vega perturbation.
*
*Revision 1.16  2004/02/17 10:56:48  emezzine
*new version for optimiser
*
*Revision 1.15  2003/12/17 09:09:05  emezzine
*correct a bug
*
*Revision 1.14  2003/12/16 13:08:40  emezzine
*Change Rtnewt
*
*Revision 1.13  2003/10/31 08:13:22  emezzine
* New version to calibrate LogProba Curve.
*
*Revision 1.12  2003/10/20 16:56:09  mab
*Manage the name space
*
*Revision 1.11  2003/10/07 13:25:27  emezzine
*Added a new argument for constructor to initilise epsilon to calculate gradient
*
*Revision 1.10  2003/07/25 15:16:52  emezzine
*  Modif of constrector of ARMQuasinewton
*
*Revision 1.9  2003/06/23 09:56:58  emezzine
*rset armquasinewton.cpp
* Use a new fct return a double
*
*Revision 1.8  2003/06/18 12:56:49  emezzine
*Corriger un bug dans Rtnewt
*
*Revision 1.7  2003/06/10 17:20:08  emezzine
*new version
*
*Revision 1.6  2003/03/24 14:44:04  emezzine
*Cut the NAG functions to compile with Unix
*
*Revision 1.5  2003/03/18 10:52:01  emezzine
*Modif de l'algo de Newton Raph si on trouve pas une racine
*
*/
//Created  25/02/03 Ezzine at 10: 25: 30

/*--------------------------------------------------------------------------*/
/*                                                                          */
/* ARM_QUASINEWTON.cpp: implementation of the ARM_QuasiNewtion class.         */
/*                                                                          */
/*--------------------------------------------------------------------------*/

#include "armfrmmodel.h"
#include "calibration.h"
#include "armquasinewton.h"



#ifndef __ARM_NO_NAMESPACE
using namespace std;
#endif

ARM_QuasiNewton::ARM_QuasiNewton(double tol,
				 long maxIter,
				 ARM_Vector* LB,
                 ARM_Vector* UB) 
{
											   
  Init();

  if (LB)
	  itsLowerBound           = (ARM_Vector *) LB->Clone();
  if (UB)
	  itsUpperBound           = (ARM_Vector *) UB->Clone();

  itsPrecision            = tol;
  itsNbMaxIter            = maxIter;
  itsVegaLevel            = 1.0e-10;
}


ARM_QuasiNewton::~ARM_QuasiNewton(void)
{
  if (itsLowerBound)
    delete itsLowerBound;
  
  if (itsUpperBound)
    delete itsUpperBound;
}


void ARM_QuasiNewton::Init(void)
{
  SetName(ARM_QUASINEWTON);
  
  itsLowerBound           = NULL;
  itsUpperBound           = NULL;
  itsPrecision            = 1.0e-3;
  itsNbMaxIter            = 0; 
  itsVegaLevel            = 1.0e-10;
}


void ARM_QuasiNewton::BitwiseCopy(const ARM_Object* ocalib)
{
  ARM_QuasiNewton* calib = (ARM_QuasiNewton *) ocalib;
  
  if (itsLowerBound)
    delete itsLowerBound;
  itsLowerBound = NULL;
  
  if (calib->itsLowerBound)
    itsLowerBound = (ARM_Vector *) calib->itsLowerBound->Clone();
  
  
  if (itsUpperBound)
    delete itsUpperBound;
  itsUpperBound = NULL;
  
  if (calib->itsUpperBound)
    itsUpperBound = (ARM_Vector *) calib->itsUpperBound->Clone();
  
  itsPrecision = calib->itsPrecision;
  itsVegaLevel = calib->itsVegaLevel;
  itsNbMaxIter = calib->itsNbMaxIter;
  
}


void ARM_QuasiNewton::Copy(const ARM_Object* calib)
{
  ARM_Calibration::Copy(calib);
  
  BitwiseCopy(calib);
}


ARM_Object* ARM_QuasiNewton::Clone(void)
{
  ARM_QuasiNewton* theClone = new ARM_QuasiNewton();
  
  theClone->Copy(this);
  
  return(theClone);
}


ARM_QuasiNewton::ARM_QuasiNewton(const ARM_QuasiNewton& calib)
{
  Init();
  
  BitwiseCopy(&calib);
}

//_________________________________________________________________________________________________

double ARM_QuasiNewton::Funcd(double x,ARM_FRMModel* model)
{
  double x0, fx0;
  x0 = x;   

  x0 += itsEpsilon;
  fx0 = model->SimpleFuncToMinimize(x0);
  
  x0 -= (2.0*itsEpsilon);
  fx0 -= model->SimpleFuncToMinimize(x0);
  
  double dfx = fx0/(2.0*itsEpsilon); 

  return dfx;
  
}


double ARM_QuasiNewton::Func(double x,ARM_FRMModel* model )
{
    double fx = model->SimpleFuncToMinimize(x);
    
    return fx;
}

ARM_Vector* ARM_QuasiNewton::Fun2m(ARM_Vector* x,size_t Idx)
{
  size_t nbParam = x->GetSize();
  
  //ARM_Vector* fvec = itsModel->SimpleFuncToMinimize(x,Idx);
  
  return (NULL);//fvec;
}


ARM_Matrix* ARM_QuasiNewton::Fun2mjac(ARM_Vector* x, size_t Idx)
{
  const double EPSILON = 1.0e-8;
  size_t nbParam = x->GetSize();
  double h, temp;
  
  ARM_Matrix* fjac= new ARM_Matrix(nbParam, nbParam);
  
  for (size_t j = 0; j < nbParam ; j++)
    {

	  ARM_Vector* Vectx0 = NULL;//itsModel->SimpleFuncToMinimize(x, Idx);
	  temp = x->Elt(j);

	  h = EPSILON*fabs(temp);
	  if (h == 0.0) h = EPSILON;

      x->Elt(j)+= h;
	  double x1 = x->Elt(j);
      
      ARM_Vector* Vectx1 = NULL;//itsModel->SimpleFuncToMinimize(x, Idx);
	  double vect1 = Vectx1->Elt(0);
	  if(nbParam==2)
		  double vect2 = Vectx1->Elt(1);

	  for (size_t i = 0; i < nbParam; i++)
	  {
		  fjac->Elt(i,j) = (Vectx1->Elt(i) - Vectx0->Elt(i))/h;
	  }
	  
	  double vect5 = fjac->Elt(0,j);//debug
	  if(nbParam==2)
		  double vect6 = fjac->Elt(1,j);
      
      x->Elt(j)= temp;
	      
      delete Vectx1;
      delete Vectx0;
    }
  
  return(fjac);
}

void ARM_QuasiNewton::Zbrak(double x1, double x2,size_t n, 
			    ARM_Vector* xb1,ARM_Vector* xb2, 
			    size_t& nroot, size_t Idx )
{
  size_t i;
  double x, fp, fc, dx;
  
  size_t nb = xb1->GetSize();
  
  nroot = 0;
  
  dx = (x2-x1)/n;
  
  x = x1;
  
  fp = Func(x, NULL);
  
  for (i = 1; i <= n; i++)
    {
      x += dx;
      fc = Func(x, NULL);
      if ((fc*fp) <= 0.0) 
        {
	  xb1->Elt(nroot) = x-dx;
	  xb2->Elt(nroot++) = x;
	  if ( nroot == nb)  
	    return;
        }
      
      fp = fc;
    }
  
  if (!nroot)
    {
      xb1->Elt(nroot) = x1;
      xb2->Elt(nroot++) = x2;
    }
  
}


double ARM_QuasiNewton::BNROptimizer(ARM_FRMModel* model)
{
    size_t j;
    double df,dx,f,rtn;

    rtn = itsInitGuess;

    for (j = 0; j < 200; j++)
    {
        f = Func(rtn,model);
        df = Funcd(rtn,model);	
        double rtn0 = rtn;
        if(fabs(df)<itsVegaLevel)      
            return rtn;
        else
        {
            dx = f/df;
            rtn -= dx;
        }

        if (rtn > 1.25*rtn0)
            rtn = 1.25*rtn0;     

        else if (rtn < 0.75*rtn0)
            rtn = 0.75*rtn0;       

        if ( fabs(dx) < itsPrecision )
        {
            Func(rtn,model);		  
            return rtn;
        }	
    }

    return rtn;
}

ARM_Vector* ARM_QuasiNewton::MultiNewton(ARM_Vector* X0, size_t Idx)
{
  size_t k,i;
  double errx,errf,xsum;//,d;
  
  size_t n = X0->GetSize();
  
  ARM_Vector* X = (ARM_Vector*) X0->Clone();
  ARM_Vector* p    = new ARM_Vector(n);
  ARM_Vector* indx = new ARM_Vector(n);
  
  for (k = 0; k < 500; k++)
    {
      ARM_Vector* fvec = Fun2m(X, Idx);
      errf = 0.0;
      xsum = 0.0;
      for (i = 0; i < n; i++)
		{
		  errf += fabs(fvec->Elt(i));
		  xsum += fabs(X->Elt(i));
		}
		  if (xsum < 1.0e-4)
		  {
              for (i = 0; i < n; i++)
              {
                  X->Elt(i) = 1.0e-03;
              }
			 delete fvec;
			 delete p;
			 delete indx;
			 
			 return X;
		  }
          if ( errf <= 1.0e-5 *itsPrecision )
		  {
              
			 delete fvec;
			 delete p;
			 delete indx;
			 
			 return X;
		  }
      
		  ARM_Matrix* fjac = Fun2mjac(X0,Idx);
		  for (i = 0; i < n; i++)
			{
			  p->Elt(i) = -fvec->Elt(i);
			  
			}

		  double a1 = fjac->Elt(0,0);
      
		  if(n == 2)
			{
			  double b1 = fjac->Elt(0,1);
			  double a2 = fjac->Elt(1,0);
			  double b2 = fjac->Elt(1,1);
			  
			  double Det = a1*b2 - a2*b1;
			  if (fabs(Det) <1.0e-7)
				{
				  delete p;
				  delete fvec;
				  delete fjac;
			  
			      return X;		
				}
			  else
				{
				  p->Elt(0) = (p->Elt(0)*b2 - p->Elt(1)*b1)/Det;
				  p->Elt(1) = (a1*p->Elt(1) - a2*p->Elt(0))/Det;
				}
			}
		  else
			{
			  p->Elt(0) =p->Elt(0)/a1;		
			}
     
		  errx=0.0;
		  for (i=0;i<n;i++)
			{
			  errx += fabs(p->Elt(i));
			  X->Elt(i) += p->Elt(i);
			  
			  if (X->Elt(i) > (1.25)*X0->Elt(i))
				{
				  X->Elt(i) = (1.25)*X0->Elt(i);
				}
			  
			  if (X->Elt(i) < (0.75)*X0->Elt(i))
				{
				  X->Elt(i) = (0.75)*X0->Elt(i);
				}
			  
			}
      
		  if ( errx <= 1.0e-5 *itsPrecision ) 
			{
			  delete p;
			  delete indx;
			  delete fvec;
			  delete fjac;
			  
			  return X;		
			}

		  delete fvec;
		  delete fjac;

		  for (i=0;i<n;i++)
			{
			  X0->Elt(i) = X->Elt(i);
			}
	}
  
  delete p;
  delete indx;
  return X;

}
ARM_Vector* ARM_QuasiNewton::NROptimizer2m(ARM_Vector* X0, size_t Idx)
{
	ARM_Vector* X = MultiNewton(X0,Idx);
	
	Fun2m(X, Idx);

	return X;

}



///___________________________________________________________________________________________________________
void lubksb(ARM_Matrix* a, ARM_Vector* indx, ARM_Vector* b)
{
	int ii=0,ip,j;
	int i;
	double sum;
	int n = a->GetNumLines();

	for (i=0 ; i<n ; i++)
	{
		ip=(int)(indx->Elt(i));

		sum=b->Elt(ip);
		b->Elt(ip) = b->Elt(i);
		if (ii!=0)
		{
			for (j = ii-1; j<i ; j++)
			{
				sum -= a->Elt(i,j)*b->Elt(j);
			}
		}
		else if(sum!=0)
		{
			ii=i+1;
		}
		b->Elt(i) = sum;
	}
	for (i=n-1;i>=1;i--)
	{
		sum = b->Elt(i);

		for (j=i+1;j<n;j++)
		{
			sum -= a->Elt(i,j)*b->Elt(j);
		}
		b->Elt(i) = sum/a->Elt(i,i);
	}
}


size_t ludcmp(ARM_Matrix* a,ARM_Vector*& indx, double& d)
{

	const double TINY = 1.0e-20; // Asmall number
	int i,imax,j,k;
	double big,dum,sum,temp;
	int n = a->GetNumLines();
		
	ARM_Vector* vv = new ARM_Vector(n);

	d=1.0;
	for (i=0;i<n;i++)
	{
		big=0.0;
		for (j=0 ; j<n ; j++)
		{
			if ((temp=fabs(a->Elt(i,j))) > big) big=temp;
		}
		if (big == 0.0) 
			return 0;

		vv->Elt(i) = 1.0/big;
	}
	for (j=0 ; j<n ; j++)
	{
		for (i=0;i<j;i++)
		{
			sum=a->Elt(i,j);
			for (k=0;k<i;k++)
			{
				sum -= a->Elt(i,k)*a->Elt(k,j);
			}
			a->Elt(i,j) = sum;
		}
		big=0.0;
		for (i=j;i<n;i++)
		{
			sum=a->Elt(i,j);
			for (k=1;k<j;k++)
			{
				sum -= a->Elt(i,k)*a->Elt(k,j);
			}
			a->Elt(i,j) = sum;
			if ( (dum=vv->Elt(i)*fabs(sum)) >= big)
			{
				big=dum;
				imax=i;
			}
		}
		if (j != imax)
		{
			for (k=0;k<n;k++)
			{
				dum=a->Elt(imax,k);
				a->Elt(imax,k) = a->Elt(j,k) ;
				a->Elt(j,k) = dum;
			}
			d = -d;
			vv->Elt(imax) = vv->Elt(j);
		}
		indx->Elt(j) = imax;
		if (a->Elt(j,j) == 0.0)
		{
			a->Elt(j,j)=TINY;
		}
		if (j != n-1)
		{
			dum=1.0/(a->Elt(j,j));
			for (i=j+1;i<n;i++)
			{
				a->Elt(i,j) *= dum;
			}
		}
	}
	delete vv;
	return 1;
}

/*!
/* Brent algorithme modified
*/


#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10

void shft(double& a, double& b, double& c, const double d)
{
    a=b,b=c,c=d;
}

template <class T>
inline const T SIGN( const T&  a, const T& b)
{
    return b>=0 ? ( a>=0? a:-a) : (a>=0? -a:a);
}

double ARM_QuasiNewton ::BrentMinimizer(const double ax, const double bx, const double cx,
                                        double tol,double& xmin)
{

    int iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=f(x);
	for (iter=1;iter<=ITMAX;iter++)
    {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a)))
        {
			xmin=x;
			return fx;
		}
		if (fabs(e) > tol1)
        {
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
			else 
            {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		}
        else
        {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=f(u);
		if (fu <= fx)
        {
			if (u >= x) a=x; else b=x;
			shft(v,w,x,u);
			shft(fv,fw,fx,fu);
		}
        else
        {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) 
            {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			}
            else if (fu <= fv || v == x || v == w)
            {
				v=u;
				fv=fu;
			}
		}
	}
    return -99999999.00;

}

#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef NRANSI
