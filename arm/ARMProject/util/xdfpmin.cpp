/* Driver for routine dfpmin */

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "linalg.h" 
#include "xdfpmin.h" 
#include "lsq.h"





long dfpmin(ARM_Vector* p,
			ARM_Matrix* hessin,
			int *iter,
			double *fret,
			void **parameters,
			ARM_Matrix* data,
            T_FUNC func)
{
  
    int check=0,i,its,j, n , calhess=0;
    double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test,test10, rien;

    n=p->GetSize();

    ARM_Vector dg(n,0.0),g(n,0.0),hdg(n,0.0),pnew(n,0.0), xi(n,0.0);
	
    fp=(*func)(data, p, parameters);
    
    NumJac(p,fp,&g,parameters,data,func);
/*	NumHess(p,hessin,parameters,data,func);
	hessin->Invert(hessin,rien);

    
    for (i=0;i<n;i++) 
    {
	    xi.Elt(i)=0.0;
     
        for (j=0;j<n;j++) 
        {
                xi.Elt(i) -= hessin->Elt(i,j)*g.Elt(j);
        }
	}*/

    for (i = 0; i < n; i++) 
    {
        xi.Elt(i) = -g.Elt(i);
        sum += p->Elt(i)*p->Elt(i);
    }
    
    stpmax=STPMX*FMAX(sqrt(sum),(double) n );
    
	for (its = 1;its <= ITMAX; its++) 
    {
        *iter=its;
    
										printf("* DFP *%d\n",*iter);
								//		printf ("  %6.2f %6.2f\n",xi.Elt(0), xi.Elt(1));
								//		printf("* DFP g*");
								//		printf ("  %6.2f %6.2f\n",g.Elt(0), g.Elt(1));
        
		lnsrch(p,fp,&g,&xi,&pnew,fret,stpmax,&check,
										parameters,data,func);
		
								//      printf("* DFP pnew*");
								//		printf ("  %6.10f %6.10f\n",pnew.Elt(0), pnew.Elt(1));
								//		char sleep=getchar();
		fp = *fret;
        
        for (i=0;i<n;i++) 
        {
            xi.Elt(i)=pnew.Elt(i)-p->Elt(i);

            p->Elt(i)=pnew.Elt(i);
        }

		
        test10=0.0;
        
        for (i=0;i<n;i++) 
        {
            temp=fabs(xi.Elt(i))/FMAX(fabs(p->Elt(i)),1.0);
        
            if ( temp > test10 )
            {
               test10=temp;
            }
        }
        
		if ( test10 < TOLX*1.0e-3 ) 
        {
			NumHess(p,hessin,parameters,data,func);

			return(5);
        }
		
		
		if (((*iter+10) % 20) == 0) 
		{
			NumHess(p,hessin,parameters,data,func);
			hessin->Invert(hessin,rien);
			
		/*	for (i=0;i<n;i++) 
			{
				for (j=0;j<n;j++) 
				{
					hessin->Elt(i,j)=-hessin->Elt(i,j);
				}
			}	*/
			calhess++;
		}

        
        
        
        for (i=0;i<n;i++) 
        {
            dg.Elt(i)=g.Elt(i);
        }

        NumJac(p,fp,&g,parameters,data,func);

        test=0.0;
        den=FMAX(*fret,1.0);
        
        for (i=0;i<n;i++) 
        {
            temp=fabs(g.Elt(i))*FMAX(fabs(p->Elt(i)),1.0)/den;

            if ( temp > test ) 
            {
               test=temp;
            }
        }
        
        if ( (test < GTOL) )
        {
			NumHess(p,hessin,parameters,data,func);
			
			return(0);
        }

if (((*iter+10) % 20) != 0) 
{
			
		      
        for (i = 0;i < n;i++) 
        {
            dg.Elt(i)=g.Elt(i)-dg.Elt(i);
        }
        
        for (i=0;i<n;i++) 
        {
            hdg.Elt(i)=0.0;

            for (j = 0;j < n; j++) 
            {
                hdg.Elt(i) += hessin->Elt(i,j)*dg.Elt(j);
            }
        }
        
        fac=0.0;
		fae=0.0;
		sumdg=0.0;
		sumxi=0.0;
        
        for (i=0;i<n;i++) 
        {
            fac += dg.Elt(i)*xi.Elt(i);
            fae += dg.Elt(i)*hdg.Elt(i);
            sumdg += SQR2(dg.Elt(i));
            sumxi += SQR2(xi.Elt(i));
        }

        if (fac*fac > EPSM*sumdg*sumxi) 
        {
            fac=1.0/fac;
            fad=1.0/fae;
            
            for (i=0;i<n;i++) 
            {
                dg.Elt(i)=fac*xi.Elt(i)-fad*hdg.Elt(i);
            }
        
            for (i=0;i<n;i++) 
            {
                for (j=0;j<n;j++) 
                {
                    hessin->Elt(i,j) += fac*xi.Elt(i)*xi.Elt(j)
                    -fad*hdg.Elt(i)*hdg.Elt(j)+fae*dg.Elt(i)*dg.Elt(j);
		         }
            }
        }
}
        

        for (i=0;i<n;i++) 
        {
            xi.Elt(i)=0.0;
        
            for (j=0;j<n;j++) 
            {
                xi.Elt(i) -= hessin->Elt(i,j)*g.Elt(j);
			}
        }
    }

    
    return(-1);
}



void lnsrch(ARM_Vector* xold,
			double fold,
			ARM_Vector* g,
            ARM_Vector* p,
			ARM_Vector* x,
            double* f,
			double stpmax,
			int *check,
			void ** parameters,
			ARM_Matrix* data,
            T_FUNC func
			)
{
    int i, n, test1;
	double a,alam,alam2 = 0.0,alamin,b,disc,f2 = 0.0,fold2 = 0.0;
    double rhs1,rhs2,slope,sum,temp,test,tmplam;
	

    n=xold->GetSize();

    for (sum=0.0,i=0;i<n;i++) 
    {
        sum += p->Elt(i)*p->Elt(i);
    }

    sum=sqrt(sum);

    if (sum > stpmax)
    {
        for (i=0;i<n;i++) 
        {
            p->Elt(i) *= stpmax/sum;
        }
    }

    for (slope=0.0,i=0;i<n;i++)
    {
        slope += g->Elt(i)*p->Elt(i);
    }

    test=0.0;
    
    for (i=0;i<n;i++) 
    {
        temp=fabs(p->Elt(i))/FMAX(fabs(xold->Elt(i)),1.0);
    
        if (temp > test) 
        {
            test=temp;
        }
    }

    alamin=TOLX/test;
    alam=1.0;
    
    for (;;) 
    {

		
		for (i=0;i<n;i++) x->Elt(i)=xold->Elt(i)+alam*p->Elt(i);
		
		*f = (*func)(data, x , parameters);
		
		printf("%6.10f %10.10f  %10.10f  %10.10f\n", *f,fold,alam,alam2);


		if ( (fabs(alam) < alamin) && (alam > 0.0) )
		{
				alam = -1.0;
				continue;
		}
		
		if ( (fabs(alam) < alamin) && (alam < 0.0) )
		{
			for(i=0 ; i<n; i++)	x->Elt(i)=xold->Elt(i);
			return;
		
		}
		
		else if (*f <= fold+ALF*alam*slope) 
		{
			if (alam2==0.0)
			{
				alam=alam*1.5;
				test1=1;
				fold=*f;
				continue;
			}

/*			if (alam2==0.0)
			{
				alam2=2.0;		
				do
				{
					alam=1.5*alam;
					fold2=*f;
					for (i=0;i<n;i++) x->Elt(i)=xold->Elt(i)+alam*p->Elt(i);
					*f = (*func)(data, x , parameters);
				}
				while (*f<fold2);
				alam2=alam;
				alam=alam2/1.5;
			
			}

			if (alam>0.0)	
					*f=dbrent(0.0, alam, alam2, 1e-5,  &alam,
						xold,
						p,
						parameters,
						data,
						func);
			else
				*f=dbrent(-2.0, alam, 0.0, 1e-5,  &alam,
						xold,
						p,
						parameters,
						data,
						func);

			for (i=0;i<n;i++) x->Elt(i)=xold->Elt(i)+alam*p->Elt(i);
			*f = (*func)(data, x , parameters);
				*/
			return;
		}
		else if (test1==1)
		{
			alam=alam/1.5;
			for (i=0;i<n;i++) x->Elt(i)=xold->Elt(i)+alam*p->Elt(i);
			*f=fold;
			return;
		}

		else 
		{
				if (alam == 1.0)
				{
					tmplam = -slope/(2.0*(*f-fold-slope));
				}
				else if ( (alam >0.0) )
				{
					rhs1 = *f-fold-alam*slope;
					rhs2=f2-fold2-alam2*slope;
					a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
					b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
            
					if (a == 0.0) 
					{
					   tmplam = -slope/(2.0*b);
					}
					else 
					{
						disc=b*b-3.0*a*slope;
                
						if (disc<0.0) 
						{
							alam=1.0e-8;
							continue;
				//            printf("\n ???>>>> Roundoff problem in Minimization: lnsrch.");
						}
						else 
						{
							tmplam=(-b+sqrt(disc))/(3.0*a);
						}
					}
                
					if (tmplam>0.5*alam)
					{
						tmplam=0.5*alam;
					}
				}
		}
		
		if (alam<0.0) 
		{
			alam = 0.1*alam;
			continue;
		}
		alam2=alam;
		f2 = *f;
		fold2=fold;
		alam=FMAX(tmplam,0.1*alam);
	}
}

double dbrent(double ax, double bx, double cx, double tol, double *xmin,
				ARM_Vector *p,
				ARM_Vector *descente,
				void ** parameters,
				ARM_Matrix* data,
				T_FUNC func)
{
	int iter,ok1,ok2,i,n;
	double a,b,d,d1,d2,du,dv,dw,dx=0.0,e=0.0;
	double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

    d = 0.0;
	
	n=p->GetSize();

	ARM_Vector nextpt(n),grad(n);

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);

	x=w=v=bx;

    for (i=0;i<n;i++) nextpt.Elt(i)=p->Elt(i)+x*descente->Elt(i);
	fw=fv=fx=(*func)(data, &nextpt , parameters);

	NumJac(&nextpt,fx,&grad,parameters,data,func);
	for (i=0;i<n;i++) dx +=grad.Elt(i)*descente->Elt(i);
	
	dw=dv=dx;

	for (iter=1;iter<=ITMAXDB;iter++)
	{
		xm=0.5*(a+b);
		tol1=tol*fabs(x)+ZEPS;
		tol2=2.0*tol1;
		if (fabs(x-xm) <= (tol2-0.5*(b-a)))
		{
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) 
		{
			d1=2.0*(b-a);
			d2=d1;
			if (dw != dx) d1=(w-x)*dx/(dx-dw);
			if (dv != dx) d2=(v-x)*dx/(dx-dv);
			u1=x+d1;
			u2=x+d2;
			ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
			ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
			olde=e;
			e=d;
			if (ok1 || ok2) 
			{
				if (ok1 && ok2)
					d=(fabs(d1) < fabs(d2) ? d1 : d2);
				else if (ok1)
					d=d1;
				else
					d=d2;
				if (fabs(d) <= fabs(0.5*olde)) 
				{
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
						d=SIGN(tol1,xm-x);
				}
				else 
				{
					d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
				}
			}
			else 
			{
				d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
			}
		} 
		else 
		{
			d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
		}
		if (fabs(d) >= tol1)
		{
			u=x+d;
		    for (i=0;i<n;i++) nextpt.Elt(i)=p->Elt(i)+u*descente->Elt(i);
			fu=(*func)(data, &nextpt , parameters);
			
		} 
		else 
		{
			u=x+SIGN(tol1,d);
		    for (i=0;i<n;i++) nextpt.Elt(i)=p->Elt(i)+u*descente->Elt(i);
			fu=(*func)(data, &nextpt , parameters);

			if (fu > fx)
			{
				*xmin=x;
				return fx;
			}
		}
	    du=0.0;
		for (i=0;i<n;i++) nextpt.Elt(i)=p->Elt(i)+u*descente->Elt(i);
		NumJac(&nextpt,fx,&grad,parameters,data,func);
		for (i=0;i<n;i++) du +=grad.Elt(i)*descente->Elt(i);

		
		if (fu <= fx)
		{
			if (u >= x) a=x; else b=x;
			MOV3(v,fv,dv, w,fw,dw)
			MOV3(w,fw,dw, x,fx,dx)
			MOV3(x,fx,dx, u,fu,du)
		}
		else
		{
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				MOV3(v,fv,dv, w,fw,dw)
				MOV3(w,fw,dw, u,fu,du)
			} 
			else if (fu < fv || v == x || v == w)
			{
				MOV3(v,fv,dv, u,fu,du)
			}
		}
	}
//	nrerror("Too many iterations in routine dbrent");
	return (9.9e99);
}


void NumJac(ARM_Vector* x,
			double f_exact,
			ARM_Vector* df,
			void ** parameters,
			ARM_Matrix* data,
			T_FUNC func )
{
    int j, n;
    double h , temp ,f_approx ;

    n=x->GetSize();
	f_exact=(*func)(data, x ,parameters);
    for (j=0 ; j<n ; j++)
    {    
        temp = x->Elt(j);
        h = EPSD * fabs(temp);

        if ( h == 0 ) 
        {
           h = EPSD;
        }

        x->Elt(j) = temp + h;
        h = x->Elt(j) - temp;
//		printf("Jacc  ");
		f_approx = (*func)(data, x ,parameters);
        x->Elt(j) = temp;
        df->Elt(j) = ( f_approx - f_exact)/h;
    }
}
    
/*
void NumJac(ARM_Vector* x,
			double f_exact,
			ARM_Vector* df,
			void ** parameters,
			ARM_Matrix* data,
			T_FUNC func )
{
    int j, n;
    double temp ,f_approx ;

    n=x->GetSize();
    f_exact = (*func)(data, x ,parameters);

	ARM_Vector dh(n,0.0);
	
	for (j=0;j<n;j++)
	{
		dh.Elt(j) =FMAX(fabs(x->Elt(j)),1e-2);
		dh.Elt(j)*=(x->Elt(j) !=0 ? x->Elt(j)/fabs(x->Elt(j)) : 1.0);
		dh.Elt(j)*=1e-8;
	}
	
	for (j=0;j<n;j++)
	{
		temp =x->Elt(j) + dh.Elt(j);
		dh.Elt(j)  =temp - x->Elt(j);
	}
    
	
	for (j=0 ; j<n ; j++)
    {    
        temp = x->Elt(j);
/*		if (j==6)  dh.Elt(j)=0.1;
		if (j==7)  dh.Elt(j)=0.4;
		if (j==8)  dh.Elt(j)=2.0;
		x->Elt(j) = temp + dh.Elt(j);
        f_approx = (*func)(data, x ,parameters);
        
		x->Elt(j) = temp;
        df->Elt(j) = ( f_approx - f_exact)/dh.Elt(j);
    }
}

*/
void NumHess(ARM_Vector* x,
			ARM_Matrix* hessien,
			void ** parameters,
			ARM_Matrix* data,
			T_FUNC func )
{
    int i,j, n;
    double temp1 , temp2, h_approx,g1_approx,g2_approx,f_exact ;
	n=x->GetSize();
	ARM_Vector dh(n,0.0);

	f_exact = (*func)(data, x ,parameters);
    
	for (j=0;j<n;j++)
	{
		dh.Elt(j) =FMAX(fabs(x->Elt(j)),1e-2);
		dh.Elt(j)*=(x->Elt(j) !=0 ? x->Elt(j)/fabs(x->Elt(j)) : 1.0);
		dh.Elt(j)*=6.0554544523933429e-6;
	}
	
	for (j=0;j<n;j++)
	{
		temp1 =x->Elt(j) + dh.Elt(j);
		dh.Elt(j)  = temp1 - x->Elt(j);
	}
    
	
	for (i=0;i<n;i++)
	{
		
		for (j=0 ; j<n ; j++)
		{    
			temp1 = x->Elt(i);
			x->Elt(i) = temp1 + dh.Elt(i);
			g1_approx = (*func)(data, x ,parameters);
			x->Elt(i) = temp1;
			
			temp2 = x->Elt(j);
			x->Elt(j) = temp2 + dh.Elt(j);
			g2_approx = (*func)(data, x ,parameters);
			x->Elt(j) = temp2;

			temp1 = x->Elt(i);
			temp2 = x->Elt(j);
			x->Elt(i) +=  dh.Elt(i);
			x->Elt(j) += dh.Elt(j);
			h_approx = (*func)(data, x ,parameters);
			x->Elt(i) = temp1;
			x->Elt(j) = temp2;
			
			hessien->Elt(i,j) = (h_approx 
									- g1_approx
									- g2_approx
									+  f_exact)/(dh.Elt(i)*dh.Elt(j));
		}
	}

}


void powell(ARM_Vector* p,
			ARM_Matrix* xi,
			double ftol,
			int *iter,
			double *fret,
			void ** parameters,
			ARM_Matrix* data,
			T_FUNC func )
{
	void linmin(ARM_Vector* p,
				ARM_Vector* xi,
				double *fret,
				void ** parameters,
				ARM_Matrix* data,
				T_FUNC func );
	int i,ibig,j;
	double del,fp,fptt,t;
	int n=p->GetSize();

	ARM_Vector pt(n),ptt(n),xit(n);
	
	*fret=(*func)(data, p ,parameters);;
	
	for (j=0;j<n;j++) pt.Elt(j)=p->Elt(j);
	for (*iter=1;;++(*iter)) {
		fp=(*fret);
		ibig=0;
		del=0.0;
		for (i=0;i<n;i++) {
			for (j=0;j<n;j++) xit.Elt(j)=xi->Elt(j,i);
			fptt=(*fret);
			linmin(p,&xit,fret,parameters,data,func);
	
			if (fabs(fptt-(*fret)) > del) {
				del=fabs(fptt-(*fret));
				ibig=i;
			}
		}
		if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) 
		{
			return;
		}
		if (*iter == ITMAX) printf("powell exceeding maximum iterations.");
		for (j=0;j<n;j++) 
		{
			ptt.Elt(j)=2.0*p->Elt(j)-pt.Elt(j);
			xit.Elt(j)=p->Elt(j)-pt.Elt(j);
			pt.Elt(j)=p->Elt(j);
		}
		fptt=(*func)(data,&ptt,parameters);
		if (fptt < fp) 
		{
			t=2.0*(fp-2.0*(*fret)+fptt)*SQR2(fp-(*fret)-del)-del*SQR2(fp-fptt);
			if (t < 0.0) 
			{
				linmin(p,&xit,fret,parameters,data,func);
		

				for (j=0;j<n;j++) 
				{
					xi->Elt(j,ibig)=xi->Elt(j,n-1);
					xi->Elt(j,n-1)=xit.Elt(j);
				}
			}
		}
	}
}


//int ncom;
//double *pcom,*xicom,(*nrfunc)(double []);

void linmin(ARM_Vector* p,
			ARM_Vector* xi,
			double *fret, 
			void ** parameters,
			ARM_Matrix* data,
			T_FUNC func )

{
	double brent(double ax, double bx, double cx,
				void ** parameters, 
				ARM_Matrix* data,   
				T_FUNC func,		
				double (*func2)(ARM_Vector*,ARM_Vector*,double,	void** ,ARM_Matrix*,T_FUNC func),
				ARM_Vector* pcom,
				ARM_Vector* xicom,
				double tol, double *xmin);

	double f1dim(ARM_Vector* pcom,
				 ARM_Vector* xicom,
				 double x,
				 void ** parameters,
				 ARM_Matrix* data,
				 T_FUNC func);

	void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
				void** parameters,
				ARM_Matrix* data,
				T_FUNC func,
				double (*func2)(ARM_Vector*,ARM_Vector*,double,	void** ,ARM_Matrix*,T_FUNC func),
				ARM_Vector* pcom,
				ARM_Vector* xicom);

	int j;
	double xx,xmin,fx,fb,fa,bx,ax;

	int n=p->GetSize();

	ARM_Vector pcom(n),xicom(n);
	
	
//	nrfunc=func;
	for (j=0;j<n;j++) {
		pcom.Elt(j)=p->Elt(j);
		xicom.Elt(j)=xi->Elt(j);
	}
	ax=0.0;
	xx=1.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,parameters,data,func,f1dim,&pcom,&xicom);
	*fret=brent(ax,xx,bx,parameters,data,func,f1dim,&pcom,&xicom,TOL,&xmin);
	for (j=0;j<n;j++) {
		xi->Elt(j) *= xmin;
		p->Elt(j) += xi->Elt(j);
	}
}

//extern int ncom;
//extern double *pcom,*xicom,(*nrfunc)(double []);

double f1dim(ARM_Vector* pcom,
			 ARM_Vector* xicom,
			 double x,
			 void ** parameters,
			 ARM_Matrix* data,
			 T_FUNC func)
{
	int j;
	double f;

	int n=pcom->GetSize();

	ARM_Vector xt(n);

	for (j=0;j<n;j++) xt.Elt(j)=pcom->Elt(j)+x*xicom->Elt(j);
	f=(*func)(data,&xt,parameters);
	return (f);
}

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
			void** parameters,
			ARM_Matrix* data,
			T_FUNC func,
			double (*func2)(ARM_Vector*,ARM_Vector*,double,	void** ,ARM_Matrix*,T_FUNC func),
			ARM_Vector* pcom,
			ARM_Vector* xicom)
{
	double ulim,u,r,q,fu,dum;

	*fa=(*func2)(pcom,xicom,*ax,parameters,data,func);
	*fb=(*func2)(pcom,xicom,*bx,parameters,data,func);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func2)(pcom,xicom,*cx,parameters,data,func);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func2)(pcom,xicom,u,parameters,data,func);
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
			fu=(*func2)(pcom,xicom,u,parameters,data,func);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func2)(pcom,xicom,u,parameters,data,func);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu, (*func2)(pcom,xicom,u,parameters,data,func) )
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func2)(pcom,xicom,u,parameters,data,func);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func2)(pcom,xicom,u,parameters,data,func);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}


double brent(double ax, double bx, double cx,
				void ** parameters, 
				ARM_Matrix* data,   
				T_FUNC func,		
				double (*func2)(ARM_Vector*,ARM_Vector*,double,	void** ,ARM_Matrix*,T_FUNC func),
				ARM_Vector* pcom,
				ARM_Vector* xicom,
				double tol, double *xmin)
{
	int iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;

    d = 0.0;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*func2)(pcom,xicom,x,parameters,data,func);
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
		fu=(*func2)(pcom,xicom,u,parameters,data,func);
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
	printf("Too many iterations in brent");
	return(999999.123456789);
}
