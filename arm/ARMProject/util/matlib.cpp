/*
 * $Log: matlib.cpp,v $
 * Revision 1.4  2003/08/05 07:57:21  jmprie
 *  ajout de GaussLegendre()
 *
 * Revision 1.3  2003/07/21 10:44:33  mab
 * suppression of ";" in :
 * if ( (stepmin<0.0) || (!_finite(stepmin)) ) ;
 * around line : 325
 *
 * Revision 1.2  2002/11/25 16:40:33  mab
 * Formatting
 *
 */


/*---------------------------------------------------------------------*/




#include "matlib.h" 
#include "armglob.h" 
#include "float.h" 





double* default_op(double* opt)
{
// Default parameters used by the optimization routines.
//    The parameters are:
//    OPTIONS(1)-
//    OPTIONS(2)-Termination tolerance for X.(Default: 1e-4).
//    OPTIONS(3)-Termination tolerance on F.(Default: 1e-4).
//    OPTIONS(4)-Termination criterion on constraint violation.(Default: 1e-6)
//    OPTIONS(5)-Algorithm: Strategy:  Not always used.
//    OPTIONS(6)-Algorithm: Optimizer: Not always used. 
//    OPTIONS(7)-Algorithm: Line Search Algorithm. (Default 0)
//    OPTIONS(8)-Function value. (Lambda in goal attainment. )
//    OPTIONS(9)-Set to 1 if you want to check user-supplied gradients
//    OPTIONS(10)-Number of Function and Constraint Evaluations.
//    OPTIONS(11)-Number of Function Gradient Evaluations.
//    OPTIONS(12)-Number of Constraint Evaluations
//    OPTIONS(13)-Number of equality constraints. 
//    OPTIONS(14)-Maximum number of iterations. (Default 100*no. of variables)
//    OPTIONS(15)-Used in goal attainment for special objectives. 
//    OPTIONS(16)-Minimum change in variables for finite difference gradients.
//    OPTIONS(17)-Maximum change in variables for finite difference gradients.
//    OPTIONS(18)- Step length. (Default 1 or less). 

    int i;


    if ( opt == NULL )
    {
       opt = new double[19];

       for (i = 0; i < 19; i++) 
       {
           opt[i]=0.0;
       }
    }
     
    if (opt[2]==0.0) opt[2]=1.0e-4;
    if (opt[3]==0.0) opt[3]=1.0e-4;
    if (opt[4]==0.0) opt[4]=1.0e-6;

    if (opt[16]==0.0) opt[16]=1.0e-8;
    if (opt[17]==0.0) opt[17]=0.1;
    if (opt[18]==0.0) opt[18]=1.0;
    
    return(opt);
}



void searchq(int* pcnt,
                double* fnew,
                double oldx,
                ARM_Vector *matl,
                ARM_Vector *matx,
                double gdold,
                double* stepsize)
{
    // SEARCHQ Line search routine for FMINU and LEASTSQ functions.
    // Performs line search procedure for unconstrained and least squares
    // optimization. Uses Quadratic Interpolation.
    // When finished pcnt returns 0.


    double newstep,tmp,tmp2;

    if ( *pcnt == 1 )
    {
       // Case 1: Next point less than initial point. 
       // Increase step-length based on last gradient evaluation

        if ( *fnew < matl->Elt(0) )
        {
           // Quadratic Extrapolation using gradient of first point and 
           // values of two other points.

           matl->Elt(1)=*fnew;
           matx->Elt(1)=*stepsize;
           newstep = -0.5*gdold*(*stepsize)*(*stepsize);
           newstep /= (*fnew-gdold*(*stepsize)-matl->Elt(0)+1.0e-16);

           if ( newstep < *stepsize )
              newstep = 1.2*(*stepsize);
            
           *stepsize = 1.2*newstep;

           *pcnt = 2;
        }
        else
        {
            // Case 2: New point greater than initial point. Decrease step-length.

            matl->Elt(2)=*fnew;
            matx->Elt(2)=*stepsize;

            // Interpolate to get stepsize

            tmp = -gdold*0.5*(*stepsize)*(*stepsize);
            tmp/=(*fnew-gdold*(*stepsize)-matl->Elt(0)+EPS);
            
            *stepsize=MAX(1.0e-8*(*stepsize),tmp);

            *pcnt=3;
        }
    }
    else if ((*pcnt==2)  && (*fnew< matl->Elt(1)) )
    {
        // Case 3: Last run was Case 1 (pcnt=2) and new point less than 
        //      both of other 2. Replace. 

            tmp=matl->Elt(2);
            tmp2=matx->Elt(2);

            matl->Elt(2)=*fnew;
            matx->Elt(2)=*stepsize;
            
            newstep=cubici2(gdold,matl,matx);
            matl->Elt(0)=matl->Elt(1);
            matx->Elt(0)=matx->Elt(1);
            matl->Elt(1)=*fnew;
            matx->Elt(1)=*stepsize;

            matl->Elt(2)=tmp;
            matx->Elt(2)=tmp2;

            *stepsize=MIN(1.0,newstep)+1.5*(*stepsize);
            
            *stepsize=MAX(1.2*newstep,1.2* (*stepsize));
    }
    else if ( (*pcnt==3) && (*fnew>=matl->Elt(0)))
    {
       // Case 4: Last run was Case 2: (pcnt=3) and new function still 
       //      greater than initial value.

        matl->Elt(1)=*fnew;
        matx->Elt(1)=*stepsize;

        if (*stepsize<1e-8)
        {
            newstep=-(*stepsize)/2.0;

            if (fabs(newstep)<EPS) 
            {
                newstep=(double)rand()/( (double) RAND_MAX) -0.5;
            }
        }
        else
        {
            matl->Elt(1)=matl->Elt(2);
            matx->Elt(1)=matx->Elt(2);
            tmp=matl->Elt(2);
            tmp2=matx->Elt(2);
            
            matl->Elt(2)=*fnew;
            matx->Elt(2)=*stepsize;

            newstep=cubici2(gdold,matl,matx);

            matl->Elt(2)=tmp;
            matx->Elt(2)=tmp2;
            
            matl->Elt(1)=*fnew;
            matx->Elt(1)=*stepsize;
        }

        matx->Elt(2)=*stepsize;

#ifdef unix         
       if (!finite(newstep)) *stepsize=*stepsize/2.0; 
#else               
       if (!_finite(newstep)) *stepsize=*stepsize/2.0; 
#endif
       
        else 
            *stepsize=newstep; 

        matl->Elt(2)=*fnew;
    }
    else if ((*pcnt==2) && (*fnew>matl->Elt(1)))
    {
       // Otherwise must have Bracketed MINmum so do quadratic interpolation.
       //  ... having just increased step.

        matx->Elt(2)=*stepsize;
        matl->Elt(2)=*fnew;

        *stepsize=cubici2(gdold,matl,matx);
        *pcnt=4;
    }
    else if ((*pcnt==3) && (*fnew<matl->Elt(0)))
    {
        // ...  having just reduced step.

        matx->Elt(1)=*stepsize;
        matl->Elt(1)=*fnew;

        *stepsize=cubici2(gdold,matl,matx);
        *pcnt=4;
    }
    else if (*pcnt==4 )
    {
        // Have just interpolated - Check to see whether function is any better 
        // - if not replace.

        *pcnt=0;
        *stepsize=fabs(*stepsize);

        // If interpolation failed use old point.

        if (*fnew>matl->Elt(1))
        {
            *fnew=matl->Elt(1);
            *stepsize=matx->Elt(1);
        }
    }
}



double cubici1(double fnew,double fold,double graddnew,
               double graddold,double stepsize)
//CUBICI1 Cubicly interpolates 2 points and gradients to estimate MINmum.
//
//    This function uses cubic interpolation and the values of two 
//    points and their gradients in order to estimate the MINmum of a 
//    a function along a line.

{
    double z,w,tmp;

#ifdef unix         
     if (!finite(fnew)) fnew=1.0e11;
#else               
     if (!_finite(fnew)) fnew=1.0e11;
#endif    
    
    z=3.0*(fold-fnew)/stepsize+graddold+graddnew;
    tmp=z*z-graddold*graddnew;
    
    if ( tmp > 0.0 )
    {
       w = sqrt(tmp);
    }
    else 
    {
       w = 0.0;
    }

    stepsize=stepsize*(z+w-graddold);

    stepsize/=(graddnew-graddold+2.0*w);

    return(stepsize);
}


 
double cubici2(double graddold, ARM_Vector* matl, ARM_Vector* matx)
{
    // CUBICI2 Cubicly interpolates 3 points and 1 gradient.
    //
    // This function uses cubic interpolation and
    // the values of 3  points and one gradient.

// FIXMEFRED: mig.vc8 (22/05/2007 15:48:54): 	register i;
	register int i;
    double stepmin,x1,root;
    ARM_Matrix A(3,3);
    ARM_Vector abd(3);



    for (i=0;i<3;i++)
        A.Elt(i,0)=1/3.0*pow(matx->Elt(i),3.0);

    for (i=0;i<3;i++)
        A.Elt(i,1)=0.5*pow(matx->Elt(i),2.0);
    
    for (i=0;i<3;i++)
        A.Elt(i,2)=1.0;

    for (i=0;i<3;i++)
        abd.Elt(i)=matl->Elt(i)-graddold*matx->Elt(i);

    
    A.QRSolve(&abd , stepmin);

    stepmin=abd.Elt(1)*abd.Elt(1)-4.0*abd.Elt(0)*graddold;
    
    root=MAX(0.0,stepmin);
    
    root=sqrt(root);
    
    x1=(-abd.Elt(1)+root)/(2.0*abd.Elt(0));
    
    if ( 2.0*abd.Elt(0)*x1+abd.Elt(1) > 0.0 )
    {
       stepmin=x1;
    }
    else
    {
       stepmin=(-abd.Elt(1)-root)/(2.0*abd.Elt(0));
    }


#ifdef unix         
     if ( (stepmin<0.0) || (!finite(stepmin)) ) 
#else               
     if ( (stepmin<0.0) || (!_finite(stepmin)) )
#endif    
     {
        for (i=0;i<3;i++)
            A.Elt(i,0)=.5*pow(matx->Elt(i),2.0);

        for (i=0;i<3;i++)
            A.Elt(i,1)=matx->Elt(i);

        for (i=0;i<3;i++)
            A.Elt(i,2)=1.0;

        for (i=0;i<3;i++)
            abd.Elt(i)=matl->Elt(i);

        A.QRSolve(&abd,stepmin);

        stepmin=fabs( -abd.Elt(1)/abd.Elt(0) );
    }

#ifdef unix         
     if (!finite(stepmin)) stepmin=matx->Elt(1)/2.0;
#else               
     if (!_finite(stepmin)) stepmin=matx->Elt(1)/2.0;
#endif    

     return(stepmin);
}    



double cubici3(double fnew,double fold,double graddnew,
               double graddold,double *stepsize)

// CUBICI3  Cubicly interpolates 2 points and gradients to find step and min.
//    This function uses cubic interpolation and the values of 
//    two points and their gradients in order estimate the MINmum of a 
//    a function along a line.

{
    double root,x1,tmp;
    double fbest;
    ARM_Matrix A(2,2);
    ARM_Vector abd(2);


#ifdef unix         
       if (!finite(fnew)) fnew=1.0e11; 
#else               
       if (!_finite(fnew)) fnew=1.0e11; 
#endif

    A.Elt(0,0)=(1/3.0)*pow((*stepsize),3.0);
    A.Elt(0,1)=0.5*(*stepsize)*(*stepsize);
    A.Elt(1,0)=(*stepsize)*(*stepsize);
    A.Elt(1,1)=(*stepsize);
    
    abd.Elt(0)=fnew-graddold*(*stepsize)-fold;
    abd.Elt(1)=graddnew-graddold;

    A.LinSolve(&abd,tmp);
    
    tmp=abd.Elt(1)*abd.Elt(1) - 4.0*abd.Elt(0)*graddold;

    root=MAX(0.0,tmp);
    root=sqrt(root);

    x1=(-abd.Elt(1)+root)/(2.0*abd.Elt(0));

    if ( 2.0*abd.Elt(0)*x1+abd.Elt(1) > 0.0 )
    {
       (*stepsize)=x1;
    }
    else
    {
       (*stepsize)=(-abd.Elt(1)-root)/(2.0*abd.Elt(0));
    }

    if ((*stepsize) < 0.0) 
       (*stepsize)=-(*stepsize);
    
    fbest=(1/3.0)*abd.Elt(0)*pow((*stepsize),3.0);
    fbest+=0.5*abd.Elt(1)*(*stepsize)*(*stepsize);
    fbest+=graddold*(*stepsize)+fold;
    

    return(fbest);
}



void updhess(ARM_Vector xnew,ARM_Vector xold,ARM_Vector gradxnew,
             ARM_Vector gradxold,ARM_Matrix *invhess,double *para,
             ARM_Vector* SD)
// UPDHESS Performs the Inverse Hessian Update.
// Returns direction of search for use with 
// unconstrained optimization problems. 

{
    int i,j,n;
    double tmp,tmp2;

    n = xnew.GetSize();

    ARM_Vector u(n),v(n);


    
    for (i = 0; i < n; i++)
    {
       u.Elt(i) = xnew.Elt(i)-xold.Elt(i);

       v.Elt(i) = gradxnew.Elt(i)-gradxold.Elt(i);
    }

    if ( para[6] == 0 )
    {
       // The BFGS Hessian Update formula:
        
        tmp = 0.0; //  tmp = (v'*u)

        for(i = 0; i < n; i++)
        {
            tmp += v.Elt(i)*u.Elt(i);
        }
        
        // tmp = u'*invhess*u
  
        tmp2=0.0;

        for(i=0;i<n;i++)
        {            
            xnew.Elt(i)=0.0;

            for(j=0;j<n;j++)
            {
                 xnew.Elt(i) += invhess->Elt(i,j)*u.Elt(j);
            }

            tmp2 +=u.Elt(i)*xnew.Elt(i);
        }
        
        // invhess + v*v'/(v'*u)

        for ( i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                invhess->Elt(i,j) += v.Elt(i)*v.Elt(j)/tmp;
            }
        }
        
        // invhess=-invhess*u*u'*invhess'/(u'*invhess*u);

        for (i = 0; i < n; i++)
        {            
            for (j = 0; j < n; j++)
            {
                invhess->Elt(i,j) -= xnew.Elt(i)*xnew.Elt(j)/tmp2;
            }
        }

        // invhess=invhess + v*v'/(v'*u)  -invhess*u*u'*invhess'/(u'*invhess*u);

        for (i = 0; i < n; i++)
        { 
            SD->Elt(i)=-gradxnew.Elt(i);
        }

        invhess->LinSolve(SD,tmp);
    }
    else if ( para[6] == 1 )
    {
    // The DFP formula
/*        a=u*u'/(u'*v);
        b=-invhess*v*v'*invhess'/(v'*invhess*v);
        invhess=invhess + a + b;
        directn=-invhess*gradxnew;*/
    }
    else if ( para[6] == 3 )
    {
    // A formula given by Gill and Murray
/*        a = 1/(v'*u);
        invhess=invhess - a*(invhess*v*u'+u*v'*invhess)+a*(1+v'*invhess*v*a)*u*u' ;
           directn=-invhess*gradxnew;*/
    }
    else if ( para[6]== 2 )
    {
       // Steepest Descent
       //directn=-gradxnew;
    }
}



void Hessien(ARM_Vector* x,
             ARM_Matrix* hessien,
             void** parameters,
             ARM_Matrix* data,
             T_FUNC func)
{
    int i,j, n;
    double temp1 , temp2, h_approx,g1_approx,g2_approx,f_exact ;

    n = x->GetSize();

    ARM_Vector dh(n,0.0);

    f_exact = (*func)(data, x, parameters);
    
    for (j = 0; j < n; j++)
    {
        dh.Elt(j) =MAX(fabs(x->Elt(j)),1e-2);
        dh.Elt(j)*=(x->Elt(j) !=0 ? x->Elt(j)/fabs(x->Elt(j)) : 1.0);
        dh.Elt(j)*=6.0554544523933429e-8;
    }
    
    for (j = 0; j < n; j++)
    {
        temp1 = x->Elt(j) + dh.Elt(j);

        dh.Elt(j) = temp1 - x->Elt(j);
    }
    
    
    for (i = 0; i < n; i++)
    {
        for (j = 0 ; j < n ; j++)
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
            
            hessien->Elt(i,j) = (h_approx-g1_approx-g2_approx
                                    +f_exact)/(dh.Elt(i)*dh.Elt(j));
        }
    }
}

/*--------------------------------------------------------------------------*/
/*  Calcul des points et ponderations pour une integration dans             */
/*  l'intervalle [-1,1] par la methode de Gauss-Legendre                    */
/*--------------------------------------------------------------------------*/
void GaussLegendre(int n, double* x, double* w)
{
	int m,j,i;
	double z1,z,pp,p3,p2,p1;

	m=(n+1)/2;
	for (i=0;i<m;++i)
    {
		z=cos(3.141592654*(i+0.75)/(n+0.5));
		do
        {
			p1=1.0;
			p2=0.0;
			for (j=0;j<n;++j)
            {
				p3=p2;
				p2=p1;
				p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
			}
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (fabs(z-z1) > EPS);
		x[i]=-z;
		x[n-1-i]=z;
		w[i]=2.0/((1.0-z*z)*pp*pp);
		w[n-1-i]=w[i];
	}
}


/*---- End of File ----*/
