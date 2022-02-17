#include <stdio.h>
#include <math.h>

#ifdef unix  
    #include<ieeefp.h>
#else 
    #include<float.h>
#endif

#include <string.h>

#include "linalg.h" 
#include "matlib.h"




#define ARM_NBITER_MAX 600 



double leastsq(ARM_Vector* p, int* CST,
               void** parameters,
               ARM_Matrix* data,
               T_FUNC func,
               double* OP)
               
{
    int nvar = p->GetSize();
    int ncst = 0;
    int l = 0;
    int nfun= *((int *) parameters[0]);
    int fin = 1; /*---- MAB : Initialized to 0 before -----*/
    int nbIter = 0;


// FIXMEFRED: mig.vc8 (22/05/2007 15:49:17): 	register i, j, k;
	register int i, j, k;
    int  PCNT = 0 , status=-1, temp = 0;
    double tmp = 0.0;
    double EstSum = 0.5, GradFactor = 0.0;

    
    double fnew = 0.0, FIRSTF = 0.0, GDOLD = 0.0, gdnew = 0.0, fbest = 0.0;
    double stepsize = 0,stepmin = 0;
    double test1 = 0.0, test2 = 0.0, test3 = 0.0;

    ncst = 0;

    for (i = 0; i < nvar;i++) 
    {
        if (CST[i]) 
           ncst++;
    }

    ARM_Vector XOUT_CST(ncst);

    nvar -= ncst;

    ARM_Vector XOUT(nvar);

    j = 0;
    k = 0;

    for (i = 0; i < nvar+ncst; i++) 
    {
        if (CST[i]) 
        {
           XOUT_CST.Elt(j) = p->Elt(i);

           j++;
        }
        else
        {
           XOUT.Elt(k) = p->Elt(i);

           k++;
        }
    }


    ARM_Vector f(nfun) ;
    ARM_Vector OLDF(nfun) ;

    ARM_Vector OLDX(XOUT);
    ARM_Vector CHG(XOUT);
    ARM_Vector OX(nvar);

    ARM_Vector gradf(nvar);
    ARM_Vector SD(nvar);


    ARM_Matrix GRAD(nvar , nfun);
    ARM_Matrix OLDG(nvar , nfun);

    ARM_Matrix TMP(nvar , nvar);

    ARM_Vector MATX(3, 0.0) , MATL(3, 0.0);

    // TMP p->Print();

    fnew = func(data,p,parameters);
    
#ifdef unix         
    if (!finite(fnew))
#else               
    if (!_finite(fnew))
#endif
    {
       // TMP printf("\nInitial parameters too distant of the solution");

       status=1;

       fin=3;
    }

    for (i = 0; i < nfun; i++)
    {
        f.Elt(i)=data->Elt(i,2);
    }

    OLDF = f;

    for (i = 0; i < nfun; i++) 
    {
        FIRSTF += f.Elt(i)  * f.Elt(i) ;
    }

    MATX.Elt(0)= FIRSTF;

    if ( OP[14] == 0.0 ) 
    {
       OP[14] = nvar*200.0;
    }

    while (( status != 1 ) && ( nbIter < ARM_NBITER_MAX ))
    {
        nbIter = nbIter+1;

        // Work Out Gradients

        for (j = 0; j < nvar; j++)
        {
            CHG.Elt(j)=MAX(fabs(XOUT.Elt(j)),1.0e-2);
        
            CHG.Elt(j)*=(XOUT.Elt(j) !=0 ? XOUT.Elt(j)/fabs(XOUT.Elt(j)) : 1.0);
            CHG.Elt(j)*=1.0e-8;
        }
        
        for (j = 0; j < nvar; j++)
        {
            tmp =XOUT.Elt(j) + CHG.Elt(j);

            CHG.Elt(j) = tmp-XOUT.Elt(j);
        }
    
        for (l = 0; l < nvar; l++)
        {    
            tmp = XOUT.Elt(l);

            XOUT.Elt(l) = tmp+CHG.Elt(l);
            
            j = 0;
            k = 0;

            for (i = 0; i < nvar+ncst; i++) 
            {
                if (CST[i]) 
                {
                   p->Elt(i) = XOUT_CST.Elt(j);    

                   j++;
                }
                else
                {
                   p->Elt(i) = XOUT.Elt(k);    

                   k++;
                }
            }

            func(data,p,parameters);
            
            for (i = 0; i < nfun; i++)
            {
                GRAD.Elt(l,i) = (data->Elt(i,2)-OLDF.Elt(i))/CHG.Elt(l);
            }

            XOUT.Elt(l) = tmp;
        }
    
        OP[10] += nvar ;
        
        for (i = 0; i < nvar; i++)
        {
            tmp = 0.0;

            for (j = 0; j < nfun; j++)
            {
                tmp += GRAD.Elt(i,j)*f.Elt(j);
            }

            gradf.Elt(i) = 2.0*tmp;
        }

        fnew = 0.0;

        for (i = 0; i < nfun; i++) 
        {
            fnew += f.Elt(i)*f.Elt(i);
        }

        //------ Initialization of Search Direction ------

        if ( status == -1 )
        {
            for (i = 0; i < nvar; i++)
            {
                for (j = 0; j < nvar; j++)
                {
                    tmp = 0.0;

                    for (k = 0; k < nfun; k++)
                    {
                        tmp += GRAD.Elt(i,k)*GRAD.Elt(j,k);
                    }
        
                    TMP.Elt(i,j) = tmp;
                }

                TMP.Elt(i,i) += GradFactor;
            }

            for (i = 0; i < nvar; i++) 
            {
                SD.Elt(i) = -gradf.Elt(i)/2.0;
            }

            TMP.LinSolve(&SD, tmp);
            
            FIRSTF = fnew;
            OLDG = GRAD;

            GDOLD = 0.0;

            for (i = 0; i < nvar; i++) 
            {
                GDOLD +=gradf.Elt(i)*SD.Elt(i);
            }
                    
            // OPTIONS(18) controls the initial starting step-size.
            // If OPTIONS(18) has been set externally then it will
            // be non-zero, otherwise set to 1.

            if ( OP[18] == 0.0 )
            {
               OP[18] = 1.0;
            }
            
            for (i = 0; i < nvar; i++) 
            {
                XOUT.Elt(i)=XOUT.Elt(i)+OP[18]*SD.Elt(i);
            }

            if ( OP[5] == 0 ) 
            {   
               // Calcul de GradFactor

               GradFactor = 0.0;

               for (i = 0; i < nfun; i++)
               {
                   tmp = 0.0;

                   for (j = 0; j < nvar; j++)
                   {
                       tmp += GRAD.Elt(j,i)*SD.Elt(j);
                   }
            
                   GradFactor += (tmp+f.Elt(i))*(tmp+f.Elt(i));                
               }

               // Cacul de SD

               for (i = 0; i < nvar; i++)
               {
                   for(j = 0; j < nvar; j++)
                   {
                      tmp = 0.0;

                      for (k = 0; k < nfun; k++)
                      {
                         tmp += GRAD.Elt(i,k)*GRAD.Elt(j,k);
                      }
            
                      TMP.Elt(i,j) = tmp;
                   }

                   TMP.Elt(i,i) +=GradFactor;
               }

                for (i = 0; i < nvar; i++) 
                {
                    SD.Elt(i)=-gradf.Elt(i)/2.0;
                }

                TMP.LinSolve(&SD,tmp);
            }

            // Calcul de EstSum

            EstSum = 0.0;

            for (i = 0; i < nfun; i++)
            {
                tmp = 0.0;

                for (j = 0; j < nvar; j++)
                {
                    tmp += GRAD.Elt(j,i)*SD.Elt(j);
                }
        
                EstSum += (tmp+f.Elt(i))*(tmp+f.Elt(i));                
            }

            for (i = 0; i < nvar; i++) 
            {
                XOUT.Elt(i) = XOUT.Elt(i)+OP[18]*SD.Elt(i);
            }

            status = 0;

            if ( OP[7] == 0 ) 
            {
               PCNT = 1;
            }
        }
        else
        {
            //-------------Direction Update------------------

            gdnew=0.0;

            for(i=0;i<nvar;i++) 
            {
               gdnew +=gradf.Elt(i)*SD.Elt(i);
            }
                
            // Case 1: New function is bigger than last and gradient w.r.t. SD -ve
            // ... interpolate. 

            if (( gdnew > 0.0 ) && ( fnew > FIRSTF ))
            {
               stepsize = cubici1(fnew,FIRSTF,gdnew,GDOLD,OP[18]);

               OP[18] = 0.9*stepsize;
            }
            else if ( fnew < FIRSTF )
            {
               //  New function less than old fun. and OK for updating 
               //        .... update and calculate new direction.
                
               stepmin = OP[18];
               fbest = cubici3(fnew,FIRSTF,gdnew,GDOLD,&stepmin);


               if ( fbest> fnew ) 
               {
                  fbest = 0.9*fnew; 
               }

               if ( gdnew < 0.0 )
               {
                  if ( stepmin < OP[18])  
                  {
                     stepmin = 2.0*OP[18]+1.0e-4;
                  }

                  OP[18] = fabs(stepmin);
               }
               else
               {
                  if (OP[18]>0.9)
                  {
                     OP[18]=MIN(fabs(stepmin),1.0);
                  }
               }
            
               // SET DIRECTION.
               // Gauss-Newton Method    

               temp=1;

               if (OP[5]==1)
               {                
                  if (OP[18]>1.0e-8)
                  {
                     // a faire

                     // TMP printf("ERROR- GN not descent direction");

                     temp=0;
                  }
                  else
                  {    
                     // TMP printf("Conditioning of Gradient Poor-Switching To LM method");

                     OP[5]=0;
                     OP[18]=fabs(OP[18]);
                  }
               }

               if (temp)
               {
                  // Levenberg_marquardt Method 
                  // N.B. EstSum is the estimated sum of squares.
                  // GradFactor is the value of lambda.
                  // Estimated Residual:

                  if ( EstSum > fbest )
                     GradFactor=GradFactor/(1.0+OP[18]);
                  else
                     GradFactor=GradFactor+(fbest-EstSum)/(OP[18]+EPS);
                    
                  for (i = 0; i < nvar; i++)
                  {
                      for (j = 0; j < nvar; j++)
                      {
                          tmp = 0.0;

                          for (k = 0; k < nfun; k++)
                          {
                              tmp +=GRAD.Elt(i,k)*GRAD.Elt(j,k);
                          }
                
                          TMP.Elt(i,j)=tmp;
                      }

                      TMP.Elt(i,i) +=GradFactor;
                  }

                  for (i = 0; i < nvar; i++) 
                  {
                      SD.Elt(i)=-gradf.Elt(i)/2.0;
                  }

                  TMP.LinSolve(&SD,tmp);

                  OP[18] = 1.0;

                  EstSum = 0.0;

                  for (i = 0; i < nfun; i++)
                  {
                      tmp = 0.0;

                      for (j = 0; j < nvar; j++)
                      {
                          tmp += GRAD.Elt(j,i)*SD.Elt(j);
                      }

                      EstSum += (tmp+f.Elt(i))*(tmp+f.Elt(i));
                  }
               }

               gdnew = 0.0;

               for (i = 0; i < nvar; i++) 
               {
                   gdnew += gradf.Elt(i)*SD.Elt(i);
               }

               OLDX = XOUT;

               // Save Variables
               FIRSTF = fnew;
               OLDG  = gradf;
               GDOLD = gdnew;    

               // If quadratic interpolation set PCNT

               if ( OP[7] == 0 )
               {
                  PCNT = 1; 
                  MATX.Elt(0) = 0.0;
                  MATX.Elt(1) = 0.0;
                  MATX.Elt(2) = 0.0;
                  MATL.Elt(0) = fnew; 
               }
            }
            else
            {
                // Halve Step-length

                if ( fabs(fnew-FIRSTF) < EPS*1.0e-3 )
                {
                   // TMP printf("\nNo improvement in search direction: Terminating");

                   status = 1;

                   fin = 2;
                }
                else
                {
                   OP[18] = OP[18]/8.0;

                   if ( OP[18] < 1.0e-8 )
                   {
                      OP[18] = -OP[18];
                   }
                }
            }

            for (i = 0; i < nvar; i++) 
            {
                XOUT.Elt(i) = OLDX.Elt(i)+OP[18]*SD.Elt(i);
            }
        }

        //----------End of Direction Update-------------------

        if ( OP[7] == 0 )
        {
           PCNT = 1;

           MATX.Elt(0)=0.0;
           MATX.Elt(1)=0.0;
           MATX.Elt(2)=0.0;
           MATL.Elt(0)=fnew;
        }

        // Check Termination 
        test1=0.0;
        test2=0.0;
        test3=0.0;

        for (i = 0; i < nvar; i++)
        {
            if ( fabs(SD.Elt(i)) > test1 ) 
            {
               test1 = fabs(SD.Elt(i));
            }
            
            if ( fabs(SD.Elt(i)*gradf.Elt(i)) > test2 ) 
            {
               test2 = fabs(SD.Elt(i)*gradf.Elt(i));
            }
            
            if ( fabs(gradf.Elt(i)) > test3 ) 
            {
               test3 = fabs(gradf.Elt(i));
            }
        }

        // TMP printf("fnew :%10.10f iter:%10.0f test1:%10.10f test2:%10.10f test3:%10.10f\n", fnew , OP[10],test1,test2,test3);

        // TMP printf ("Parametre :");

        for (i = 0; i < nvar; i++) 
        {
            printf("   v[%d]:%10.10f",i,XOUT.Elt(i) );
        }

     /* TMP
        printf ("\nResidus :");

        for (i = 0; i < nfun; i++) 
            printf(" - R[%d]:%10.10f",i,f.Elt(i) );
        
     */

        if (( test1 < OP[2] ) && ( test2 < OP[3] ) 
            && ( test3 < 10.0*(OP[3]+OP[2])))
        {
           // TMP printf("\nOptimization Terminated Successfully");

           status = 1;

           fin = 0;
        }
        else if ( OP[10] > OP[14] )
        {
           // TMP printf("\nMaximum number of iterations has been exceeded");

           status = 1;

           fin = 1;
        }
        else
        {
           // Line search using mixed polynomial interpolation and extrapolation.

           if ( PCNT != 0 )
           {
              // <= used in case when no improvement found.

              tmp = 0.0;

              for (i = 0; i < nfun; i++) 
              {
                  tmp += OLDF.Elt(i)*OLDF.Elt(i);
              }

              while ( PCNT > 0 )
              {
                  j = 0;

                  k = 0;

                  for (i = 0; i < nvar+ncst; i++) 
                  {
                      if (CST[i]) 
                      {
                         p->Elt(i) = XOUT_CST.Elt(j);    

                         j++;
                      }
                      else
                      {
                         p->Elt(i) = XOUT.Elt(k);    

                         k++;
                      }
                  }

                  fnew = func(data,p,parameters);

                  for (i = 0; i < nfun; i++) 
                  {
                      f.Elt(i) = data->Elt(i,2);
                  }

                  OP[10] = OP[10]+1;

                  if ( OP[10] >= OP[14] ) 
                  {
                     // TMP printf("\nMaxmum number of iterations has been exceeded");

                     status = 1;

                     fin = 3;

                     PCNT = 0;
                  }

                  fnew = 0.0;

                  for (i = 0; i < nfun; i++) 
                  {
                      fnew += f.Elt(i)*f.Elt(i);
                  }

#ifdef unix         
                    if (!finite(fnew)) fnew=1.0e11;
#else               
                    if (!_finite(fnew)) fnew=1.0e11;
#endif

                    if ( fnew <= tmp )
                    {
                        OX = XOUT;

                        OLDF = f;

                        tmp = 0.0;

                        for (i = 0; i < nfun; i++) 
                        {
                            tmp += OLDF.Elt(i)*OLDF.Elt(i);
                        }
                    }
                    
                    searchq(&PCNT,&fnew,0.0,&MATL,&MATX,GDOLD,&OP[18]);
                    // TMP searchq(&PCNT,&fnew,OLDX,&MATL,&MATX,GDOLD,&OP[18]);

                    for (i = 0; i < nvar; i++) 
                    {
                        XOUT.Elt(i) = OLDX.Elt(i)+OP[18]*SD.Elt(i);
                    }

                    if ( fabs(fnew-FIRSTF) < EPS*1.0e-3 ) 
                    {
                       PCNT = 0 ; 
                    }

                    if ( OP[10] >= OP[14] ) 
                    {
                        // TMP printf("\nMaxmum number of iterations has been exceeded");

                        status=1;

                        fin = 3;

                        PCNT = 0;
                    }
                }

                XOUT = OX;

                f = OLDF;
            }
            else
            {
               j = 0;
               k = 0;

               for (i=0;i<nvar+ncst;i++) 
               {
                   if (CST[i]) 
                   {
                      p->Elt(i)=XOUT_CST.Elt(j);    

                     j++;
                   }
                   else
                   {
                      p->Elt(i)=XOUT.Elt(k);    

                      k++;
                   }
               }

               func(data, &XOUT, parameters);

               for (i = 0; i < nfun; i++)
               {
                   f.Elt(i)=data->Elt(i,2);
               }

               OP[10] = OP[10]+1;
            }
        }
    }

    OP[8] = fnew;

    j = 0;
    k = 0;

    for (i = 0; i < nvar+ncst; i++) 
    {
        if (CST[i]) 
        {
           p->Elt(i) = XOUT_CST.Elt(j);    

           j++;
        }
        else
        {
           p->Elt(i)=OLDX.Elt(k);    

           k++;
        }
    }
    
    return(fin);

    // fin=0 => Optimization Terminated Successfully
    // fin=1 => Maxmum number of iterations has been exceeded
    // fin=2 => No improvement in search direction: Terminating
    // fin=3 => Initial parameters too distant of the solution
}


/*---- End of File ----*/
