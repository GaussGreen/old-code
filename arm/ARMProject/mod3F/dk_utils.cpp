/*===========================================================
  Name    : hw_vfdk_LDHD_lattice.cpp
  Owner   : DK
  Created : DK  
  Comment : Utility functions.
  Time    : 16:38 11/27/2002
=============================================================*/
#include <iostream>
#include <strstream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifndef WIN32
#include <sys/param.h>
#endif 
#define ACC 1.0e-04
#define JMAX 5000

#include "dk_utils.h"

#include "DKMaille.h"
#include "DKMaille2D.h"
double g_s;
double g_tau;

#define TOKYO
#define DEBUG_ROTATION(a,i,j,k,l) g=a.at(i,j);h=a.at(k,l);a.at(i,j)=g-g_s*(h+g*g_tau); a.at(k,l)=h+g_s*(g-h*g_tau);

#define M_PI 3.14159265358979

double norcoef_new[6] = {
                            1.330274429 ,-1.821255978 ,1.781477937 ,-.356563782 , .31938153, .0 };
double Poly5_new(double dX, double dCoeffs[6])
{
    return dCoeffs[5]
           +(dCoeffs[4]*dX)
           +(dCoeffs[3]*pow(dX,2))
           +(dCoeffs[2]*pow(dX,3))
           +(dCoeffs[1]*pow(dX,4))
           +(dCoeffs[0]*pow(dX,5));
}

int ordering_dk_maille(double target_point, int *nearest_node,const DKMaille<double> &time_array, int array_length)
{
    int k=1;
    int zflag=0;

    if(target_point>time_array[array_length])
    {
        *nearest_node=array_length;
        return 0;
    }

    if(target_point<time_array[1])
    {
        *nearest_node=1;
        return 0;
    }

    if(target_point==time_array[array_length])
    {
        *nearest_node=array_length;
        return 1;
    }

    while((k<=array_length-1)&&(zflag==0))
    {
        if((target_point>=time_array[k])&&(target_point<time_array[k+1]))
        {
            *nearest_node=k;
            zflag=1;
        }
        k+=1;
    }
    return 1;
}

int ordering2(double target_point, int *nearest_node,const DKMaille<double> &time_array, int array_length)
{
    int k=0;
    int zflag=0;

    if(target_point>time_array[array_length-1])
    {
        *nearest_node=array_length-1;
        return 0;
    }

    if(target_point<time_array[0])
    {
        *nearest_node=0;
        return 0;
    }

    if(target_point==time_array[array_length-1])
    {
        *nearest_node=array_length-1;
        return 1;
    }

    while((k<=array_length-2)&&(zflag==0))
    {
        if((target_point>=time_array[k])&&(target_point<time_array[k+1]))
        {
            *nearest_node=k;
            zflag=1;
        }
        k+=1;
    }
    return 1;
}


int ordering(double target_point, int *nearest_node,const DKMaille<double> &time_array)
{
    int k=0;
    int zflag=0;

    int array_length=time_array.entries();

    if(target_point>time_array[array_length-1])
    {
        *nearest_node=array_length-1;
        return 0;
    }

    if(target_point<time_array[0])
    {
        *nearest_node=0;
        return 0;
    }

    if(target_point==time_array[array_length-1])
    {
        *nearest_node=array_length-1;
        return 1;
    }

    while((k<array_length-1)&&(zflag==0))
    {
        if((target_point>=time_array[k])&&(target_point<time_array[k+1]))
        {
            *nearest_node=k;
            zflag=1;
        }
        k+=1;
    }
    return 1;
}


double interpollinear_dk_maille(int nearest_node, double target_point, const DKMaille<double> &time_array,const DKMaille<double> &rate_array)
{
    double a,b;
    double target_rate;

    a=(time_array[nearest_node+1]-target_point)/(time_array[nearest_node+1]-time_array[nearest_node]);
    b=(target_point-time_array[nearest_node])/(time_array[nearest_node+1]-time_array[nearest_node]);
    target_rate=a*rate_array[nearest_node]+b*rate_array[nearest_node+1];

    return target_rate;
}


double interpol_dk_maille(int nearest_node, double target_point, const DKMaille<double> &time_array, const DKMaille<double> &rate_array, int array_length)
{
    double term1,term2,term3,term4;
    double target_rate;
    int i,s;
    s=array_length;

    if((nearest_node==1)||(nearest_node==2))
    {
        term1=rate_array[1]*((target_point-time_array[2])*(target_point-time_array[3])*(target_point-time_array[4]))
              /((time_array[1]-time_array[2])*(time_array[1]-time_array[3])
                *(time_array[1]-time_array[4]));
        term2=rate_array[2]*((target_point-time_array[1])*(target_point-time_array[3])*(target_point-time_array[4]))
              /((time_array[2]-time_array[1])*(time_array[2]-time_array[3])
                *(time_array[2]-time_array[4]));
        term3=rate_array[3]*((target_point-time_array[1])*(target_point-time_array[2])*(target_point-time_array[4]))
              /((time_array[3]-time_array[1])*(time_array[3]-time_array[2])
                *(time_array[3]-time_array[4]));
        term4=rate_array[4]*((target_point-time_array[1])*(target_point-time_array[2])*(target_point-time_array[3]))
              /((time_array[4]-time_array[1])*(time_array[4]-time_array[2])
                *(time_array[4]-time_array[3]));
    }
    else if((nearest_node==s-1))
    {
        term1=rate_array[s-3]*((target_point-time_array[s-2])*(target_point-time_array[s-1])*(target_point-time_array[s]))
              /((time_array[s-3]-time_array[s-2])*(time_array[s-3]-time_array[s-1])
                *(time_array[s-3]-time_array[s]));
        term2=rate_array[s-2]*((target_point-time_array[s-3])*(target_point-time_array[s-1])*(target_point-time_array[s]))
              /((time_array[s-2]-time_array[s-3])*(time_array[s-2]-time_array[s-1])
                *(time_array[s-2]-time_array[s]));
        term3=rate_array[s-1]*((target_point-time_array[s-3])*(target_point-time_array[s-2])*(target_point-time_array[s]))
              /((time_array[s-1]-time_array[s-3])*(time_array[s-1]-time_array[s-2])
                *(time_array[s-1]-time_array[s]));
        term4=rate_array[s]*((target_point-time_array[s-3])*(target_point-time_array[s-2])*(target_point-time_array[s-1]))
              /((time_array[s]-time_array[s-3])*(time_array[s]-time_array[s-2])
                *(time_array[s]-time_array[s-1]));
    }
    else if(nearest_node==s)
    {
        target_rate=rate_array[s];
        return target_rate;
    }
    else
    {
        i=nearest_node;
        term1=rate_array[i-1]*((target_point-time_array[i])*(target_point-time_array[i+1])*(target_point-time_array[i+2]))
              /((time_array[i-1]-time_array[i])*(time_array[i-1]-time_array[i+1])
                *(time_array[i-1]-time_array[i+2]));
        term2=rate_array[i]*((target_point-time_array[i-1])*(target_point-time_array[i+1])*(target_point-time_array[i+2]))
              /((time_array[i]-time_array[i-1])*(time_array[i]-time_array[i+1])
                *(time_array[i]-time_array[i+2]));
        term3=rate_array[i+1]*((target_point-time_array[i-1])*(target_point-time_array[i])*(target_point-time_array[i+2]))
              /((time_array[i+1]-time_array[i-1])*(time_array[i+1]-time_array[i])
                *(time_array[i+1]-time_array[i+2]));
        term4=rate_array[i+2]*((target_point-time_array[i-1])*(target_point-time_array[i])*(target_point-time_array[i+1]))
              /((time_array[i+2]-time_array[i-1])*(time_array[i+2]-time_array[i])
                *(time_array[i+2]-time_array[i+1]));
    }
    target_rate=term1+term2+term3+term4;
    return target_rate;
}


double interpol2(int nearest_node, double target_point, const DKMaille<double> &time_array, const DKMaille<double> &rate_array)
{
    double term1,term2,term3,term4;
    double target_rate;
    int i,s;
    int array_length=time_array.entries();
    s=array_length;

    if((nearest_node==0)||(nearest_node==1))
    {
        term1=rate_array[0]*((target_point-time_array[1])*(target_point-time_array[2])*(target_point-time_array[3]))
              /((time_array[0]-time_array[1])*(time_array[0]-time_array[2])
                *(time_array[0]-time_array[3]));
        term2=rate_array[1]*((target_point-time_array[0])*(target_point-time_array[2])*(target_point-time_array[3]))
              /((time_array[1]-time_array[0])*(time_array[1]-time_array[2])
                *(time_array[1]-time_array[3]));
        term3=rate_array[2]*((target_point-time_array[0])*(target_point-time_array[1])*(target_point-time_array[3]))
              /((time_array[2]-time_array[0])*(time_array[2]-time_array[1])
                *(time_array[2]-time_array[3]));
        term4=rate_array[3]*((target_point-time_array[0])*(target_point-time_array[1])*(target_point-time_array[2]))
              /((time_array[3]-time_array[0])*(time_array[3]-time_array[1])
                *(time_array[3]-time_array[2]));
    }
    else if((nearest_node==s-2))
    {
        term1=rate_array[s-4]*((target_point-time_array[s-3])*(target_point-time_array[s-2])*(target_point-time_array[s-1]))
              /((time_array[s-4]-time_array[s-3])*(time_array[s-4]-time_array[s-2])
                *(time_array[s-4]-time_array[s-1]));
        term2=rate_array[s-3]*((target_point-time_array[s-4])*(target_point-time_array[s-2])*(target_point-time_array[s-1]))
              /((time_array[s-3]-time_array[s-4])*(time_array[s-3]-time_array[s-2])
                *(time_array[s-3]-time_array[s-1]));
        term3=rate_array[s-2]*((target_point-time_array[s-4])*(target_point-time_array[s-3])*(target_point-time_array[s-1]))
              /((time_array[s-2]-time_array[s-4])*(time_array[s-2]-time_array[s-3])
                *(time_array[s-2]-time_array[s-1]));
        term4=rate_array[s-1]*((target_point-time_array[s-4])*(target_point-time_array[s-3])*(target_point-time_array[s-2]))
              /((time_array[s-1]-time_array[s-4])*(time_array[s-1]-time_array[s-3])
                *(time_array[s-1]-time_array[s-2]));
    }
    else if(nearest_node==s-1)
    {
        target_rate=rate_array[s-1];
        return target_rate;
    }
    else
    {
        i=nearest_node;
        term1=rate_array[i-1]*((target_point-time_array[i])*(target_point-time_array[i+1])*(target_point-time_array[i+2]))
              /((time_array[i-1]-time_array[i])*(time_array[i-1]-time_array[i+1])
                *(time_array[i-1]-time_array[i+2]));
        term2=rate_array[i]*((target_point-time_array[i-1])*(target_point-time_array[i+1])*(target_point-time_array[i+2]))
              /((time_array[i]-time_array[i-1])*(time_array[i]-time_array[i+1])
                *(time_array[i]-time_array[i+2]));
        term3=rate_array[i+1]*((target_point-time_array[i-1])*(target_point-time_array[i])*(target_point-time_array[i+2]))
              /((time_array[i+1]-time_array[i-1])*(time_array[i+1]-time_array[i])
                *(time_array[i+1]-time_array[i+2]));
        term4=rate_array[i+2]*((target_point-time_array[i-1])*(target_point-time_array[i])*(target_point-time_array[i+1]))
              /((time_array[i+2]-time_array[i-1])*(time_array[i+2]-time_array[i])
                *(time_array[i+2]-time_array[i+1]));
    }
    target_rate=term1+term2+term3+term4;
    return target_rate;
}




int spline_dk_maille(double yp1,double yp2,DKMaille<double> y,DKMaille<double> time_array,DKMaille<double> rate_array, int array_length)
{
    register int i;
    register int j;
    register int k;
    double p,qn,sig,un;
    DKMaille<double> u(array_length+1);
    for(j=1;j<=array_length;j++) y[j]=0.;
    for(j=1;j<=array_length;j++) u[j]=0.;

    if(yp1>0.99e30)
        y[1]=u[1]=0.00;
    else
    {
        y[1]=0.0;
        u[1]=(3./(time_array[2]-time_array[1]))*((rate_array[2]-rate_array[1])/(time_array[2]-time_array[1])-yp1);
    }

    for(i=2;i<=array_length-1;i++)
    {
        sig=(time_array[i]-time_array[i-1])/(time_array[i+1]-time_array[i-1]);
        p=sig*y[i-1]+2.;
        y[i]=(sig-1.)/p;
        u[i]=(rate_array[i+1]-rate_array[i])/(time_array[i+1]-time_array[i])-(rate_array[i]-rate_array[i-1])/(time_array[i]-time_array[i-1]);
        u[i]=(6.*u[i]/(time_array[i+1]-time_array[i-1])-sig*u[i-1])/p;
    }

    if(yp2>0.99e30)
        qn=un=0.001;
    else
    {
        qn=0.0;
        un=(3.0/(time_array[array_length]-time_array[array_length-1]))*(yp2-(rate_array[array_length]-rate_array[array_length-1])/(time_array[array_length]-time_array[array_length-1]));
    }

    y[array_length]=(un-qn*u[array_length-1])/(qn*y[array_length-1]+1.);
    for(k=array_length-1;k>=1;k--)
        y[k]=y[k]*y[k+1]+u[k];
    return 1;
}
double cubicspline_dk_maille(int nearest_node, double target_point,const DKMaille<double> &time_array,const DKMaille<double> &rate_array, int array_length)
{
    double a,b,c,d;
    double target_rate;
    double y1,y2;
    int j=0;
    DKMaille<double> y(array_length+1);
    a=(time_array[nearest_node+1]-target_point)/(time_array[nearest_node+1]-time_array[nearest_node]);
    b=1.-a;
    c=(1./6.)*(a*a*a-a)*pow((time_array[nearest_node+1]-time_array[nearest_node]),2.);
    d=(1./6.)*(b*b*b-b)*pow((time_array[nearest_node+1]-time_array[nearest_node]),2.);
    y1=(rate_array[2]-rate_array[1])/(time_array[2]-time_array[1]);
    y2=(rate_array[array_length]-rate_array[array_length-1])/(time_array[array_length]-time_array[array_length-1]);
    spline_dk_maille(y1,y2,y,time_array,rate_array,array_length);
    target_rate=a*rate_array[nearest_node]+b*rate_array[nearest_node+1]+c*y[nearest_node]+d*y[nearest_node+1];
    return target_rate;
}

int spline2(double yp1,double yp2,DKMaille<double> &y,DKMaille<double> time_array,DKMaille<double> rate_array)
{
    register int i;
    register int j;
    register int k;
    double p,qn,sig,un;
    int array_length=time_array.entries();
    DKMaille<double> u(array_length-1);
    for(j=0;j<array_length;j++) y.at(j)=0.;
    for(j=0;j<array_length;j++) u.at(j)=0.;

    if(yp1>0.99e30)
        y.at(0)=u.at(0)=0.00;
    else
    {
        y.at(0)=0.0;
        u.at(0)=(3./(time_array.at(1)-time_array.at(0)))*((rate_array.at(1)-rate_array.at(0))/(time_array.at(1)-time_array.at(0))-yp1);
    }

    for(i=1;i<array_length-1;i++)
    {
        sig=(time_array.at(i)-time_array.at(i-1))/(time_array.at(i+1)-time_array.at(i-1));
        p=sig*y.at(i-1)+2.;
        y.at(i)=(sig-1.)/p;
        u.at(i)=(rate_array.at(i+1)-rate_array.at(i))/(time_array.at(i+1)-time_array.at(i))-(rate_array.at(i)-rate_array.at(i-1))/(time_array.at(i)-time_array.at(i-1));
        u.at(i)=(6.*u.at(i)/(time_array.at(i+1)-time_array.at(i-1))-sig*u.at(i-1))/p;
    }

    if(yp2>0.99e30)
        qn=un=0.001;
    else
    {
        qn=0.0;
        un=(3.0/(time_array.at(array_length-1)-time_array.at(array_length-2)))*(yp2-(rate_array.at(array_length-1)-
                rate_array.at(array_length-2))/(time_array.at(array_length-1)-time_array.at(array_length-2)));
    }

    y.at(array_length-1)=(un-qn*u.at(array_length-2))/(qn*y.at(array_length-2)+1.);
    for(k=array_length-2;k>=0;k--)
        y.at(k)=y.at(k)*y.at(k+1)+u.at(k);
    return 1;
}

double cubicspline2(int nearest_node, double target_point,const DKMaille<double> &time_array,const DKMaille<double> &rate_array)
{
    int array_length=time_array.entries();
    double a,b,c,d;
    double target_rate;
    double y1,y2;
    int j=0;
    DKMaille<double> y(array_length);
    a=(time_array[nearest_node+1]-target_point)/(time_array[nearest_node+1]-time_array[nearest_node]);
    b=1.-a;
    c=(1./6.)*(a*a*a-a)*pow((time_array[nearest_node+1]-time_array[nearest_node]),2.);
    d=(1./6.)*(b*b*b-b)*pow((time_array[nearest_node+1]-time_array[nearest_node]),2.);
    y1=(rate_array[1]-rate_array[0])/(time_array[1]-time_array[0]);
    y2=(rate_array[array_length-1]-rate_array[array_length-2])/(time_array[array_length-1]-time_array[array_length-2]);
    spline2(y1,y2,y,time_array,rate_array);
    target_rate=a*rate_array[nearest_node]+b*rate_array[nearest_node+1]+c*y[nearest_node]+d*y[nearest_node+1];
    return target_rate;
}

double rateinterpolation_dk_maille(int method,
                                   double target_point,
                                   const DKMaille<double> &time_array,
                                   const DKMaille<double> &rate_array,
                                   int array_length)
{
    register int i;
    int nearest_node;
    double target_rate;

    if(target_point<0.)
    {
        throw("ERROR: Negative Maturity?");
        return 0.0;
    }

    else if((target_point<time_array[1])&&(target_point>=0.))
    {
        target_rate=rate_array[1];
        return target_rate;
    }

    else if(target_point>time_array[array_length])
    {
        target_rate=rate_array[array_length];
        return target_rate;
    }

    for(i=1;i<=array_length;i++)
    {
        if(target_point==time_array[i])
        {
            target_rate=rate_array[i];
            return target_rate;
        }
    }

    ordering_dk_maille(target_point,&nearest_node,time_array,array_length);
    if(method==2) target_rate=interpollinear_dk_maille(nearest_node,target_point,time_array,rate_array);
    if(method==1) target_rate=cubicspline_dk_maille(nearest_node,target_point,time_array,rate_array,array_length);
    if(method==0) target_rate=interpol_dk_maille(nearest_node,target_point,time_array,rate_array,array_length);
    return target_rate;
}

double interpolation2(int method,
                      double target_point,
                      const DKMaille<double> &time_array,
                      const DKMaille<double> &rate_array)
{

    int array_length=time_array.entries();
    if(time_array.entries()!=rate_array.entries()) throw("Inconsistent arrays");

    register int i;
    int nearest_node;
    double target_rate;

    if((target_point<time_array[0]))
    {
        target_rate=rate_array[0];
        return target_rate;
    }

    else if(target_point>time_array[array_length-1])
    {
        target_rate=rate_array[array_length-1];
        return target_rate;
    }

    for(i=0;i<array_length;i++)
    {
        if(target_point==time_array[i])
        {
            target_rate=rate_array[i];
            return target_rate;
        }
    }

    ordering2(target_point,&nearest_node,time_array,array_length);
    if(method==2) target_rate=interpollinear_dk_maille(nearest_node,target_point,time_array,rate_array);
    if(method==1) target_rate=cubicspline2(nearest_node,target_point,time_array,rate_array);
    if(method==0) target_rate=interpol2(nearest_node,target_point,time_array,rate_array);
    return target_rate;
}

void AppendDates(const DKMaille<double> &dNoticeDates,
                 const DKMaille<double> &dSwaptionDates,
                 const DKMaille<double> &dBasisSwaptionDates,
                 DKMaille<double> &dSliceDates)
{

    for(unsigned int uiM=0;uiM<dNoticeDates.entries();uiM++)
    {
        if(dNoticeDates.at(uiM)-dSwaptionDates.at(0)<-1.e-9) dSliceDates.insert(dNoticeDates.at(uiM));
    }
    for(uiM=0;uiM<dSwaptionDates.entries();uiM++)
    {
        dSliceDates.insert(dSwaptionDates.at(uiM));
    }
    for(uiM=0;uiM<dBasisSwaptionDates.entries();uiM++)
    {
        if(dBasisSwaptionDates.at(uiM)!=dSwaptionDates.at(uiM+1)) throw("Unauthorised funding leg in analytical bootstrapping");
    }
}

void AppendFXDates(const DKMaille<double> &dNoticeDates,
                   double dOptionMaturityDate,
                   DKMaille<double> &dSliceDates)
{

    for(unsigned int uiM=0;uiM<dNoticeDates.entries();uiM++)
    {
        if(dNoticeDates.at(uiM)-dOptionMaturityDate<-1.e-9) dSliceDates.insert(dNoticeDates.at(uiM));
    }
    dSliceDates.insert(dOptionMaturityDate);

}


double LinearInterpolation(double target_point,const DKMaille<double> &x,const DKMaille<double> &y)
{
    register int i;
    int nearest_node;
    double target_rate;
    if ( target_point < x[0] )
    {
        target_rate = y[0];
        return target_rate;
    }
    else if ( target_point > x[x.entries()-1] )
    {
        target_rate = y[x.entries()-1];
        return target_rate;
    }
    else
    {
        for (i=0; i<x.entries(); i++)
        {
            if ( target_point == x[i] )
            {
                target_rate = y[i];
                return target_rate;
            }
        }
        ordering(target_point, &nearest_node, x);
        target_rate = interpollinear_dk_maille(nearest_node, target_point, x, y);
        return target_rate;
    }
}

// LatticeJacobiTransformation picked up from numerical recipes
// utilised in DK 5F equity smile model
void LatticeJacobiTransformation(DKMaille2D<double> &a,
                                 int n,
                                 DKMaille<double> &d,
                                 DKMaille2D<double> &v,
                                 int *nrot)
{
    int j,iq,ip,i;
    double tresh,theta,t,sm,h,g,c;
    // Globals g_tau,g_s

    DKMaille<double> b(n);
    DKMaille<double> z(n);

    for(ip=0;ip<n;ip++)
    {
        for(iq=0;iq<n;iq++) v.at(ip,iq)=0.;
        v.at(ip,ip)=1.;
    }

    for(ip=0;ip<n;ip++)
    {
        b.at(ip)=d.at(ip)=a.at(ip,ip);
        z.at(ip)=0.;
    }

    *nrot=0;

    for(i=1;i<51;i++)
    {
        sm=0.;
        for(ip=0;ip<n-1;ip++) // Sum off diagonal elements
        {
            for(iq=ip+1;iq<n;iq++) sm+=fabs(a.at(ip,iq));
        }
        if(sm==0.)
        {
            return;
        }

        if(i<4)
            tresh=0.2*sm/(9);
        else
            tresh=0.0;

        for(ip=0;ip<n-1;ip++)
        {
            for(iq=ip+1;iq<n;iq++)
            {
                g=100*fabs(a.at(ip,iq));
                if(i>4 && (double)(fabs(d.at(ip))+g) == (double)fabs(d.at(ip))
                        && (double)(fabs(d.at(iq))+g) == (double)fabs(d.at(iq)) )
                    a.at(ip,iq)=0.;
                else if (fabs(a.at(ip,iq))>tresh)
                {
                    h=d.at(iq)-d.at(ip);
                    if((double)(fabs(h)+g) == (double)fabs(h))
                        t=a.at(ip,iq)/h;
                    else
                    {
                        theta=0.5*h/a.at(ip,iq);
                        t=1./(fabs(theta)+sqrt(1.+theta*theta));
                        if(theta<0.0) t=-t;
                    }
                    c=1./sqrt(1.+t*t);
                    g_s=t*c;
                    g_tau=g_s/(1.+c);
                    h=t*a.at(ip,iq);
                    z.at(ip) -= h;
                    z.at(iq) += h;
                    d.at(ip) -= h;
                    d.at(iq) += h;
                    a.at(ip,iq) = 0.0;
                    for(j=0;j<=ip-1;j++)
                    {
                        DEBUG_ROTATION(a,j,ip,j,iq)
                    }
                    for(j=ip+1;j<=iq-1;j++)
                    {
                        DEBUG_ROTATION(a,ip,j,j,iq)
                    }
                    for(j=iq+1;j<n;j++)
                    {
                        DEBUG_ROTATION(a,ip,j,iq,j)
                    }
                    for(j=0;j<n;j++)
                    {
                        DEBUG_ROTATION(v,j,ip,j,iq)
                    }
                    ++(*nrot);
                }
            }
        }
        for(ip=0;ip<n;ip++)
        {
            b.at(ip) += z.at(ip);
            d.at(ip) = b.at(ip);
            z.at(ip) = 0.0;
        }
    }
    throw("Too many Jacobi iterations");
} // LatticeJacobiTransformation()

void EigenSort(DKMaille<double> &d,
               DKMaille2D<double> &v,
               int n)
{
    int k,j,i;
    double p;
    for(i=0;i<n-1;i++)
    {
        p=d.at(k=i);
        for(j=i+1;j<=n-1;j++)
            if(d.at(j)>=p) p=d.at(k=j);
        if(k!=i)
        {
            d.at(k)=d.at(i);
            d.at(i)=p;
            for(j=0;j<=n-1;j++)
            {
                p=v.at(j,i);
                v.at(j,i)=v.at(j,k);
                v.at(j,k)=p;
            }
        }
    }
}

void multiple(double dCoupon, DKMaille<double> *Payoff)
{
    static unsigned int ii;
    for(ii=0;ii<Payoff->entries();ii++)
        Payoff->at(ii)*=dCoupon;
}

DKMaille<double> Reverse(DKMaille<double> &Input)
{
    unsigned int uiSize=Input.entries();
    static DKMaille<double> dRet(uiSize);
    dRet.resize(uiSize);
    for(unsigned int ii=0;ii<uiSize;ii++)
    {
        dRet.at(ii)=Input.at(uiSize-ii-1);
    }
    return dRet;
}

double dItoIntegral_HW1F_TD(double dMeanFieldDecay,
                            double dTerminalDate,
                            const DKMaille<double> &dDatesStrip,
                            const DKMaille<double> &dSigmaStrip)
{
    int i=0;
    double dSum=0.;
    while(dDatesStrip[i]<dTerminalDate-0.00001)
    {
        dSum+=(dSigmaStrip[i]/dMeanFieldDecay)*(
                  (dDatesStrip[i+1]-dDatesStrip[i])
                  -(exp(-dMeanFieldDecay*(dTerminalDate-dDatesStrip[i+1]))
                    -exp(-dMeanFieldDecay*(dTerminalDate-dDatesStrip[i])))/(dMeanFieldDecay)
              );
        i++;
    }
    return dSum;
} // dItoIntegral_HW1F_TD

double dItoCrossIntegral_HW1F_TD(double dMeanFieldDecay,
                                 double dTerminalDate,
                                 const DKMaille<double> &dDatesStrip,
                                 const DKMaille<double> &dSigmaStrip,
                                 const DKMaille<double> &dFXSigmaStrip)
{
    int i=0;
    double dSum=0.;
    while(dDatesStrip[i]<dTerminalDate-0.00001)
    {
        dSum+=dFXSigmaStrip[i]*(dSigmaStrip[i]/dMeanFieldDecay)*(
                  (dDatesStrip[i+1]-dDatesStrip[i])
                  -(exp(-dMeanFieldDecay*(dTerminalDate-dDatesStrip[i+1]))
                    -exp(-dMeanFieldDecay*(dTerminalDate-dDatesStrip[i])))/(dMeanFieldDecay)
              );
        i++;
    }
    return dSum;
} // dItoIntegral_HW1F_TD

double dKacDeterminant2F_TD(double dMeanFieldDecay, double dMeanFieldDecay2,
                            double dObservationDate, double dTerminalDate,
                            double dSigma1, double dSigma2, double dDate1,
                            double dDate2)
{
    double dTerm1 = (dTerminalDate-dObservationDate);
    double dTerm2 = -(1. /dMeanFieldDecay)*(exp(-dMeanFieldDecay*(dDate1-dTerminalDate))
                                            -exp(-dMeanFieldDecay*(dDate1-dObservationDate)));
    double dTerm3 = -(1. /dMeanFieldDecay2)*(exp(-dMeanFieldDecay2*(dDate2-dTerminalDate))
                    -exp(-dMeanFieldDecay2*(dDate2-dObservationDate)));
    double dTerm4 = +(1./(dMeanFieldDecay+dMeanFieldDecay2))
                    *(
                        exp(-dMeanFieldDecay*(dDate1-dTerminalDate)-dMeanFieldDecay2*(dDate2-dTerminalDate))
                        -exp(-dMeanFieldDecay*(dDate1-dObservationDate)-dMeanFieldDecay2*(dDate2-dObservationDate))
                    );
    return (dTerm1+dTerm2+dTerm3+dTerm4)*dSigma1*dSigma2;
} // double dKacDeterminant2F_TD

double dKacDeterminant2F_TD_DK(double dMeanFieldDecay,
                               double dMeanFieldDecay2,
                               double dObservationDate,
                               double dTerminalDate,
                               double dDate1,
                               double dDate2,
                               const DKMaille<double> &dSlices,
                               const DKMaille<double> &dShortRateOnSlice,
                               const DKMaille<double> &dShortRateOnSlice2)
{
    // Slices run through dObservationDate,dTerminalDate and dDate1,dDate2 by construction
    // Get slice indices and short rates
    unsigned int ui1=0; unsigned int ui2=0; unsigned int ui3=0; unsigned int ui4=0;
    for(unsigned int ui=0;ui<dSlices.entries();ui++)
    {
        if(fabs(dSlices.at(ui)-dObservationDate)<1.e-9) ui1=ui;
        if(fabs(dSlices.at(ui)-dTerminalDate)<1.e-9) ui2=ui;
        if(fabs(dSlices.at(ui)-dDate1)<1.e-9) ui3=ui;
        if(fabs(dSlices.at(ui)-dDate2)<1.e-9) ui4=ui;
    }
    double dVolatility=0.;
    double dVolatilityTest=0.;
    double dTerm1=0; double dTerm2=0.;
    for(unsigned int uiK=ui1;uiK<ui2;uiK++)
    {
        dTerm1=0.;
        dTerm2=0.;
        for(unsigned int uiM=uiK;uiM<ui3;uiM++)
        {
            dTerm1+=dShortRateOnSlice.at(uiM)*(exp(-dMeanFieldDecay*(dSlices.at(uiM)-dSlices.at(uiK)))-exp(-dMeanFieldDecay*(dSlices.at(uiM+1)-dSlices.at(uiK))));
        }
        for(unsigned int uiL=uiK;uiL<ui4;uiL++)
        {
            dTerm2+=dShortRateOnSlice2.at(uiL)*(exp(-dMeanFieldDecay2*(dSlices.at(uiL)-dSlices.at(uiK)))-exp(-dMeanFieldDecay2*(dSlices.at(uiL+1)-dSlices.at(uiK))));
        }
        dVolatility+=(dSlices.at(uiK+1)-dSlices.at(uiK))*dTerm1*dTerm2;
    }
    return dVolatility;
} // double dKacDeterminant2F_TD_DK



double isConstant(DKMaille2D<double> &dAbsVol)
{
    // return dAbsVol if the surface is constant, return 0 if not
    unsigned int uC=0;
    double dPrevious=0.;
    for(unsigned int ui=0;ui<dAbsVol.rows();ui++)
    {
        for(unsigned int uk=0;uk<dAbsVol.columns();uk++)
        {
            double dCurrent=dAbsVol.at(ui,uk);
            if(!(ui==0&&uk==0))
            {
                if(dPrevious-dAbsVol.at(ui,uk)!=0.) return 0.;
            }
            dPrevious=dAbsVol.at(ui,uk);
        }
    }
    return dPrevious;
} // double isConstant

double dItoIntegral_HW2F_TD(double dMeanFieldDecay,
                            double dMeanFieldDecay2,
                            double dObservationDate,
                            double dTerminalDate,
                            const DKMaille<double> &dDatesStrip,
                            const DKMaille<double> &dSigmaStrip,
                            const DKMaille<double> &dSigmaStrip2)
{
    int i=0;
    double dSum=0.;

    while(dDatesStrip[i]<dTerminalDate-0.00001)
    {
        double dDate1=dTerminalDate;
        double dDate2=dTerminalDate;
        dSum+=dKacDeterminant2F_TD(dMeanFieldDecay,dMeanFieldDecay2,dDatesStrip[i],dDatesStrip[i+1],
                                   dSigmaStrip[i],dSigmaStrip2[i],dDate1,dDate2);
        i++;
    }
    return dSum;
} // dItoIntegral_HW2F_TD

DKMaille < double >&dCheckDates(DKMaille < double >&dStrip,
                                double dStartDate,
                                double dEndDate)
{
    static DKMaille < double >dIntegrationStrip(0);

    dIntegrationStrip.resize(0);
    // Modification pour ne pas rejeter une borne quand c'est seulement un probleme de precision des doubles
    // A revoir de maniere plus globale
    if (dStartDate < dStrip.at(0) || dEndDate > dStrip.at(dStrip.entries() - 1) + 1.e-10)
        throw("Integration Dates Are Out Of Bounds");
    if (dStartDate > dStrip.at(dStrip.entries() - 1) || dEndDate < dStrip.at(0))
        throw("Integration Dates Are Out Of Bounds");
    if (dStartDate > dEndDate)
        throw("Start Date after End Date");
    unsigned int ui = 0;

    dIntegrationStrip.insert(dStartDate);
    while (ui < dStrip.entries() && dStrip.at(ui) < dEndDate)
    {
        if (dStrip.at(ui) > dStartDate && dStrip.at(ui) < dEndDate)
            dIntegrationStrip.insert(dStrip.at(ui));
        ui++;
    }
    dIntegrationStrip.insert(dEndDate);
    return dIntegrationStrip;
}

DKMaille<double> &dNewVols(DKMaille<double> &dStrip,DKMaille<double> &dVols,DKMaille<double> &dNewStrip)
{
    static DKMaille<double> dRet;
    dRet.resize(dNewStrip.entries());
    for(unsigned int ui=0;ui<dNewStrip.entries()-1;ui++)
    {
        for(unsigned int um=0;um<dStrip.entries()-1;um++)
        {
            if(dNewStrip.at(ui)>=dStrip.at(um)&&dNewStrip.at(ui)<dStrip.at(um+1))
                dRet.at(ui)=dVols.at(um);
        }
    }
    dRet.at(dNewStrip.entries()-1)=dRet.at(dNewStrip.entries()-2);
    return dRet;
}

double SpotFXVolIntegral(DKMaille<double> &dStrip,
                         DKMaille<double> &dSpotFXVol,
                         double dStartDate,
                         double dEndDate)
{
    double dSum=0.;

    // Checks on dates
    DKMaille<double> dIntegrationStrip=dCheckDates(dStrip,dStartDate,dEndDate);
    DKMaille<double> dIntegrationStripVols=dNewVols(dStrip,dSpotFXVol,dIntegrationStrip);

    for(unsigned ui=0;ui<dIntegrationStrip.entries()-1;ui++)
    {
        dSum+=(dIntegrationStripVols[ui]*dIntegrationStripVols[ui])*(dIntegrationStrip[ui+1]-dIntegrationStrip[ui]);
    }
    return dSum;

} // SpotFXVolIntegral

double IRStdDevIntegral(DKMaille<double> &dStrip,
                        DKMaille<double> &dStripIR1StdDev,
                        DKMaille<double> &dStripIR2StdDev,
                        double dMeanReversion1,
                        double dMeanReversion2,
                        double dStartDate,
                        double dEndDate,
                        double dForwardMaturityDate)
{

    if(dEndDate>dForwardMaturityDate)
        throw("Integration dates not correct");

    double dSum=0.;

    // Checks on dates
    DKMaille<double> dIntegrationStrip=dCheckDates(dStrip,dStartDate,dEndDate);
    DKMaille<double> dIntegrationStripVols1=dNewVols(dStrip,dStripIR1StdDev,dIntegrationStrip);
    DKMaille<double> dIntegrationStripVols2=dNewVols(dStrip,dStripIR2StdDev,dIntegrationStrip);


    for(unsigned ui=0;ui<dIntegrationStrip.entries()-1;ui++)
    {
        dSum+=(1./(dMeanReversion1*dMeanReversion2))*dKacDeterminant2F_TD(dMeanReversion1,
                dMeanReversion2,
                dIntegrationStrip.at(ui),
                dIntegrationStrip.at(ui+1),
                dIntegrationStripVols1[ui],
                dIntegrationStripVols2[ui],
                dForwardMaturityDate,
                dForwardMaturityDate);

    }
    return dSum;
} // IRStdDevIntegral

double IRStdDevIntegral_DK(DKMaille<double> &dStrip,
                           DKMaille<double> &dStripIR1StdDev,
                           DKMaille<double> &dStripIR2StdDev,
                           double dMeanReversion1,
                           double dMeanReversion2,
                           double dStartDate,
                           double dEndDate,
                           double dForwardMaturityDate,
                           const DKMaille<double> &shortRateDates,
                           const DKMaille<double> &shortRates,
                           const DKMaille<double> &shortRates2)
{

    if(dEndDate>dForwardMaturityDate)
        throw("Integration dates not correct");

    double dSum=0.;

    // Checks on dates
    DKMaille<double> dIntegrationStrip=dCheckDates(dStrip,dStartDate,dEndDate);
    DKMaille<double> dIntegrationStripVols1=dNewVols(dStrip,dStripIR1StdDev,dIntegrationStrip);
    DKMaille<double> dIntegrationStripVols2=dNewVols(dStrip,dStripIR2StdDev,dIntegrationStrip);
    // Integration slicing has to include notice dates, and bond maturity dates
    DKMaille<double> dDatesForSlices;
    DKMaille<double> dSlices;
    DKMaille<double> dShortRateOnSlices;
    DKMaille<double> dShortRateOnSlices2;
    AppendFXDates(dIntegrationStrip,dEndDate,dDatesForSlices);
    GetXCCYSlices(shortRateDates,shortRates,shortRates2,dDatesForSlices,dSlices,dShortRateOnSlices,dShortRateOnSlices2,0.5);

    for(unsigned ui=0;ui<dIntegrationStrip.entries()-1;ui++)
    {
        dSum+=(1./(dMeanReversion1*dMeanReversion2))*dKacDeterminant2F_TD_DK(dMeanReversion1,
                dMeanReversion2,
                dIntegrationStrip.at(ui),
                dIntegrationStrip.at(ui+1),
                dForwardMaturityDate,
                dForwardMaturityDate,
                dSlices,
                dShortRateOnSlices,
                dShortRateOnSlices2)*dIntegrationStripVols1[ui]*dIntegrationStripVols2[ui];

    }
    return dSum;
} // IRStdDevIntegral_DK






double SpotFXVolIRStdDevCrossIntegral(DKMaille<double> &dStrip,
                                      DKMaille<double> &dSpotFXVol,
                                      DKMaille<double> &dStripIRStdDev,
                                      double dMeanReversion,
                                      double dStartDate,
                                      double dEndDate,
                                      double dForwardMaturityDate)
{

    if(dEndDate>dForwardMaturityDate)
        throw("Integration dates not correct");

    double dSum=0.;

    // Checks on dates
    DKMaille<double> dIntegrationStrip=dCheckDates(dStrip,dStartDate,dEndDate);
    DKMaille<double> dIntegrationStripVols1=dNewVols(dStrip,dSpotFXVol,dIntegrationStrip);
    DKMaille<double> dIntegrationStripVols2=dNewVols(dStrip,dStripIRStdDev,dIntegrationStrip);

    for(unsigned ui=0;ui<dIntegrationStrip.entries()-1;ui++)
    {
        dSum+=dIntegrationStripVols1[ui]*(dIntegrationStripVols2[ui]/dMeanReversion)*(
                  (dIntegrationStrip[ui+1]-dIntegrationStrip[ui])
                  -(exp(-dMeanReversion*(dForwardMaturityDate-dIntegrationStrip[ui+1]))
                    -exp(-dMeanReversion*(dForwardMaturityDate-dIntegrationStrip[ui])))/(dMeanReversion)
              );
    }
    return dSum;
} // SpotFXVolIRStdDevCrossIntegral

double dKacDeterminant1F_TD_DK(double dMeanFieldDecay,
                               double dObservationDate,
                               double dTerminalDate,
                               double dDate1,
                               const DKMaille<double> &dSlices,
                               const DKMaille<double> &dShortRateOnSlice)
{
    // Slices run through dObservationDate,dTerminalDate and dDate1 by construction
    // Get slice indices and short rates
    unsigned int ui1=0; unsigned int ui2=0; unsigned int ui3=0;
    for(unsigned int ui=0;ui<dSlices.entries();ui++)
    {
        if(fabs(dSlices.at(ui)-dObservationDate)<1.e-9) ui1=ui;
        if(fabs(dSlices.at(ui)-dTerminalDate)<1.e-9) ui2=ui;
        if(fabs(dSlices.at(ui)-dDate1)<1.e-9) ui3=ui;
    }
    double dVolatility=0.;
    double dTerm1=0; double dTerm2=0.;
    for(unsigned int uiK=ui1;uiK<ui2;uiK++)
    {
        dTerm1=0.;
        for(unsigned int uiM=uiK;uiM<ui3;uiM++)
        {
            dTerm1+=dShortRateOnSlice.at(uiM)*(exp(-dMeanFieldDecay*(dSlices.at(uiM)-dSlices.at(uiK)))-exp(-dMeanFieldDecay*(dSlices.at(uiM+1)-dSlices.at(uiK))));
        }
        dVolatility+=(dSlices.at(uiK+1)-dSlices.at(uiK))*dTerm1;
    }
    return dVolatility;
} // double dKacDeterminant1F_TD_DK

double SpotFXVolIRStdDevCrossIntegral_DK(DKMaille<double> &dStrip,
        DKMaille<double> &dSpotFXVol,
        DKMaille<double> &dStripIRStdDev,
        double dMeanReversion,
        double dStartDate,
        double dEndDate,
        double dForwardMaturityDate,
        const DKMaille<double> &shortRateDates,
        const DKMaille<double> &shortRates)
{

    if(dEndDate>dForwardMaturityDate)
        throw("Integration dates not correct");

    double dSum=0.;

    // Checks on dates
    DKMaille<double> dIntegrationStrip=dCheckDates(dStrip,dStartDate,dEndDate);
    DKMaille<double> dIntegrationStripVols1=dNewVols(dStrip,dSpotFXVol,dIntegrationStrip);
    DKMaille<double> dIntegrationStripVols2=dNewVols(dStrip,dStripIRStdDev,dIntegrationStrip);
    // Integration slicing has to include notice dates, and bond maturity dates
    DKMaille<double> dDatesForSlices;
    DKMaille<double> dSlices;
    DKMaille<double> dShortRateOnSlices;
    AppendFXDates(dIntegrationStrip,dEndDate,dDatesForSlices);
    GetSlices(shortRateDates,shortRates,dDatesForSlices,dSlices,dShortRateOnSlices,0.5);

    for(unsigned ui=0;ui<dIntegrationStrip.entries()-1;ui++)
    {
        dSum+=dIntegrationStripVols1[ui]*(dIntegrationStripVols2[ui]/dMeanReversion)*
              dKacDeterminant1F_TD_DK(dMeanReversion,
                                      dIntegrationStrip.at(ui),
                                      dIntegrationStrip.at(ui+1),
                                      dForwardMaturityDate,
                                      dSlices,
                                      dShortRateOnSlices);
    }
    return dSum;
} // SpotFXVolIRStdDevCrossIntegral_DK



double SpotFXVolIRStdDevQuantoIntegral(DKMaille<double> &dStrip,
                                       DKMaille<double> &dSpotFXVol,
                                       DKMaille<double> &dStripIRStdDev,
                                       double dMeanReversion,
                                       double dStartDate,
                                       double dEndDate)
{

    if(dStartDate<dStrip.at(0))
        throw("Integration dates not correct");
    double dSum=0.;

    // Checks on dates
    DKMaille<double> dIntegrationStrip=dCheckDates(dStrip,dStartDate,dEndDate);
    DKMaille<double> dIntegrationStripVols1=dNewVols(dStrip,dSpotFXVol,dIntegrationStrip);
    DKMaille<double> dIntegrationStripVols2=dNewVols(dStrip,dStripIRStdDev,dIntegrationStrip);

    for(unsigned ui=0;ui<dIntegrationStrip.entries()-1;ui++)
    {
        dSum+=dIntegrationStripVols1[ui]*(dIntegrationStripVols2[ui]/dMeanReversion)*
              (exp(-dMeanReversion*dIntegrationStrip[ui])-exp(-dMeanReversion*dIntegrationStrip[ui+1]));
    }
    return dSum;
} // SpotFXVolIRStdDevQuantoIntegral

double GetForwardFXVol(double dSpotDate,
                       double dStartDate,
                       double dEndDate,
                       double dForwardMaturityDate,
                       DKMaille<double> dStrip,
                       DKMaille<double> &dStripDomesticStdDev,
                       DKMaille<double> &dStripForeignStdDev,
                       DKMaille<double> &dStripSpotFXVol,
                       double dDomesticMeanReversion,
                       double dForeignMeanReversion,
                       double dCorrelationSpotFXDomestic,
                       double dCorrelationSpotFXForeign,
                       double dCorrelationForeignDomestic)
{


    // Adjust to ModelDates
    if(dSpotDate!=0.)
    {
        for(unsigned int ui=0;ui<dStrip.entries();ui++)
        {
            dStrip.at(ui)=(dStrip.at(ui)-dSpotDate)/365.;
        }
        dStartDate=(dStartDate-dSpotDate)/365.;
        dEndDate=(dEndDate-dSpotDate)/365.;
        dForwardMaturityDate=(dForwardMaturityDate-dSpotDate)/365.;
    }
    else
    {
        // do nothing
    }
    double dTerm1=SpotFXVolIntegral(dStrip,dStripSpotFXVol,dStartDate,dEndDate);
    double dTerm2=IRStdDevIntegral(dStrip,dStripForeignStdDev,dStripForeignStdDev,dForeignMeanReversion,dForeignMeanReversion,dStartDate,dEndDate,dForwardMaturityDate);
    double dTerm3=IRStdDevIntegral(dStrip,dStripDomesticStdDev,dStripDomesticStdDev,dDomesticMeanReversion,dDomesticMeanReversion,dStartDate,dEndDate,dForwardMaturityDate);
    double dTerm4=-2.*dCorrelationSpotFXForeign*
                  SpotFXVolIRStdDevCrossIntegral(dStrip,dStripSpotFXVol,dStripForeignStdDev,dForeignMeanReversion,dStartDate,dEndDate,dForwardMaturityDate);
    double dTerm5=+2.*dCorrelationSpotFXDomestic*
                  SpotFXVolIRStdDevCrossIntegral(dStrip,dStripSpotFXVol,dStripDomesticStdDev,dDomesticMeanReversion,dStartDate,dEndDate,dForwardMaturityDate);

    double dTerm6=-2.*dCorrelationForeignDomestic*
                  IRStdDevIntegral(dStrip,dStripDomesticStdDev,dStripForeignStdDev,dDomesticMeanReversion,dForeignMeanReversion,dStartDate,dEndDate,dForwardMaturityDate);
    double dResult=sqrt((dTerm1+dTerm2+dTerm3+dTerm4+dTerm5+dTerm6)/(dEndDate-dStartDate));
    return dResult;

}

void GetShortRateOnSlice(const DKMaille<double> &shortRateDates,
                         const DKMaille<double> &shortRates,
                         const DKMaille<double> &dSlices,
                         DKMaille<double> &dShortRateOnSlices,
                         DKMaille<double> &dShortRateKOnSlices)
{
    for(unsigned ui=0;ui<dSlices.entries();ui++)
    {
        double dShortRate=rateinterpolation_dk_maille(2,dSlices.at(ui),shortRateDates,shortRates,shortRateDates.entries()-1);
        dShortRateOnSlices.insert(dShortRate);
        if(dShortRate>0.002)
            dShortRateKOnSlices.insert(dShortRate);
        else
            dShortRateKOnSlices.insert(0.002);
    }
}


void GetSlices(const DKMaille<double> &shortRateDates,
               const DKMaille<double> &shortRates,
               const DKMaille<double> &dInputDates,
               DKMaille<double> &dSlices,
               DKMaille<double> &dShortRateOnSlice,
               double dt)
{
    dSlices.insert(dInputDates.at(0));
    unsigned int uiSlice=0;
    double dTestSlice=dSlices.at(0)+dt;
    for(unsigned int ui=0;ui<dInputDates.entries()-1;ui++)
    {
        while(dInputDates.at(ui+1)-dTestSlice>1.e-9)
        {
            dSlices.insert(dTestSlice);
            uiSlice++;
            dTestSlice=dSlices.at(uiSlice)+dt;
        }
        dSlices.insert(dInputDates.at(ui+1));
        uiSlice++;
        dTestSlice=dSlices.at(uiSlice)+dt;
    }
    for(ui=0;ui<dSlices.entries();ui++)
    {
        dShortRateOnSlice.insert(rateinterpolation_dk_maille(2,dSlices.at(ui),shortRateDates,shortRates,shortRateDates.entries()-1));
    }
}


void GetXCCYSlices(const DKMaille<double> &shortRateDates,
                   const DKMaille<double> &shortRates,
                   const DKMaille<double> &shortRatesForeign,
                   const DKMaille<double> &dInputDates,
                   DKMaille<double> &dSlices,
                   DKMaille<double> &dShortRateOnSlice,
                   DKMaille<double> &dShortRateForeignOnSlice,
                   double dt)
{
    dSlices.insert(dInputDates.at(0));
    unsigned int uiSlice=0;
    double dTestSlice=dSlices.at(0)+dt;
    for(unsigned int ui=0;ui<dInputDates.entries()-1;ui++)
    {
        while(dInputDates.at(ui+1)-dTestSlice>1.e-9)
        {
            dSlices.insert(dTestSlice);
            uiSlice++;
            dTestSlice=dSlices.at(uiSlice)+dt;
        }
        dSlices.insert(dInputDates.at(ui+1));
        uiSlice++;
        dTestSlice=dSlices.at(uiSlice)+dt;
    }
    for(ui=0;ui<dSlices.entries();ui++)
    {
        dShortRateOnSlice.insert(rateinterpolation_dk_maille(2,dSlices.at(ui),shortRateDates,shortRates,shortRateDates.entries()-1));
        dShortRateForeignOnSlice.insert(rateinterpolation_dk_maille(2,dSlices.at(ui),shortRateDates,shortRatesForeign,shortRateDates.entries()-1));
    }
}



double GetForwardFXVol_DK(double dSpotDate,
                          double dStartDate,
                          double dEndDate,
                          double dForwardMaturityDate,
                          DKMaille<double> dStrip,
                          DKMaille<double> &dStripDomesticStdDev,
                          DKMaille<double> &dStripForeignStdDev,
                          DKMaille<double> &dStripSpotFXVol,
                          double dDomesticMeanReversion,
                          double dForeignMeanReversion,
                          double dCorrelationSpotFXDomestic,
                          double dCorrelationSpotFXForeign,
                          double dCorrelationForeignDomestic,
                          const DKMaille<double> &shortRateDates,
                          const DKMaille<double> &shortRates,
                          const DKMaille<double> &shortRatesForeign)
{
    // Adjust to ModelDates
    if(dSpotDate!=0.)
    {
        for(unsigned int ui=0;ui<dStrip.entries();ui++)
        {
            dStrip.at(ui)=(dStrip.at(ui)-dSpotDate)/365.;
        }
        dStartDate=(dStartDate-dSpotDate)/365.;
        dEndDate=(dEndDate-dSpotDate)/365.;
        dForwardMaturityDate=(dForwardMaturityDate-dSpotDate)/365.;
    }
    else
    {
        // do nothing
    }

    double dTerm1=SpotFXVolIntegral(dStrip,dStripSpotFXVol,dStartDate,dEndDate);
    double dTerm2=IRStdDevIntegral_DK(dStrip,dStripForeignStdDev,dStripForeignStdDev,dForeignMeanReversion,dForeignMeanReversion,dStartDate,dEndDate,dForwardMaturityDate,
                                      shortRateDates,shortRatesForeign,shortRatesForeign);
    double dTerm3=IRStdDevIntegral_DK(dStrip,dStripDomesticStdDev,dStripDomesticStdDev,dDomesticMeanReversion,dDomesticMeanReversion,dStartDate,dEndDate,dForwardMaturityDate,
                                      shortRateDates,shortRates,shortRates);
    double dTerm4=-2.*dCorrelationSpotFXForeign*
                  SpotFXVolIRStdDevCrossIntegral_DK(dStrip,dStripSpotFXVol,dStripForeignStdDev,dForeignMeanReversion,dStartDate,dEndDate,dForwardMaturityDate,
                                                    shortRateDates,shortRatesForeign);
    double dTerm5=+2.*dCorrelationSpotFXDomestic*
                  SpotFXVolIRStdDevCrossIntegral_DK(dStrip,dStripSpotFXVol,dStripDomesticStdDev,dDomesticMeanReversion,dStartDate,dEndDate,dForwardMaturityDate,
                                                    shortRateDates,shortRates);
    double dTerm6=-2.*dCorrelationForeignDomestic*
                  IRStdDevIntegral_DK(dStrip,dStripDomesticStdDev,dStripForeignStdDev,dDomesticMeanReversion,dForeignMeanReversion,dStartDate,dEndDate,dForwardMaturityDate,
                                      shortRateDates,shortRates,shortRatesForeign);
    double dResult=sqrt((dTerm1+dTerm2+dTerm3+dTerm4+dTerm5+dTerm6)/(dEndDate-dStartDate));
    return dResult;
}






DKMaille<double> CreateDateStrip(double dSpotDate,
                                 DKMaille<double> &dDomesticStrip,
                                 DKMaille<double> &dForeignStrip,
                                 DKMaille<double> &dSpotFXStrip)
{

    if((dSpotDate==0.&&dDomesticStrip.at(0)!=0.)||(dSpotDate==0.&&dForeignStrip.at(0)!=0.)||(dSpotDate==0.&&dSpotFXStrip.at(0)!=0.))
        throw("Strip dates do not start at 0. Not allowed !!!");
    static DKMaille<double> dStrip;
    dStrip.resize(0);
    unsigned int uiSize=0;
    bool bAlreadyPresent=false;
    // Create the strip points
    if(dDomesticStrip.at(0)!=0.||dSpotDate==0.)
    {
        dStrip.insert(dSpotDate);
        uiSize+=1;
    }
    for(unsigned int ui=0;ui<dDomesticStrip.entries();ui++)
    {
        bAlreadyPresent=false;
        if(dDomesticStrip.at(ui)==dSpotDate) bAlreadyPresent=true;
        if(bAlreadyPresent==false)
        {
            dStrip.insert(dDomesticStrip.at(ui));
            uiSize+=1;
        }
    }
    for(unsigned int uiK=0;uiK<dForeignStrip.entries();uiK++)
    {
        bAlreadyPresent=false;
        if(dForeignStrip.at(uiK)==dSpotDate) bAlreadyPresent=true;
        for(unsigned int uiM=0;uiM<dDomesticStrip.entries();uiM++)
        {
            if(dForeignStrip.at(uiK)==dDomesticStrip.at(uiM)) bAlreadyPresent=true;
        }
        if(bAlreadyPresent==false)
        {
            dStrip.insert(dForeignStrip.at(uiK));
            uiSize+=1;
        }
    }
    for(uiK=0;uiK<dSpotFXStrip.entries();uiK++)
    {
        bAlreadyPresent=false;
        if(dSpotFXStrip.at(uiK)==dSpotDate) bAlreadyPresent=true;
        for(unsigned int uiM=0;uiM<dDomesticStrip.entries();uiM++)
        {
            if(dSpotFXStrip.at(uiK)==dDomesticStrip.at(uiM)) bAlreadyPresent=true;
        }
        for(uiM=0;uiM<dForeignStrip.entries();uiM++)
        {
            if(dSpotFXStrip.at(uiK)==dForeignStrip.at(uiM)) bAlreadyPresent=true;
        }
        if(bAlreadyPresent==false)
        {
            dStrip.insert(dSpotFXStrip.at(uiK));
            uiSize+=1;
        }
    }
    // Sort dates
    DKMaille<double> dStrip_(uiSize);
    for(unsigned int uiM=0;uiM<uiSize;uiM++)
    {
        dStrip_.at(uiM)=dStrip.at(uiM);
    }
    dStrip.resize(uiSize);
    for(uiM=0;uiM<uiSize;uiM++)
    {
        dStrip.at(uiM)=dStrip_.at(uiM);
    }
    Sort(dStrip,dStrip.entries());
    dStrip=Reverse(dStrip);
    return dStrip;
}

void Sort(DKMaille<double> &d,
          int n)
{
    int k,j,i;
    double p;
    for(i=0;i<n-1;i++)
    {
        p=d.at(k=i);
        for(j=i+1;j<=n-1;j++)
            if(d.at(j)>=p) p=d.at(k=j);
        if(k!=i)
        {
            d.at(k)=d.at(i);
            d.at(i)=p;
        }
    }
}
DKMaille<double> dKacVector(unsigned int uiSize, double dSwapRate,
                            DKMaille<double> dPeriods, DKMaille<double> dDiscountFactors)
{
    DKMaille<double> dVector(uiSize);
    unsigned int uiIndex;
    dVector[0]=dDiscountFactors[0];
    dVector[uiSize-1]=-(1.+dPeriods[uiSize-1]*dSwapRate)*dDiscountFactors[uiSize-1];
    if(uiSize>2)
    {
        for(uiIndex=1; uiIndex<uiSize-1; uiIndex++)
        {
            dVector[uiIndex]=-dPeriods[uiIndex]*dSwapRate*dDiscountFactors[uiIndex];
        }
    }
    return dVector;
}

double dKacDeterminant(double dMeanFieldDecay, double dObservationDate, double dTerminalDate,
                       double dSigmaDate1, double dSigmaDate2, double dDate1, double dDate2)
{
    double dTerm1 = (dTerminalDate-dObservationDate);

    double dTerm2 = -(1. /dMeanFieldDecay)*(exp(-dMeanFieldDecay*(dDate1-dTerminalDate))
                                            -exp(-dMeanFieldDecay*(dDate1-dObservationDate)));

    double dTerm3 = -(1. /dMeanFieldDecay)*(exp(-dMeanFieldDecay*(dDate2-dTerminalDate))
                                            -exp(-dMeanFieldDecay*(dDate2-dObservationDate)));

    double dTerm4 = +(0.5/dMeanFieldDecay)*(exp(-dMeanFieldDecay*(dDate1+dDate2-2.*dTerminalDate))
                                            -exp(-dMeanFieldDecay*(dDate1+dDate2-2.*dObservationDate)));

    return (dTerm1+dTerm2+dTerm3+dTerm4)*dSigmaDate1*dSigmaDate2;
}


int normalDK(double *res,double xx)
{
    try
    {
        double y;
        double dd;
        if (xx < 0.0)
        {
            dd = xx*-1.0;
            normalDK(res,dd);
            *res = 1.0 - *res;
            return 0;
        }
        y = (1 + 0.231641910 * xx);
        if (y==0.0) return -1;
        y = 1. / y;
        *res = 1. - Poly5_new(y,norcoef_new) * exp(- xx * xx / 2.) / sqrt (2. * M_PI);
        return 0;

    }
    catch(...)
    {
  #ifdef TOKYO
        // Never caused problems in London !!!???!!!!????
        //if (AbsoluteVariate > 50)
        if (xx > 35)
        {
            if (xx < 0.0)  *res = 0;
            else           *res = 1;

            return 0;
        }
  #endif // TOKYO

        // Dimitri will be _SO_ cross when he sees that I've made a change in this file :-)
        std::cerr << "res=" << *res << ", xx=" << xx << std::endl;
        throw("ERROR: Exception thrown in normalDK.");
    }
}

double ATMpricevol_DK(double tte, double inputvol) {
    double probability=0.;
    double argument=0.5*inputvol*sqrt(tte);
    normalDK(&probability,argument);
    return 2.*probability-1.;
}

double ATMpricevol_derivative_DK(double tte, double inputvol) {
    double probability=0.;
    double argument=0.5*inputvol*sqrt(tte);
    double derivative=sqrt(tte);
    return derivative*(1./sqrt(2.*M_PI))*exp(-0.5*argument*argument);
}

double dBlackScholes(double dTimeToExpiry,
                     double dVol,
                     double dStrike,
                     double dForward,
                     double dCallPut)
{
    double d1=(log(dForward/dStrike)+0.5*dVol*dVol*dTimeToExpiry)/(dVol*sqrt(dTimeToExpiry));
    double d2=(log(dForward/dStrike)-0.5*dVol*dVol*dTimeToExpiry)/(dVol*sqrt(dTimeToExpiry));
    double dND1=0.;
    double dND2=0.;
    normalDK(&dND1,d1);
    normalDK(&dND2,d2);
    double dCall=dForward*dND1-dStrike*dND2;
    double dResult;
    if(dCallPut==0.) dResult=dCall;
    if(dCallPut==1.) dResult=dCall-dForward+dStrike;
    return dResult;
}


// BiLinear Interpolation routine
double BiLinearInterpolation(DKMaille<double> &dXArray,
                             DKMaille<double> &dYArray,
                             DKMaille<double> &dZArray,
                             double dX,
                             double dY)
{
    // If out of range throw message
    if( dX<dXArray[0] || dX>dXArray[1] || dY>dYArray[0] || dY<dYArray[1] )
        throw("ERROR: Out of Range in BiLinearInterpolation");

    // We only treat here the case of a square grid
    double dt = ( dX - dXArray[0] ) / ( dXArray[1] - dXArray[0] );
    double du = ( dY - dYArray[0] ) / ( dYArray[1] - dYArray[0] );
    return dZArray[0]*(1.-dt)*(1.-du)+dZArray[1]*(dt)*(1.-du)+dZArray[2]*(dt)*(du)+dZArray[3]*(1.-dt)*(du);
} // BiLinearInterpolation




double LinearInterpolation(double dMidPos,
                           double dFirstPos,
                           double dSecondPos,
                           double dFirstVal,
                           double dSecondVal)
{
    // Interpolate dDistance from dFirstPos to dSecondPos
    double dGapBetween = dSecondPos - dFirstPos;
    double  dFirstScalar = (dSecondPos - dMidPos)/dGapBetween;
    double dSecondScalar = (dMidPos - dFirstPos)/dGapBetween;
    return (dFirstVal * dFirstScalar) + (dSecondVal * dSecondScalar);
} // LinearInterpolation(...)




double InterpolateMatrix(const DKMaille2D<double> &matrix,
                         double dX,
                         double dY,
                         DKMaille<double>  &xAxis,
                         DKMaille<double>  &yAxis)
{
    unsigned uiAbsX1Pos=1; // absolute X and Y positions, in case dX and/or dY don't hit the matrix
    unsigned uiAbsY1Pos=1;
    unsigned uiAbsX2Pos=2;
    unsigned uiAbsY2Pos=2;

    // If we don't supply xAxis and/or yAxis parameters we assume that they are
    // the same size as the rows or columns of the matrix and numbered in integers
    // from zero.
    if (0 == xAxis.entries())
    {
        unsigned uiLimit = matrix.rows();
        xAxis.resize(uiLimit);
        for (unsigned ii=0;ii<uiLimit;ii++)
        {
            xAxis[ii] = ii;
        }
    }

    if (0 == yAxis.entries())
    {
        unsigned uiLimit = matrix.columns();
        yAxis.resize(uiLimit);
        for (unsigned ii=0;ii<uiLimit;ii++)
        {
            yAxis[ii] = ii;
        }
    }
    // Done setting up defaults

    // Identify whether it's in the matrix
    if (dX >= xAxis[xAxis.entries()-1])
    {
        uiAbsX1Pos = xAxis.entries()-1;
        uiAbsX2Pos = uiAbsX1Pos;
    }
    else if((dX <= xAxis[0])
            || (1 == xAxis.entries()))
    {
        uiAbsX1Pos = 0;
        uiAbsX2Pos = 0;
    } // if dX in matrix
    else
    {
        // Find the positions along the xAxis
        for (unsigned ii=0; ii<xAxis.entries()-1; ii++)
        {
            if (dX == xAxis[ii+1])
            {
                uiAbsX1Pos = ii+1;
                uiAbsX2Pos = uiAbsX1Pos;
                break;
            }
            else if ((dX > xAxis[ii])
                     && (dX < xAxis[ii+1]))
            {
                uiAbsX1Pos = ii;
                uiAbsX2Pos = ii+1;
                break;
            }
        } // for
    } // finding the xAxis positions



    if (dY >= yAxis[yAxis.entries()-1])
    {
        uiAbsY1Pos = yAxis.entries()-1;
        uiAbsY2Pos = uiAbsY1Pos;
    }
    else if((dY <= yAxis[0])
            || (1 == yAxis.entries()))
    {
        uiAbsY1Pos = 0;
        uiAbsY2Pos = 0;
    } // if dY in matrix
    else
    {
        // Find the positions along the yAxis
        for (unsigned ii=0; ii<yAxis.entries()-1; ii++)
        {
            if (dY == yAxis[ii+1])
            {
                uiAbsY1Pos = ii+1;
                uiAbsY2Pos = uiAbsY1Pos;
                break;
            }
            else if ((dY > yAxis[ii])
                     && (dY < yAxis[ii+1]))
            {
                uiAbsY1Pos = ii;
                uiAbsY2Pos = ii+1;
                break;
            }
        } // for
    } // finding the yAxis positions

    // If it's off one of the corners then just return the resulting position
    if ((uiAbsX1Pos == uiAbsX2Pos)
            && (uiAbsY1Pos == uiAbsY2Pos))
    {
        return matrix.at(uiAbsX1Pos, uiAbsY1Pos);
    }
    // Checks for linear interpolation on one side
    else if (uiAbsX1Pos == uiAbsX2Pos)
    {
        // Just need to interpolate yAxis
        return LinearInterpolation(dY,
                                   yAxis[uiAbsY1Pos],
                                   yAxis[uiAbsY2Pos],
                                   matrix.at(uiAbsX1Pos, uiAbsY1Pos),
                                   matrix.at(uiAbsX1Pos, uiAbsY2Pos));
    }
    else if (uiAbsY1Pos == uiAbsY2Pos)
    {
        // Just need to interpolate xAxis
        return LinearInterpolation(dX,
                                   xAxis[uiAbsX1Pos],
                                   xAxis[uiAbsX2Pos],
                                   matrix.at(uiAbsX1Pos, uiAbsY1Pos),
                                   matrix.at(uiAbsX2Pos, uiAbsY1Pos));
    } // if it's off the edge of the matrix
    else
    {
        // If we've reached this point then the required value must be within the matrix.
        // So we need to do bilinear interpolation.
        DKMaille<double> dXArray(2);
        DKMaille<double> dYArray(2);
        DKMaille<double> dZArray(4);

        dXArray[0] = uiAbsX1Pos;
        dXArray[1] = uiAbsX2Pos;
        dYArray[1] = uiAbsY1Pos;
        dYArray[0] = uiAbsY2Pos;
        dZArray[0] = matrix.at(uiAbsX1Pos, uiAbsY2Pos);
        dZArray[1] = matrix.at(uiAbsX2Pos, uiAbsY2Pos);
        dZArray[2] = matrix.at(uiAbsX2Pos, uiAbsY1Pos);
        dZArray[3] = matrix.at(uiAbsX1Pos, uiAbsY1Pos);

        // Scale dX and dY so that they fall in the correct position between the absolute positions
        // in the matrix.
        double dXAbsPos = uiAbsX1Pos+(dX-xAxis[uiAbsX1Pos])/(xAxis[uiAbsX2Pos]-xAxis[uiAbsX1Pos]);
        double dYAbsPos = uiAbsY1Pos+(dY-yAxis[uiAbsY1Pos])/(yAxis[uiAbsY2Pos]-yAxis[uiAbsY1Pos]);

        return BiLinearInterpolation(dXArray,
                                     dYArray,
                                     dZArray,
                                     dXAbsPos,
                                     dYAbsPos);
    }

    return 0;
} // InterpolateMatrix(...)



void CreateVanillaSwap(double dSpotDate,
                       double dXX,
                       double dYY,
                       DKMaille<double> &dSwaptionDates,
                       DKMaille<double> &dAccrualPeriods,
                       DKMaille<double> &dBasisSwaptionDates,
                       DKMaille<double> &dBasisAccrualPeriods,
                       DKMaille<double> &dBasis,
                       double dNoticePeriod,
                       DKMaille<double> &dCurveDates,
                       DKMaille<double> &dCurve,
                       DKMaille<double> &dCurveAdjusted)
{

    DKMaille<double> dZCRates(dCurveDates.entries());
    DKMaille<double> dZCAdjRates(dCurveDates.entries());

    for(unsigned ui=0;ui<dCurveDates.entries();ui++)
    {
        if(dCurve[ui]==1.)
        {
            dZCRates[ui]=-log(dCurve[ui+1])/((dCurveDates[ui+1]-dSpotDate)/365.);
            dZCAdjRates[ui]=-log(dCurveAdjusted[ui+1])/((dCurveDates[ui+1]-dSpotDate)/365.);
        }
        else
        {
            dZCRates[ui]=-log(dCurve[ui])/((dCurveDates[ui]-dSpotDate)/365.);
            dZCAdjRates[ui]=-log(dCurveAdjusted[ui])/((dCurveDates[ui]-dSpotDate)/365.);
        }
    }

    // find the size of the arrays
    // insert here bank-specific functionality to calculate number of payment dates for fixed and floating leg
    // currently defaults to semi-annual and I am still waiting

    double dFixedPaymentFrequency=2.;
    double dFloatPaymentFrequency=2.;

    // This is to make sure we always have real tenors
    int i=1;
    double dIncrement=0.5;
    double dLow=0.4;
    double dHigh=0.5;
    while(i>0&&i<=41)
    {
        if(dLow<dXX&&dXX<dHigh)
        {
            dXX=dHigh;
            i=-1;
        }
        dLow=dLow+dIncrement;
        dHigh=dHigh+dIncrement;
        i++;
    }

    unsigned int uiStart=0; unsigned int uiEnd=-1;  unsigned um=0;
    unsigned int uiSize= (unsigned int)(dFixedPaymentFrequency *(dXX)+1);
    unsigned int uiSizeFloat= (unsigned int)(dFloatPaymentFrequency *(dXX));



    dSwaptionDates.resize(uiSize);
    dAccrualPeriods.resize(uiSize);
    // Artificially set basis data
    dBasisSwaptionDates.resize(uiSizeFloat);
    dBasisAccrualPeriods.resize(uiSizeFloat);
    dBasis.resize(uiSizeFloat);
    for(ui=0;ui<uiSize;ui++)
    {
        // Insert here CDC functionality to calculate 1. adjusted payment dates and 2. accrual periods for the
        // fixed and the floating leg based on the currency and its corresponding plain vanilla swap/swaption convention
        // currently defaults to ACT/365.
        // currently defaults to no business day adjustment
        dSwaptionDates.at(ui)=dSpotDate+dNoticePeriod*365.+dYY*365.+(double)(ui)*365./dFixedPaymentFrequency;
        if(ui==0) dAccrualPeriods.at(ui)=0.0;
        else dAccrualPeriods.at(ui)=(dSwaptionDates.at(ui)-dSwaptionDates.at(ui-1))/365.;
        if(ui!=0)
        {
            dBasisSwaptionDates.at(ui-1)=dSwaptionDates.at(ui);
            dBasisAccrualPeriods.at(ui-1)=(dSwaptionDates.at(ui)-dSwaptionDates.at(ui-1))/360.;
            double dDate1=dSwaptionDates.at(ui-1);
            double dDate2=dSwaptionDates.at(ui);
            double dBasisDiscountAtDate1 = exp(-rateinterpolation_dk_maille(2,dDate1,dCurveDates,dZCRates,dCurveDates.entries()-1)*((dDate1-dSpotDate)/365.))
                                           /exp(-rateinterpolation_dk_maille(2,dDate1,dCurveDates,dZCAdjRates,dCurveDates.entries()-1)*((dDate1-dSpotDate)/365.));
            double dBasisDiscountAtDate2 = exp(-rateinterpolation_dk_maille(2,dDate2,dCurveDates,dZCRates,dCurveDates.entries()-1)*((dDate2-dSpotDate)/365.))
                                           /exp(-rateinterpolation_dk_maille(2,dDate2,dCurveDates,dZCAdjRates,dCurveDates.entries()-1)*((dDate2-dSpotDate)/365.));
            dBasis.at(ui-1)=(dBasisDiscountAtDate1/dBasisDiscountAtDate2-1.)/dBasisAccrualPeriods.at(ui-1);
        }
    }
}

void FromDiscountToZero(double dSpotDate,
                        DKMaille<double> &dZCRates,
                        DKMaille<double> &dZCAdjRates,
                        const DKMaille<double> &dDiscountDates,
                        const DKMaille<double> &dDiscountRates,
                        const DKMaille<double> &dAdjustedDiscountRates)
{
    // Move to zero curve
    for(unsigned int ui=0;ui<dDiscountDates.entries();ui++)
    {
        if(dDiscountRates[ui]==1.)
        {
            dZCRates.at(ui)=-log(dDiscountRates[ui+1])/((dDiscountDates[ui+1]-dSpotDate)/365.);
            dZCAdjRates.at(ui)=-log(dAdjustedDiscountRates[ui+1])/((dDiscountDates[ui+1]-dSpotDate)/365.);
        }
        else
        {
            dZCRates.at(ui)=-log(dDiscountRates[ui])/((dDiscountDates[ui]-dSpotDate)/365.);
            dZCAdjRates.at(ui)=-log(dAdjustedDiscountRates[ui])/((dDiscountDates[ui]-dSpotDate)/365.);
        }
    }
} // FromDiscountToZero()

void AggregateVolStrip(DKMaille<double> dStrip,DKMaille<double> dInputStrip,DKMaille<double> dVolStrip,DKMaille<double> &dStripOutput)
{
    if(dInputStrip.at(0)!=dStrip.at(0))
        throw("Error input Dates in vol strip conversion. Spot dates not consistent");
    for(unsigned int ui=0;ui<dStrip.entries();ui++)
    {
        for(unsigned int uiM=0;uiM<dInputStrip.entries()-1;uiM++)
        {
            if((dStrip.at(ui)>=dInputStrip.at(uiM))&&(dStrip.at(ui)<=dInputStrip.at(uiM+1)))
                dStripOutput.at(ui)=dVolStrip.at(uiM);
        }
        // extend beyond input size
        if(dStrip.at(ui)>=dInputStrip.at(dInputStrip.entries()-1))
            dStripOutput.at(ui)=dVolStrip.at(dVolStrip.entries()-1);
    }
} // AggregateVolStrip()

double signature(double x)
{
    double dXX=0.;
    if(x>=0.) dXX=1.;
    if(x<0.) dXX=-1.;
    return dXX;
}


double ArgumentFunction(double dX, double dY, double da, double db, double drho)
{
    double daprime = da/sqrt(2.*(1.-drho*drho));
    double dbprime = db/sqrt(2.*(1.-drho*drho));
    return exp(daprime*(2.*dX-daprime)+dbprime*(2.*dY-dbprime)+2.*drho*(dX-daprime)*(dY-dbprime));
} // ArgumentFunction

double CoreCumulativeBivariate(double da, double db, double drho)
{
    double dConstant = sqrt(1.-drho*drho) / 3.141592653589790 ;
    DKMaille<double> A(4); A[0]=0.3253030; A[1]=0.4211071; A[2]=0.1334425; A[3]=0.006374323;
    DKMaille<double> B(4); B[0]=0.1337764; B[1]=0.6243247; B[2]=1.3425378; B[3]=2.2626645;
    double dSum=0.;
    for(unsigned int uI=0;uI<A.entries();uI++)
    {
        for(unsigned int uL=0;uL<B.entries();uL++)
        {
            dSum+=dConstant*A[uI]*A[uL]*ArgumentFunction(B[uI],B[uL],da,db,drho);
        }
    }
    return dSum;
} // CoreCumulativeBivariate

double CumulativeNormalBivariateIntegral(double da, double db, double drho)
{
    double dX;
    double dNa=0.;
    double dNb=0.;
    normalDK(&dNa,da);
    normalDK(&dNb,db);
    double dTerma=0.;
    double dTermb=0.;
    if(da*db*drho<=0.)
    {
        if((da<=0.)&&(db<=0.)&&(drho<=0.)) dX=CoreCumulativeBivariate(da,db,drho);
        if((da<=0.)&&(db>=0.)&&(drho>=0.)) dX=dNa-CoreCumulativeBivariate(da,-db,-drho);
        if((da>=0.)&&(db<=0.)&&(drho>=0.)) dX=dNb-CoreCumulativeBivariate(-da,db,-drho);
        if((da>=0.)&&(db>=0.)&&(drho<=0.)) dX=dNa+dNb-1.+CoreCumulativeBivariate(-da,-db,drho);
    }

    else
    {
        double drhoab = (drho*da-db)*signature(da)/sqrt(da*da-2.*drho*da*db+db*db);
        double drhoba = (drho*db-da)*signature(db)/sqrt(da*da-2.*drho*da*db+db*db);
        double dDelta = 0.25*(1.-signature(da)*signature(db));
        if((da<0.0)&&(drhoab<0.0)) dTerma=CoreCumulativeBivariate(da,0.0,drhoab);
        if((da<0.0)&&(drhoab>0.0)) dTerma=dNa-CoreCumulativeBivariate(da,0.0,-drhoab);
        if((da>0.0)&&(drhoab>0.0)) dTerma=0.5-CoreCumulativeBivariate(-da,0.0,-drhoab);
        if((da>0.0)&&(drhoab<0.0)) dTerma=dNa+0.5-1.0+CoreCumulativeBivariate(-da,0.0,drhoab);
        if((db<0.0)&&(drhoba<0.0)) dTermb=CoreCumulativeBivariate(db,0.0,drhoba);
        if((db<0.0)&&(drhoba>0.0)) dTermb=dNb-CoreCumulativeBivariate(db,0.0,-drhoba);
        if((db>0.0)&&(drhoba>0.0)) dTermb=0.5-CoreCumulativeBivariate(-db,0.0,-drhoba);
        if((db>0.0)&&(drhoba<0.0)) dTermb=dNb+0.5-1.0+CoreCumulativeBivariate(-db,0.0,drhoba);
        dX=dTerma+dTermb-dDelta;
    }
    return dX;
} // CumulativeNormalBivariateIntegral
