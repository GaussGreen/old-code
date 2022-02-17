/*===========================================================
  Name    : dk_multifactor_skewtree.cpp
  Owner   : DK
  Created : DK  
  Comment : Don't touch !
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

#include "DKMaille.h"
#include "DKMaille2D.h"
#include "dk_utils.h"
#include "dk_multifactor_skewtree.h"

#define M_PI 3.14159265358979

static int smilemodel=0;
static int normal_lognormal_map=1;

int Round(double payf)
{
	int pay=0;
	pay=(int)payf;  // cast RDJ 31/7/00
	if((payf-pay)>=0.5) pay+=1;
	else pay=pay;
	return pay;
}
double fonction(int n,double poly[],double logg) {
double f_of_logg=0.;
int j=0;
f_of_logg=poly[n];
for(j=n-1;j>=0;j--) f_of_logg=f_of_logg*logg+poly[j];
return f_of_logg;
}

double fonction_derivee(int n,double poly[],double logg) {
double f_of_loggp=0.;
int j=0;
f_of_loggp=(double)n*poly[n];
for(j=n-1;j>=1;j--) f_of_loggp=f_of_loggp*logg+(double)j*poly[j];
return f_of_loggp;
}

double integrale(double x1,double x2,int n,double poly[])
{
double res=0.;
int i=0;
for(i=0;i<=n;i++) res+=(poly[i]/(double)(i+1))*(pow(x2,i+1)-pow(x1,i+1));
return res;
}

void find_poly(double logg[],double underlying[],int n,double *poly)
{

/* The below are the variables needed for determining the polynomial */
int k,j,i;
double phi,ff,b;
double s[4];

/* The function that calculates the polynomial's coefficients */
for(i=0;i<=n;i++) s[i]=poly[i]=0.;
s[n]=-logg[0];
for(i=1;i<=n;i++)	{
for(j=n-i;j<=n-1;j++) s[j]-=logg[i]*s[j+1];
s[n]-=logg[i];
			}

for(j=0;j<=n;j++)	{
phi=n+1;
for(k=n;k>=1;k--) phi=k*s[k]+logg[j]*phi;
ff=underlying[j]/phi;
b=1.;
for(k=n;k>=0;k--) 		{
poly[k]+=b*ff;
b=s[k]+logg[j]*b;
				}
			}
}      


void find_limits(double cutoff,int howmany,double underlyingfunction[],int *index1,int *index2) 
{
int i=howmany;
int out=0;


while((i>=1)&&(out==0)) {

if(underlyingfunction[i]<underlyingfunction[i-1])
if(cutoff<underlyingfunction[i-1]&&cutoff>=underlyingfunction[i]) {
*index1=i;
*index2=i-1;
out=1;
}

if(underlyingfunction[i]>underlyingfunction[i-1])
if(cutoff>=underlyingfunction[i-1]&&cutoff<underlyingfunction[i]) {
*index1=i;
*index2=i-1;
out=1;
}

i-=1;

}

}






int GetSub(DKMaille<double> pdArray, int iLimit,double dCutOff, int iNumPoints)
{
  // The limit of the index of the matrix is iLimit for an array 1..iLimit
  // 2D array - take only a 1D component...
  bool bFound = false;
  int ii = 1;

  if (pdArray[1] > dCutOff)
  {
    // We're starting above zero so need to see it hit zero
    for (ii=1; ((ii <= iLimit) && (pdArray[ii] > dCutOff)); ii++);

    if ((ii <= iLimit) && (pdArray[ii] <= dCutOff))
      bFound = true;
  }
  else
  {
    // Starting -ve and going up
    for (ii=1; ((ii <= iLimit) && (pdArray[ii] < dCutOff)); ii++);

    if ((ii <= iLimit) && (pdArray[ii] >= dCutOff))
      bFound = true;
  }

  if (bFound)
  {
    // Now ii denotes the point either just before or just after the zero point.
    // Therefore ii-1 is the other side of zero.
    // Check we're not close to the limits of the array
//    if ((ii < iLimit) && (ii > 2))
	int iPtsBeforeAxis = div(iNumPoints, 2).quot;
	if ((ii < iLimit) && (ii > iPtsBeforeAxis))
//      return ii-2; // in the middle
      return ii-iPtsBeforeAxis; // in the middle
//    else if (3 > ii)
    else if ((iPtsBeforeAxis+1) > ii)
      return 1; // at the start
    else
//      return iLimit-3;
		return iLimit-(iNumPoints-1);
    }
  else
  {
    // Failed to find it
    return -1;
  }
} // GetSub(...)


double smooching(int totalnodes,
						     DKMaille<double> *optiontree,
            		 int howmany,
						     DKMaille<double> g_tree,
						     DKMaille<double> array_intrinsic,
						     double cutoff,
                 unsigned int uiNumFactors)
{
double polynomial[4];
int ii=0;
int iii=0;
double limit1=0.;
double limit2=0.;
double ff=0.;
double dff=0.;
double dxx=0.;
double comeback=0.;
int jjj=100;
double accuracy=1.e-13;
int count=0;
double middle=0.;
int side=0;
int i1=0;
int i2=0;
double spacing=0.;
double treelogg[4];
double underlyingfunction[4];
int correctionindex=0;
double correction=0.;
int loop=0;
bool bDoSmoothing=true;			

memset(polynomial,0, 4*sizeof (double) );
memset(treelogg,0, 4*sizeof (double) );
memset(underlyingfunction,0, 4*sizeof (double) );

int iFirstIndex = GetSub(array_intrinsic,totalnodes,cutoff, howmany+1);

if (-1 == iFirstIndex)
	bDoSmoothing=false;

if(bDoSmoothing==true) 
{
				for(loop=0;loop<=howmany;loop++) 
        {
					treelogg[loop]=g_tree[iFirstIndex+loop];
					underlyingfunction[loop]=array_intrinsic[iFirstIndex+loop];
				}

spacing=fabs(treelogg[0]-treelogg[1]);
find_limits(cutoff,howmany,underlyingfunction,&i1,&i2);
find_poly(treelogg,underlyingfunction,howmany,polynomial);

/***** Newton-Raphson *****/

comeback=0.5*(treelogg[i1]+treelogg[i2]);
ii=1;
while(ii<=jjj)
{
ff=fonction(howmany,polynomial,comeback);
dff=fonction_derivee(howmany,polynomial,comeback);
dxx=ff/dff;
comeback-=dxx;
count+=1;
if(fabs(dxx)<accuracy) ii=jjj;
ii+=1;
}

middle=0.5*(treelogg[i1]+treelogg[i2]);
if(comeback<middle) side=0;
if(comeback>middle) side=1;
if(fabs(comeback-middle)<1.e-13) side=3;

if (side==3) {
correction=0.;
correctionindex=i1;
}
if(side==0) {
correction=fabs(integrale(comeback,middle,howmany,polynomial))/spacing;
correctionindex=i1;
}
if(side==1) {
correction=fabs(integrale(middle,comeback,howmany,polynomial))/spacing;
correctionindex=i2;
}

for(iii=0;iii<=howmany;iii++) 
	{
	if(iii==correctionindex) 
    optiontree->at(iFirstIndex+iii)+=correction/(double)uiNumFactors;
	}			
}
return 1.0;
}


double binomial_coeff(int i,int j)
{
	int k,l,m;
	double result,use;
	i=i-1;
	j=j-1;
	k=i-j;

	result=1.0;
	if(j>k)
	{
		m=k;
		for(l=i;l>=j+1;l--)
		{
			if(m>1) use=((double)l) / ((double)m);
			else use=l;
			result*=use;
			m--;
		}
	}
	else 
	{
		m=j;
		for(l=i;l>=k+1;l--)
		{
			if(m>1)	use=((double)l) / ((double)m);
			else use=l;
			result *= use;
			m--;
		}
	}
return(result) ;
} 

int orderingintree(DKMaille<double> time_array, int array_size, double target_time, int *nearest_node)
{
	register int i=0;
  register int k=0;
	int flag=0;
  if(target_time<time_array[0]) 
	{
		*nearest_node=0;
		return 0;
	}
	if(target_time>time_array[array_size]) 
	{
		*nearest_node=array_size+1;              
		return 0;
	}
	for(i=0;i<=array_size;i++)
	{
		if(target_time==time_array[i])
		{
			*nearest_node=i+1;
			return 0;	
		}
	}
	while((k<array_size)&&(flag==0)) 
	{
		if((target_time>time_array[k])&&(target_time<time_array[k+1]))
		{
			*nearest_node=k+1;
			if(target_time>=time_array[k]+(time_array[k+1]-time_array[k])/2.) *nearest_node=k+2;
			flag=1;
		}
		k+=1;
	}      
	return 0;
}

int calage(DKMaille<double> time_array, 
           int array_size, 
           int *nearest_node2, 
           int *nearest_node1, 
           double target_time1, 
           double target_time2)
{
	int pos1=0,pos2=0;
	orderingintree(time_array,array_size,target_time1,&pos1);
	orderingintree(time_array,array_size,target_time2,&pos2);
	*nearest_node1=pos1;
	*nearest_node2=pos2;
	return 0;
}

