
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


#include "dk_smoothing.h"
#include "DKMaille.h"
#include "DKMaille2D.h"



int iGetAdjustmentIndex(int ii, double dAdjustment, double dDifferenceInOptionValue)
{
if(dDifferenceInOptionValue==0.) return ii+1;
else return ii;
}



double LinearSmoothing(int buy_sell, double uNow, double uNext, double &optionValue)
{  
  double nextAdj = 0.0;

  double adj;
  double width;
  double hight;
  double xm;

  if (fabs(uNow)+fabs(uNext) > 0.0 && uNow*uNext <= 0.0)
  {
    xm = - uNow / (uNext - uNow);
    hight = (uNow + uNext)/2.0;
    width = 0.5 - xm;
    adj = 0.5 * width * hight;

    if ((0==buy_sell && uNow > 0.0) || (1==buy_sell && uNow < 0.0))
    {
      adj *= -1.0;
    }
  
    if (xm < 0.5)
    {
      optionValue += adj;
    }
    else
    {
      nextAdj = adj;
    }
  }

  return nextAdj;
} // LinearSmoothing(...)




//=========================================================================================================
bool SolveCubic(double &x, double *y, double &a, double &b, double &c, double &d)
{
  // Find the following equation, for given x[i], y[i], (i=0,1,2,3)
  // Then, solve the equation. (the solution is "x"))
  //
  // y[i] = a*x[i]^3 + b*x[i]^2 + c*x[i] + d
  //
  // x[0] = -1;
  // x[1] =  0;
  // x[2] =  1;
  // x[3] =  2;
  //
  // We know that the solution is between x[1] and x[2], i.e. between 0.0 and 1.0

  // Therefore,
  //
  // Inverse of [-1 1 -1 1]
  //            [ 0 0  0 1]
  //            [ 1 1  1 1]
  //            [ 8 4  2 1] gives the a,b,c,d

  double d6 = 1.0/6.0;
  double d3 = 1.0/3.0;
  double d2 = 1.0/2.0;

  a =  -d6*y[0] + d2*y[1] - d2*y[2] + d6*y[3];
  b =   d2*y[0] -    y[1] + d2*y[2];
  c =  -d3*y[0] - d2*y[1] +    y[2] - d6*y[3];
  d =                y[1];

  double A=0.0;
  double B=0.0;
  double xx, xxx;
  double newX;

  x = 0.0;
  int loopLimit = 10000;
  for (int i=0; i<loopLimit; ++i)
  {
    if (i == loopLimit-1)
    {
      return false;
    }
    
    xx = x*x;
    xxx= xx*x;

    A = 3*a*xx + 2*b*x + c;
    B =(a*xxx + b*xx + c*x + d) - A*x;

    if (fabs(A) < 1.0e-300)
    {
      x = x+1.001;// new assumption
      continue;
    }    
    newX = -B/A;

    if (fabs(newX-x) < 0.00001)
    {
      x = newX;
      break;
    }
    else
    {
      x = newX;
    }
  }

  
  if (x < 0.0 || 1.0 < x)
  {
    return false;
  }    
  return true;
} // SolveCubic(...)







//=========================================================================================================
double CubicSmoothing(int buy_sell, double uBefore, double uNow, double uNext, double uAfterNext, double &optionValue)
{
  double nextAdj = 0.0;

  double adj;
  
  if (fabs(uNow)+fabs(uNext) > 0.0 && uNow*uNext <= 0.0)
  {
    double y[4];
    y[0] = uBefore;
    y[1] = uNow;
    y[2] = uNext;
    y[3] = uAfterNext;

    double a,b,c,d;// where y = a*x*x*x + b*x*x + c*x + d
    double xm;
    if (!SolveCubic(xm, y, a,b,c,d))
    {
      // cubic failed
      // try linear
      nextAdj = LinearSmoothing(buy_sell, uNow, uNext, optionValue);
      return nextAdj;

      throw("Error : payoff.cpp : SolveCubic : No solution found");// problem;
    }

    // Integrate the area
    {
      double x = 0.5;
      double xx  = x*x;
      double xxx = x*xx;
      double xxxx= x*xxx;

      adj  = a*xxxx/4.0 + b*xxx/3.0 + c*xx/2.0 + d*x;

      x   = xm;
      xx  = x*x;
      xxx = x*xx;
      xxxx= x*xxx;

      adj -= a*xxxx/4.0 + b*xxx/3.0 + c*xx/2.0 + d*x;
    }

    if (0==buy_sell)
    {
      adj =  fabs(adj);
    }
    else
    {
      adj = -fabs(adj);
    }

    if (xm <= 0.5)
    {
      optionValue += adj;
    }
    else
    {
      nextAdj = adj;
    }
  }

  return nextAdj;
} // CubicSmoothing(...)


void SmoothingAlgorithm(int &iAdjustmentIndex,double &Adjustment,DKMaille<double> &dFunctionToSmooth,DKMaille<double> &dFunctionToAdjust)
{
	if(dFunctionToSmooth.entries()!=dFunctionToAdjust.entries())
		throw("Size misfit between functions in smoothing");
  unsigned uiNumNodesInVector = dFunctionToSmooth.entries();

  double dAdjustment=0.;

  double dOptionValue;
  double dOptionValueBefore;

  Adjustment=0.;

  for (int ii=0; ii < (int)uiNumNodesInVector; ii++)
  {
    dFunctionToAdjust.at(ii)+=dAdjustment;

    dAdjustment = 0.;

    if ((ii+1) == (int)uiNumNodesInVector)
    {
      // Do nothing to the top node 
    }
    else if ((0 <= (ii-1)) && ((ii+2) < (int)uiNumNodesInVector))
    {
      // Cubic smoothing on the nodes in the middle
      dOptionValue = dFunctionToAdjust.at(ii);
      dOptionValueBefore=dOptionValue;
      dAdjustment = CubicSmoothing(0,
                                  dFunctionToSmooth.at(ii-1),
                                  dFunctionToSmooth.at(ii),
                                  dFunctionToSmooth.at(ii+1),
                                  dFunctionToSmooth.at(ii+2),
                                  dOptionValue);

      // The adjustment may be added to this node, in which case it is included in
      // dOptionValue, alternatively it may need to be added to the next node up, in
      // which case it is included in dAdjustment and added next time round the loop.
      dFunctionToAdjust.at(ii)=dOptionValue;
      if((dAdjustment!=0.)||(dOptionValue-dOptionValueBefore)!=0.) 
      {
        iAdjustmentIndex=iGetAdjustmentIndex(ii,dAdjustment,-dOptionValueBefore+dOptionValue);
        if(ii==iAdjustmentIndex) Adjustment=-dOptionValueBefore+dOptionValue;
				else Adjustment=dAdjustment;
      }                                              
    }
    else if ((0 == ii) || ((int)uiNumNodesInVector > (ii+1)))
    {
      // Linear smoothing at the next to bottom node and the next to top node
      dOptionValue = dFunctionToAdjust.at(ii);
      dOptionValueBefore=dOptionValue;
      dAdjustment = LinearSmoothing(0,
																		dFunctionToSmooth.at(ii),
																		dFunctionToSmooth.at(ii+1),
																		dOptionValue);

      // See comment above about adding values.
      dFunctionToAdjust.at(ii)=dOptionValue;
      if((dAdjustment!=0.)||(dOptionValue-dOptionValueBefore)!=0.) 
      {
        iAdjustmentIndex=iGetAdjustmentIndex(ii,dAdjustment,-dOptionValueBefore+dOptionValue);
        if(ii==iAdjustmentIndex) Adjustment=-dOptionValueBefore+dOptionValue;
				else Adjustment=dAdjustment;
      } // if we made an adjustment
    } // if...else  chooses cubic or linear smoothing
  } // for
} // SmoothingAlgorithm(...)



