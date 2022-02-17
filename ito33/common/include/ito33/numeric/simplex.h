///////////////////////////////////////////////////////////////////////////////
// File:             ito33/numeric/simplex.h                                 //
// Purpose:          Simplex                                                 // 
// Author:           Ito33                                                   //
// Created:          2004/12/16                                              //
// RCS-ID:           $Id: simplex.h,v 1.4 2005/12/02 14:41:00 wang Exp $                                                  //
// Copyright (C) Trilemma LLP  2000 -                                        //
///////////////////////////////////////////////////////////////////////////////

/**
   @file  ito33/numeric/simplex.h

   @brief Simplex method numerical recipies
  
   It is rather a slow method
 */
#ifndef _ITO33_NUMERIC_SIMPLEX_H_
#define _ITO33_NUMERIC_SIMPLEX_H_

#include "ito33/beforestd.h"
#ifdef SIMPLEXDEBUG
  #include <iostream>
#endif
#include <cmath>
#include <vector>
#include "ito33/afterstd.h"

#include "ito33/useexception.h"

#include "ito33/numeric/exception.h"
#include "ito33/numeric/numericerror.h"

extern const ito33::Error ITO33_NEG_TOL;

namespace ito33
{

namespace numeric
{

class Simplex
{
public:

  /**
     Ctor  
     
     @param dTolerance is optional
     @param nNbMaxFuncCall is optional, maximum number of function call allowed
     @param dParam1LowerBound, dParam1UpperBound of the first parameter
     @param dParam2LowerBound, dParam2UpperBound of the second parameter
   */
  Simplex(double dParam1LowerBound,
         double dParam1UpperBound,
         double dParam2LowerBound,
         double dParam2UpperBound,
         double dTolerance = 1e-6, 
         size_t nNbMaxFuncCall = 500)
            : m_dParam1LowerBound(dParam1LowerBound),
              m_dParam1UpperBound(dParam1UpperBound),
              m_dParam2LowerBound(dParam2LowerBound),
              m_dParam2UpperBound(dParam2UpperBound),
              m_dTolerance(dTolerance),
              m_nNbMaxFuncCall(nNbMaxFuncCall), 
              m_nNbFuncCall(0)             
  {
    
    if ( dTolerance <= 0.0 )
      throw EXCEPTION_MSG
            (
              ITO33_NEG_TOL,
              TRANS("Negative tolerance in Newton solver.")
            );
            
  }

  /**
     @param func is a template function that must take 
     dParam1,dParam2,dPrice1,dPrice2,dF,dPrice1Error,dPrice2Error
     7 parameters:
          dParam1 current value of first variable
          dParam2 current value of second variable
          &dPrice1 first cds price
          &dPrice2 second cds price
          &dF      value of the objective function
          &dPrice1Error difference between the market price and first cds
          &dPrice2Error difference between the market price and second cds


     @param dInitial1 is the first  guess
     @param dInitial2 is the second guess
   */
  template<typename T>
  ito33::numeric::NumericError operator()(T& func, double &dParam1, 
                                          double &dParam2)
  {
    //reset values  
    m_nNbFuncCall = 0;
        
    std::vector< std::vector<double> > pP;//contains the simplex points
    std::vector<double> pSum(2);
    std::vector<double> pY(3); //value at each points
    
   size_t nIdx,nIhi,nIlo,nInhi, nMpts = 3,nDim = 2;
   double dRtol,dYsave,dYtry;

    double dPrice1       = 0.;
    double dPrice2       = 0.;
    double dPrice1Error  = 0.;
    double dPrice2Error  = 0.;
    double dF            = 0.;

   //initialization
   //we are defining a starting point
   //to be of size n+1 where n is
   //the number of unknows
   std::vector<double> P0(2);
   P0[0] = dParam1;
   P0[1] = dParam2;

   func(P0[0],P0[1],dPrice1,dPrice2,dF,dPrice1Error,dPrice2Error);
   m_nNbFuncCall++;
   pY[0] = dF;  

   //P1
   std::vector<double> P1(2);
   P1[0] = m_dParam1LowerBound;
   P1[1] = m_dParam2LowerBound;

   func(P1[0],P1[1],dPrice1,dPrice2,dF,dPrice1Error,dPrice2Error);
   m_nNbFuncCall++;
   pY[1] = dF;

   //P2
   std::vector<double> P2(2);
   P2[0] = m_dParam1UpperBound;
   P2[1] = m_dParam2UpperBound;
   
   func(P2[0],P2[1],dPrice1,dPrice2,dF,dPrice1Error,dPrice2Error);
   m_nNbFuncCall++;
   pY[2]  = dF;  

   //store the starting point
   pP.push_back(P0);
   pP.push_back(P1);
   pP.push_back(P2);

   //compute sum and check constraints are not violated
   pSum[0] = pP[0][0] + pP[1][0] + pP[2][0];
   pSum[1] = pP[0][1] + pP[1][1]+ pP[2][1];

   CheckConstraints(pSum[0],pSum[1]);

  
  for (;;) 
  {
    nIlo = 0;

    //First we must determine which point is the highest (worst), next-highest, and lowest
    //(best), by looping over the points in the simplex.

    nIhi = pY[0]>pY[1] ? (nInhi=1,0) : (nInhi=0,1); 

    for ( nIdx = 0; nIdx < nMpts; nIdx++) 
    {

      if ( pY[nIdx] <= pY[nIlo]) 
        nIlo = nIdx;

      if ( pY[nIdx] > pY[nIhi]) 
      {
        nInhi = nIhi;
        nIhi  = nIdx;
      } 
      else if ( pY[nIdx] > pY[nInhi] && nIdx != nIhi) 
        nInhi = nIdx;
   
    } //end for loop

  dRtol = 2.0*fabs(pY[nIhi]-pY[nIlo])/(fabs(pY[nIhi])+fabs(pY[nIlo])+1.e-10);
  //Compute the fractional range from highest to lowest and return if satisfactory.

  if (dRtol < m_dTolerance) //If returning, put best point and value in slot 1.   
  { 
    Swap(pY[0],pY[nIlo]);

    for ( nIdx = 0; nIdx < nDim;nIdx++) 
      Swap(pP[0][nIdx],pP[nIlo][nIdx]);

# ifdef SIMPLEXDEBUG   
    std::cout << y[0] << " " << pY[1] << " " << pY[2] << std::endl;
    std::cout << "alpha = " << pP[0][0] << " beta = " << pP[0][1] << std::endl;
    std::cout << "nFunction call " << m_nNbFuncCall << std::endl;
#endif

    dParam1 = pP[0][0];
    dParam2 = pP[0][1];

    if ( pY[0] < m_dTolerance )
      return ITO33_CONVERGED;

    return ITO33_NOT_CONVERGED;
  }//end if

   
  if ( m_nNbFuncCall >= m_nNbMaxFuncCall )
  {
	  dParam1 = pP[0][0];
    dParam2 = pP[0][1];

    return ITO33_TOO_MANY_FUNCTION_CALL;
  } //end if

  m_nNbFuncCall += 2;

 //   Begin a new iteration. First extrapolate by a factor -1 through the face of the simplex
//across from the high point, i.e., reflect the simplex from the high point.

  dYtry = Amotry(func,pP,pY,pSum,nIhi,-1.0);

  if ( dYtry <= pY[nIlo])
  {
  //Gives a result better than the best point, so try an additional extrapolation by a
  //factor 2.
   dYtry = Amotry(func,pP,pY,pSum,nIhi,2.0);

  }
  else if ( dYtry >= pY[nInhi]) 
  {
  //The reflected point is worse than the second-highest, so look for an intermediate
  //lower point, i.e., do a one-dimensional contraction.
    dYsave = pY[nIhi];
    dYtry  = Amotry(func,pP,pY,pSum,nIhi,0.5);

    if ( dYtry >= dYsave) 
    { 
      //Can’t seem to get rid of that high point. Better
      //contract around the lowest (best) point. 
      for ( nIdx = 0; nIdx< nMpts; nIdx++) 
      {
        if ( nIdx != nIlo) 
        {     

          pP[nIdx][0]=pSum[0]=0.5*(pP[nIdx][0] + pP[nIlo][0]);
          pP[nIdx][1]=pSum[1]=0.5*(pP[nIdx][1] + pP[nIlo][1]);

          CheckConstraints(pP[nIdx][0],pP[nIdx][1]);
          CheckConstraints(pSum[0],pSum[1]);
      
          func(pSum[0],pSum[1],dPrice1,dPrice2,dF,dPrice1Error,dPrice2Error);
          
          pY[nIdx] = dF;
        }//end if

      }//end for


      m_nNbFuncCall += nDim; 
      //Keep track of function evaluations.
      //Recompute psum.
      
      pSum[0] = pP[0][0] + pP[1][0] + pP[2][0]; 
      pSum[1] = pP[0][1] + pP[1][1]+ pP[2][1];

      CheckConstraints(pSum[0],pSum[1]);

    }

  } 
  else 
     m_nNbFuncCall--; //Correct the evaluation count.

  } //Go back for the test of doneness and the next
    //iteration. 

} //Simplex


template<typename T>
double Amotry(T& func,std::vector< std::vector<double> > &pP, 
              std::vector<double> &pY, std::vector<double> &pSum,
              size_t nIhi, double dFac)
//Extrapolates by a factor fac through the face of the simplex across from the high point, tries
//it, and replaces the high point if the new point is better.
{
   double dFac1,dFac2,dYtry;
   double dF = 0.0;
   std::vector<double> pPtry(2);
 
   dFac1 = (1.0 - dFac)/2.;
   dFac2 = dFac1 - dFac;

   pPtry[0] = pSum[0]*dFac1-pP[nIhi][0]*dFac2;
   pPtry[1] = pSum[1]*dFac1-pP[nIhi][1]*dFac2;

   CheckConstraints(pPtry[0],pPtry[1]);

   double dPrice1      = 0.;
   double dPrice2      = 0.;
   double dPrice1Error = 0.;
   double dPrice2Error = 0.;

   func(pPtry[0],pPtry[1],dPrice1,dPrice2,dF,dPrice1Error,dPrice2Error);

   dYtry = dF; //Evaluate the function at the trial point.

    
  if ( dYtry < pY[nIhi])   
  { //If it’s better than the highest, then replace the highest.
      pY[nIhi] = dYtry;
  
        pSum[0] +=  pPtry[0]-pP[nIhi][0];
        pSum[1] +=  pPtry[0]-pP[nIhi][1];
        pP[nIhi][0] = pPtry[0];
        pP[nIhi][1] = pPtry[1];
   
    }//end if

  return dYtry;
}//HazardRateSpotComponentPowerCalibrator::amotry


  /**
    Get the number of function calls
  */
  size_t GetNbFuncCall() const 
  { 
    return m_nNbFuncCall; 
  }


private: 

  double m_dTolerance;
  size_t m_nNbMaxFuncCall;
  size_t m_nNbFuncCall;
  double m_dParam1LowerBound;
  double m_dParam1UpperBound;
  double m_dParam2LowerBound;
  double m_dParam2UpperBound;

  //function to swap two double values
  void Swap(double &dVal1,double &dVal2)
  {
    double dTmp = dVal1;
    dVal1 = dVal2;
    dVal2 = dTmp;
  }//end swap


  //void check constraints
void CheckConstraints(double &dParam1,double &dParam2)
{
  if ( dParam1 > m_dParam1UpperBound )
    dParam1 = m_dParam1UpperBound;
     
  if ( dParam1 < m_dParam1LowerBound )
    dParam1 = m_dParam1LowerBound;
   
  if ( dParam2 > m_dParam2UpperBound )
    dParam2 = m_dParam2UpperBound;

  if ( dParam2 < m_dParam2LowerBound )
    dParam2 = m_dParam2LowerBound;

  }//checkConstraints

};


}  // namespace numeric 

}  // namespace ito33


#endif // #ifndef _ITO33_NUMERIC_SIMPLEX_H_
