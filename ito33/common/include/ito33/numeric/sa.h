///////////////////////////////////////////////////////////////////////////////
// File:             ito33/numeric/sa.h                                     
// Purpose:          Simulated annealing                                    
// Author:           ITO 33                                                   
// Created:          April 20, 2005                                            
// RCS-ID:           $Id: sa.h,v 1.4 2006/03/23 09:27:14 yann Exp $                                                  
// Copyright (C) Trilemma LLP  2005 -                                        
///////////////////////////////////////////////////////////////////////////////

/**
   @file  ito33/numeric/sa.h

   @brief Simulated Annealing. Minimization
   method, that is slow but usually is used to
   avoid being stuck in local mnimums
 
 */
#ifndef _ITO33_NUMERIC_SIMULATEDANNEALING_H_
#define _ITO33_NUMERIC_SIMULATEDANNEALING_H_

#include "ito33/beforestd.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <time.h>
#include "ito33/afterstd.h"

#include "ito33/useexception.h"

#include "ito33/numeric/exception.h"
#include "ito33/numeric/numericerror.h"

extern const ito33::Error ITO33_NEG_TOL, ITO33_BAD_PARAM;

namespace ito33
{

namespace numeric
{
/// Declaration of the SA class
class SA
{
public:

  /**
     Ctor  
     
     @param dTolerance is optional
     @param pdParamLowerBound, lower bound vector
     @param pdParamUpperBound, upper bound vector
     @param nNStep Set number of iterations before altering step    
     @param nNTemp Number of iterations before varying temperature
   */
  SA(std::vector<double> pdParamLowerBound,
           std::vector<double> pdParamUpperBound,
           double dTolerance = 1.e-8)
           : m_pdParamLowerBound(pdParamLowerBound),
             m_pdParamUpperBound(pdParamUpperBound),
             m_dTolerance(dTolerance),
             m_dCoolingRate(.85),
             m_dStartingTemperature(150000.),
             m_nMaxTempReductions(500),
             m_nNStep(20),
             m_nNEpsilon(6),
             m_bHasOutput(false),
             m_dStepMax(.5)
  {   

    if ( dTolerance <= 0.0 )
      throw EXCEPTION_MSG
            (
              ITO33_NEG_TOL,
              TRANS("Negative tolerance in simulating annealing solver.")
            );
    
    if ( pdParamLowerBound.size() != pdParamUpperBound.size() )
       throw EXCEPTION_MSG
            (
              ITO33_BAD_PARAM,
              TRANS("Lower Bound and Upper bound constraints \n"
                "do not have the same size.")
            );

     size_t nIdx;
     for ( nIdx = 0 ; nIdx < pdParamUpperBound.size(); nIdx++)
     {
       if ( pdParamLowerBound[nIdx] > pdParamUpperBound[nIdx] )
           throw EXCEPTION_MSG
            (
              ITO33_BAD_PARAM,
              TRANS("Error lower bound greater than upper bound.")
            );

     }

     m_nSize = pdParamUpperBound.size();

     m_nNTemp = 100 > m_nSize? 5*m_nSize : 100;
            
  }

  /**
    Set the temperature cooling rate
    recommended value by the paper of Corona, Marchesi, Martini and Ridella
    m_dCoolingRate = .85

    @param dCoolingRate cooling rate
  */
  void SetCoolingRate(double dCoolingRate)
  {
    if ( dCoolingRate <= 0.0 || dCoolingRate > 1.0)
      throw EXCEPTION_MSG
            (
              ITO33_BAD_PARAM,
              TRANS("Cooling rate negative or too large.")
            );

    m_dCoolingRate = dCoolingRate;
  }

  /**
    Get the cooling rate

    @return cooling rate
  */
  double GetCoolingRate() const
  {
    return m_dCoolingRate;
  }

  /**
    Set the starting temperature

    @param dStartingTemperature starting temperature
  */
  void SetStartingTemperature(double dStartingTemperature)
  {
    if ( dStartingTemperature < 1.e5 || dStartingTemperature > 1.e10)
      throw EXCEPTION_MSG
            (
              ITO33_BAD_PARAM,
              TRANS("Starting temperature too small or too large.")
            );

    m_dStartingTemperature = dStartingTemperature;
  }

  /**
    Get the starting temperature

    @return the starting temperature
  */
  double GetStartingTemperature() const
  {
    return m_dStartingTemperature;
  }

  /**
   Set the maximum number of temperature reduction
   before forced termination

   @param nMaxTempReductions
  */
  void SetMaxTemperatureReduction(size_t nMaxTempReductions)
  {
    if ( nMaxTempReductions < 100 || nMaxTempReductions > 1000)
      throw EXCEPTION_MSG
            (
              ITO33_BAD_PARAM,
              TRANS("Invalid number of temperature reduction.")
            );

    m_nMaxTempReductions = nMaxTempReductions;
  }

  /**
   Get the maximum number of temperature reduction
   before forced termination

   @return nMaxTempReductions
  */
  size_t GetMaxTemperatureReduction() const
  {
    return m_nMaxTempReductions;
  }


  /**
    Set the number of steps before altering direction

    @param nNStep
  */
  void SetNumberOfStepsBeforeAlteringDirection(size_t nNStep)
  {
    if ( nNStep < 1 || nNStep > 1000)
      throw EXCEPTION_MSG
            (
              ITO33_BAD_PARAM,
              TRANS("Invalid number of steps before altering direction.")
            );

    m_nNStep = nNStep;
  }

   /**
    Get the number of steps before altering direction

    @return nNStep
  */
  size_t GetNumberOfStepsBeforeAlteringDirection() const
  {
    return m_nNStep;
  }
  
  /**
    Number of iterations before varying temperature.  
    Recommned: nNtemp = 100 > 5*nSize? 100: 5*nSize;

    @param nNTemp
  */
  void SetNumberOfStepsBeforeAlteringTemperature(size_t nNTemp)
  {

    if ( nNTemp < 1)
      throw EXCEPTION_MSG
       (
         ITO33_BAD_PARAM,
         TRANS("Invalid number of steps before altering temperature.")
       );


    m_nNTemp = nNTemp;
  }

  /**
    Get the number of iterations before varying temperature.  

    @return nNTemp
  */
  size_t SetNumberOfStepsBeforeAlteringTemperature() const
  {
    return m_nNTemp;
  }

  /**
   AT THE END OF A TEMPERATURE LOOP: If the current energy is within
   epsilon of the previous Nepsilon energies (from the end of the previous
   temperature loops) AND the current temperature is within epsilon of the
   optimal energy found so far, terminate. This is carried out below.
   
    @param nNEpsilon
  */
  void SetNEpsilon(size_t nNEpsilon)
  {
    if ( nNEpsilon < 3)
      throw EXCEPTION_MSG
       (
         ITO33_BAD_PARAM,
         TRANS("Invalid number of steps before altering temperature.")
       );

    m_nNEpsilon = nNEpsilon;
  }

  /**
   Get nNEpsilon

   @return nNEpsilon
  */
  size_t GetNEpsilon() const
  {
    return m_nNEpsilon;
  }

  /**
    Set the output to true
  */
  void SetOutputToTrue()
  {
    m_bHasOutput = true;
  }

  /**
   Set the maximum relative change in the search direction

   @param dStepMax
  */
  void SetMaximumRelativeChangeInSearchDirection(double dStepMax)
  {
    m_dStepMax = dStepMax;
  }

  /**
   Get the maximum relative change in the search direction

   @return dStepMax
  */
  double GetMaximumRelativeChangeInSearchDirection() const
  {
    return m_dStepMax;
  }

  /**
     @param func is a template function that must take 
     
          pdParam current value of variables
          &dF     return value of the objective function

   */
  template<typename T>
    ito33::numeric::NumericError operator()(T& func, std::vector<double> &pdParam) 
  { 
    
    if ( pdParam.size() != m_nSize )    
      throw EXCEPTION_MSG
       (  
         ITO33_BAD_PARAM,
         TRANS("Constraint vector size different from number of variables.")
        );
    
    size_t nIdx;

    // Seed the random-number generator with current time so that
    // the random numbers will be different every time 
    srand( (unsigned)time( NULL ) );


    std::vector<double> pdParamNew(m_nSize);
    std::vector<double> pdParamOpt(m_nSize);
    std::vector<double> pdParamDir(m_nSize);
    double dFNew;
    double dFOpt;
    double dF;
    double dTemp = m_dStartingTemperature;

    //Initialization
    for (nIdx = 0; nIdx < m_nSize ; nIdx++)
      pdParamDir[nIdx]  = m_pdParamUpperBound[nIdx]- m_pdParamLowerBound[nIdx];

    func(pdParam,dF);
    m_nNbFunctionEvaluations = 1;

    pdParamOpt = pdParam;
    dFOpt      = dF;
    dFNew      = dF;

   //  Set the parameter used to adjust the stepsize.
   double C = 2.0;

   // Allocate array to store final temperature at each temperature.
   std::vector<double> Toptvec(m_nMaxTempReductions);


   /// Initialize counters.
   size_t nCounterK = 1;
   size_t nCounterJ = 1;
   size_t nCounterM = 1;
   size_t nNumAcceptances = 0;
   bool bConv = false;
   bool bNoConvergence = true;
   size_t nRestart = 10;

   // Do SA.

   while ( ( nCounterK <= m_nMaxTempReductions) && (bNoConvergence) )
   {     
     nCounterM = 1;

      // Periodically restart by using the best point so far. 
      if ( ( nCounterK % nRestart)==0 )
      {
        pdParam = pdParamOpt;
      }

      while ( nCounterM <= m_nNTemp)
      {
         // Over a loop for a fixed stepsize, record the number of acceptances.
         nNumAcceptances = 0;
         nCounterJ = 1;

         while ( nCounterJ <= m_nNStep)
         {
            
           // Generate a new point. 
           size_t nCounterInnerLoop = 0;
           
           do
           {
            
             for (nIdx = 0 ; nIdx < m_nSize ; nIdx++)
             {
               pdParamNew[nIdx] = pdParam[nIdx]
                + (2.0*rand()/double(RAND_MAX)-1.0)*pdParamDir[nIdx]; 
             }

             nCounterInnerLoop++;

             if ( nCounterInnerLoop > 10000 )
             {
               pdParam = pdParamOpt;
               nCounterInnerLoop = 0;
        
               for (nIdx = 0; nIdx < m_nSize ; nIdx++)
                 pdParamDir[nIdx]  = m_pdParamUpperBound[nIdx]- m_pdParamLowerBound[nIdx];
             }

           } while ( !CheckConstraints(pdParamNew) );

           
           func(pdParamNew, dFNew);

           m_nNbFunctionEvaluations++;

            // Update according to the Metropolis criterion.
            if ( BoltzProb(dF, dFNew, dTemp) > double(rand())/double(RAND_MAX) )
            {
               pdParam = pdParamNew;
               dF = dFNew;
               nNumAcceptances++;
            }  

            // If E is the best so far, update accordingly. 
            if ( dFNew < dFOpt)
            { 
              pdParamOpt = pdParamNew;
               dFOpt = dFNew;
            }

            nCounterJ++;

         } //end j-loop

    
         // Update the step vector. 
         for (nIdx = 0 ; nIdx < m_nSize; nIdx++)
         {
            double dOldDir = pdParamDir[nIdx];

            if (nNumAcceptances > 0.6*m_nNStep)
            {
               pdParamDir[nIdx]= pdParamDir[nIdx]*
                   (1.0 + C*(nNumAcceptances/double(m_nNStep) - 0.6)/0.4);
            }
            else 
            {
               pdParamDir[nIdx]= pdParamDir[nIdx]/
                  (1.0 + C*(0.4-nNumAcceptances/double(m_nNStep))/0.4);
            }

            if ( pdParamDir[nIdx]/ dOldDir - 1 > m_dStepMax)
               pdParamDir[nIdx] = dOldDir;

         }
         
 
         if ( m_bHasOutput )
         {
           std::cout << "Number of steps before altering Temperature: " 
             << m_nNTemp - nCounterM << std::endl;
           std::cout << "Number of Accepted point: " 
             << nNumAcceptances << std::endl;
           std::cout << "dF: " << dF << " dFOpt " << dFOpt 
             << " dTemp: " << dTemp << std::endl;
         }

         nCounterM++;

      } //end m-loop

      // Store the final value for this temperature.
      Toptvec[nCounterK-1] = dF;

      //Check for convergence
      bConv = false;
      
      if ( nCounterK > m_nNEpsilon)
      {
         bConv = true;
         
         for ( nIdx = 1; nIdx <= m_nNEpsilon; nIdx++)
         {
            bConv = bConv
              && (fabs(dF-Toptvec[nCounterK - nIdx - 1]) < m_dTolerance);
         }

         bConv = bConv && (fabs( dF - dFOpt) < m_dTolerance);
      }

      bNoConvergence = !(bConv);

      //Reduce the temperature.
      dTemp = m_dCoolingRate*dTemp;

      nCounterK++;

      if ( m_bHasOutput )
      {
        std::cout << std::endl;
        std::cout << "Current Temperature: " << dTemp 
          << " Best optimal found so far: " << dFOpt << std::endl;
        std::cout << std::endl;
      }

   } //end k-loop


   if ( !(bNoConvergence) )
   {
      return ITO33_NO_ERROR;
   }

   return ITO33_NOT_CONVERGED;
  
  
} //end SA

/**
  Get the number of function evaluation

  @return the number of function evaluation

*/
  size_t GetNbEvaluationFunction() const
  {
    return m_nNbFunctionEvaluations;
  }
  

private: 

  ///tolerance
  double m_dTolerance;

  ///lower bound vector
  std::vector<double> m_pdParamLowerBound;

  ///upper bound vector
  std::vector<double> m_pdParamUpperBound;

  ///Number of function evaluation
  size_t m_nNbFunctionEvaluations;

  ///cooling rate
  double m_dCoolingRate;

  ///Starting Temperature
  double m_dStartingTemperature;

  ///Maximum number of temperature reduction
  size_t m_nMaxTempReductions;
    
  ///Set number of iterations before altering step.
  size_t m_nNStep;

  ///Number of iterations before varying temperature.  
  size_t m_nNTemp;

  ///Number of unknows
  size_t m_nSize;
  
  /// AT THE END OF A TEMPERATURE LOOP: If the current energy is within
  ///   epsilon of the previous Nepsilon energies (from the end of the previous
  ///   temperature loops) AND the current temperature is within epsilon of the
  ///   optimal energy found so far, terminate. This is carried out below.
  size_t m_nNEpsilon;

  ///maximum relative change in the search direction
  double m_dStepMax;

  ///Set falg to indicate output
  bool m_bHasOutput;

  /**
    Check whether or not the contraints have been violated

   @return true/false if constraints are not violated/violated
  */
  bool CheckConstraints(const std::vector<double> &pdParam)
  {    
    size_t nIdx;
    for ( nIdx = 0 ; nIdx < m_nSize ; nIdx++)
    {
      if ( pdParam[nIdx] > m_pdParamUpperBound[nIdx] ||
           pdParam[nIdx] < m_pdParamLowerBound[nIdx] )
         return false;
    }    
    
    return true;
  }

  /**
    Probability that we accept point x. If E2<E1 then the probability
    is 1. If E2>E1, the probability is calculated using a Boltzmann-type
   distribution. This is the Metropolis criterion.

   @param dE1 old energy value
   @param dE2 new energy value
   @param dT current Temperature
  */
  double BoltzProb(double dE1, double dE2, double dT)
  {
    if ( dE2 <= dE1)
      return 1.0;
 
    return exp(-(dE2 - dE1) / dT);
  } 


};


}  // namespace numeric 

}  // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_SA_H_

