///////////////////////////////////////////////////////////////////////////////
// File:             ito33/numeric/asa.h                                     
// Purpose:          Adaptive simulated annealing                                    
// Author:           ITO 33                                                   
// Created:          May 17, 2005                                            
// RCS-ID:           $Id: asa.h,v 1.5 2006/03/23 09:27:14 yann Exp $                                                  
// Copyright (C) Trilemma LLP  2005 -                                        
///////////////////////////////////////////////////////////////////////////////

/**
   @file  ito33/numeric/asa.h

   @brief Adaptive Simulated Annealing. Minimization
   method, supposidely faster than regular simulated
   annealing inspired by the Lester Ingber paper
   www.ingber.com
 
 */
#ifndef _ITO33_NUMERIC_ADAPTIVE_SIMULATED_ANNEALING_H_
#define _ITO33_NUMERIC_ADAPTIVE_SIMULATED_ANNEALING_H_


#include "ito33/beforestd.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <time.h>
#include <fstream>
#include <stdio.h>
#include "ito33/afterstd.h"

#include "ito33/useexception.h"
#include "ito33/constants.h"

#include "ito33/numeric/exception.h"
#include "ito33/numeric/numericerror.h"
#include "ito33/numeric/predicatedouble.h"

/// Declaration of error code
extern const ito33::Error ITO33_NEG_TOL, ITO33_BAD_PARAM;

// normally there should be no "using" in header but this file is strictly
// internal and its usage is very limited ...
using namespace std;

namespace ito33
{

namespace numeric
{

/**
   Adaptive simulated annealing class.

   Method for finding a global minimum of a nonlinear function of
   several variables.
*/
class ASA
{
public:

  /**
     Ctor  
     
     @param pdLowerBound lower bound vector
     @param pdUpperBound upper bound vector
     @param dTolerance is optional
   */
  ASA(vector<double> pdLowerBound, vector<double> pdUpperBound, double dTolerance = 1.e-7)
           : m_pdLowerBound(pdLowerBound), m_pdUpperBound(pdUpperBound),
             m_dTolerance(dTolerance),
             m_dTemperatureRatioScale(1.e-5),
             m_dCostParameterScaleRatio(1.),
             m_nTemperatureAnnealScale(100),
             m_bHasOutput(false),
             m_bHasOutputDetails(false),
             m_nReannealingFrequency(10),
             m_nGeneratedFrequency(10000),
             m_nNumberSample(5),
             m_nLimitAcceptances(10000),
             m_nLimitGenerated(10000),
             m_nMaximumCostRepeat(5),
             m_dDoubleTolerance(1.e-18)
  {   

    if ( dTolerance <= 0.0 )
      throw EXCEPTION_MSG
            (
              ITO33_NEG_TOL,
              TRANS("Negative tolerance in simulating annealing solver.")
            );
    
    if ( pdLowerBound.size() != pdUpperBound.size() )
       throw EXCEPTION_MSG
            (
              ITO33_BAD_PARAM,
              TRANS("Lower Bound and Upper bound constraints \n"
                "do not have the same size.")
            );

     size_t nIdx;
     for ( nIdx = 0 ; nIdx < pdUpperBound.size(); nIdx++)
     {
       if ( pdLowerBound[nIdx] > pdUpperBound[nIdx] )
           throw EXCEPTION_MSG
            (
              ITO33_BAD_PARAM,
              TRANS("Error lower bound greater than upper bound.")
            );
     }

     m_nNumberParameters = pdUpperBound.size();
     m_dNumberParameters = double(m_nNumberParameters);
      
     //compute the difference between the bounds
     m_pdRangeBound.resize(m_nNumberParameters);

     for ( nIdx = 0; nIdx < m_nNumberParameters; nIdx++)
       m_pdRangeBound[nIdx] = pdUpperBound[nIdx] - pdLowerBound[nIdx];
    
 
  }//end constructor


  /**
    Set the output to true

    @param fileName name of the output file
  */
  void SetOutput(const char *fileName)
  {
    m_bHasOutput = true;
    outputFile.open(fileName);

    if ( !outputFile.is_open() )
      throw EXCEPTION_MSG
        (
           ITO33_BAD_PARAM,
           TRANS("Error opening debug output file.")
        );

    outputFile.setf(ios::scientific);
    outputFile.setf(ios::fixed);
    outputFile.setf(ios::showpoint);
    outputFile.precision(7);
    outputFile.width(14);
  }

  /**
   Set the detailed output level to true
   if the output has not been set previously
   then it is set automatically to write to
   the file asa.out
  */
  void SetOutputDetailsTrue()
  {
    m_bHasOutputDetails = true;

    if ( !m_bHasOutput )
      SetOutput("asa.out");
  }

  /**
    Default value is 1.0
    for quick decreasing temperature is worth experimenting
    and decreasing that value

    @param dCostParameterScaleRatio
  */
  void SetCostParameterScaleRatio(double dCostParameterScaleRatio)
  {
    m_dCostParameterScaleRatio = dCostParameterScaleRatio;
  }

  /**
   Set how often the temperature is reannealed
   default value is 100. Note that for fast decreasing
   temperature it is sometimes good to decrease
   that number see Rosenbrock example.

   @param nReannealingFrequency reannealing frequency
  */
  void SetReannealingFrequency(long int nReannealingFrequency)
  {
    m_nReannealingFrequency = nReannealingFrequency;
  }

  /// Get the generating frequency
  void SetGeneratingFrequency(long int nGeneratedFrequency)
  {
    m_nGeneratedFrequency = nGeneratedFrequency;
  }

  /**
   Set the number of sample point to set the initial cost
   temperature default 5
 
   @param nNumberSample number of sample point
  */
  void SetSamplingNumber(long int nNumberSample)
  {
    m_nNumberSample = nNumberSample;
  }

  /**
    Set the temperature ratio scale default is 1.e-5
    usually fine

    @param dTemperatureRatioScale temperature scale ratio
  */
  void SetTemperatureRatioScale(double dTemperatureRatioScale)
  {
    m_dTemperatureRatioScale = dTemperatureRatioScale;
  }

  /**
    Set the maximum number of accepted point before
    ending code

    @param nLimitAcceptances maximum number of accepted limit
  */
  void SetLimitAcceptance(long int nLimitAcceptances)
  {
    m_nLimitAcceptances = nLimitAcceptances;
  }

  /**
    Set the maximum limit of generated number of points
    before terminating

    @param nLimitGenerated maximum number of points generated
  */
  void SetLimitGenerated(long int nLimitGenerated)
  {
    m_nLimitGenerated = nLimitGenerated;
  }

  /**
    Set the maximum number of repeat when
    current objective function is less
    than tolerance and the abs difference
    between the last saved state and
    the best state is below the tolerance value
    as well
  */
  void SetMaximumCostRepeat(long int nMaximumCostRepeat)
  {
     m_nMaximumCostRepeat = nMaximumCostRepeat;
  }
  /**
    Perform Adaptive simulated annealing
     
     @param func is a template function that must take 
     @param pdParameter array containing first guess
  */
  template<typename T> 
   NumericError operator()(T& func, vector<double> &pdParameter) 
  { 
    m_nNbFunctionEvaluations = 0;

    if ( pdParameter.size() != m_nNumberParameters )    
      throw EXCEPTION_MSG
       (  
         ITO33_BAD_PARAM,
         TRANS("Constraint vector size different from number of variables.")
        );
    
    size_t nIdx;
    for (nIdx = 0; nIdx < m_nNumberParameters; nIdx++)
    {
      if ( pdParameter[nIdx] > m_pdUpperBound[nIdx] || 
           pdParameter[nIdx] < m_pdLowerBound[nIdx] )
         throw EXCEPTION_MSG
          (  
            ITO33_BAD_PARAM,
            TRANS("Error initial guess violates boundaries.")
          );
    }

    rand_seed = 696969;
    resettable_randflt(&rand_seed, 1);

    AsaState currentGeneratedState(m_nNumberParameters);
    AsaState lastSavedState(m_nNumberParameters);
    AsaState bestGeneratedState(m_nNumberParameters);

    currentGeneratedState.pdParameter = pdParameter;
    lastSavedState.pdParameter        = pdParameter;
    bestGeneratedState.pdParameter    = pdParameter;

    vector<double> pdDerivative(m_nNumberParameters, 0.0);
    double dMaxDerivative = 0.0; 
    double dAcceptedToGeneratedRatio = 0.0;

    // temperature parameters 
    double dCurrentCostTemperature = 0.0;
    double dInitialCostTemperature = 0.0;
    double dTemperatureScaleCost   = 0.0;
    vector<double> pdTemperatureScaleParameters(m_nNumberParameters, 0.0);
    vector<double> pdCurrentUserParameterTemp(m_nNumberParameters, 1.0);
    vector<double> pdInitialUserParameterTemp(m_nNumberParameters, 1.0);

    // counts of generated states and acceptances 
    vector<long int> pnIndexParameterGenerations(m_nNumberParameters, 1);
    long int nIndexCostRepeat         = 0;
    long int nNumberGenerated         = 0;
    long int nRecentNumberGenerated   = 0;
    long int nNumberAccepted          = 0;
    long int nRecentNumberAcceptances = 0;
    long int nIndexCostAcceptances    = 0;
    long int nReannealingCounter      = 0;
    long int nNumberAcceptanceSaved   = 0;
    
    outputFile <<"Initial Random Seed = "<< rand_seed << endl;
    outputFile << "number_parameters =" << m_nNumberParameters << endl;

    //create initial temperature by averaging a few runs
    //for the objective function
    double dF = 0.0;
    for (nIdx = 0; nIdx < size_t(m_nNumberSample); nIdx++) 
    {
      GenerateNewState(currentGeneratedState, lastSavedState,
                       pdCurrentUserParameterTemp);
      
      func(currentGeneratedState.pdParameter, currentGeneratedState.dCost);
      
      dF += fabs(currentGeneratedState.dCost);
    }

    dF = dF/double(m_nNumberSample);

    dInitialCostTemperature = dF;
    dCurrentCostTemperature = dF;

    // set all parameters to the initial parameter values 
    currentGeneratedState.pdParameter = pdParameter;
    lastSavedState.pdParameter        = pdParameter;
    bestGeneratedState.pdParameter    = pdParameter;

    // let asa generate valid initial parameters 
    GenerateNewState(currentGeneratedState, lastSavedState,
                       pdCurrentUserParameterTemp);
 
    func(currentGeneratedState.pdParameter, currentGeneratedState.dCost);
    m_nNbFunctionEvaluations = 1;
                            
    //set all states to the last one generated 
    bestGeneratedState.pdParameter = currentGeneratedState.pdParameter;
    lastSavedState.pdParameter     = currentGeneratedState.pdParameter;

    // set all costs to the last one generated 
    bestGeneratedState.dCost = currentGeneratedState.dCost;
    lastSavedState.dCost     = currentGeneratedState.dCost;

    dAcceptedToGeneratedRatio = 1.0;

    PrintState(dCurrentCostTemperature, 
               pdCurrentUserParameterTemp, 
               pdDerivative,
               dAcceptedToGeneratedRatio,
               nNumberAccepted, 
               nIndexCostAcceptances, 
               nNumberGenerated, 
               lastSavedState, 
               bestGeneratedState);


    while ( nNumberAccepted <= m_nLimitAcceptances &&  
            nNumberGenerated <= m_nLimitGenerated ) 
    {

      double tmp_var_db1 = -log(m_dTemperatureRatioScale);
      double tmp_var_db2 =  log( double(m_nTemperatureAnnealScale) );

      for (nIdx = 0; nIdx < m_nNumberParameters; nIdx++)
      {
        pdTemperatureScaleParameters[nIdx] = tmp_var_db1 
          * exp(-tmp_var_db2/m_dNumberParameters);
      }

      
     dTemperatureScaleCost = tmp_var_db1 * 
      exp( -tmp_var_db2/m_dNumberParameters)*m_dCostParameterScaleRatio;


    //Calculate new temperature T[i]_k
    for ( nIdx = 0 ; nIdx < m_nNumberParameters; nIdx++)
    {
      //skip parameters with too small range
      if ( m_pdRangeBound[nIdx] < m_dDoubleTolerance )
        continue;

      double log_new_temperature_ratio = -pdTemperatureScaleParameters[nIdx] *
        pow ( double(pnIndexParameterGenerations[nIdx]), 1./m_dNumberParameters);

      pdCurrentUserParameterTemp[nIdx] = pdInitialUserParameterTemp[nIdx]
        * exp(log_new_temperature_ratio);


      // check for too small a parameter temperature 
      if ( pdCurrentUserParameterTemp[nIdx] < m_dDoubleTolerance) 
      {
        pdParameter = bestGeneratedState.pdParameter;
        return  ITO33_TEMPERATURE_TOO_SMALL;
      }
    }

    // calculate new cost temperature T_k
    double log_new_temperature_ratio = -dTemperatureScaleCost * 
      pow( double (nIndexCostAcceptances) , 1./m_dNumberParameters );

    dCurrentCostTemperature = 
      dInitialCostTemperature * exp(log_new_temperature_ratio);

    // check for too small a cost temperature
    if ( dCurrentCostTemperature < m_dDoubleTolerance)
    {
      pdParameter = bestGeneratedState.pdParameter;
      return  ITO33_TEMPERATURE_TOO_SMALL;
    }

    // generate a new valid set of parameters     
    GenerateNewState(currentGeneratedState, lastSavedState,
                     pdCurrentUserParameterTemp);
 
    func(currentGeneratedState.pdParameter, currentGeneratedState.dCost);
    m_nNbFunctionEvaluations++;


    //decide to accept/reject the new state 
    AcceptNewState(dCurrentCostTemperature,
                   nRecentNumberAcceptances, 
                   nNumberAccepted,
                   nIndexCostAcceptances,
                   nNumberAcceptanceSaved,
                   nRecentNumberGenerated,
                   nNumberGenerated,
                   pnIndexParameterGenerations,
                   currentGeneratedState, 
                   lastSavedState);


    // calculate the ratio of acceptances to generated states
    dAcceptedToGeneratedRatio = double(nRecentNumberAcceptances + 1) /
                                double(nRecentNumberGenerated   + 1);
    
    
    // check for new minimum
    if ( currentGeneratedState.dCost < bestGeneratedState.dCost) 
    {     
      nRecentNumberAcceptances = 0;
      nRecentNumberGenerated   = 0;
      nIndexCostRepeat         = 0;

      // copy the current state into the best_generated state
      bestGeneratedState.dCost = currentGeneratedState.dCost;
      bestGeneratedState.pdParameter.swap(currentGeneratedState.pdParameter);

      if ( m_bHasOutput )
      {
        outputFile << "best...->cost=" << bestGeneratedState.dCost 
                   << ", number_accepted=" << nNumberAccepted 
                   << ", number_generated=" << nNumberGenerated << endl;
   
        if ( m_bHasOutputDetails )
         PrintState(dCurrentCostTemperature, 
                    pdCurrentUserParameterTemp, 
                    pdDerivative,
                    dAcceptedToGeneratedRatio,
                    nNumberAccepted, 
                    nIndexCostAcceptances, 
                    nNumberGenerated, 
                    lastSavedState, 
                    bestGeneratedState);
                  
      }

    }    


    //try to see if we have THE new minimum, abemus minimum
    if ( fabs(lastSavedState.dCost-bestGeneratedState.dCost) 
            < m_dDoubleTolerance   && 
            fabs(bestGeneratedState.dCost) < m_dTolerance )
    {
         nIndexCostRepeat++;

        if ( nIndexCostRepeat >= m_nMaximumCostRepeat )
        {
          pdParameter = bestGeneratedState.pdParameter;
          return ITO33_NO_ERROR;
         }

     }


    //periodic reannealing
    if ( ( (nNumberAccepted % m_nReannealingFrequency ) == 0 &&
            nNumberAcceptanceSaved == nNumberAccepted ) 
          || (nNumberGenerated % m_nGeneratedFrequency) == 0 )
    {

      nReannealingCounter++;

      ComputeDerivatives(func, pdDerivative, dMaxDerivative , currentGeneratedState,
                          bestGeneratedState);

      Reanneal(pdDerivative, dMaxDerivative,
               dCurrentCostTemperature, dInitialCostTemperature,
               dTemperatureScaleCost, pdCurrentUserParameterTemp,
               pdInitialUserParameterTemp, pdTemperatureScaleParameters, 
               nIndexCostAcceptances, pnIndexParameterGenerations,
               lastSavedState, bestGeneratedState);
 
      if ( m_bHasOutput )
        PrintState(dCurrentCostTemperature, 
                  pdCurrentUserParameterTemp, 
                  pdDerivative,
                  dAcceptedToGeneratedRatio,
                  nNumberAccepted, 
                  nIndexCostAcceptances, 
                  nNumberGenerated, 
                  lastSavedState, 
                 bestGeneratedState);

    }//reannealing

  } //end while loop
  
    pdParameter = bestGeneratedState.pdParameter;

    return ITO33_NO_ERROR;
  

  } //end ASA


  ///**
  //  Get the number of function evaluation

  //  @return the number of function evaluation
  //*/
  long int GetNbEvaluationFunction() const
  {
    return m_nNbFunctionEvaluations;
  }
  

private: 

  ///tolerance
  double m_dTolerance;

  ///lower bound vector
  vector<double> m_pdLowerBound;

  ///upper bound vector
  vector<double> m_pdUpperBound;

  ///range vector
  vector<double> m_pdRangeBound;

  ///Number of function evaluation
  long int m_nNbFunctionEvaluations;

  ///Number of unknows
  size_t m_nNumberParameters;

  ///double representation of unknowns
  double m_dNumberParameters;

  ///Set flag to indicate output
  bool m_bHasOutput;

  ///Set flag to indicate detailed output
  bool m_bHasOutputDetails;

  ///output file
  ofstream outputFile;

  ///how often the temperature are reannealed based on the derivatives
  long int m_nReannealingFrequency;

  ///how often reannealing is called based on the number of generated points
  long int m_nGeneratedFrequency;

  ///number of sampling point to get the first temperature
  long int m_nNumberSample;

  ///Scale guide to the expected cost temperature of convergence 
  ///within a small range of the global minimum
  double m_dTemperatureRatioScale;

  double m_dCostParameterScaleRatio;

  ///This scale is a guide to achieve the expected cost temperature
  long int m_nTemperatureAnnealScale;

  ///Maximum number of point accepted
  long int m_nLimitAcceptances;

  ///Maximum number of point generated
  long int m_nLimitGenerated;

  ///tolerance double
  double m_dDoubleTolerance;

  ///Number of repeated small state
  long int m_nMaximumCostRepeat;

  ///Random Number seed
  long int rand_seed;

  /// The state of the system in terms of parameters and function value 
  class AsaState 
  {
 
  public:
    double dCost;                //current obj function cost value
    vector<double> pdParameter;  //parameters value
   
  
    AsaState(size_t nSize)  
    {
      dCost = 0.0;
      pdParameter.clear();
      pdParameter.resize(nSize, 0.0); 
    }

    NO_COPY_CLASS(AsaState);
  };


  /**
    Generate a scalar according to ingber special distribution

    @param dTemperature current temperature
    
    @return new computed random scalar
  */
  double GenerateAsaState(double dTemperature)
  {
   double dX, dY, dZ;
   
   dX = randflt(&rand_seed);
   dY = dX < .5 ? -1. : 1.;
   dZ = dY * dTemperature*( pow( (1.+ 1./dTemperature), fabs(2.*dX-1.) ) - 1.);
   
   return dZ;
  }
  

  /**
   Generate a new set of valid solution for the optimization problem

   @param currentGeneratedState current generated state
   @param lastSavedState last saved state
   @param pdCurrentUserParameterTemp current parameter temperature
  */
  
  void GenerateNewState(AsaState &currentGeneratedState,
                        const AsaState &lastSavedState,
                        const vector<double> &pdCurrentUserParameterTemp)  
  {   

    // generate a new value for each parameter 
    size_t nIdx;
    size_t nFailMaxCounter = 10;
    for ( nIdx = 0 ; nIdx < m_nNumberParameters; nIdx++)  
    { 
      // ignore parameters that have too small a range 
      if (fabs(m_pdRangeBound[nIdx]) < m_dDoubleTolerance)
        continue;

      double dTemperature  = pdCurrentUserParameterTemp[nIdx];
      double dParameter    = lastSavedState.pdParameter[nIdx];

      // generate a new state x within the parameter bounds 
      double dX = 0.;
      size_t nFailCounter = 0;

      for (;;)
      {
        nFailCounter++;

        dX = dParameter + GenerateAsaState(dTemperature)*m_pdRangeBound[nIdx];

        // exit the loop if within its valid parameter range 
        if ( dX <= m_pdUpperBound[nIdx] && dX >= m_pdLowerBound[nIdx])
          break;

        if (nFailCounter >= nFailMaxCounter)
        {
          dX = dParameter;
          break;
        }

       
     }

    // save the newly generated value 
    currentGeneratedState.pdParameter[nIdx] = dX;

    }//loop over parameters

  } // GenerateNewState


  /**
   Accept or reject a newly created state

   @param dCurrentCostTemperature current globa temperature
   @param nRecentNumberAcceptances current number of accepted point
   @param nNumberAccepted total number of accepted point
   @param nIndexCostAcceptances k
   @param nRecentNumberGenerated current number of point generated
   @param nNumberGenerated total number of point generated
   @param pnIndexParameterGenerations k_i
   @param currentGeneratedState new found state
   @param lastSavedState last saved state
  */
  void AcceptNewState(double dCurrentCostTemperature,
                      long int &nRecentNumberAcceptances,
                      long int &nNumberAccepted,
                      long int &nIndexCostAcceptances,
                      long int &nNumberAcceptancesSaved,
                      long int &nRecentNumberGenerated,
                      long int &nNumberGenerated,
                      vector<long int> &pnIndexParameterGenerations,
                      const AsaState &currentGeneratedState,
                      AsaState &lastSavedState)
  {
    // update accepted and generated count
    nNumberAcceptancesSaved++;
    nRecentNumberGenerated++;
    nNumberGenerated++;
    
    size_t nIdx;
    for (nIdx = 0 ; nIdx < m_nNumberParameters; nIdx++)
    {
      if ( m_pdRangeBound[nIdx] > m_dDoubleTolerance)
        pnIndexParameterGenerations[nIdx]++;
    }

    // effective cost function for testing acceptance criteria,
    //   calculate the cost difference and divide by the temperature 
    double dDeltaCost = (currentGeneratedState.dCost - lastSavedState.dCost)
                      / (dCurrentCostTemperature + m_dDoubleTolerance);

    double dProbTest = min( 1., exp( -dDeltaCost) );
   
    // accept/reject the new state 
    if ( dProbTest >= randflt(&rand_seed) ) 
    {
      lastSavedState.dCost       = currentGeneratedState.dCost;
      lastSavedState.pdParameter =  currentGeneratedState.pdParameter;

     // update acceptance counts 
     nRecentNumberAcceptances++;
     nNumberAccepted++;
     nNumberAcceptancesSaved = nNumberAccepted;
     nIndexCostAcceptances++;
    }
  } //Accept New State

  /**
    Debug output function

    @param dCurrentCostTemperature current temperature cost
    @param pdCurrentUserParameterTemp temperature of each parameters
    @param pdDerivatives current derivatives
    @param dAcceptedToGeneratedRatio point accepted versus generated
    @param nNumberAccepted number of point accepted
    @param nIndexCostAcceptance k
    @param nNumberGenerated number of point generate
    @param lastSavedState last best point found
    @param bestGeneratedState best point found so far

  */
 
  void PrintState(double dCurrentCostTemperature, 
                  const vector<double> &pdCurrentUserParameterTemp, 
                  const vector<double> &pdDerivatives,
                  double dAcceptedToGeneratedRatio, 
                  long int nNumberAccepted, 
                  long int nIndexCostAcceptance, 
                  long int nNumberGenerated, 
                  const AsaState &lastSavedState, 
                  const AsaState  &bestGeneratedState) 
  {
    outputFile << endl << endl;
    
    outputFile << "index_cost_acceptances = " << nIndexCostAcceptance << 
                ", current_cost_temperature = " << dCurrentCostTemperature << endl;
    
    outputFile << "accepted_to_generated_ratio = " 
               << dAcceptedToGeneratedRatio << endl;

    outputFile << "number_generated = " << nNumberGenerated <<
                ", number_accepted  = " << nNumberAccepted << endl;;
   
    outputFile <<"best...->cost = " << bestGeneratedState.dCost 
               <<", last...->cost = " << lastSavedState.dCost << endl;

    outputFile <<"# best...->parameter \t current_parameter_temp \t derivative" << endl;

    size_t nIdx;
    for ( nIdx = 0 ; nIdx < m_nNumberParameters; nIdx++)
    {
      outputFile << nIdx <<  "\t" << bestGeneratedState.pdParameter[nIdx]
                 << "\t\t" << pdCurrentUserParameterTemp[nIdx] 
                 << "\t\t" <<pdDerivatives[nIdx]<<endl;
    }

    outputFile << endl;

  } //PrintState

  
  /**
    Compute sensitivity of parameters with temperature
    also recent the current generated state to the best saved state

    @param func addresse of objective function
    @param pdDerivatives array of derivative to compute
    @param currentGeneratedState current generated state
    @param bestSavedState best state found so far

  */
  template<typename T> 
    void ComputeDerivatives(T &func, vector<double> &pdDerivative, 
                           double &dMaxDerivative,
                           AsaState &currentGeneratedState,
                           const AsaState &bestGeneratedState )
  { 
    double dRecentBestCost = bestGeneratedState.dCost;

    currentGeneratedState.pdParameter = bestGeneratedState.pdParameter;
    currentGeneratedState.dCost       = dRecentBestCost;

    // calculate tangents 
    size_t nIdx;
    for (nIdx = 0; nIdx < m_nNumberParameters; nIdx++)
    {
        
      if ( m_pdRangeBound[nIdx] < m_dDoubleTolerance)        
      {
        pdDerivative[nIdx] = .0;
        continue;
      }
           
      // save the v_th parameter and delta_parameter 
      double dParameter      = bestGeneratedState.pdParameter[nIdx];
      double dDeltaParameter = 1.e-3;

      double dParameterOffset = (1. + dDeltaParameter) * dParameter;

      if ( dParameterOffset > m_pdUpperBound[nIdx] ||
           dParameterOffset < m_pdLowerBound[nIdx]) 
      {
        dDeltaParameter  = -dDeltaParameter;
        dParameterOffset = (1. + dDeltaParameter) * dParameter;
      }

      // generate the first sample point 
      currentGeneratedState.pdParameter[nIdx] = dParameterOffset;

      func(currentGeneratedState.pdParameter, currentGeneratedState.dCost);

      // restore the parameter state 
      currentGeneratedState.pdParameter[nIdx] = dParameter;

      // calculate the numerical derivative *
      pdDerivative[nIdx] = (currentGeneratedState.dCost - dRecentBestCost) / 
                           (dDeltaParameter * dParameter + m_dDoubleTolerance);

    } //end loop over point


    // find the maximum |tangent| from all tangents 
    dMaxDerivative = 0.;

    for( nIdx = 0 ; nIdx < m_nNumberParameters; nIdx++) 
    {
       if ( m_pdRangeBound[nIdx] < m_dDoubleTolerance)        
        continue;

      // find the maximum |tangent| (from all tangents) 
      if ( fabs(pdDerivative[nIdx]) > dMaxDerivative ) 
        dMaxDerivative = fabs( pdDerivative[nIdx] ); 
    }
   
    // restore the best cost function value 
    currentGeneratedState.dCost = dRecentBestCost;

  }//end ComputeDerivatives


  /**
   Reanneal the different parameters 
   depending on the sensitivity 

   @param nReannealCounter number of reannealing
   @param pdDerivative vector of computed derivatives
   @param dCurrentCostTemperature T
   @param dInitialCostTemperature To
   @param dTemperatureScaleCost   C
   @param pdCurrentUserParameterTemp T_i
   @param pdInitialUserScaleParameters T_oi
   @param pdTemperatureScaleParameters C_i
   @param nIndexCostAcceptances k
   @param pnIndexParameterGeneration k_i
   @param lastSavedState last saved state
   @param bestGeneratedState best generated state
  */
  void Reanneal(const vector<double> &pdDerivative, 
                double dMaxDerivative, 
                const double &dCurrentCostTemperature, 
                double &dInitialCostTemperature, 
                const double &dTemperatureScaleCost, 
                const vector<double> &pdCurrentUserParameterTemp,
                vector<double> &pdInitialUserParameterTemp, 
                const vector<double> &pdTemperatureScaleParameters,
                long int &nIndexCostAcceptances, 
                vector<long int> &pnIndexParameterGenerations,
                const AsaState &lastSavedState, 
                const AsaState &bestGeneratedState)
  {
  
    size_t nIdx;
    for (nIdx = 0; nIdx < m_nNumberParameters; nIdx++)
    {
      // use the temp double to prevent overflow 
      double dTmp = 0.0;

      // skip parameters with too small range or integer parameters 
      if ( m_pdRangeBound[nIdx] < m_dDoubleTolerance )
        continue;
      
      // ignore parameters with too small tangents 
      if ( fabs( pdDerivative[nIdx] ) < m_dDoubleTolerance)
        continue;

      // reset the index of parameter generations appropriately 
      double dNewTemperature = fabs ( 
        pdCurrentUserParameterTemp[nIdx]* dMaxDerivative/pdDerivative[nIdx] );

      if ( dNewTemperature < pdInitialUserParameterTemp[nIdx] ) 
      {
        double log_init_cur_temp_ratio =
          fabs( log( ( m_dDoubleTolerance  + pdInitialUserParameterTemp[nIdx] )
                   / ( m_dDoubleTolerance + dNewTemperature) ) );

        dTmp =  m_dDoubleTolerance + pow( log_init_cur_temp_ratio
                   / pdTemperatureScaleParameters[nIdx], m_dNumberParameters);
      } 
      else 
      {
        dTmp = 1.0;
      }

      // Reset index_parameter_generations if index reset too large,
      //   and also reset the initial_user_parameter_temp, to achieve
      //   the same new temperature. 
      while ( dTmp > 50000. ) 
      {
        double log_new_temperature_ratio = -pdTemperatureScaleParameters[nIdx] 
          * pow (dTmp, 1./m_dNumberParameters);

        dNewTemperature = pdInitialUserParameterTemp[nIdx] *
          exp( log_new_temperature_ratio );
        
        dTmp /= 10.;

        double temperature_rescale_power = 1./pow(10., 1./m_dNumberParameters);

        pdInitialUserParameterTemp[nIdx] = dNewTemperature * 
          pow ( pdInitialUserParameterTemp[nIdx] /
                dNewTemperature, temperature_rescale_power);
      }
      
      // restore from temporary double 
      pnIndexParameterGenerations[nIdx] = long(dTmp);
    }

    double cost_best = bestGeneratedState.dCost;
    double cost_last = lastSavedState.dCost;

    // (re)set the initial cost_temperature 
    double dTmp = max( fabs(cost_last), fabs(cost_best) );
    dTmp = max( dTmp, fabs(cost_best - cost_last) );
    dTmp = max( m_dDoubleTolerance, dTmp);
    dInitialCostTemperature = min( dInitialCostTemperature, dTmp);

    dTmp = double(nIndexCostAcceptances);

    double dTmp1 = max(fabs(cost_last - cost_best), dCurrentCostTemperature);
    dTmp1 = max( m_dDoubleTolerance, dTmp1);
    dTmp1 = min( dTmp1, dInitialCostTemperature);


    if ( dCurrentCostTemperature > dTmp1 ) 
    {
      double tmp_var_db3 =
        fabs( log( (m_dDoubleTolerance + dInitialCostTemperature) / dTmp1 ) );

      dTmp = m_dDoubleTolerance +
        pow( tmp_var_db3 / dTemperatureScaleCost, m_dNumberParameters);

    } 
    else 
    {
      double log_init_cur_temp_ratio =
        fabs( log( ( m_dDoubleTolerance + dInitialCostTemperature) /
                   ( m_dDoubleTolerance + dCurrentCostTemperature) ) );

      dTmp = m_dDoubleTolerance 
        + pow( log_init_cur_temp_ratio/ dTemperatureScaleCost, m_dNumberParameters);

    }

    // reset index_cost_temperature if index reset too large 
    while (dTmp > 50000 ) 
    {
      double log_new_temperature_ratio = -dTemperatureScaleCost 
        * pow( dTmp,1. / m_dNumberParameters);

      double new_temperature = 
        dInitialCostTemperature * exp(log_new_temperature_ratio);

      dTmp /= 10.;

      double temperature_rescale_power = 
        1. / pow ( 10. , 1./ m_dNumberParameters);

      dInitialCostTemperature = new_temperature * 
        pow(dInitialCostTemperature / new_temperature, temperature_rescale_power);
    }
    
    
    nIndexCostAcceptances = long(dTmp);

 }

 
 
  /* This RNG is a modified algorithm of that presented in
   * %A K. Binder
   * %A D. Stauffer
   * %T A simple introduction to Monte Carlo simulations and some
   *    specialized topics
   * %B Applications of the Monte Carlo Method in statistical physics
   * %E K. Binder
   * %I Springer-Verlag
   * %C Berlin
   * %D 1985
   * %P 1-36
   * where it is stated that such algorithms have been found to be
   * quite satisfactory in many statistical physics applications. */
   double myrand (long int *rand_seed)
  {                   
    *rand_seed = long( (25173 * (*rand_seed) + 13849) % 65536);
    return ((double) (*rand_seed) / 65536.0); 
  }




  double randflt (long int * rand_seed)
  {
    return (resettable_randflt (rand_seed, 0));
  }

  double resettable_randflt (long int * rand_seed, int reset) 
  /* shuffles random numbers in random_array[SHUFFLE] array */
  {
    double rranf;
    unsigned kranf;
    int n;
    static int initial_flag = 0;
    long int initial_seed;
    static double random_array[256];  /* random variables */


    if (*rand_seed < 0)
      *rand_seed = -*rand_seed;

    if ((initial_flag == 0) || reset) 
    {
      initial_seed = *rand_seed;

      for (n = 0; n < 256; ++n)
        random_array[n] = myrand (&initial_seed);

      initial_flag = 1;

      for (n = 0; n < 1000; ++n)  /* warm up random generator */
        rranf = randflt (&initial_seed);

      rranf = randflt (rand_seed);

      return (rranf);
    }

    kranf = (unsigned) (myrand (rand_seed) * 256) % 256;
    rranf = *(random_array + kranf);
    *(random_array + kranf) = myrand (rand_seed);

    return (rranf);
  }


  NO_COPY_CLASS(ASA);

};


  
}  // namespace numeric 

}  // namespace ito33

#endif // _ITO33_NUMERIC_ADAPTIVE_SIMULATED_ANNEALING_H_

