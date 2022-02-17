//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Doxygen.hpp
//
//   Description : Doxygen documentation for the EDR library
//
//   Date        : Sep 2004
//
//
//----------------------------------------------------------------------------

#ifndef DOXYGEN_HPP
#define DOXYGEN_HPP

#include "edginc/config.hpp"
#include <vector>

DRLIB_BEGIN_NAMESPACE

//! Monte carlo path generator for scalar lognormal process
/** A lognormal simulation of a single equity price with constant
    drift and volatility. The Stochastic Differential Equation for
    this model is 
    \f$ dS_t = \mu S_t dt + \sigma S_t dW_t \f$ where 
    \f$ \mu \f$ is the equity drift, 
    \f$ \sigma \f$ is the equity volatility and 
    \f$ W_t \f$ is a standard Brownian Motion. */
class TOOLKIT_DLL MCPathGenLNDoxygen {
public:
    //! Class constructor
    /** The inputs are: 
        \param spot is the equity price at simulation start, 
        \param dirft is the continuously compounded equity growth rate, 
        \param vol is the annualized equity volatility and \param dates is 
        the requested set of simulatiod dates */
    MCPathGenLNDoxygen(double spot, double drift, double vol, const vector<double>& dates);

    //! Virtual destructor
    virtual ~MCPathGenLNDoxygen() {}

    //! Generates the lognormal path
    /** Draws a set of standard normal variables using \fn standardNormal() 
        and diffuses the equity path according to the lognormal SDE */
    void generatePath();

    //! Gives access to path member variable
    const vector<double>& getPath() const;
private:
    //! Draws a standard normal random variable
    /** Based on gasdev taken from Numerical Recipes in C p.289.
        It uses random rejection of uniform random numbers. 
        \sa unifrom() */
    static double standardNormal();
    
    //! Draws a uniform random variable
    static double uniform();

    // Mandatory fields
    double spot;            //!< Equity spot price at simulation start \f$ S_0 \f$
    double drift;           //!< Equity drift i.e. interest rate - dividend (continuously compounded)
    double vol;             //!< Equity volatility
    vector<double> dates;   //!< Simulation dates
    
    // Transient fields
    vector<double> path;    //!< Equity spot price path
};

DRLIB_END_NAMESPACE

#endif
