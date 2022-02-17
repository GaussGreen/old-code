//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : FindRootBrent.hpp - from CMLib FindRoot.h
//
//   Description : Finds roots of functions using Brent methodology
//
//   Author      : CMLib
//
//   Date        : Dec 17, 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_FINDROOTBRENT_HPP
#define EDR_FINDROOTBRENT_HPP

#include "edginc/FindRoot.hpp"

DRLIB_BEGIN_NAMESPACE

/** For finding roots of functions using Brent methodology - see the
    solve method below */
class FindRootBrent: public FindRoot{
public:
    class UTIL_DLL Function {
    public:
        // function sets isWithinDomain to true if argument is within domain
        virtual double operator()(double arg,
                                  bool&  isWithinDomain) = 0;
        
#ifdef CM_STORE_SOLVER_POINTS
        vector<std::pair<double,double> > m_point; // debug output
#endif
    };
    
    /** To be able to pass anything, even C function pointers, to
        the Brent root finder, the following template is used
        by the root finder to wrap the functionoid that's passed in. */
    template<class F>
    class Proxy : public Function {
    public:
        Proxy(F& function): m_function(function) {}
        
        double operator()(double arg,
                          bool&  isWithinDomain) {
            isWithinDomain = true;
            double y = m_function(arg);
#ifdef CM_STORE_SOLVER_POINTS
            m_point.push_back(make_pair(arg,y));
#endif
            return y;
        }
        
    private:
        F& m_function;
    };


    /** This is the version of the Brent root finder which takes
        any functionoid, even C function pointers. The functionoid must
        have the signature double f( double ). */
    template<class F> static Point solve(
        F&     function,                        // (I) Function to be solved
        double guess,                       // (I) Initial guess
        double lowerBound                   // (I) Result cannot be smaller
        = LOWER_BOUND,                                                
        double upperBound                   // (I) Result cannot be larger
        = UPPER_BOUND,                                        
        double xAccuracy                    // (I) X accuracy tolerance
        = X_TOLERANCE,                                        
        double fAccuracy                    // (I) Function accuracy tolerance
        = F_TOLERANCE,        
        int    maxNumOfIterations           // (I) Maximum number of iterations
        = MAX_NUM_OF_ITERATIONS,
        double initialXStepSize             // (I) Size of step in x 
        = INITIAL_X_STEP_SIZE) {
        
        Proxy<F> proxy(function); // Wrap functionoid into base proxy.
        
        return solveForFunction(static_cast<Function&>(proxy),
                                guess,        
                                lowerBound, upperBound,
                                xAccuracy, fAccuracy,
                                maxNumOfIterations,
                                initialXStepSize);
    }
    //// Was named just solve() but MSVC fell over 
    static UTIL_DLL Point solveForFunction(
        Function& function,                  // (I) Function to be solved
        double    guess,                    // (I) Initial guess
        double    lowerBound,               // (I) Result cannot be smaller
        double    upperBound,               // (I) Result cannot be larger
        double    xAccuracy,                // (I) X accuracy tolerance
        double    fAccuracy,                // (I) Function accuracy tolerance
        int       maxNumOfIterations,       // (I) Maximum number of iterations
        double    initialXStepSize);        // (I) Derivative, defaults to 0

private:
    class Finder;
    FindRootBrent(); // do not create instances of this class
};

DRLIB_END_NAMESPACE

#endif
