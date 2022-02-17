//----------------------------------------------------------------------------
//
//   Group       : QR Credit
//
//   Filename    : OptimizerSpellucci.hpp
//
//   Description : Wrapper class for Spellucci's C code for the DONLP2 minimization algorithm
//                 This wrapper supports only unconstrained optimization
//
//   Date        : 31 August 2006
//
//----------------------------------------------------------------------------

#ifndef EDR_OPTIMIZERSPELLUCCI_HPP
#define EDR_OPTIMIZERSPELLUCCI_HPP

#include <valarray>
#include <vector>
#include <map>
#include "Optimizer.hpp"
/*
#include "JAKOB_matrix.h"
#include "JAKOB_linalg.h"
#include "JAKOB_norm.h"
#include "JAKOB_gradient.h"
*/

DRLIB_BEGIN_NAMESPACE

// Forward declaration of QSpellucci
class QSpellucci;
typedef smartPtr<QSpellucci> QSpellucciSP;
typedef smartConstPtr<QSpellucci> QSpellucciConstSP;
#ifndef QLIB_OPTIMIZER_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<QSpellucci>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<QSpellucci>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<QSpellucci>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<QSpellucci>);
#endif

/** Wrapper around JS's Spellucci class */
class UTIL_DLL QSpellucci : public OptimizerND
{
public:
    static CClassConstSP const TYPE;

// foundation
    ~QSpellucci(); 

		/*! 
    Default construction

    \param del0     Parameter used to determine binding inequality constaints 
                    in the initial phase of the %optimization. A constraint is 
                    considered binding if
                    \f[ g_i(x) / \max \{ 1, \|\nabla g_i(x) \| \} 
                    \leq del0 \f]
    \param tau0     Bound determine how much the unscaled penalty-term may
                    derivate from zero. The algorithm assumes that within the
                    region described by
                    \f[ \sum_{i=1}^{N_H} |h_i(x)| - 
                    \sum_{i=1}^{N_G} \min\{ 0, g_i(x) \} \leq tau0 \f]
                    all functions may be evaluated safely.
    \param taubnd   Amount by which bounds may be violated if numerical 
                    differentation is used.
    \param epsdif   The expected relative precision of the gradient evaluation.
    \param epsfcn   The expected relative precision of the object function 
                    evaluation.
    \param difftype The numerical differention schema: Valid difftypes 
                    are 1, 2, and 3. Numerical differentation uses \f$n\f$, 
                    \f$2n\f$, and \f$6n\f$ function evaluations for each 
                    gradient if difftype is 1, 2, 3 respectively.
                    - 1: Ordinary forward difference quotient with step size 
                        \f[ 0.1 epsfcn^{\frac{1}{2}}.\f]
                    - 2: Symmetric difference quotient with step size 
                        \f[ 0.1 epsfcn^{\frac{1}{3}}.\f]
                    - 3: Sixth order approximation based on Richardson 
                        extrapolation with step size
                        \f[ 0.1 epsfcn^{\frac{1}{7}}.\f]
		\param iterma     The maximal number of main loop iterations.
		\param epsx       One of the termination criteria. 
                    \f[ \|\nabla L(x,\mu,\lambda) \| 
                    \leq epsx(1 + \|\nabla f(x) \|), \f]
                    with \f$ L \f$ being the Lagrange function corresponding to 
                    the constrained %optimization problem.
                    This is not the only termination criterion but the others
                    are purely internal.
		*/
		QSpellucci(double del0     = 0.2e0,  
			 	       double tau0     = 1.e0,
							 double taubnd   = 1.0,
               double epsdif   = 1.e-16,
               double epsfcn   = 1.e-16,
               int    difftype = 1,
               int    iterma   = 1000,
               double epsx     = 1.e-8);

    void minimize(const MFunctionND&  func,
                  const CDoubleArray& guess,        // initial guess
                  CDoubleArray&       x) const;     // result

    void minimize(const MFunctionND&  func,
                  const CDoubleArray& xguess,        // initial guess
                  const CStringArray& ids,           // identifiers
                  CDoubleArray&       x) const;      // result

private:
// types
  //! Forward declaration of Spellucci (implemented in the .cpp file)
  class Spellucci;

  static void load(CClassSP& clazz);

  static IObject* defaultCtor();
// data 
  //! Optimizer 
  Spellucci* mp_spellucci;
}; // QSpellucci


DRLIB_END_NAMESPACE

#endif



