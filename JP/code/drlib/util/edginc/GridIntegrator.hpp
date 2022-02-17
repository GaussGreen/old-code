//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : GridIntegrator.hpp
//
//   Description : 
//
//   Author      : Jakob Sidenius
//
//   Date        : 5 October 2006
//
//----------------------------------------------------------------------------
#include "edginc/Object.hpp"
#include "edginc/Function.hpp"
#include <cmath>
#include <map>
#include <set>

#ifndef EDR_GRIDINTEGRATOR_HPP
#define EDR_GRIDINTEGRATOR_HPP

DRLIB_BEGIN_NAMESPACE

using std::map;
using std::set;

#define TRY
#define CATCH(arg1)



/*! 
  \brief A GridIntegrator is a quadrature rule

  This type holds the definition of a quadrature rule and
  can apply this to anything which allows function call syntax with argument type double and 
  return type double or valarray<double>.

  Different algorithm for setting nodes and weights are implemented using the 'named ctor' idiom,
  ie, static member functions which return a GridIntegrator. 
*/
class UTIL_DLL GridIntegrator
{
public:
// foundation
  //! Trivial default ctor (should not be allowed - blame QLib design...)
  GridIntegrator() {}

  //! Construct from nodes and weights
  GridIntegrator(const map<double,double>& nodesandweights);

// named 'ctors'
  //! Construct a GridIntegrator representing Simpson quadrature with n nodes over [a,b]
  static GridIntegrator Simpson(int n,
                                double a,
                                double b);

  /*! Construct a GridIntegrator representing Simpson quadrature with n nodes over [a,b].
      This overload supports functions with discontinuities by computing a correction to
      the "naive" simpson rule. There is one additive correction for each discontinuity,
      each term being proportional to the size of the discontinuity. 

      The size of the discontinuity at a point x is found as the difference between 
      integrand values at x+endpointFuzz and x-endpointFuzz.

      Note that the nodes of the standard Simpson are unaffected by the location and 
      number of discountinuities.

      This functionality is intended for use for functions with discontinuities.   

      Note that 'discSet' entries outside of [a,b] are ignored.
  */
  static GridIntegrator Simpson(int n,
                                double a,
                                double b,
                                const set<double>& discSet,
                                double endpointFuzz = 1.0e-10);

  /*! Construct a GridIntegrator representing Simpson quadrature with n nodes over [a,b].
      This ctor splits the interval [a,b] into disjoint subintervals with 
      endpoints in 'discSet' and uses Simpson on each subinterval. 
      
      Note that endpoints of subintervals are actually shifted into the interval by 'endpointFuzz'
      such that a (small) positive value of this parameter serves to avoid valuation at these
      points.

      This functionality is intended for use for functions with singularities.   

      Note that 'discSet' entries outside of [a,b] are ignored.
  */
/*  static GridIntegrator SimpsonWithSplit(int n,
                                         double a,
                                         double b,
                                         const set<double>& singSet,
                                         double endpointFuzz = 1.0e-10);
*/


// inspection
  //! Return nodes and weights
  const map<double,double>& NodesAndWeights() const { return m_nodesandweights; }
  //! Return integral of function
  //! Use this for functions returning double
  template<typename Func>
  double operator()(const Func& f) const;
  //! Use this overload for functions returning a vector<double>
  //! Integral of vector function will be written to 'res', which MUST 
  //! be initialized to zero and have the correct size on input.
  template<typename Func>
  void operator()(const Func& f, vector<double>& res) const;
  //! Use this overload for functions represented by MFunctionND
  //! Integral of vector function will be added component-wise to 'res'
  //! Note that the range of f must be one-dimensional for now. 
  void operator()(const MFunctionND& f, vector<double>& res) const;

// modification
  //! If a factor in the integrand is available separately, this function can be used to apply
  //! it to the quadrature weights before doing the quadrature. 
  //! Especially relevant for integration of vector functions with a common factor for all components.
  template<typename W>
  void ScaleWeightsByWeightingFunction(const W weightFunction);

  //! This method 'un-applies' a factor set in ScaleWeightsByWeightingFunction()
  template<typename W>
  void DivideWeightsByWeightingFunction(const W weightFunction);

private:
// data
  // Nodes and weights
  map<double,double> m_nodesandweights;
}; // GridIntegrator


template<typename W>
void GridIntegrator::ScaleWeightsByWeightingFunction(const W weightFunction)
{
  map<double,double>::iterator iNode = m_nodesandweights.begin();
  for (; iNode != m_nodesandweights.end(); ++iNode)
    iNode->second *= weightFunction(iNode->first);
}

template<typename W>
void GridIntegrator::DivideWeightsByWeightingFunction(const W weightFunction)
{
  map<double,double>::iterator iNode = m_nodesandweights.begin();
  for (; iNode != m_nodesandweights.end(); ++iNode)
    iNode->second /= weightFunction(iNode->first);
}

inline 
void GridIntegrator::operator()(const MFunctionND& f, vector<double>& res) const
{
TRY
#ifdef _DEBUG
  double ws = 0.0;
  map<double,double>::const_iterator iNode = m_nodesandweights.begin();
  for (; iNode != m_nodesandweights.end(); ++iNode)
    ws += iNode->second;
#endif
  // Container to hold argument
  CDoubleArray arg(1);
  // Container to hold node value
  CDoubleArray fp(res.size());
  // Loop over nodes  
  map<double,double>::const_iterator iNode = m_nodesandweights.begin();
  for (; iNode != m_nodesandweights.end(); ++iNode)
  {
    arg[0] = iNode->first;
    // Get new point value
    f(arg, fp);
    // Check
/*    if (fp.size() != res.size())
      throw Logical("Size mismatch."); */
    // Loop over vectors
    vector<double>::const_iterator iF = fp.begin();
    vector<double>::iterator iVal = res.begin();
    for (; iVal != res.end(); ++iVal, ++iF)
      *iVal += iNode->second * *iF;
  }
CATCH(GridIntegrator::operator())
}


template<typename Func>
void GridIntegrator::operator()(const Func& f, vector<double>& res) const
{
TRY
#ifdef _DEBUG
  double ws = 0.0;
  map<double,double>::const_iterator iNode = m_nodesandweights.begin();
  for (; iNode != m_nodesandweights.end(); ++iNode)
    ws += iNode->second;
#endif
  // Loop over nodes  
  map<double,double>::const_iterator iNode = m_nodesandweights.begin();
  for (; iNode != m_nodesandweights.end(); ++iNode)
  {
    // Get new point value
    const vector<double>& fp = f(iNode->first);
    // Check
/*    if (fp.size() != res.size())
      throw Logical("Size mismatch."); */
    // Loop over vectors
    vector<double>::const_iterator iF = fp.begin();
    vector<double>::iterator iVal = res.begin();
    for (; iVal != res.end(); ++iVal, ++iF)
      *iVal += iNode->second * *iF;
  }
CATCH(GridIntegrator::operator())
}

template<typename Func>
double GridIntegrator::operator()(const Func& f) const
{
  double res = 0.0;

  map<double,double>::const_iterator iNode = m_nodesandweights.begin();
  for (; iNode != m_nodesandweights.end(); ++iNode)
    res += iNode->second * f(iNode->first);
  return res;
}

/*!
  \brief Trapezoidal rule.

  This routine computes the nth stage of refinement of an extended trapezoidal rule. 
  When called with n=1, the routine returns the crudest estimate of \f$\int_a^bf(x)dx\f$. 
  Subsequent calls with n=2,3,... (in that sequential order) will improve the accuracy 
  by adding 2n-2 additional interior points.

  \par Template parameters
  \param Func     Type of function to be integrated.
  \param AreaType Type of integral value.

  \par Function arguments
  \param f  Function to be integrated.
  \param a  Lower integration limit.
  \param b  Upper integration limit.
  \param n  Order of refinement.  
*/
template<class Func, typename AreaType>
AreaType TrapezoidalRule(const Func& f, 
                         typename Func::ArgumentType a, 
                         typename Func::ArgumentType b, 
                         unsigned int n, 
                         AreaType& s)
{
  if (n == 1) 
    return (s = 0.5 * (b-a) * (f(a) + f(b)));
  else 
  {
    unsigned int it,j;
    for (it = 1, j = 1; j < n - 1; ++j) 
      it <<= 1;
    // This is the spacing of the points to be added.
    typename Func::ArgumentType del = (b - a) / it; 
    typename Func::ArgumentType x = a + 0.5 * del;
    typename Func::ReturnType sum = 0.0;
    for (j = 1; j <= it; ++j, x += del) 
      sum += f(x);
    // This replaces s by its refined value.
    s = 0.5 * (s + (b - a) * sum / it); 
    return s;
  }
}




DRLIB_END_NAMESPACE

#endif

