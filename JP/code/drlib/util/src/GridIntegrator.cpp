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
//   Todo        : Exceptions
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GridIntegrator.hpp"

DRLIB_BEGIN_NAMESPACE

GridIntegrator::GridIntegrator(const map<double,double>& nodesandweights)
TRY
: m_nodesandweights(nodesandweights)
{
  // Check sizes
/*  if (m_nodes.size() != m_weights.size())
    throw Constructor("Unequal numbers of quadrature nodes and weights.");
  if (m_nodes.size() == 0)
    throw Constructor("No quadrature nodes given."); */
}
CATCH(GridIntegrator::GridIntegrator)


GridIntegrator 
GridIntegrator::Simpson(const int n,
                        const double a,
                        const double b)
TRY
{
/*  if (b <= a)
    throw Argument("Upper limit less than lower limit.");
  if (n < 2)
    throw Argument("Less than two points requested for Simpson quadrature."); */
  // Round n up to nearest odd number
  const int N = n + 1 - (n % 2); 
  const double h = (b - a) / (N - 1);

  // Set nodes and weights
  map<double,double> nodesandweights;
  nodesandweights[a] = h / 3.0;
  double node = a;
  for (unsigned int i = 1; i < N-1; ++i)
  {
    node += h;
    double weight = h / 3.0;
    weight *= ((i+1) % 2) ? 4.0 : 2.0; 
    nodesandweights[node] = weight;
  }  
  nodesandweights[b] = h / 3.0;

#ifdef _DEBUG
  double ws = 0.0;
  map<double,double>::const_iterator iNode = m_nodesandweights.begin();
  for (; iNode != m_nodesandweights.end(); ++iNode)
    ws += iNode->second;
#endif

  return GridIntegrator(nodesandweights);
}
CATCH(GridIntegrator::Simpson)


GridIntegrator 
GridIntegrator::Simpson(const int n,
                        const double a,
                        const double b,
                        const set<double>& discSet,
                        const double endpointFuzz)
{
  // Standard Simpson
  GridIntegrator ss = GridIntegrator::Simpson(n, a, b);
  map<double,double> nw = ss.NodesAndWeights();

  // Get distance between 'small' nodes
  map<double,double>::const_iterator iD_second = nw.begin();
  ++iD_second;
  const double h = iD_second->first - nw.begin()->first;

  // Now insert two nodes per discontinuity (one on either side) and corresponding weights
  // Find first discontinuity >= a
  set<double>::const_iterator iD_begin = discSet.lower_bound(a);
  // Find first discontinuity > b
  set<double>::const_iterator iD_end = discSet.upper_bound(b);
  for (set<double>::const_iterator iD = iD_begin; iD != iD_end; ++iD)
  {
    // Discard discontinuities at end points
    if (*iD - a < endpointFuzz || b - *iD < endpointFuzz)
      continue;

    // Find index of Simpson node immediately to the left of discontinuity
    const double d = floor(*iD / h);
    // Compute weight
    double w = h * (d + 1.0/3.0) - *iD;

    // Insert nodes and weights
    nw[*iD - endpointFuzz] = -w;
    nw[*iD + endpointFuzz] = w;
  }

  return GridIntegrator(nw);
}

/*
GridIntegrator 
GridIntegrator::SimpsonWithSplit(const int n,
                                 const double a,
                                 const double b,
                                 const set<double>& discSet,
                                 const double endpointFuzz)
TRY
{
/*  if (b <= a)
    throw Argument("Upper limit less than lower limit.");
  if (n < 2)
    throw Argument("Less than two points requested for Simpson quadrature."); 
  // Find first discontinuity >= a
  set<double>::const_iterator iD_begin = discSet.lower_bound(a);
  // Find first discontinuity > b
  set<double>::const_iterator iD_end = discSet.upper_bound(b);
  // Construct set of discontinuities in [a,b]
  set<double> effDiscSet(iD_begin,iD_end);
  // Insert a and b; then the set defines a disjoint cover of [a,b] by discontinuity-free intervals
  effDiscSet.insert(a);
  effDiscSet.insert(b);

  set<double> nodes;
  vector<double> weights;
  weights.reserve(n + effDiscSet.size());

  // Loop over subintervals to generate nodes and weights
  set<double>::const_iterator iD = effDiscSet.begin();
  set<double>::const_iterator iD_ = iD++;
  for (; iD != effDiscSet.end(); ++iD, ++iD_)
  {
    const int nn = static_cast<int>(ceil(n * (*iD - *iD_) / (b - a))); 
    if (*iD - *iD_ <= 2.0 * endpointFuzz)
      // Ignore "zero-length" intervals
      continue;
    else if (nn < 2)
    {
      // Always have at least two points
      nodes.insert(*iD_ + endpointFuzz);
      nodes.insert(*iD - endpointFuzz);
      weights.insert(weights.end(), 2, 0.5 * (*iD - *iD_));
    }
    else
    {
      // Normal case
      GridIntegrator singleInterval = GridIntegrator::Simpson(nn, *iD_ + endpointFuzz, *iD - endpointFuzz);
      nodes.insert(singleInterval.m_nodes.begin(), singleInterval.m_nodes.end());
      weights.insert(weights.end(), singleInterval.m_weights.begin(), singleInterval.m_weights.end());
    }
  }

  return GridIntegrator(nodes, weights);
}
CATCH(GridIntegrator::SimpsonWithSplit)
*/
  
  
DRLIB_END_NAMESPACE
  
  
  
  
  
  
  
