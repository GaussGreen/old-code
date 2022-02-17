/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/domain_fixedspacemesh.h
// Purpose:     Domain class with same space mesh for all times 
// Author:      WANG Xuewen
// Created:     2004/05/03
// RCS-ID:      $Id: domain_fixedspacemesh.h,v 1.10 2005/09/14 14:39:32 nabil Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/domain_fixedspacemesh.h
    @brief Numeric Domain class defines a fixed space mesh. 
 */

#ifndef _ITO33_NUMERIC_DOMAIN_FIXEDSPACEMESH_H_
#define _ITO33_NUMERIC_DOMAIN_FIXEDSPACEMESH_H_

#include <math.h>

#include "ito33/constants.h"
#include "ito33/numeric/domain.h"

namespace ito33
{

namespace numeric
{

/**
   Domain class with same space mesh for all times. 
 */
class DomainFixedSpaceMesh : public Domain
{
 
public:

  /**
     Ctor constructs a trivial fixed space mesh from a single point.

     This can be used for problems that the resulted prices depend only
     on time.

     @param dS the single point at the space mesh 

     @todo Should we construct trivial space mesh with two points so that
           the surface that the values depend only on time can still be 
           outputed/plotted as a 'real' surface by the usual tools?
   */
  DomainFixedSpaceMesh(double dS) : m_pdS(1)
  {
    m_pdS[0] = dS;
  }

  /**
     Ctor constructs the fixed space mesh from a double pointer and 
     uses the fixed space mesh for default end user space mesh.

     @param pdS the space points
     @param nNbS the size of the space mesh
   */
  DomainFixedSpaceMesh(const double *pdS, size_t nNbS)
    : m_pdS(pdS, pdS + nNbS)
  {
  }

  /**
     Get the space mesh
     
     @return a reference to the internal space mesh container 
   */
  const Spots& GetSpaceMesh() const
  {
    return m_pdS;
  }

  /**
     Gets the first space mesh at the given time index

     @param nIdxT the index of the time at which the space mesh is required

     @return the space mesh at current time
   */
  const Spots& GetFirstSpaceMeshAt(size_t nIdxT) const;
 
  /**
     Gets the last space mesh at the given time index

     @param nIdxT the index of the time at which the space mesh is required

     @return the space mesh at current time
   */
  const Spots& GetLastSpaceMeshAt(size_t nIdxT) const;

  /**
     Gets the output space mesh at the given time index

     @param nIdxT the index of the time at which the space mesh is required

     @return the space mesh at current time
   */
  const Spots& GetOutputSpaceMeshAt(size_t nIdxT) const;

  /**
     Add a time point to the numerical surface

     @param dTime A time point to be added to the surface
   */
  void AddTime(double dTime)
  {
    numeric::Domain::AddTime(dTime);
  }

private:

  /// the space mesh points
  Spots m_pdS;
   
}; // class DomainFixedSpaceMesh

inline const Domain::Spots& 
DomainFixedSpaceMesh::GetFirstSpaceMeshAt(size_t /* nIdxT */) const
{
  return m_pdS;
}

inline const Domain::Spots& 
DomainFixedSpaceMesh::GetLastSpaceMeshAt(size_t /* nIdxT */) const
{
  return m_pdS;
}

inline const Domain::Spots& 
DomainFixedSpaceMesh::GetOutputSpaceMeshAt(size_t /* nIdxT */) const
{
  return m_pdS;
}


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_DOMAIN_FIXEDSPACEMESH_H_
