/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/fixedspacemeshmanager.h
// Purpose:     fixed space mesh manager for PDE problems
// Author:      Wang
// Created:     2004/02/11
// RCS-ID:      $Id: fixedspacemeshmanager.h,v 1.5 2004/10/05 09:13:39 pedro Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_PRICING_FIXEDSPACEMESHMANAGER_H_
#define _ITO33_PRICING_FIXEDSPACEMESHMANAGER_H_

#include "ito33/array.h"

namespace ito33
{

namespace pricing
{


/// Base class for fixed space management
class FixedSpaceMeshManager
{
public:
  
  /// Dummy virtual ctor
  FixedSpaceMeshManager() { }

  /// Dummy virtual dtor 
  virtual ~FixedSpaceMeshManager() { }

  size_t GetNbS() const { return m_nNbS; }
  
  /// Get a pointer to the fixed space mesh
  const double *GetS() const 
  {
    return m_pdS.Get();
  }

  /// Get a pointer to the fixed space mesh
  const double *GetLogS() const 
  {
    return m_pdLogS.Get();
  }


protected:
 
  /// Construct the space mesh
  virtual void ConstructSpaceMesh() = 0;

  size_t m_nNbS;

  Array<double> m_pdS;

  Array<double> m_pdLogS;


private:

  NO_COPY_CLASS(FixedSpaceMeshManager);

}; // class FixedSpaceMeshManager


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_FIXEDSPACEMESHMANAGER_H_
