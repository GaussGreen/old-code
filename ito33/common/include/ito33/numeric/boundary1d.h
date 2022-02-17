/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/boundary1d.h
// Purpose:     a class for the boundary condition of a 1D PDE
// Author:      WANG Xuewen
// Created:     2003/08/22
// RCS-ID:      $Id: boundary1d.h,v 1.3 2004/10/05 09:13:38 pedro Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////


#ifndef _ITO33_NUMERIC_BOUNDARY1D_H_
#define _ITO33_NUMERIC_BOUNDARY1D_H_

#include "ito33/debug.h"
#include "ito33/numeric/boundaryconditiontype.h"

namespace ito33{

namespace numeric
{

/**
  @brief Class representing boundary conditions for 1D PDE problems

  Provides functions for getting and setting the left and right
  boundary conditions for one dimensional PDE problems. Current
  boundary condition types include Dirchlet and 
  linear (gamma = U_{SS} = 0).

*/
class Boundary1D
{
public:

  Boundary1D() : m_BCTypeLeft(BCType_Gamma),
                 m_BCTypeRight(BCType_Gamma),
                 m_dValueLeft(0), 
                 m_dValueRight(0)
  {
  }

  void SetLeft(BoundaryConditionType leftBC, double dValue = 0.0)
  {
    m_BCTypeLeft = leftBC;
    m_dValueLeft = dValue;
  }
  
  void SetRight(BoundaryConditionType rightBC, double dValue = 0.0)
  {
    m_BCTypeRight = rightBC;
    m_dValueRight = dValue;
  }

  BoundaryConditionType GetLeftType() const
  {
    return m_BCTypeLeft;
  }

  double GetLeftValue() const
  {
    ASSERT_MSG(m_BCTypeLeft == BCType_Dirichlet,
               "Can't get the value for a non Dirichlet boundary condition!"); 

    return m_dValueLeft;
  }

  BoundaryConditionType GetRightType() const
  {
    return m_BCTypeRight;
  }

  double GetRightValue() const
  {
    ASSERT_MSG(m_BCTypeRight == BCType_Dirichlet,
               "Can't get the value for a non Dirichlet boundary condition!"); 

    return m_dValueRight;
  }

private:
  BoundaryConditionType
    m_BCTypeLeft,
    m_BCTypeRight;
  double
    m_dValueLeft, 
    m_dValueRight;
};

} // namespace numeric

} // namespace ito33

#endif // #ifdef _ITO33_NUMERIC_BOUNDARY1D_H_

