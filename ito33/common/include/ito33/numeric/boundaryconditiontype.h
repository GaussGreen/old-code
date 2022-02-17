/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/boundaryconditiontype.h
// Purpose:     enum of boundary condition type
// Author:      WANG Xuewen
// Created:     2003/07/25
// RCS-ID:      $Id: boundaryconditiontype.h,v 1.3 2004/10/05 09:13:38 pedro Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////


#ifndef _ITO33_NUMERIC_BOUNDARYCONDITIONTYPE_H_
#define _ITO33_NUMERIC_BOUNDARYCONDITIONTYPE_H_

namespace ito33{

namespace numeric{

enum BoundaryConditionType
{
  BCType_Dirichlet,
  BCType_Gamma,
  BCType_Neumann,

  BCType_Max
};

} // namespace numeric

} // namespace ito33

#endif // _ITO33_NUMERIC_BOUNDARYCONDITIONTYPE_H_

