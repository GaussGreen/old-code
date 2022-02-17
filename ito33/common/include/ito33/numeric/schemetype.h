/******************************************************************************
 * File.........: ito33/schemetype.h
 * Purpose......: enum of the time discretization scheme type 
 * Author.......: WANG Xuewen 
 * Created......: 3/19/2003
 * RCS-ID.......: $Id: schemetype.h,v 1.4 2004/10/05 09:13:38 pedro Exp $
 * Copyright....: (c) 2003 Trilemma LLP
 ******************************************************************************/

/**
  \file ito33/numeric/schemetype.h
  \brief enum of the time discretization scheme type 
*/

#ifndef _ITO33_NUMERIC_SCHEMETYPE_H_
#define _ITO33_NUMERIC_SCHEMETYPE_H_

namespace ito33
{

namespace numeric
{

enum SchemeType
{
  SchemeType_Implicit = 0,
  SchemeType_CrankNicolson = 1,
  SchemeType_ThreeLevel = 2,
  SchemeType_Explicit = 3,
 
  SchemeType_Max
};

} // namespace numeric

} // namespace ito33


#endif // _ITO33_NUMERIC_SCHEMETYPE_H_
