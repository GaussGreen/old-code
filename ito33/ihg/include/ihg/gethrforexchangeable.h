/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/gethrforexchangeable.h
// Purpose:     Get the hr of the issuer of an exchangeable
// Created:     2006/05/24
// RCS-ID:      $Id: gethrforexchangeable.h,v 1.3 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _IHG_GETHRFOREXCHANGEABLE_H_
#define _IHG_GETHRFOREXCHANGEABLE_H_

#include "ito33/finance/issuer.h"
#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/convertiblelike.h"

#include "ito33/ihg/hazardratetimeonly.h"

extern const ito33::finance::BondError
    ITO33_EXCHANGEABLE_NO_ISSUER_DEFAULT;

namespace ito33
{

namespace ihg
{

inline shared_ptr<HazardRateTimeOnly>
GetHRForExchangeable(const finance::ConvertibleLike& cb)
{
  ASSERT( cb.IsExchangeable() );

  shared_ptr<finance::Issuer> pIssuer( cb.GetIssuer() );

  CHECK_COND( pIssuer && pIssuer->HasDefaultIntensity(),
              ITO33_EXCHANGEABLE_NO_ISSUER_DEFAULT );

  shared_ptr<HazardRateTimeOnly>
    hr( new HazardRateTimeOnly(pIssuer->GetDatesOfDefaultIntensity(),
                               pIssuer->GetValuesOfDefaultIntensity()) );

  return hr;
}

} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_GETHRFOREXCHANGEABLE_H_
