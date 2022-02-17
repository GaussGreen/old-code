/////////////////////////////////////////////////////////////////////////////
// Name:        src/finance/objectloader.cpp
// Purpose:     implementation of finance::ObjectLoader for loading financial
//              objects from datastore
// Created:     2006/05/09
// RCS-ID:      $Id: objectloader.cpp,v 1.1 2006/08/24 10:03:59 zhang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/useexception.h"

#include "ito33/finance/objectloader.h"

#include "ito33/finance/ratedata.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/moneymarket.h"
#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/reset.h"
#include "ito33/finance/bondlike/attachedwarrantconvertiblebond.h"

#include "ito33/datastore/reader.h"
#include "ito33/datastore/typeconv.h"
#include "ito33/datastore/yieldcurveloader.h"


extern const ito33::Error
  ITO33_BAD_DATE;

namespace ito33
{

using datastore::Reader;

namespace finance
{

// ============================================================================
// implementation
// ============================================================================
ObjectLoader::ObjectLoader(const char* dsn)
            : m_sDSN(dsn),
              m_pReader( Reader::Get(dsn).release() ),
              m_valuationDate(ito33::Date::Today())
{
}

void ObjectLoader::SetValuationDate(const Date& valuationDate)
{
  CHECK_COND ( valuationDate.IsValid(), ITO33_BAD_DATE );

  m_valuationDate = valuationDate;
}

Date ObjectLoader::GetValuationDate() const
{
  return m_valuationDate;
}

shared_ptr<RateData> ObjectLoader::GetAllRateData() const
{
  datastore::YieldCurveLoader ycLoader(*m_pReader, m_valuationDate);

  shared_ptr<RateData> pRateData(ycLoader.GetRateData());

  const RateData::YCElements& elements(pRateData->GetAll());
  for ( RateData::YCElements::const_iterator iter = elements.begin(),
        end = elements.end();
        iter != end;
        ++iter )
  {
    iter->second->GetYieldCurve()->SetReferenceDate(m_valuationDate);
  }

  return pRateData;
}

void ObjectLoader::InvalidateCache()
{
  return m_pReader->InvalidateCache();
}

// macro for implementing ObjectLoader::GetXXX(id) where XXX is a class name under
// finance namespace
#define IMPLEMENT_LOAD(Type)                                                  \
  shared_ptr< Type> ObjectLoader::Get ##Type(ObjectIdType id) const           \
{                                                                             \
  shared_ptr<datastore:: Type> pDS(m_pReader->Get ##Type(id));                \
                                                                              \
  ASSERT_MSG( pDS, "Reader::Get() should not return a null pointer" );        \
                                                                              \
  return pDS->GetCoreObject();                                                \
}

// macro for implementing ObjectLoader::GetXXX(id) where XXX is a general derivative
#define IMPLEMENT_LOAD_DERIVATIVE(Type) IMPLEMENT_LOAD(Type)

// macro for implementing ObjectLoader::GetXXX(id) where XXX is a CB type
// in particular XXX is derived from CBBase
#define IMPLEMENT_LOAD_CB(Type)                                               \
  shared_ptr< Type> ObjectLoader::Get ##Type(ObjectIdType id) const           \
{                                                                             \
  shared_ptr<datastore::ConvertibleBond>                                      \
    pDS(m_pReader->GetConvertibleBond(id));                                   \
                                                                              \
  ASSERT_MSG( pDS, "Reader::Get() should not return a null pointer" );        \
                                                                              \
  /* FIXME : throw if casted pointer is NULL */                               \
  return dynamic_pointer_cast<Type>(pDS->GetCoreObject());                    \
}

IMPLEMENT_LOAD(Equity)

IMPLEMENT_LOAD_DERIVATIVE(Option)
IMPLEMENT_LOAD_DERIVATIVE(CDS)
IMPLEMENT_LOAD_DERIVATIVE(EDS)

IMPLEMENT_LOAD_CB(ConvertibleBond)
IMPLEMENT_LOAD_CB(Reset)
IMPLEMENT_LOAD_CB(AttachedWarrantConvertibleBond)

} // namespace finance

} // namespace ito33
