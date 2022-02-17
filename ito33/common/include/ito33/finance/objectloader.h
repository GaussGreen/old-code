/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/objectloader.h
// Purpose:     declaration of ObjectLoader for loading from datastore
// Created:     2006/05/09
// RCS-ID:      $Id: objectloader.h,v 1.1 2006/08/24 10:03:59 zhang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/objectloader.h
    @brief declaration of ObjectLoader for loading finance objects
           from datastore
 */

#ifndef _ITO33_FINANCE_OBJECTLOADER_H_
#define _ITO33_FINANCE_OBJECTLOADER_H_

#include "ito33/common.h"

#include "ito33/sharedptr.h"
#include "ito33/string.h"
#include "ito33/date.h"

#include "ito33/finance/objectloaderdlldecl.h"

namespace ito33
{

namespace datastore
{
  class Reader;
}

namespace finance
{
  class ITO33_DLLDECL RateData;
  class ITO33_DLLDECL Equity;
  class ITO33_DLLDECL Option;
  class ITO33_DLLDECL CDS;
  class ITO33_DLLDECL EDS;
  class ITO33_DLLDECL ConvertibleBond;
  class ITO33_DLLDECL Reset;
  class ITO33_DLLDECL AttachedWarrantConvertibleBond;

/**
    ObjectLoader loads finance objects from datastore.
 */
class ITO33_OBJECTLOADER_DLLDECL ObjectLoader
{
public:

  /// type of the object id
  typedef long ObjectIdType;

  /**
      Constructor takes the data source name of the underlying datastore.

      May throw if it's the first time it's called and the initialization
      fails.

      @param dsn to use for datastore plugin initialization
   */
  ObjectLoader(const char* dsn = "");

  // Default dtor is ok.
  
  /**
      Valuation date.

      @param valuationDate valuation date
   */
  void SetValuationDate(const Date& valuationDate);

  /**
      Valuation date.

      @return Valuation date.
   */
  Date GetValuationDate() const;

  /**
      Loads all rate data from the datastore.

      @return Rate data pointer
   */
  shared_ptr<RateData> GetAllRateData() const;
  
  /**
      Loads an equity from the datastore.

      Throws if the id doesn't correspond to an equity.

      @return an equity pointer
   */
  shared_ptr<Equity> GetEquity(ObjectIdType id) const;

  /**
      Loads an option from the datastore.

      Throws if the id doesn't correspond to an option.

      @return an option pointer
   */
  shared_ptr<Option> GetOption(ObjectIdType id) const;

  /**
      Loads a CDS from the datastore.

      Throws if the id doesn't correspond to a CDS.

      @return a CDS pointer
   */
  shared_ptr<CDS> GetCDS(ObjectIdType id) const;

  /**
      Loads an EDS from the datastore.

      Throws if the id doesn't correspond to an EDS.

      @return an EDS pointer
   */
  shared_ptr<EDS> GetEDS(ObjectIdType id) const;

  /**
      Loads a convertible bond from the datastore.

      Throws if the id doesn't correspond to a ConvertibleBond.

      @return a convertible bond pointer
   */
  shared_ptr<ConvertibleBond> GetConvertibleBond(ObjectIdType id) const;

  /**
      Loads a reset convertible bond from the datastore.

      Throws if the id doesn't correspond to a Reset.

      @return a reset convertible bond pointer
   */
  shared_ptr<Reset> GetReset(ObjectIdType id) const;

  /**
      Loads a convertible bond with share dependent conversion from the 
      datastore.

      Throws if the id doesn't correspond to a AttachedWarrantConvertibleBond.

      @return a pointer to a convertible bond with share dependent conversion 
   */
  shared_ptr<AttachedWarrantConvertibleBond> 
  GetAttachedWarrantConvertibleBond(ObjectIdType id) const;

  /**
      Invalidates the datastore cache entirely, in case that there is a 
      datastore change.
   */
  void InvalidateCache();

protected:
  
  /// pointer to the datastore reader
  shared_ptr<datastore::Reader> m_pReader;

private:

  std::string m_sDSN;

  Date m_valuationDate;

  NO_COPY_CLASS(ObjectLoader);
};

} // namespace datastore

} // namespace ito33

#endif // _ITO33_FINANCE_OBJECTLOADER_H_
