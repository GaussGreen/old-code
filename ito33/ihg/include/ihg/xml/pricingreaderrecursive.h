/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/xml/pricingreaderrecursive.h
// Purpose:     reading pricing recursively XML files
// Author:      ITO33 Canada
// Created:     April 4, 2005
// RCS-ID:      $Id: pricingreaderrecursive.h,v 1.6 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/xml/pricingreaderrecursive.h
    @brief Classes for reading pricing XML-files generated in testsuite/bondlike/xmlfiles
 */

#ifndef _ITO33_IHG_XML_PRICINGREADERRECURSIVE_H_
#define _ITO33_IHG_XML_PRICINGREADERRECURSIVE_H_

#include "ito33/sharedptr.h"

#include "ito33/xml/reader.h"

namespace ito33
{

namespace ihg
{

  class TheoreticalModel;

namespace XML
{

/**
    This class reads pricing XML documents that were
    created for ihg/test/testsuite/bondlike/xmlfiles and provides 
    easy access to the data in them.

    Right now XML data can only be loaded from a file but we could also load it
    from a string if needed. Unfortunately, xmlwrapp doesn't take std::istream
    as input...
 */
class PricingReaderRecursive : public ito33::XML::Reader
{
public:
  /**
    Load XML from the given file.

    If an error occurs while loading date, i.e. file is not found, couldn't be
    read or doesn't contain a valid XML document with the IHG roto tag, an
    exception is thrown.

    @param filename the name of the file to read XML from
   */
  PricingReaderRecursive(const char *filename)
                       : ito33::XML::Reader(filename)
  {
  }

  PricingReaderRecursive(const char *data, size_t len)
                       : ito33::XML::Reader(data, len)
  {
  }

  // default dtor OK
   ~PricingReaderRecursive() {}

  /**
     Fill in the passed in SessionData object with the data from XML.

     If any data needed for SessionData initialization is missing, an 
     exception is thrown.

     @sa GetSessionDataFromNode() in include/ito33/xml/finance/sessiondata.h.
 
     @return ptr to sessiondata to fill in with information we read from node
   */
  shared_ptr<finance::SessionData> ReadSessionData() const; // throws

  /**
     Be notified about all derivatives found in the XML document.

     For each subsection of the "<derivatives>" tag, the corresponding method
     of visitor object will be called.

     This method may throw if a DerivativeVisitor method throws.
   */
  void ReadDerivatives(finance::DerivativeVisitor& visitor) const;
  
  /**
     Read the first derivative found in the XML document. We don't care about 
     the type of the derivative.
   */
  void ReadDerivative(shared_ptr<finance::Derivative>& pDerivative) const;

  /**
     Fill in the passed in TheoreticalModel object with model parameters.

     @param pModel the object to fill in
   */
  void ReadTheoreticalModel(shared_ptr<TheoreticalModel>& pModel) const;


private:

  NO_COPY_CLASS(PricingReaderRecursive);
};

} // namespace XML

} // namespace ihg

} // namespace ito33

#endif // _ITO33_IHG_XML_PRICINGREADERRECURSIVE_H_
