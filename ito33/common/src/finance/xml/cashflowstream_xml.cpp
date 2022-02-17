/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/cashflowstream_xml.cpp
// Purpose:     Restore cashflowstream object from XML document
// Author:      Yann d'halluin
// Created:     2004-06-25
// RCS-ID:      $Id: cashflowstream_xml.cpp,v 1.18 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"
#include "ito33/date.h"
#include "ito33/useexception.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/frequency.h"
#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cashflowstream_general.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/read_frequency.h"
#include "ito33/xml/finance/read_daycountconvention.h"
#include "ito33/xml/finance/cashflowstream_all.h"
#include "ito33/xml/finance/common.h"

namespace ito33
{
namespace XML
{
 
/**
    Exception thrown when the expected node is not found in XML document.
 */
class MissingCashFlowStream : public Exception
{
public:
  /**
      Constructs the exception thrown when an expected XML node is not found.

      This exception is thrown when cash flow stream data is missed.

      @param node the parent node
      @param name the name of child note we didn't find
      @param filename the name of the file where the exception occured
                      (usually just __FILE__)
      @param line the line number where it occured (usually __LINE__)
      @param function the name of the function where the exception occured
                      (__FUNCTION__ is unfortunately not yet supported by all
                      compilers so this argument is left empty for them)
   */
  MissingCashFlowStream(const char *filename,
                       size_t line,
                       const char *function)
    : Exception(ITO33_BAD_DATA,
                "Cash flow stream data is missed.",
                filename,
                line,
                function)
  {
  }
};

#define THROW_MISSING_CASH_FLOW_STREAM_EXCEPTION                      \
  throw ito33::XML::MissingCashFlowStream(__FILE__, __LINE__, __FUNCTION__)

shared_ptr<finance::CashFlowStreamUniform>
GetCashFlowStreamUniformFromNode(const xml::node& node);

shared_ptr<finance::CashFlowStreamGeneral>
GetCashFlowStreamGeneralFromNode(const xml::node& node);

/*
  Since CashFlowStream has only two sub classes. we don't need to use visitor
  and just do if/else
  */
shared_ptr<finance::CashFlowStream> 
GetCashFlowStreamInNode(const xml::node& node)
{
  shared_ptr<finance::CashFlowStream> pCashFlowStream;

  xml::node::const_iterator i;
  if ( ( i = node.find(XML_TAG_CASHFLOWSTREAMUNIFORM_ROOT ) ) != node.end() )
    pCashFlowStream = GetCashFlowStreamUniformFromNode(*i);
  else if( (i = node.find(XML_TAG_CASHFLOWSTREAMGENERAL_ROOT)) != node.end() )
    pCashFlowStream = GetCashFlowStreamGeneralFromNode(*i);
  else
    THROW_MISSING_CASH_FLOW_STREAM_EXCEPTION;

  return pCashFlowStream;
}


shared_ptr<finance::CashFlowStreamUniform> 
GetCashFlowStreamUniformInNode(const xml::node& node)
{
  shared_ptr<finance::CashFlowStreamUniform> pCashFlowStream;

  xml::node::const_iterator i;
  if ( ( i = node.find(XML_TAG_CASHFLOWSTREAMUNIFORM_ROOT ) ) != node.end() )
    pCashFlowStream = GetCashFlowStreamUniformFromNode(*i);
  else
  {
    typedef MissingNodeException Exception;
    throw EXCEPTION_MSG(node, XML_TAG_CASHFLOWSTREAMUNIFORM_ROOT);
  }

  return pCashFlowStream;
}


shared_ptr<finance::CashFlowStreamGeneral> 
GetCashFlowStreamGeneralInNode(const xml::node& node)
{
  shared_ptr<finance::CashFlowStreamGeneral> pCashFlowStream;

  xml::node::const_iterator i;
  if ( ( i = node.find(XML_TAG_CASHFLOWSTREAMGENERAL_ROOT ) ) != node.end() )
    pCashFlowStream = GetCashFlowStreamGeneralFromNode(*i);
  else
  {
    typedef MissingNodeException Exception;
    throw EXCEPTION_MSG(node, XML_TAG_CASHFLOWSTREAMGENERAL_ROOT);
  }

  return pCashFlowStream;
}

/////
shared_ptr<finance::CashFlowStreamUniform>
GetCashFlowStreamUniformFromNode(const xml::node& node)
{
  finance::LastPaymentType lastPaymentType = finance::LastPaymentType_Max;

  Date contractingDate;
  Date first;
  Date preLast;
  Date last;
  double annualAmt;
  Date::DayCountConvention dcm;
  finance::Frequency freq;
  
  contractingDate 
        = GetDateFromName(node, XML_TAG_CASHFLOWSTREAM_CONTRACTINGDATE);
    
  freq = GetFrequencyFromName(node, XML_TAG_PAYMENTFREQUENCY);

  dcm = GetDayCountConventionFromName(node, XML_TAG_DAYCOUNTCONVENTION); 
  
  annualAmt = GetDoubleFromName(node,
                                XML_TAG_CASHFLOWSTREAMUNIFORM_ANNUALAMOUNT); 
  xml::node::const_iterator iter;

  // check if the cash flow is adjusted by presence of cash flows
  std::vector<Date> pDates;
  std::vector<double> pdRates;

  for(iter = node.begin(); iter != node.end(); iter++)
  {
    if( strcmp(iter->get_name(), XML_TAG_CASHFLOWSTREAMGENERAL_CASHFLOW) == 0)
    {
      pDates.push_back(
        GetDateFromName(*iter, XML_TAG_FINANCE_DATE));
      pdRates.push_back(
        GetDoubleFromName(*iter, XML_TAG_FINANCE_RATE));
    }
  }

  // Adjusted
  if ( !pDates.empty() )
    return shared_ptr<finance::CashFlowStreamUniform>(
             new finance::CashFlowStreamUniform
               (contractingDate, pDates, pdRates, annualAmt, dcm, freq));
      
  iter = node.find(XML_TAG_CASHFLOWSTREAMUNIFORM_DURATION);

  // Defined by giving a duration in months
  if ( iter != node.end() )
  {
    size_t nMonths = GetLongFromNode(*iter);

    return shared_ptr<finance::CashFlowStreamUniform>(
             new finance::CashFlowStreamUniform
               (contractingDate, nMonths, freq, annualAmt, dcm));
  }

  first = GetDateFromName(node,
                          XML_TAG_CASHFLOWSTREAMUNIFORM_FIRSTDATE);
  
  xml::node::const_iterator iterSub;
  if ( (iterSub = node.find(XML_TAG_LASTPAYMENTTYPE) )
        != node.end() )
  {
    lastPaymentType = GetEnumFromNode
                      (
                        *iterSub,
                        SIZEOF(g_LastPaymentType),
                        g_LastPaymentType
                      );
  }

  last = GetDateFromName(node, XML_TAG_CASHFLOWSTREAMUNIFORM_LASTDATE);
  
  if (lastPaymentType != finance::LastPaymentType_Max )
    return shared_ptr<finance::CashFlowStreamUniform>(
             new finance::CashFlowStreamUniform(
                          contractingDate,first,last,
                          annualAmt,dcm,freq,lastPaymentType ));
  else
    return shared_ptr<finance::CashFlowStreamUniform>(
             new finance::CashFlowStreamUniform(
                                contractingDate,first,last,
                                annualAmt,dcm,freq ));
}

shared_ptr<finance::CashFlowStreamGeneral>
GetCashFlowStreamGeneralFromNode(const xml::node& node)
{
  Date contractingDate 
        = GetDateFromName(node, XML_TAG_CASHFLOWSTREAM_CONTRACTINGDATE);

  xml::node::const_iterator iter;

  std::vector<Date> pDates;
  std::vector<double> pdRates;

  for(iter = node.begin(); iter != node.end(); iter++)
  {
    if( strcmp(iter->get_name(), XML_TAG_CASHFLOWSTREAMGENERAL_CASHFLOW) == 0)
    {
      pDates.push_back(
        GetDateFromName(*iter, XML_TAG_FINANCE_DATE));
      pdRates.push_back(
        GetDoubleFromName(*iter, XML_TAG_FINANCE_RATE));
    }
  }
  
  finance::Frequency 
    freq = GetFrequencyFromName(node, XML_TAG_PAYMENTFREQUENCY);

  Date::DayCountConvention 
    dcm = GetDayCountConventionFromName(node, XML_TAG_DAYCOUNTCONVENTION);

  return shared_ptr<finance::CashFlowStreamGeneral>(
          new finance::CashFlowStreamGeneral(contractingDate, pDates, pdRates, 
    dcm, freq));
}


} // namespace XML

} // namespace ito33

