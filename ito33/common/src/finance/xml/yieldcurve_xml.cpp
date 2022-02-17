/////////////////////////////////////////////////////////////////////////////
// Name:        yieldcurve_xml.cpp
// Purpose:     Restore YieldCurve object from XML document
// Author:      Vadim Zeitlin
// Created:     2004-05-04
// RCS-ID:      $Id: yieldcurve_xml.cpp,v 1.21 2006/08/23 10:39:04 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/finance/yieldcurve.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/yieldcurve_annuallycompounded.h"
#include "ito33/finance/yieldcurve_swap.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/yieldcurve.h"
#include "ito33/xml/finance/read_daycountconvention.h"
#include "ito33/xml/finance/read_timeunit.h"
#include "ito33/xml/finance/read_frequency.h"

// define the field of the Factory declared in ito33/xml/finance/yieldcurve.h
ITO33_IMPLEMENT_THE_FACTORY(ito33::YieldCurveFactory);

ITO33_FORCE_LINK_THIS_MODULE(yieldcurve_xml);

namespace ito33
{

namespace XML
{

void GetYieldCurveLegPart(const xml::node& node,
                     double& dRate,
                     size_t& nMaturityDuration,
                     finance::TimeUnit& maturityUnit
                     )
{  
  dRate = GetDoubleFromName(node, XML_TAG_SWAPCURVELEG_RATE);
  nMaturityDuration = GetLongFromName
                      (node, XML_TAG_SWAPCURVELEG_MATURITY_DURATION);
  maturityUnit = GetTimeUnitFromName(node, XML_TAG_SWAPCURVELEG_MATURITY_UNIT);
}

/// restores CashRate object from given node known as 'root' of CashRate
finance::CashRate GetCashRateFromNode(const xml::node& node)
{
  double dRate;
  size_t nMaturityDuration;
  finance::TimeUnit maturityUnit;
  GetYieldCurveLegPart(node, dRate, nMaturityDuration, maturityUnit);
  return finance::CashRate(dRate, nMaturityDuration, maturityUnit);
}

/// restores SwapRate object from given node known as 'root' of SwapRate
finance::SwapRate GetSwapRateFromNode(const xml::node& node)
{
  double dRate;
  size_t nMaturityDuration;
  finance::TimeUnit maturityUnit;
  finance::Frequency freq;

  GetYieldCurveLegPart(node, dRate, nMaturityDuration, maturityUnit);
  freq = GetFrequencyFromName(node, XML_TAG_PAYMENTFREQUENCY);

  return finance::SwapRate(dRate, nMaturityDuration, maturityUnit, freq);
}

/**
    Restores a flat YieldCurve object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the session tag in DOM tree
    @return the new option object to be deleted by the caller
 */
static
finance::YieldCurve *ReadYieldCurveFlat(const xml::node *pNode)
{
  const xml::node& node = *pNode;

  double dYieldRate = GetDoubleFromName(node, XML_TAG_FINANCE_FLAT);

  return new finance::YieldCurveFlat( dYieldRate );
}

ITO33_DEFINE_YIELDCURVE_READER(XML_TAG_YIELDCURVEFLAT_ROOT, YieldCurveFlat);


/**
    Restores an annually compounded YieldCurve object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the session tag in DOM tree
    @return the new option object to be deleted by the caller
 */
static
finance::YieldCurve *ReadYieldCurveAnnuallyCompounded(const xml::node *pNode)
{
 const xml::node& node = *pNode;
 Date referenceDate = GetDateFromName(node,XML_TAG_YIELDCURVE_REFERENCE_DATE);

 finance::YieldCurveAnnuallyCompounded *yc;
 yc =  new finance::YieldCurveAnnuallyCompounded(referenceDate);

  for (  xml::node::const_iterator j = node.begin(); j != node.end(); ++j )
  {
    if ( strcmp(j->get_name(), XML_TAG_YIELDCURVE_LEG) == 0 )
    {
      yc->AddLeg
         (
          static_cast<int>(GetLongFromNode
                           (GetNodeByName(*j, XML_TAG_YIELDCURVE_LEG_DAY))),
          GetDoubleFromNode(GetNodeByName(*j, XML_TAG_FINANCE_VALUE))
         );
    }
  }

  yc->Validate();

 return yc;
}

ITO33_DEFINE_YIELDCURVE_READER(XML_TAG_YIELDCURVEANNUALLYCOMPOUNDED_ROOT,
                               YieldCurveAnnuallyCompounded);

/**
    Restores an swap Curve object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the session tag in DOM tree
    @return the new option object to be deleted by the caller
 */
static
finance::YieldCurve *ReadYieldCurveSwap(const xml::node *pNode)
{
 const xml::node& node = *pNode;
 Date referenceDate = GetDateFromName(node,XML_TAG_YIELDCURVE_REFERENCE_DATE);

 finance::YieldCurveSwap *yc =  new finance::YieldCurveSwap(referenceDate);


 finance::CashRates cashRates;
 finance::SwapRates swapRates;

 cashRates.SetBasis( GetDayCountConventionFromName(node,
                              XML_TAG_YIELDCURVE_SWAP_CASHBASIS) );
 swapRates.SetBasis( GetDayCountConventionFromName(node,
                              XML_TAG_YIELDCURVE_SWAP_SWAPBASIS) );
 
  for (  xml::node::const_iterator j = node.begin(); j != node.end(); ++j )
  {
    if ( strcmp(j->get_name(), XML_TAG_CASHRATE_ROOT) == 0 )
    {
      finance::CashRate rate(GetCashRateFromNode(*j));
      cashRates.AddLeg(rate.GetRate(),
                      rate.GetMaturityDuration(),
                      rate.GetMaturityUnit());

    }
    else if ( strcmp(j->get_name(), XML_TAG_SWAPRATE_ROOT) == 0 )
    {
      finance::SwapRate rate(GetSwapRateFromNode(*j));
      swapRates.AddLeg(rate.GetRate(),
                      rate.GetMaturityDuration(),
                      rate.GetMaturityUnit(),
                      rate.GetPaymentFrequency());

    }
      //end if

  } //end for

  yc->SetSwapRates(swapRates);
  yc->SetCashRates(cashRates);
  yc->Validate();

 return yc;
}

ITO33_DEFINE_YIELDCURVE_READER(XML_TAG_YIELDCURVESWAP_ROOT,
                               YieldCurveSwap);

shared_ptr<finance::YieldCurve>
GetYieldCurveFromNode(const xml::node& node)
{
  shared_ptr<finance::YieldCurve> pYc;

  // we have to iterate here, as the first subnode is not real yield curve node
  xml::node::const_iterator i = node.begin();
  for(; i != node.end() && !pYc; i++)
      pYc = shared_ptr<finance::YieldCurve>
                    (YieldCurveFactory::Create(i->get_name(), &(*i)));

  if(!pYc)
  {
    // name is not exactly message but we can still (ab)use this macro here
    typedef TypeMismatchException Exception;
    throw EXCEPTION_MSG(ITO33_BAD_DATA,
                        String::Printf(TRANS("Nothing found in Node \"%s\"."),
                                       node.get_name())
                        );
  }

  return pYc;
}

} // namespace XML

} // namespace ito33
