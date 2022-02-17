/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/xml/hazardrate.h
// Purpose:     Names of elements and attributes used in XML for HazardRate
// Author:      Yann d'Halluin
// Created:     2004-05-19
// RCS-ID:      $Id: hazardrate.h,v 1.11 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/xml/hazardrate.h
    @brief Contains the names of the elements used in the XML description of a
           hazardrate. 
 */

#ifndef _IHG_XML_HAZARDRATE_H_
#define _IHG_XML_HAZARDRATE_H_

#include "ito33/sharedptr.h"
#include "ito33/factory.h"
#include "ito33/string.h"
#include "ito33/enum_values_names.h"

#define XML_TAG_HAZARDRATE_ROOT              "hazard_rate"
#define XML_TAG_HAZARDRATEFLAT_ROOT          "hazard_rate_flat"
#define XML_TAG_HAZARDRATETIMEONLY_ROOT      "hazard_rate_time_only"
#define XML_TAG_HAZARDRATECOMBO_ROOT         "hazard_rate_combo"
#define XML_TAG_HAZARDRATEDECAY_ROOT         "hazard_rate_decay"
#define XML_TAG_HAZARDRATELINEAR_ROOT        "hazard_rate_linear"
#define XML_TAG_HAZARDRATEPOWER_ROOT         "hazard_rate_power"
#define XML_TAG_HAZARDRATETIMECOMPONENT_ROOT "hazard_rate_time_component"
#define XML_TAG_HAZARDRATETIMEONLY_CDSCURVE_ROOT  "hazard_rate_time_only_by_cds_curve"

#define XML_TAG_HAZARDRATE_FLAT                 "flat"
#define XML_TAG_HAZARDRATE_TIMEONLY             "timeonly"
#define XML_TAG_HAZARDRATE_DECAY                "decay"
#define XML_TAG_HAZARDRATE_LINEAR               "linear"
#define XML_TAG_HAZARDRATE_POWER                "power"
#define XML_TAG_HAZARDRATE_ALPHA                "alpha"
#define XML_TAG_HAZARDRATE_S0                   "S0"
#define XML_TAG_HAZARDRATE_SLOPE                "slope"
#define XML_TAG_HAZARDRATE_B                    "b"
#define XML_TAG_HAZARDRATE_BETA                 "beta"
#define XML_TAG_HAZARDRATE_SPOT                 "spot"


#define XML_TAG_HAZARDRATE_TIMECOMPONENT        "time_component"

namespace xml 
{ 
  class node; 
}

namespace ito33
{

namespace ihg
{
  class HazardRate;
}

/**
    HazardRateFactory allows creation of hazard rate by name.

    This factory allows creating objects of type ihg::HazardRate using
    string and creator functions taking xml::node as parameter.
 */
typedef Factory<std::string, ihg::HazardRate, xml::node> HazardRateFactory;

/**
    Macro to be used to define a factory for a specific hazardrate-derived
    class.

    This macro supposes that a function ReadHazardRate() taking an xml::node
    and returning a new HazardRate is defined.

    @param name of the XML tag from which hazardrate is read
    @param hazardrate the concrete type to be created for this key
 */
#define ITO33_DEFINE_HAZARDRATE_READER(name, hazardrate)                      \
    /* use unnamed namespace to hide the declarations inside it */            \
    namespace                                                                 \
    {                                                                         \
        static HazardRateFactory                                              \
        FactoryFor ## hazardrate (name, Read ## hazardrate);                  \
    } struct Dummy /* just to force semicolon after the macro */


} // namespace ito33


namespace xml { class node; }

namespace ito33
{

namespace ihg
{

class HazardRate;


namespace XML
{

/**
    Restore a hazard rate object from XML. Otherwise, a NULL pointer
    is returned.

    @param node the root tag in DOM tree where we should find hazard rate
    @return the new hazard rate object
 */
 shared_ptr<HazardRate> ReadHazardRate(const xml::node &node);

} // namespace XML

} // namespace ihg

} // namespace ito33

#endif // _ITO33_XML_HAZARDRATE_H_
