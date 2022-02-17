/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/xml/parametrization.h
// Purpose:     Names of elements and attributes in XML for parametrizations
// Author:      ITO33
// Created:     2004/11/25
// RCS-ID:      $Id: parametrization.h,v 1.11 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/xml/parametrization.h
    @brief Contains the names of the elements used in the XML description of a
           parametrization. 
 */

#ifndef _IHG_XML_PARAMETRIZATION_H_
#define _IHG_XML_PARAMETRIZATION_H_

#include "ito33/factory.h"
#include "ito33/string.h"
#include "ito33/sharedptr.h"

#define XML_TAG_PARAMETRIZATION_HRWITHTIMECOMPONENT_ROOT             "parametrization_hrwithtimecomponent"
#define XML_TAG_PARAMETRIZATION_HRWITHSPOTCOMPONENTPOWER_ROOT        "parametrization_hrwithspotcomponentpower"

#define XML_TAG_PARAMETRIZATION_VOLFLATHRWITHTIMECOMPONENT_ROOT      "parametrization_volflathrwithtimecomponent"
#define XML_TAG_PARAMETRIZATION_VOLFLATHRWITHSPOTCOMPONENTPOWER_ROOT "parametrization_volflathrwithspotcomponentpower"
#define XML_TAG_PARAMETRIZATION_VOLFLATHRPOWER_ROOT                  "parametrization_volflathrpower"

#define XML_TAG_PARAMETRIZATION_VOLPOWER_ROOT                        "parametrization_volpower"
#define XML_TAG_PARAMETRIZATION_VOLPOWERHRPOWER_ROOT                 "parametrization_volpowerhrpower"
#define XML_TAG_PARAMETRIZATION_VOLPOWERHRWITHTIMECOMPONENT_ROOT     "parametrization_volpowerhrwithtimecomponent"

#define XML_TAG_PARAMETRIZATION_VOLTANH_ROOT                         "parametrization_voltanh"
#define XML_TAG_PARAMETRIZATION_VOLTANHHRWITHSPOTCOMPONENTPOWER_ROOT "parametrization_voltanh_hrwithspotcomponentpower"
#define XML_TAG_PARAMETRIZATION_VOLTANHHRWITHTIMECOMPONENT_ROOT      "parametrization_voltanh_hrwithtimecomponent"
#define XML_TAG_PARAMETRIZATION_VOLTANHHRPOWER_ROOT                  "parametrization_voltanhhrpower"

#define XML_TAG_PARAMETRIZATION_VOLWITHTIMECOMPONENT_ROOT            "parametrization_volwithtimecomponent"

#define XML_TAG_PARAMETRIZATION_VOLTIMEONLYHRWITHTIMECOMPONENT_ROOT  "parametrization_voltimeonlyhrwithtimecomponent"



namespace xml { class node; }

namespace ito33
{

namespace ihg
{

class Parametrization;
class ParametrizationHRWithTimeComponent;
class ParametrizationHRWithSpotComponentPower;
class ParametrizationVolFlatHRWithTimeComponent;
class ParametrizationVolFlatHRWithSpotComponentPower;
class ParametrizationVolFlatHRPower;
class ParametrizationVolPower;
class ParametrizationVolPowerHRPower;
class ParametrizationVolPowerHRWithTimeComponent;
class ParametrizationVolTanh;
class ParametrizationVolTanhHRWithSpotComponentPower;
class ParametrizationVolTanhHRWithTimeComponent;
class ParametrizationVolTanhHRPower;
class ParametrizationVolWithTimeComponent;

/**
    ParametrizationFactory allows creation of Parametrizations by name.

    This factory allows creating objects of type ihg::Parametrization using
    string and creator functions taking xml::node as parameter.

    Because the creator in Factory specifies that the third parameter is passed
    as pointer(so can be optional), this requires xml::node be passed as 
    pointer to the specific ReadParametrization function(weired). 
 */
typedef Factory<std::string, Parametrization, xml::node>
          ParametrizationFactory;

/**
    Macro to be used to define a factory for a specific parametrization
    class.

    This macro supposes that a function ReadParametrization() taking an
    xml::node and returning a new Parametrization is defined.

    @param name of the XML tag from which Parametrization is read
    @param parametrization the concrete type to be created for this key
 */
#define ITO33_DEFINE_PARAMETRIZATION_READER(name, parametrization)            \
    /* use unnamed namespace to hide the declarations inside it */            \
    namespace                                                                 \
    {                                                                         \
        static ParametrizationFactory                                         \
        FactoryFor ## parametrization (name, Read ## parametrization);        \
    } struct Dummy /* just to force semicolon after the macro */


namespace XML
{


/**
    Read a ParametrizationHRWithTimeComponent from given node
    if the node has right tag.
 */
bool Restore(const xml::node &node,
             shared_ptr<ParametrizationHRWithTimeComponent>& pParam);

/**
    Read a ParametrizationHRWithSpotComponentPower from given node
    if the node has right tag.
 */
bool Restore(const xml::node &node,
             shared_ptr<ParametrizationHRWithSpotComponentPower>& pParam);

/**
    Read a ParametrizationVolFlatHRWithTimeComponent from given node
    if the node has right tag.
 */
bool Restore(const xml::node &node,
             shared_ptr<ParametrizationVolFlatHRWithTimeComponent>& pParam);

/**
    Read a ParametrizationVolFlatHRWithSpotComponentPower from given node
    if the node has right tag.
 */
bool Restore(const xml::node &node,
             shared_ptr<ParametrizationVolFlatHRWithSpotComponentPower>& pParam);

/**
    Read a ParametrizationVolFlatHRPower from given node if the node has 
    right tag.
 */
bool Restore(const xml::node &node,
             shared_ptr<ParametrizationVolFlatHRPower>& pParam);

/**
    Read a ParametrizationVolPower from given node
    if the node has right tag.
 */
bool Restore(const xml::node &node,
             shared_ptr<ParametrizationVolPower>& pParam);

/**
    Read a ParametrizationVolPowerHRPower from given node
    if the node has right tag.
 */
bool Restore(const xml::node &node,
             shared_ptr<ParametrizationVolPowerHRPower>& pParam);

/**
    Read a ParametrizationVolPowerHRWithTimeComponent from given node
    if the node has right tag.
 */
bool Restore(const xml::node &node,
             shared_ptr<ParametrizationVolPowerHRWithTimeComponent>& pParam);

/**
    Read a ParametrizationVolTanh from given node
    if the node has right tag.
 */
bool Restore(const xml::node &node,
             shared_ptr<ParametrizationVolTanh>& pParam);

/**
    Read a ParametrizationVolTanhHRWithSpotComponentPower from given node
    if the node has right tag.
 */
bool Restore(const xml::node &node,
             shared_ptr<ParametrizationVolTanhHRWithSpotComponentPower>& pParam);

/**
    Read a ParametrizationVolTanhHRWithTimeComponent from given node
    if the node has right tag.
 */
bool Restore(const xml::node &node,
             shared_ptr<ParametrizationVolTanhHRWithTimeComponent>& pParam);

/**
    Read a ParametrizationVolTanhHRPower from given node
    if the node has right tag.
 */
bool Restore(const xml::node &node,
             shared_ptr<ParametrizationVolTanhHRPower>& pParam);

/**
    Read a ParametrizationVolwithTimeComponent from given node
    if the node has right tag.
 */
bool Restore(const xml::node &node,
             shared_ptr<ParametrizationVolWithTimeComponent>& pParam);

} // namespace xml

} // namespace ihg

} // namespace ito33

#endif // _IHG_XML_PARAMETRIZATION_H_

