/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/xml/spotcomponent.h
// Purpose:     Names of elements and attributes used in XML for SpotComponent
// Author:      David Pooley
// Created:     2004-06-18
// RCS-ID:      $Id: spotcomponent.h,v 1.6 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/xml/spotcomponent.h
    @brief Contains the names of the elements used in the XML description of a
           spotcomponent. 
 */

#ifndef _IHG_XML_SPOTCOMPONENT_H_
#define _IHG_XML_SPOTCOMPONENT_H_

#include "ito33/sharedptr.h"
#include "ito33/factory.h"

#define XML_TAG_SPOTCOMPONENT_ROOT      "spot_component"

#define XML_TAG_SPOTCOMPONENT_POWER_ROOT      "spot_component_power"
#define XML_TAG_SPOTCOMPONENT_BETA      "beta"
#define XML_TAG_SPOTCOMPONENT_S0        "S0"


namespace xml { class node; }

namespace ito33
{

namespace ihg
{
class SpotComponent;
}

/**
    SpotComponentFactory allows creation of hazard rate by name.

    This factory allows creating objects of type ihg::SpotComponent using
    string and creator functions taking xml::node as parameter.
 */
  typedef Factory<std::string, ihg::SpotComponent, xml::node>
                                                        SpotComponentFactory;

/**
    Macro to be used to define a factory for a specific spotcomponent-derived
    class.

    This macro supposes that a function ReadSpotComponent() taking an xml::node
    and returning a new SpotComponent is defined.

    @param name of the XML tag from which spotcomponent is read
    @param derivative the concrete type to be created for this key
 */
#define ITO33_DEFINE_SPOTCOMPONENT_READER(name, spotcomponent)                \
    /* use unnamed namespace to hide the declarations inside it */            \
    namespace                                                                 \
    {                                                                         \
        static SpotComponentFactory                                           \
        FactoryFor ## spotcomponent (name, Read ## spotcomponent);            \
    } struct Dummy /* just to force semicolon after the macro */


namespace ihg
{


namespace XML
{

/**
    Restores a spotcomponent object from XML. Otherwise, a NULL pointer
    is returned.

    @param node the SpotComponent tag in DOM tree
    @return the new spotcomponent object
 */
 shared_ptr<SpotComponent> ReadSpotComponent(const xml::node &node);

/**
    Gets a spotcomponent object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the root tag where we should find "spot_component" tag
    @return the spotcomponent object
 */
 shared_ptr<SpotComponent> GetSpotComponentFromNode(const xml::node &node);


} // namespace XML

} // namespace ihg

} // namespace ito33

#endif // _IHG_XML_SPOTCOMPONENT_H_

