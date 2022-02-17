/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/xml/spotcomponent_xml.cpp
// Purpose:     Restore spot component object from XML document
// Author:      Yann d'Halluin
// Created:     2004-05-18
// RCS-ID:      $Id: spotcomponent_xml.cpp,v 1.6 2006/08/20 09:31:05 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////


#include "ito33/useexception.h"
#include "ito33/vector.h"

#include "ito33/ihg/spotcomponent.h"
#include "ito33/ihg/hrspotcomponentpower.h"

#include "ihg/xml/spotcomponent.h"

#include "ito33/xml/read.h"

extern const ito33::Error ITO33_UNEXPECTED;

using namespace ito33;
using namespace ito33::XML;

const extern ito33::Error ITO33_BAD_DATA;

ITO33_IMPLEMENT_THE_FACTORY(ito33::SpotComponentFactory);

namespace ito33
{
namespace ihg
{
namespace XML
{


ito33::shared_ptr<ito33::ihg::SpotComponent>
ReadSpotComponent(const xml::node &node)
{
  shared_ptr<ihg::SpotComponent> spotcomponent;
  xml::node::const_iterator i;
  for ( i = node.begin(); i != node.end(); ++i ) 
  {
    spotcomponent = make_ptr( SpotComponentFactory::Create
                              (i->get_name(), &(*i)) );
    if ( spotcomponent )
      break;
  }
  return spotcomponent;
}


ito33::shared_ptr<ito33::ihg::SpotComponent>
GetSpotComponentFromNode(const xml::node &node)
{
  xml::node nodeRoot(GetNodeByName(node, XML_TAG_SPOTCOMPONENT_ROOT));

  shared_ptr<ihg::SpotComponent> spotcomponent(ReadSpotComponent(nodeRoot));
  if ( !spotcomponent )
  {
      typedef TypeMismatchException Exception;

      throw EXCEPTION_MSG
            (
              ITO33_BAD_DATA,
              String::Printf
              (
                TRANS("Node \"%s\" doesn't contain valid spot component data."),
                XML_TAG_SPOTCOMPONENT_ROOT
              )
            );
  }

  return spotcomponent;
}

}
}
}

namespace ito33
{

//------------------------------- power ---------------------------------------
static ito33::ihg::SpotComponent*
ReadHRSpotComponentPower(const xml::node *pNode) 
{
  return 
    new ito33::ihg::HRSpotComponentPower
            (
              GetDoubleFromName(*pNode, XML_TAG_SPOTCOMPONENT_BETA),
              GetDoubleFromName(*pNode, XML_TAG_SPOTCOMPONENT_S0)
            );
}



ITO33_DEFINE_SPOTCOMPONENT_READER(XML_TAG_SPOTCOMPONENT_POWER_ROOT, HRSpotComponentPower);

}
