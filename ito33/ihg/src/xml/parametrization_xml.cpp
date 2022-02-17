/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/xml/parametrization_xml.cpp
// Purpose:     Restore parametrization objects from XML document
// Author:      ITO33
// Created:     2004/11/25
// RCS-ID:      $Id: parametrization_xml.cpp,v 1.14 2006/08/20 09:31:05 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/*
              --   README   --

  How to add a new parametrization, for example ParametrizationXXX?

  1. Go to "FIRST PART", add code looks like

      static Parametrization*
      ReadParametrizationXXX(const xml::node *pNode);

      IMPLEMENT_RESTORE_PARAM                     \
          (XML_TAG_PARAMETRIZATION_XXX_ROOT, \
            ParametrizationXXX)

  2. Go to "SECOND PART", add implementation of the function
     ReadParametrizationXXX()

  */
#include "ito33/useexception.h"
#include "ito33/sharedptr.h"
#include "ito33/link.h"


#include "ito33/ihg/parametrization_hrwithtimecomponent.h"
#include "ito33/ihg/parametrization_hrwithspotcomponentpower.h"
#include "ito33/ihg/parametrization_volflat_hrwithtimecomponent.h"
#include "ito33/ihg/parametrization_volflat_hrwithspotcomponentpower.h"
#include "ito33/ihg/parametrization_volflat_hrpower.h"
#include "ito33/ihg/parametrization_volpower.h"
#include "ito33/ihg/parametrization_volpower_hrpower.h"
#include "ito33/ihg/parametrization_volpower_hrwithtimecomponent.h"
#include "ito33/ihg/parametrization_voltanh.h"
#include "ito33/ihg/parametrization_voltanh_hrwithtimecomponent.h"
#include "ito33/ihg/parametrization_voltanh_hrwithspotcomponentpower.h"
#include "ito33/ihg/parametrization_voltanh_hrpower.h"
#include "ito33/ihg/parametrization_volwithtimecomponent.h"
#include "ito33/ihg/parametrization_voltimeonly_hrwithtimecomponent.h"

#include "ihg/xml/parametrization.h"
#include "ihg/xml/spotcomponent.h"
#include "ihg/xml/volatility.h"
#include "ihg/xml/hazardrate.h"

#include "ito33/xml/read.h"

extern const ito33::Error ITO33_UNEXPECTED;

using namespace ito33;
using namespace ito33::XML;


// define the field of the Factory declared in ihg/xml/parametrization.h
ITO33_IMPLEMENT_THE_FACTORY(ito33::ihg::ParametrizationFactory);

ITO33_FORCE_LINK_THIS_MODULE(ihg_parametrization_xml);

/*
  macro to 
  1. define function Resotre(const xml::node&, shared_ptr<ParamName>&)
  2. ITO33_DEFINE_PARAMETRIZATION_READER for tag ParamName 
  */
#define IMPLEMENT_RESTORE_PARAM(tag, name)                                   \
  bool Restore(const xml::node &node, shared_ptr< name >& pParam)            \
  {                                                                          \
    if ( strcmp(node.get_name(), tag) != 0 )                                 \
      return false;                                                          \
                                                                             \
    pParam = shared_ptr< name >                                              \
                  (static_cast< name *> (Read ## name (&node)) );            \
    return true;                                                             \
  }                                                                          \
  ITO33_DEFINE_PARAMETRIZATION_READER(tag, name);                            \


namespace ito33
{

namespace ihg
{

namespace XML
{

/**************************** FIRST PART *************************************/

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static Parametrization*
ReadParametrizationHRWithTimeComponent(const xml::node *pNode);

IMPLEMENT_RESTORE_PARAM                     \
     (XML_TAG_PARAMETRIZATION_HRWITHTIMECOMPONENT_ROOT, \
      ParametrizationHRWithTimeComponent);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static Parametrization*
ReadParametrizationHRWithSpotComponentPower(const xml::node *pNode);

IMPLEMENT_RESTORE_PARAM                          \
     (XML_TAG_PARAMETRIZATION_HRWITHSPOTCOMPONENTPOWER_ROOT, \
      ParametrizationHRWithSpotComponentPower);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static Parametrization*
ReadParametrizationVolFlatHRWithTimeComponent(const xml::node *pNode);

IMPLEMENT_RESTORE_PARAM                          \
     (XML_TAG_PARAMETRIZATION_VOLFLATHRWITHTIMECOMPONENT_ROOT, \
      ParametrizationVolFlatHRWithTimeComponent);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static Parametrization*
ReadParametrizationVolFlatHRWithSpotComponentPower(const xml::node *pNode);

IMPLEMENT_RESTORE_PARAM                          \
     (XML_TAG_PARAMETRIZATION_VOLFLATHRWITHSPOTCOMPONENTPOWER_ROOT, \
      ParametrizationVolFlatHRWithSpotComponentPower);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static Parametrization*
ReadParametrizationVolFlatHRPower(const xml::node *pNode);

IMPLEMENT_RESTORE_PARAM                          \
     (XML_TAG_PARAMETRIZATION_VOLFLATHRPOWER_ROOT, \
      ParametrizationVolFlatHRPower);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static Parametrization*
ReadParametrizationVolPower(const xml::node *pNode);

IMPLEMENT_RESTORE_PARAM                          \
     (XML_TAG_PARAMETRIZATION_VOLPOWER_ROOT,                 \
      ParametrizationVolPower);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static Parametrization*
ReadParametrizationVolPowerHRPower(const xml::node *pNode);

IMPLEMENT_RESTORE_PARAM                          \
     (XML_TAG_PARAMETRIZATION_VOLPOWERHRPOWER_ROOT,          \
      ParametrizationVolPowerHRPower);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static Parametrization*
ReadParametrizationVolPowerHRWithTimeComponent(const xml::node *pNode);

IMPLEMENT_RESTORE_PARAM                             \
     (XML_TAG_PARAMETRIZATION_VOLPOWERHRWITHTIMECOMPONENT_ROOT, \
      ParametrizationVolPowerHRWithTimeComponent);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static Parametrization*
ReadParametrizationVolTanh(const xml::node *pNode);

IMPLEMENT_RESTORE_PARAM                         \
     (XML_TAG_PARAMETRIZATION_VOLTANH_ROOT,                 \
      ParametrizationVolTanh);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static Parametrization*
ReadParametrizationVolTanhHRWithSpotComponentPower(const xml::node *pNode);

IMPLEMENT_RESTORE_PARAM                         \
     (XML_TAG_PARAMETRIZATION_VOLTANHHRWITHSPOTCOMPONENTPOWER_ROOT,                 \
      ParametrizationVolTanhHRWithSpotComponentPower);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static Parametrization*
ReadParametrizationVolTanhHRWithTimeComponent(const xml::node *pNode);

IMPLEMENT_RESTORE_PARAM                         \
     (XML_TAG_PARAMETRIZATION_VOLTANHHRWITHTIMECOMPONENT_ROOT,                 \
      ParametrizationVolTanhHRWithTimeComponent);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static Parametrization*
ReadParametrizationVolTanhHRPower(const xml::node *pNode);

IMPLEMENT_RESTORE_PARAM                         \
     (XML_TAG_PARAMETRIZATION_VOLTANHHRPOWER_ROOT,                 \
      ParametrizationVolTanhHRPower);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static Parametrization*
ReadParametrizationVolWithTimeComponent(const xml::node *pNode);

IMPLEMENT_RESTORE_PARAM                         \
     (XML_TAG_PARAMETRIZATION_VOLWITHTIMECOMPONENT_ROOT,                 \
      ParametrizationVolWithTimeComponent);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static Parametrization*
ReadParametrizationVolTimeOnlyHRWithTimeComponent(const xml::node *pNode);

IMPLEMENT_RESTORE_PARAM                         \
     (XML_TAG_PARAMETRIZATION_VOLTIMEONLYHRWITHTIMECOMPONENT_ROOT,                 \
      ParametrizationVolTimeOnlyHRWithTimeComponent);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// remove local macro IMPLEMENT_RESTORE_PARAM
#undef IMPLEMENT_RESTORE_PARAM


/**************************** SECOND PART **************************************

/**
    Restore a ParametrizationHRWithTimeComponent object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node containing the parametrization data tag in DOM tree
    @return the new parametrization object
 */
/*static */Parametrization*
ReadParametrizationHRWithTimeComponent(const xml::node *pNode)
{
  // Create the return object
  ParametrizationHRWithTimeComponent* 
    pParametrization = new ParametrizationHRWithTimeComponent();

  // This parametrization can accept a volatility, spot component.  
  // Loop over the contents of this node looking for data
  xml::node::const_iterator i;

  // Read the spot component, if any
  shared_ptr<ihg::SpotComponent> pSpotComponent;
  if ( ( i = pNode->find(XML_TAG_SPOTCOMPONENT_ROOT) ) != pNode->end() )
    pSpotComponent = ReadSpotComponent(*i);

  // Read the volatility, if any
  shared_ptr<ihg::Volatility> pVolatility = ReadVolatility(*pNode);

  if ( pVolatility )
    pParametrization->SetVolatility(pVolatility);

  if ( pSpotComponent )
  {
    pParametrization->SetSpotComponent(pSpotComponent);
  }

  return pParametrization;
}

/**
    Restore a ParametrizationHRWithSpotComponentPower object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node containing the parametrization data tag in DOM tree
    @return the new parametrization object
 */
Parametrization*
ReadParametrizationHRWithSpotComponentPower(const xml::node *pNode)
{
  // This parametrization is constructed with a volatility. It does not
  // accept any other data

  // Read the volatility
  shared_ptr<Volatility> pVolatility = ReadVolatility(*pNode);

  // Verify that a volatility was read
  if ( !pVolatility )
  {
    throw ito33::EXCEPTION_MSG
            (
              ITO33_UNEXPECTED,
              TRANS("Volatility undefined when reading " \
                    "ParametrizationHRWithSpotComponentPower.")
            );
  }

  // Create and return the parametrization
  ParametrizationHRWithSpotComponentPower* 
    pParametrization = new ParametrizationHRWithSpotComponentPower(pVolatility);

  return pParametrization;
}


/**
    Restore a ParametrizationVolFlatHRWithTimeComponent object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node containing the parametrization data tag in DOM tree
    @return the new parametrization object
 */
Parametrization*
ReadParametrizationVolFlatHRWithTimeComponent(const xml::node *pNode)
{
  // Create the return object
  ParametrizationVolFlatHRWithTimeComponent* 
    pParametrization = new ParametrizationVolFlatHRWithTimeComponent();

  // This parametrization can accept a spot component.  
  // Loop over the contents of this node looking for data
  xml::node::const_iterator i;

  // Read the spot component, if any
  shared_ptr<ihg::SpotComponent> pSpotComponent;
  if ( ( i = pNode->find(XML_TAG_SPOTCOMPONENT_ROOT) ) != pNode->end() )
    pSpotComponent = ReadSpotComponent(*i);

  if ( pSpotComponent )
  {
    pParametrization->SetSpotComponent(pSpotComponent);
  }

  return pParametrization;

}


/**
    Restore a ParametrizationVolFlatHRWithSpotComponentPower object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node containing the parametrization data tag in DOM tree
    @return the new parametrization object
 */
Parametrization*
ReadParametrizationVolFlatHRWithSpotComponentPower(const xml::node* /* pNode */)
{

  // Create the return object
  ParametrizationVolFlatHRWithSpotComponentPower* 
    pParametrization = new ParametrizationVolFlatHRWithSpotComponentPower();

  // Nothing to read.  No input data in this class
  return pParametrization;

}

/**
    Restore a ParametrizationVolFlatHRPower object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node containing the parametrization data tag in DOM tree
    @return the new parametrization object
 */
Parametrization*
ReadParametrizationVolFlatHRPower(const xml::node* /* pNode */)
{

  // Create the return object
  ParametrizationVolFlatHRPower* 
    pParametrization = new ParametrizationVolFlatHRPower();

  // Nothing to read.  No input data in this class
  return pParametrization;

}


/**
    Read a ParametrizationVolPower from XML.

    @param pNode the parametrization node in DOM tree
    @return pointer to new parametrization object
 */
Parametrization*
ReadParametrizationVolPower(const xml::node *pNode)
{

  // This parametrization is constructed with a hazard rate. It does not
  // accept any other data

  // Read the hazard rate
  shared_ptr<HazardRate> pHazardRate = ReadHazardRate(*pNode);

  // Verify that a hazard rate was read
  if ( !pHazardRate )
  {
    throw ito33::EXCEPTION_MSG
            (
              ITO33_UNEXPECTED,
              TRANS("Hazard rate undefined when reading " \
                    "ParametrizationVolPower.")
            );
  }

  // Create and return the parametrization
  ParametrizationVolPower* 
    pParametrization = new ParametrizationVolPower(pHazardRate);

  return pParametrization;

}

/**
    Restore a ParametrizationVolPowerHRPower object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node containing the parametrization data tag in DOM tree
    @return the new parametrization object
 */
Parametrization*
ReadParametrizationVolPowerHRPower(const xml::node* /* pNode */)
{

  // Nothing to read.  No input data in this class
  // Create the return object
  return new ParametrizationVolPowerHRPower();
}



/**
    Restore a ParametrizationVolPowerHRWithTimeComponent object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node containing the parametrization data tag in DOM tree
    @return the new parametrization object
 */
Parametrization*
ReadParametrizationVolPowerHRWithTimeComponent(const xml::node *pNode)
{

  // Create the return object
  ParametrizationVolPowerHRWithTimeComponent* 
    pParametrization = new ParametrizationVolPowerHRWithTimeComponent();

  // This parametrization can accept a spot component.  
  // Loop over the contents of this node looking for data
  xml::node::const_iterator i;

  // Read the spot component, if any
  shared_ptr<ihg::SpotComponent> pSpotComponent;
  if ( ( i = pNode->find(XML_TAG_SPOTCOMPONENT_ROOT) ) != pNode->end() )
    pSpotComponent = ReadSpotComponent(*i);

  if ( pSpotComponent )
    pParametrization->SetSpotComponent(pSpotComponent);

  return pParametrization;

}


/**
    Read a ParametrizationVolTanh from XML.

    @param pNode the parametrization node in DOM tree
    @return pointer to new parametrization object
 */
Parametrization*
ReadParametrizationVolTanh(const xml::node *pNode)
{

  // This parametrization is constructed with a hazard rate. It does not
  // accept any other data

  // Read the hazard rate
  shared_ptr<HazardRate> pHazardRate = ReadHazardRate(*pNode);

  // Verify that a hazard rate was read
  if ( !pHazardRate )
  {
    throw ito33::EXCEPTION_MSG
            (
              ITO33_UNEXPECTED,
              TRANS("Hazard rate undefined when reading " \
                    "ParametrizationVolPower.")
            );
  }

  // Create and return the parametrization
  return new ParametrizationVolTanh(pHazardRate);
}

/**
    Restore a ParametrizationVolTanhHRWithSpotComponentPower object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node containing the parametrization data tag in DOM tree
    @return the new parametrization object
 */
Parametrization*
ReadParametrizationVolTanhHRWithSpotComponentPower(const xml::node* /* pNode */)
{

  // Create the return object
  ParametrizationVolTanhHRWithSpotComponentPower* 
    pParametrization = new ParametrizationVolTanhHRWithSpotComponentPower();

  // Nothing to read.  No input data in this class
  return pParametrization;

}

/**
    Restore a ParametrizationVolTanhHRPower object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node containing the parametrization data tag in DOM tree
    @return the new parametrization object
 */
Parametrization*
ReadParametrizationVolTanhHRPower(const xml::node* /* pNode */)
{

  // Create the return object
  ParametrizationVolTanhHRPower* 
    pParametrization = new ParametrizationVolTanhHRPower();

  // Nothing to read.  No input data in this class
  return pParametrization;

}

/**
    Restore a ParametrizationVolTanhHRWithTimeComponent object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node containing the parametrization data tag in DOM tree
    @return the new parametrization object
 */
Parametrization*
ReadParametrizationVolTanhHRWithTimeComponent(const xml::node *pNode)
{

  // Create the return object
  ParametrizationVolTanhHRWithTimeComponent* 
    pParametrization = new ParametrizationVolTanhHRWithTimeComponent();

  // This parametrization can accept a spot component.
  // Loop over the contents of this node looking for data
  xml::node::const_iterator i;

  // Read the spot component, if any
  shared_ptr<ihg::SpotComponent> pSpotComponent;
  if ( ( i = pNode->find(XML_TAG_SPOTCOMPONENT_ROOT) ) != pNode->end() )
    pSpotComponent = ReadSpotComponent(*i);

  if ( pSpotComponent )
    pParametrization->SetSpotComponent(pSpotComponent);

  return pParametrization;
}

/**
    Read a ParametrizationVolWithTimeComponent from XML.

    @param pNode the parametrization node in DOM tree
    @return pointer to new parametrization object
 */
Parametrization*
ReadParametrizationVolWithTimeComponent(const xml::node *pNode)
{
  // This parametrization is constructed with a hazard rate. It does not
  // accept any other data

  // Read the hazard rate
  shared_ptr<HazardRate> pHazardRate = ReadHazardRate(*pNode);

  // Verify that a hazard rate was read
  if ( !pHazardRate )
  {
    throw ito33::EXCEPTION_MSG
            (
              ITO33_UNEXPECTED,
              TRANS("Hazard rate undefined when reading " \
                    "ParametrizationVolPower.")
            );
  }

  return new ParametrizationVolWithTimeComponent(pHazardRate);
}

/**
    Restore a ParametrizationVolTimeOnlyHRWithTimeComponent object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node containing the parametrization data tag in DOM tree
    @return the new parametrization object
 */
Parametrization*
ReadParametrizationVolTimeOnlyHRWithTimeComponent(const xml::node *pNode)
{

  // Create the return object
  ParametrizationVolTimeOnlyHRWithTimeComponent* 
    pParametrization = new ParametrizationVolTimeOnlyHRWithTimeComponent();

  // This parametrization can accept a spot component.  
  // Loop over the contents of this node looking for data
  xml::node::const_iterator i;

  // Read the spot component, if any
  shared_ptr<ihg::SpotComponent> pSpotComponent;
  if ( ( i = pNode->find(XML_TAG_SPOTCOMPONENT_ROOT) ) != pNode->end() )
    pSpotComponent = ReadSpotComponent(*i);

  if ( pSpotComponent )
    pParametrization->SetSpotComponent(pSpotComponent);

  return pParametrization;
}


} // namespace xml

} // namespace ihg

} // namespace ito33
