/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/xml/hedgereader.cpp
// Purpose:     code for reading pricing HG XML files
// Created:     2005/06/10
// RCS-ID:      $Id: hedgereader.cpp,v 1.5 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/finance/derivatives.h"

#include "hg/xml/common.h"
#include "hg/xml/hedgereader.h"

#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/derivative.h"

#include "ito33/xml/read.h"

#include "ito33/hg/theoreticalmodel.h"

namespace ito33
{

namespace hg
{


shared_ptr<finance::Derivative> 
XML::HedgeReader::ReadTarget() const
{
  shared_ptr<finance::Derivative> pDerivative;

  xml::node mainNode = GetMainNode();

  xml::node::const_iterator 
    pNodeDerivative = mainNode.find(XML_TAG_HG_TARGET);

  // If the target tag is not found, return empty derivative pointer
  if ( pNodeDerivative == mainNode.end() )
    return pDerivative;

  shared_ptr<finance::SessionData> pSessionData( ReadSessionData() );

  const xml::node::const_iterator end = pNodeDerivative->end();
  for ( xml::node::const_iterator i = pNodeDerivative->begin(); i != end; ++i )
  {
    shared_ptr<finance::Derivative>
      derivative( DerivativeFactory::Create(i->get_name(), &(*i)) );

    if ( derivative )
    {
      derivative->SetSessionData(pSessionData);

      pDerivative = derivative;

      return pDerivative;
    }
  }

  return pDerivative;
}

shared_ptr<finance::Derivatives>
XML::HedgeReader::ReadHedgeInstruments() const
{
  shared_ptr<finance::Derivatives> pDerivatives;

  xml::node mainNode = GetMainNode();

  xml::node::const_iterator 
    pNodeDerivatives = mainNode.find(XML_TAG_HG_HEDGEINSTRUMENTS);

  // If the hedging instrument tag is not found, return empty derivatives 
  // list. Possibly hedge only with the underlying
  if ( pNodeDerivatives == mainNode.end() )
    return pDerivatives;

  pDerivatives = shared_ptr<finance::Derivatives>(new finance::Derivatives);

  shared_ptr<finance::SessionData> pSessionData( ReadSessionData() );

  const xml::node::const_iterator end = pNodeDerivatives->end();
  for ( xml::node::const_iterator i = pNodeDerivatives->begin(); i != end; ++i )
  {
    shared_ptr<finance::Derivative>
      pDerivative( DerivativeFactory::Create(i->get_name(), &(*i)) );

    if ( pDerivative )
    {
      pDerivative->SetSessionData(pSessionData);

      pDerivatives->Add(pDerivative);
    }
    //else: ignore unknown derivatives
  }

  return pDerivatives;
}

bool
XML::HedgeReader::ReadComputeHERO() const
{
  xml::node mainNode = GetMainNode();

  xml::node::const_iterator 
    pNodeComputeHERO = mainNode.find(XML_TAG_HG_COMPUTEHERO);

  // If the compute hero tag is not found, return false by default
  if ( pNodeComputeHERO == mainNode.end() )
    return false;

  return ito33::XML::GetBoolFromNode(*pNodeComputeHERO);
}

} // namespace hg

} // namespace ito33
