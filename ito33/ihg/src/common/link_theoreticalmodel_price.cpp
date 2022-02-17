/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/link_theoreticalmodel_price.cpp
// Purpose:     force linking of theoreticalmodel_price*.cpp
// Author:      Wang
// Created:     2004/10/25
// RCS-ID:      $Id: link_theoreticalmodel_price.cpp,v 1.15 2006/08/10 23:29:23 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// this file can be included in the main project when linking of all 
// theoreticalmodel_price*.cpp is required.

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(IHGPriceAsianOption);

ITO33_FORCE_LINK_MODULE(IHGPriceOption);

ITO33_FORCE_LINK_MODULE(IHGPriceCDS);

ITO33_FORCE_LINK_MODULE(IHGPriceEDS);

ITO33_FORCE_LINK_MODULE(IHGPriceOneTouch);

ITO33_FORCE_LINK_MODULE(IHGPriceVarianceSwap);

ITO33_FORCE_LINK_MODULE(IHGPriceParBond);

ITO33_FORCE_LINK_MODULE(IHGPriceBond);

ITO33_FORCE_LINK_MODULE(IHGPriceAttachedWarrantConvertibleBond);

ITO33_FORCE_LINK_MODULE(IHGPriceCB);

ITO33_FORCE_LINK_MODULE(IHGPriceReset);

ITO33_FORCE_LINK_MODULE(IHGPricePEPSLike);

ITO33_FORCE_LINK_MODULE(IHGPriceGeneralizedPEPSLike);

ITO33_FORCE_LINK_MODULE(IHGPricePERCSLike);

ITO33_FORCE_LINK_MODULE(IHGPriceCBOption);
