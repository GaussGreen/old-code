/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/mesh/optionspacemesh.h
// Purpose:     do space mesh for options
// Author:      (z)
// Created:     03/10/22
// RCS-ID:      $Id: optionspacemesh.h,v 1.9 2006/03/22 13:10:22 yann Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_PRICING_GRIDS_OPTIONSPACEMESH_H_
#define _ITO33_PRICING_GRIDS_OPTIONSPACEMESH_H_

#include "ito33/beforestd.h"
#include <vector>
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/debug.h"
#include "ito33/common.h"
#include "ito33/pricing/optionliketype.h"


namespace ito33
{

namespace numeric
{

namespace mesh
{

/// Declaration of the MeshForLogS class
class MeshForLogS
{
public:
  MeshForLogS() : m_dSMin(-8), m_dSMax(8), m_dRho(1.1),
    m_dDxMax(0.5), m_dMaxDxDiff(0.15) {}

  /**
    Set stretch value

    @param dRho (input) the stretch value

    REQUIRE: Stretch value must be >= 1
    */
  void SetStretch(double dRho)
  {
    ASSERT_MSG(dRho >= 1,
      "SetStretch(): the stretch value must be greater than 1!");

    m_dRho = dRho;
  }

  
  /**
    Set the maximum tolerated mesh size at interested region

    @param dDx (input) the maximum tolerated mesh size

    REQUIRE: input value must be greater than 0
    */
  void SetMaxDxDiff(double dDx)
  {
    ASSERT_MSG(dDx <= 0,
      "SetMaxDxDiff(): the input value must be greater than 0!");

    m_dMaxDxDiff = dDx;
  }

protected:
  double
    /// Left limit point
    m_dSMin,
    /// Right limit point
    m_dSMax,
    /// stretch value
    m_dRho,
    /// the max tolerated mesh size at infinity
    m_dDxMax,
    /// the min tolerated mesh size at interested region
    m_dMaxDxDiff;
};

/// Declaration of the OptionSpaceMesh
class OptionSpaceMesh : public MeshForLogS
{
public:
  OptionSpaceMesh() :
      m_dHorizonMin(2), m_dHorizonMax(4),
      m_nMinNumber(200), m_nMaxNumber(6000),
      m_dDiff(0), m_dConv(0),
      m_dExtraHorizon(1.0),
      m_optiontype(ito33::pricing::Option_Call) {}

  void Build(size_t nNbNodes, double dLogSpot, std::vector<double> &vecGrid);

  /**
    Set the base diffusion size of the PDE

    @param dDiffusion (input) diffusion size
    */
  void SetDiffusionSize(double dDiffusion)
  {
    ASSERT_MSG(dDiffusion >= 0,
      "The diffusion size must zero or positive.");

    m_dDiff = dDiffusion; 

    // Limit the range of diffusion.  Small/large values can be obtained
    // for some of the volatility types (e.g. power vol at strike of zero
    // for put-call parity calculations is huge).
    if (m_dDiff < 0.01)
      m_dDiff = 0.01;
    else if (m_dDiff > 5.0)
      m_dDiff = 5.0;
  }

  
  /**
    Set the base convection size of the PDE

    @param dConvection (input) convection size
    */
  void SetConvectionSize(double dConvection)
  {
    m_dConv = dConvection;
  }


  /**
    Set option type

    @param optionType (input) the type of option
    */
  void SetOptionType(pricing::OptionLikeType optionType)
  {
    m_optiontype = optionType; 
  }

  /**
    Set extra refinement when spot is far from strike

    @param dExtraHorizon Multiplication to horizon
    */
  void SetExtraHorizon(double dExtraHorizon)
  {
    m_dExtraHorizon = dExtraHorizon;
    if (dExtraHorizon > 6.0)
      m_dExtraHorizon = 6.0;
    else if (dExtraHorizon < 1.2)
      m_dExtraHorizon = 1.0;
  }

protected:

  void BuildWithBarriers(size_t nNbNodes, double dLogSpot, std::vector<double> &vecGrid);
  void BuildWithoutBarriers(size_t nNbNodes, double dLogSpot, std::vector<double> &vecGrid);

  const double
    m_dHorizonMin,
    m_dHorizonMax;

  const size_t
    m_nMinNumber,
    m_nMaxNumber;

  double
    /// base diffusion size
    m_dDiff;
  double
    /// estimate length of the convection zone
    m_dConv;

  double
    /// Extra refinement in case the spot is far from the strike
    m_dExtraHorizon;

  ito33::pricing::OptionLikeType m_optiontype;

  double
    m_dUpperBarrier,
    m_dLowerBarrier;

  NO_COPY_CLASS(OptionSpaceMesh);
};


} // namespace mesh
} // namespace numeric
} // namespace ito33

#endif
