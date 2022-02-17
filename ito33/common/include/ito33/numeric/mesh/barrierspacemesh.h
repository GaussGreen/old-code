/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/mesh/barrierspacemesh.h
// Purpose:     Space mesh for contracts with barriers (eds, onetouch, etc.)
// Created:     2005/07/18
// RCS-ID:      $Id: barrierspacemesh.h,v 1.1 2005/07/21 20:52:03 dave Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_NUMERIC_MESH_BARRIERSPACEMESH_H_
#define _ITO33_NUMERIC_MESH_BARRIERSPACEMESH_H_

#include "ito33/beforestd.h"
#include <vector>
#include "ito33/afterstd.h"

#include "ito33/debug.h"
#include "ito33/common.h"

namespace ito33
{

namespace numeric
{

namespace mesh
{


class BarrierSpaceMesh
{

public:
  
  BarrierSpaceMesh()
    : m_dLowerBarrier(-1.0), m_dUpperBarrier(-1.0), 
      m_dDiffusion(0.01), m_dConvection(0.2),
      m_dSpecialPoint(0.0), m_bSpecialRefine(false)
  {}

  /** 
    Set the lower barrier.

    @param dLeftBarrier The position of the left barrier.
  */
  void SetLowerBarrier(double dLowerBarrier)
  {
    m_dLowerBarrier = dLowerBarrier;
  }

  /** 
    Set the upper barrier.

    @param dUpperBarrier The position of the upper barrier.
  */
  void SetUpperBarrier(double dUpperBarrier)
  {
    m_dUpperBarrier = dUpperBarrier;
  }

  /**
    Build the mesh

    @param nNbNodes The number of nodes to use (approximate)
    @param dSpot The current spot price
    @param vecGrid The final grid (output)
   */
  void Build(size_t nNbNodes, double dSpot, std::vector<double> &vecGrid);

  /**
    Set the base diffusion size of the PDE

    @param dDiffusion (input) diffusion size
    */
  void SetDiffusionSize(double dDiffusion)
  {
    ASSERT_MSG(dDiffusion >= 0,
      "The diffusion size must zero or positive.");

    m_dDiffusion = dDiffusion; 

    if (m_dDiffusion < 0.01)
      m_dDiffusion = 0.01;
  }

  /**
    Set the base convection size of the PDE

    @param dConvection (input) convection size
    */
  void SetConvectionSize(double dConvection)
  {
    m_dConvection = dConvection;
  }

  /**
    Set a special point to be in the mesh

    @param dPoint the position of the special point
    @param bRefine whether or not to refine around the special point
   */
  void SetSpecialPoint(double dPoint, bool bRefine)
  {
    m_dSpecialPoint = dPoint;
    m_bSpecialRefine = bRefine;
  }

protected:

  /// base diffusion size
  double m_dDiffusion;

  /// estimate length of the convection zone
  double m_dConvection;

  /// The upper barrier position
  double m_dUpperBarrier;

  /// The lower barrier position
  double m_dLowerBarrier;

  /// Special point to include in mesh
  double m_dSpecialPoint;

  /// Whether or not to refine at the special point
  bool m_bSpecialRefine;

  NO_COPY_CLASS(BarrierSpaceMesh);
};


} // namespace mesh
} // namespace numeric
} // namespace ito33

#endif  // _ITO33_NUMERIC_MESH_BARRIERSPACEMESH_H_
