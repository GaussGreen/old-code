#include "ito33/beforestd.h"
#include <vector>
//#include <iostream>
#include "ito33/afterstd.h"
#include "ito33/numeric/mesh/genraf.h"
#include "ito33/numeric/mesh/barrierspacemesh.h"

using ito33::numeric::mesh::GenRaf;

namespace ito33
{

namespace numeric
{

namespace mesh
{

void BarrierSpaceMesh::Build(size_t nNbNodes, double dSpot,
                             std::vector<double> &pdMesh)
{
 
  // If the barriers are not defined, use sensible bounds.
  // Since this is internal code, assume barriers have been properly set.
  // Do not start at zero to be consistent with option meshes.
  // Refine around barriers if they are defined.
  double dMinSpot = 1.e-6;
  bool bRefineMinSpot = false;
  if (m_dLowerBarrier > 0.0)
  {
    dMinSpot = m_dLowerBarrier;
    bRefineMinSpot = true;
  }

  double dMaxSpot;
  if (dSpot > m_dLowerBarrier)
    dMaxSpot = dSpot * 100.0;
  else
    dMaxSpot = m_dLowerBarrier * 100.0;
  bool bRefineMaxSpot = false;
  if (m_dUpperBarrier > 0.0)
  {
    dMaxSpot = m_dUpperBarrier;
    bRefineMaxSpot = true;
  }

  bool bHasMidSpot = false;
  double dMidSpot = (dMaxSpot + dMinSpot)/2.0;
  if (m_dSpecialPoint > dMinSpot && m_dSpecialPoint < dMaxSpot)
  {
    dMidSpot = m_dSpecialPoint;
    bHasMidSpot = true;
  }

  // Brief testing shows that generating the mesh in log space gives
  // better results.
  dMinSpot = log(dMinSpot);  
  dMaxSpot = log(dMaxSpot);
  dMidSpot = log(dMidSpot);

  // Only refine the special point if requested, and if it is far enough 
  // away from the barriers.
  bool bRefineMidSpot = m_bSpecialRefine && bHasMidSpot;
  if (bRefineMidSpot)
  {
    double dMinSpotCheck = dMinSpot + m_dDiffusion + m_dConvection;
    double dMaxSpotCheck = dMaxSpot - m_dDiffusion - m_dConvection;

    if ( dMidSpot < dMinSpotCheck || dMidSpot > dMaxSpotCheck )
      bRefineMidSpot = false;
  }

  // Some default params for genraf
  double dDelta = (dMaxSpot - dMinSpot) / nNbNodes;
  double dDeltaMin = dDelta / 100.0;
  double dDeltaMax = (dMaxSpot - dMinSpot) / 20.0;
  
  // Let the stretch by 'x' be for 150 nodes. Then scale so that the
  // appropriate number of nodes are generated for larger meshes
  double dStretchTarget = 1.15;
  double dStretch = exp( log(dStretchTarget) / (nNbNodes / 150.0 ) );

  // Generate in a temp mesh first.  Copy later.
  int iNbSGuess;
  int iNbS;
  int iNbSTmp;
  std::vector<double> pdTmpMesh;

  // There are three possible special points (min, mid, max), and each one
  // can either be refined or not.  This gives 8 possible combinations.
  // While we could simply generate a grid for each section, splitting
  // them up into different cases allows for different mesh generation
  // parameters in each case.  Although tedious, hopefully the code is
  // easy to follow and debug.
  if ( bRefineMinSpot == true && 
       bRefineMidSpot == false && 
       bRefineMaxSpot == false )
  {
    // Refine at left and expand to right.
    dStretchTarget = 1.05;
    dStretch = exp( log(dStretchTarget) / (nNbNodes / 150.0 ) );

    iNbSGuess = GenRafSize(dMaxSpot-dMinSpot, dDeltaMin, dDeltaMax, dStretch);
    pdTmpMesh.resize(iNbSGuess);

    if (bHasMidSpot == false)
      numeric::mesh::GenRaf(dMinSpot, dMaxSpot, dDeltaMin, dDeltaMax, 
                            dStretch, 1, &pdTmpMesh[0], iNbS);
    else
    {
      numeric::mesh::GenRaf(dMinSpot, dMidSpot, dDeltaMin, dDeltaMax, 
                            dStretch, 1, &pdTmpMesh[0], iNbSTmp);
      iNbS = iNbSTmp - 1;

      double dSpace = pdTmpMesh[iNbS] - pdTmpMesh[iNbS-1];
      numeric::mesh::GenRaf(dMidSpot, dMaxSpot, dSpace, dDeltaMax, 
                            dStretch, 1, &pdTmpMesh[iNbS], iNbSTmp);
      iNbS += iNbSTmp;
    }

  }
  else if ( bRefineMinSpot == false && 
            bRefineMidSpot == false && 
            bRefineMaxSpot == true )
  {
    // Refine at right and expand to left. 
    dStretchTarget = 1.04;
    dStretch = exp( log(dStretchTarget) / (nNbNodes / 150.0 ) );

    iNbSGuess = GenRafSize(dMaxSpot-dMinSpot, dDeltaMin, dDeltaMax, dStretch);
    pdTmpMesh.resize(iNbSGuess);
    
    if (bHasMidSpot == false)
      numeric::mesh::GenRaf(dMinSpot, dMaxSpot, dDeltaMin, dDeltaMax, 
                            dStretch, 0, &pdTmpMesh[0], iNbS);
    else
    {
      numeric::mesh::GenRaf(dMidSpot, dMaxSpot, dDeltaMin, dDeltaMax, 
                            dStretch, 1, &pdTmpMesh[0], iNbSTmp);

      double dSpace = pdTmpMesh[1] - pdTmpMesh[0];

      numeric::mesh::GenRaf(dMinSpot, dMidSpot, dSpace, dDeltaMax, 
                            dStretch, 1, &pdTmpMesh[0], iNbSTmp);
      iNbS = iNbSTmp - 1;

      numeric::mesh::GenRaf(dMidSpot, dMaxSpot, dDeltaMin, dDeltaMax, 
                            dStretch, 1, &pdTmpMesh[iNbS], iNbSTmp);
      iNbS += iNbSTmp;
    }
  }
  else if ( bRefineMinSpot == true && 
            bRefineMidSpot == false && 
            bRefineMaxSpot == true )
  {
    // Refine at left and right barriers, but not in the middle.
    // dMidSpot should be defined (defaults to halfway between min and max)
    iNbSGuess = GenRafSize(dMidSpot-dMinSpot, dDeltaMin, dDeltaMax, dStretch);
    iNbSGuess += GenRafSize(dMaxSpot-dMidSpot, dDeltaMin, dDeltaMax, dStretch);

    pdTmpMesh.resize(iNbSGuess);

    numeric::mesh::GenRaf(dMinSpot, dMidSpot, dDeltaMin, dDeltaMax, 
                          dStretch, 0, &pdTmpMesh[0], iNbSTmp);
    iNbS = iNbSTmp - 1;
    
    numeric::mesh::GenRaf(dMidSpot, dMaxSpot, dDeltaMin, dDeltaMax, 
                          dStretch, 1, &pdTmpMesh[iNbS-1], iNbSTmp);
    iNbS += iNbSTmp;

  }
  else if ( bRefineMinSpot == true && 
            bRefineMidSpot == true && 
            bRefineMaxSpot == false )
  {
    // Refine at left barrier, and the spot
    // The mesh will be in three parts...min <-> midLeft, midLeft <-> mid, 
    // mid <-> max
    double dMidLeft = (dMinSpot + dMidSpot)/2.0;

    iNbSGuess =  GenRafSize(dMidLeft-dMinSpot, dDeltaMin, dDeltaMax, dStretch);
    iNbSGuess += GenRafSize(dMidSpot-dMidLeft, dDeltaMin, dDeltaMax, dStretch);
    iNbSGuess += GenRafSize(dMaxSpot-dMidSpot, dDeltaMin, dDeltaMax, dStretch);

    pdTmpMesh.resize(iNbSGuess);

    // Generate the mesh in three stages
    numeric::mesh::GenRaf(dMinSpot, dMidLeft, dDeltaMin, dDeltaMax, 
                          dStretch, 1, &pdTmpMesh[0], iNbSTmp);
    iNbS = iNbSTmp - 1;

    numeric::mesh::GenRaf(dMidLeft, dMidSpot, dDeltaMin, dDeltaMax, 
                          dStretch, 0, &pdTmpMesh[iNbS], iNbSTmp);
    iNbS += iNbSTmp - 1;

    numeric::mesh::GenRaf(dMidSpot, dMaxSpot, dDeltaMin, dDeltaMax, 
                          dStretch, 1, &pdTmpMesh[iNbS], iNbSTmp);
    iNbS += iNbSTmp;

  }
  else if ( bRefineMinSpot == false && 
            bRefineMidSpot == true && 
            bRefineMaxSpot == true )
  {
    // Refine at right barrier, and the spot

    // The mesh will be in three parts...min <-> midSpot, midSpot <-> midRight, 
    // midRight <-> max
    double dMidRight = (dMidSpot + dMaxSpot)/2.0;
    dStretchTarget = 1.1;
    dStretch = exp( log(dStretchTarget) / (nNbNodes / 150.0 ) );

    iNbSGuess =  GenRafSize(dMidSpot-dMinSpot, dDeltaMin,
                            dDeltaMax, dStretch);
    iNbSGuess += GenRafSize(dMidRight-dMidSpot, dDeltaMin,
                            dDeltaMax, dStretch);
    iNbSGuess += GenRafSize(dMaxSpot-dMidRight, dDeltaMin/20.0,
                            dDeltaMax, dStretch);

    pdTmpMesh.resize(iNbSGuess);

    // Generate the mesh in three stages
    numeric::mesh::GenRaf(dMinSpot, dMidSpot, dDeltaMin, dDeltaMax, 
                          dStretch, 0, &pdTmpMesh[0], iNbSTmp);
    iNbS = iNbSTmp - 1;

    numeric::mesh::GenRaf(dMidSpot, dMidRight, dDeltaMin, dDeltaMax, 
                          dStretch, 1, &pdTmpMesh[iNbS], iNbSTmp);
    iNbS += iNbSTmp - 1;

    numeric::mesh::GenRaf(dMidRight, dMaxSpot, dDeltaMin/20.0, dDeltaMax, 
                          dStretch, 0, &pdTmpMesh[iNbS], iNbSTmp);
    iNbS += iNbSTmp;

  }
  else
  {
    // The mesh will be in four parts...min <-> midLeft, midLeft <-> mid, 
    // mid <-> midRight, midRight <-> max
    double dMidLeft  = (dMinSpot + dMidSpot)/2.0;
    double dMidRight = (dMidSpot + dMaxSpot)/2.0;

    iNbSGuess  = GenRafSize(dMidLeft-dMinSpot, dDeltaMin, dDeltaMax, dStretch);
    iNbSGuess += GenRafSize(dMidSpot-dMidLeft, dDeltaMin, dDeltaMax, dStretch);
    iNbSGuess += GenRafSize(dMidRight-dMidSpot, dDeltaMin,dDeltaMax, dStretch);
    iNbSGuess += GenRafSize(dMaxSpot-dMidRight, dDeltaMin,dDeltaMax, dStretch);

    pdTmpMesh.resize(iNbSGuess);

    // Generate the mesh in four stages
    numeric::mesh::GenRaf(dMinSpot, dMidLeft, dDeltaMin, dDeltaMax, 
                          dStretch, 0, &pdTmpMesh[0], iNbSTmp);
    iNbS = iNbSTmp - 1;

    numeric::mesh::GenRaf(dMidLeft, dMidSpot, dDeltaMin, dDeltaMax, 
                          dStretch, 1, &pdTmpMesh[iNbS], iNbSTmp);
    iNbS += iNbSTmp - 1;

    numeric::mesh::GenRaf(dMidSpot, dMidRight, dDeltaMin, dDeltaMax, 
                          dStretch, 0, &pdTmpMesh[iNbS], iNbSTmp);
    iNbS += iNbSTmp - 1;

    numeric::mesh::GenRaf(dMidRight, dMaxSpot, dDeltaMin, dDeltaMax, 
                          dStretch, 1, &pdTmpMesh[iNbS], iNbSTmp);
    iNbS += iNbSTmp;

  }

  // Copy to the real grid
  pdMesh.resize(iNbS);
  for (int nIdx = 0; nIdx < iNbS; nIdx++)
  {
    //pdMesh[nIdx] = pdTmpMesh[nIdx];
    pdMesh[nIdx] = exp(pdTmpMesh[nIdx]);
    //std::cout << nIdx << " " << pdMesh[nIdx] << std::endl;
  }

  // Smoothing
  //for (int nIdx = 1; nIdx < iNbS - 1; nIdx++)
  //{
  //  pdMesh[nIdx] = (pdMesh[nIdx - 1] + pdMesh[nIdx + 1]) / 2.0;
  //}

  //std::cout << "requested = " << nNbNodes << ", generated = " << iNbS << std::endl;

}

} // namespace mesh

} // namespace numeric

} // namespace ito33
