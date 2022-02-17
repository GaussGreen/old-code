#include "ito33/beforestd.h"
#include <vector>
#include "ito33/afterstd.h"
#include "ito33/numeric/mesh/genraf.h"
#include "ito33/numeric/mesh/optionspacemesh.h"

#include "ito33/pricing/optionliketype.h"

using ito33::numeric::mesh::GenRaf;
using ito33::numeric::mesh::OptionSpaceMesh;


void OptionSpaceMesh::Build(size_t nNbNodes, double dLogSpot,
                            std::vector<double> &vecMeshLogS)
{
  /*
  The idea is to build this kind of mesh

                    A
                    *----------\
                                 \  B
                                   \*___________
                                                 \  C              0
                                                   \*______________*

  on the left side of 0. 
  - A is an approximation of -infty, where the mesh size can be large, let it 
    be m_dDxMax
  - [B,C] is a zone which play some role in the accuracy of the result of 
    calculation, but we don't want put too many mesh point here. Let the mesh
    size at B be 2 * dDx
  - [C, 0] is the zone that we want to put the most points, the mesh size at C
    and 0 is dDx.
  
  The right side is somehow symetric.

  How to calculate A, B, C

  - A is defined by m_dSMin
  - B is defined dHorizon * m_dDiff, where m_dDiff is the size of diffusion
    dHorizon is calculated in this way : its values is between m_dHorizonMin
    and m_dHorizonMax, and is a function of the number of nodes. When nNbNodes
    is greater than m_nMaxNumber(6000), we use m_dHorizonMax, when it is less
    than m_nMinNumber(200), we use m_dHorizonMin.
  - C is somehow m_dHorizonMin * m_dDiff.
  
  */

  double
    dHorizon,
    vdLeftSpot,
    vdRightSpot,
    dStrikeLeftG,
    dStrikeRightG;

  int
    viShift = 1;

  double
    dDxDiff;

  double
    dTmp;

  if(nNbNodes < m_nMinNumber)
    dHorizon = m_dHorizonMin;
  else if(nNbNodes > m_nMaxNumber)
    dHorizon = m_dHorizonMax;
  else
    dHorizon = double(nNbNodes - m_nMinNumber) / double(m_nMaxNumber - m_nMinNumber)
    * (m_dHorizonMax - m_dHorizonMin) + m_dHorizonMin;

  vdRightSpot = m_dDiff * m_dHorizonMin;
  vdLeftSpot = -vdRightSpot;

  // Add extra points in the right side of the grid when the spot is much
  // larger than the strike.  Do not add if the diffusion is already
  // quite large
  if (m_dExtraHorizon > 1.0 && m_dDiff < 3.0)
    vdRightSpot *= m_dExtraHorizon;
  

  // limit the convection size as the user may give ridiculous hazard rate.
  const double MAXCONVSIZE = 4;
  if(m_dConv > MAXCONVSIZE)
    m_dConv = MAXCONVSIZE;
  else if(m_dConv < -MAXCONVSIZE)
    m_dConv = -MAXCONVSIZE;

  if( m_dConv > 0)
    vdRightSpot += m_dConv;
  else
    vdLeftSpot += m_dConv;

  if(vdLeftSpot > dLogSpot - 0.5 * m_dDiff)
    vdLeftSpot = dLogSpot - m_dDiff;

  dStrikeLeftG = vdLeftSpot - (dHorizon - m_dHorizonMin) * m_dDiff;
  dStrikeRightG = vdRightSpot + (dHorizon - m_dHorizonMin) * m_dDiff;

  dDxDiff = (vdRightSpot - vdLeftSpot) / (double)nNbNodes;

  if(dDxDiff > m_dMaxDxDiff)
    dDxDiff = m_dMaxDxDiff;

  if(-vdLeftSpot < dDxDiff * 6)
    vdLeftSpot = - dDxDiff * 6;
  else if(dLogSpot - vdLeftSpot < 6. * dDxDiff)
    vdLeftSpot = dLogSpot - 6. * dDxDiff;
  if(vdRightSpot < dDxDiff * 6)
    vdRightSpot = dDxDiff * 6;

  double
    dDxPetit = 0.2 *
              dDxDiff,
    dRhoSpot = 1.1;
  double
    dRhoTmp = 0.;

  int I, J;

  J = nNbNodes;
  J += (int)( 2. * log(dDxDiff / dDxPetit) / log(dRhoSpot)
    + (vdRightSpot - vdLeftSpot) / dDxDiff);
  if( dStrikeRightG - dStrikeLeftG - vdRightSpot + vdLeftSpot > 0)
  {
    J += (int) 
      (
        ( dStrikeRightG - dStrikeLeftG - vdRightSpot + vdLeftSpot )
        * 0.5 / dDxDiff
      );
  }
  if((m_dSMax - m_dSMin - (vdRightSpot - vdLeftSpot)) > 0)
    J += (int) ( (m_dSMax - m_dSMin - (vdRightSpot - vdLeftSpot)) / m_dDxMax + 2. *
    log(m_dDxMax / dDxDiff) / log(2.) );
  
  int nLogSpots = J;
  double *pdLogSpots = new double [nLogSpots];
  nNbNodes = 1;
    
  if(viShift)
  {
    dRhoTmp = 2;
    pdLogSpots[0] = vdLeftSpot;
    double dDxTmp = dDxDiff;
    I = 0;
    while(pdLogSpots[I] > dStrikeLeftG)
    {
      pdLogSpots[I + 1] = pdLogSpots[I] - dDxTmp;
      I++;
      if(dDxTmp < dDxDiff * 2.)
      {
        dDxTmp *= dRhoSpot;
        if(dDxTmp > dDxDiff * 2.)
          dDxTmp = dDxDiff * 2.;
      }
    }
    if(dDxTmp > m_dDxMax)
      dDxTmp = m_dDxMax;
    while ( I < nLogSpots && pdLogSpots[I] > m_dSMin + m_dDxMax * 0.5
             )
    {
      pdLogSpots[I + 1] = pdLogSpots[I] - dDxTmp;
      I++;
      if(dDxTmp < m_dDxMax)
      {
        dDxTmp *= dRhoTmp;
        if(dDxTmp > m_dDxMax)
          dDxTmp = m_dDxMax;
      }
    }

    if(I == nLogSpots)
      throw "Critical error: the mesh size has not been well evaluated!";

    nNbNodes = I + 1;
    for(I = 0, J = nNbNodes - 1; I < J; I++, J--)
    {
      dTmp = pdLogSpots[I];
      pdLogSpots[I] = pdLogSpots[J];
      pdLogSpots[J] = dTmp;
    }
  }
  else
  {
    GenRaf( m_dSMin, -8, dDxDiff, m_dDxMax, 1.1, 0, pdLogSpots + nNbNodes - 1, J);
    nNbNodes += J - 1;
    GenRaf( -8, vdLeftSpot, dDxDiff, dDxDiff, 1.1, 0, pdLogSpots + nNbNodes - 1, J);
    nNbNodes += J - 1;
  }
  
//  cout << "Nb points a gauche : " << nNbNodes << endl;

//  if(iAmericanEuropean)
  {
    if(m_optiontype == pricing::Option_Put)
    {
      if(dLogSpot > 4 * dDxDiff)
      {
      GenRaf( vdLeftSpot, -dDxDiff, dDxDiff, dDxDiff, dRhoSpot, 0, pdLogSpots + nNbNodes - 1, J);
      nNbNodes += J - 1;
      pdLogSpots[nNbNodes++] = 0;
      pdLogSpots[nNbNodes++] = dDxDiff;
      GenRaf(dDxDiff, vdRightSpot, dDxDiff, dDxDiff, dRhoSpot, 1, pdLogSpots + nNbNodes - 1, J);
      nNbNodes += J - 1;
      }
      else
      {
      // on garantit la continute de resultat en fonction de spot initial
      dTmp = (dLogSpot < - 4.2 * dDxDiff) ? dLogSpot : -4.2 * dDxDiff;
      GenRaf( vdLeftSpot, dTmp - dDxPetit, dDxPetit, dDxDiff, dRhoSpot, 0, pdLogSpots + nNbNodes - 1, J);
      nNbNodes += J - 1;
      pdLogSpots[nNbNodes++] = dTmp;
      pdLogSpots[nNbNodes++] = dTmp + dDxPetit;
      dRhoTmp = (dDxPetit - (-dDxDiff - (dTmp - dDxPetit))) / (dDxPetit * (dDxDiff / dDxPetit) - (-dDxDiff - (dTmp - dDxPetit)));
      if(dRhoTmp < dRhoSpot)
        dRhoTmp = dRhoSpot;
      GenRaf(dTmp + dDxPetit, -dDxDiff, dDxPetit, dDxDiff, dRhoTmp, 1, pdLogSpots + nNbNodes - 1, J);
      nNbNodes += J - 1;
      pdLogSpots[nNbNodes++] = 0;
      pdLogSpots[nNbNodes++] = dDxDiff;
      GenRaf(dDxDiff, vdRightSpot, dDxDiff, dDxDiff, dRhoSpot, 1, pdLogSpots + nNbNodes - 1, J);
      nNbNodes += J - 1;
      }
    }
    else
    {
      // on garantit la continute de resultat en fonction de spot initial
      //dTmp = (dLogSpot > 4.2 * dDxDiff) ? dLogSpot : 4.2 * dDxDiff;
      GenRaf( vdLeftSpot, -dDxDiff, dDxDiff, dDxDiff, dRhoSpot, 0, pdLogSpots + nNbNodes - 1, J);
      nNbNodes += J - 1;
      pdLogSpots[nNbNodes++] = 0;
      pdLogSpots[nNbNodes++] = dDxDiff;
      GenRaf( dDxDiff, vdRightSpot, dDxDiff, dDxDiff, dRhoTmp, 0, pdLogSpots + nNbNodes - 1, J);
      nNbNodes += J - 1;
    }
  }
  /*
  else
  {
    iToto = 3;
    pdS[0] = vdLeftSpot;
    pdS[1] = 0;
    pdS[2] = vdRightSpot;

    // OPTION: on peut mettre en commentaire ceci, 
    //         mais il me semble qu'il donne pas beaucoup plus de precision
    dDxPetit = dDxDiff;

    GenMaillC(2, dDxPetit, dDxDiff, dRhoSpot, pdS, iToto,
      pdLogSpots + nNbNodes - 1, J);
    nNbNodes += J - 1;
  }
  */

  if(viShift)
  {
    dRhoTmp = 2;
    double dDxTmp = dDxDiff;
    I = nNbNodes - 1;
    while(pdLogSpots[I] < dStrikeRightG)
    {
      pdLogSpots[I + 1] = pdLogSpots[I] + dDxTmp;
      I++;
      if(dDxTmp < 2. * dDxDiff)
      {
        dDxTmp *= dRhoSpot;
        if(dDxTmp > 2. * dDxDiff)
          dDxTmp = 2. * dDxDiff;
      }
    }
    while(pdLogSpots[I] < m_dSMax - m_dDxMax * 0.5)
    {
      pdLogSpots[I + 1] = pdLogSpots[I] + dDxTmp;
      I++;
      if(dDxTmp < m_dDxMax)
      {
        dDxTmp *= dRhoTmp;
        if(dDxTmp > m_dDxMax)
          dDxTmp = m_dDxMax;
      }
    }
    nNbNodes = I + 1;
  }
  else
  {
    GenRaf( -8, m_dSMin, dDxDiff, m_dDxMax, 1.3, 1, pdLogSpots + nNbNodes - 1, J);
    nNbNodes += J - 1;
    GenRaf( vdLeftSpot, -8, dDxDiff, dDxDiff, 2, 1, pdLogSpots + nNbNodes - 1, J);
    nNbNodes += J - 1;
  }

  vecMeshLogS.resize(nNbNodes);
  for(size_t nI = 0; nI < nNbNodes; nI++)
    vecMeshLogS[nI] = pdLogSpots[nI];

  delete [] pdLogSpots;
}
