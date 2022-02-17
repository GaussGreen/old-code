/**************************************************************************
 * file      : src/numeric/mesh/GenRaf.cpp
 * Purpose   : "genraf" mesh generator
 * Author    : ZHANG Yunzhi
 * RCS-ID    : $Id: genraf.cpp,v 1.4 2004/10/11 14:14:16 wang Exp $
 * Copyright : (c) 1999 - 2003 Trilemma LLP all rights reserved
 **************************************************************************/

#include "ito33/beforestd.h"
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/numeric/mesh/genraf.h"

namespace ito33
{

namespace numeric
{

namespace mesh
{

size_t GenMeshSize(int iNRaf, double dDeltax1, double dDeltax, double dRho, 
                   double *pdS, int iDimS)
{
  int
    I;

  // On calcule le min de DeltaS
  double
    dDeltaS = pdS[iDimS - 1] - pdS[0],
    dLong = dDeltaS;

  for(I = 1; I < iDimS; I++) 
    if(pdS[I] - pdS[I - 1] < dDeltaS - 1.0e-10)
      dDeltaS = pdS[I] - pdS[I - 1];
    
  if(dDeltax1 > dDeltaS / (2.0 * iNRaf) + 1.0e-10)
    dDeltax1 = dDeltaS / (2.0 * iNRaf);

  return 2 * (size_t) ceil (2.0 * iDimS * (log(dDeltax / dDeltax1) / log(dRho))
                            + ceil(dLong / dDeltax) + 2
                            + iDimS);
}


void GenRaf(double dA, double dB, double dDeltax1, double dDeltax, double dRho,
            int iRafA,
            double *pdX, int &iDimX) 
{
  int
    I,
    J;

  double
    dTemp,
    dC = dDeltax / dDeltax1;
  const double
    dEpsilon = 1.0e-6;

  double
    dXfinal;
  int
    iSgn = 1;
  if(iRafA)
    dXfinal = dA;
  else
  {
    iSgn = - 1;
    dXfinal = dB;
  }
    
  double
    dDeltaxI = dDeltax1;
  int
    iNpt = 1;

  do
  {
    dXfinal += iSgn * dDeltaxI;
    iNpt++;
    dDeltaxI *= dRho;
    if(dDeltaxI > dDeltax)
      dDeltaxI = dDeltax;
  }
  while((dXfinal - dA) * (dXfinal - dB) < -(dB - dA) * dDeltax1 * dEpsilon);
    
  if(iRafA) { 

    dDeltax1 = dDeltax1 * (dB - dA) / (dXfinal - dA);
  }
  else {
    
    dDeltax1 = dDeltax1 * (dB - dA) / (dB - dXfinal);
    }
  dDeltax = dC * dDeltax1;

  iDimX = iNpt;
  pdX[0] = iRafA ? dA : dB;
  dDeltaxI = dDeltax1;

  for(I = 1; I < iNpt - 1; I++)
  {
    pdX[I] = pdX[I - 1] + iSgn * dDeltaxI;
    dDeltaxI *= dRho;
    if(dDeltaxI > dDeltax)
      dDeltaxI = dDeltax;
  }
   
  if(iRafA)
    pdX[iNpt - 1] = dB;
  else
  { 
    pdX[iNpt - 1] = dA;

    for (I = 0, J = iNpt - 1; I < J; I++, J--)
    {
      dTemp = pdX[I];
      pdX[I] = pdX[J];
      pdX[J] = dTemp;
    }
  }
}


void GenMesh(int iNraf, double dDeltax1, double dDeltax, double dRho, 
             double *pdS, int iDimS, 
             double *pdGrille, int &iDimGrille)
{
  int
    I,
    iLastG = 0,
    iUsed = iDimGrille;
  double
    dDxLoc = dDeltax1;

  if(iDimS == 2)
  {
    GenRaf(pdS[0], pdS[1], dDeltax, dDeltax, dRho, 0, 
           pdGrille + iLastG, iDimGrille);
    return;
  }

  if(dDxLoc > (pdS[1] - pdS[0]) / iNraf)
    dDxLoc = (pdS[1] - pdS[0]) / iNraf;
  if(dDxLoc > (pdS[2] - pdS[1]) / (2.0 * iNraf))
    dDxLoc = (pdS[2] - pdS[1]) / (2.0 * iNraf);
  
  GenRaf(pdS[0], pdS[1], dDxLoc, dDeltax, dRho, 0, pdGrille + iLastG, iUsed);
  iLastG += iUsed - 1;


  for (I = 1; I < iDimS - 2; I++)
  {

    iUsed = iDimGrille - iLastG;
    GenRaf(pdS[I], (pdS[I + 1] + pdS[I]) * 0.5, dDxLoc, dDeltax, dRho, 1, pdGrille +
      iLastG, iUsed);
  
    iLastG += iUsed - 1;

    dDxLoc = dDeltax1;
    if(dDxLoc > (pdS[I + 1] - pdS[I]) / (2.0 * iNraf))
      dDxLoc = (pdS[I + 1] - pdS[I]) / (2.0 * iNraf);
    if(I + 2 == iDimS-1)
      if(dDxLoc > (pdS[I + 2] - pdS[I + 1]) / iNraf)
        dDxLoc = (pdS[I + 2] - pdS[I + 1]) / iNraf;
      else
        if(dDxLoc > (pdS[I + 2] - pdS[I + 1]) / (2.0 * iNraf))
          dDxLoc = (pdS[I + 2] - pdS[I + 1]) / (2.0 * iNraf);
    
    iUsed = iDimGrille-iLastG;
    GenRaf((pdS[I + 1] + pdS[I]) * 0.5, pdS[I + 1], dDxLoc, dDeltax, dRho, 0,
           pdGrille + iLastG, iUsed);
    iLastG += iUsed - 1;
  }

  iUsed = iDimGrille - iLastG;
  GenRaf(pdS[iDimS - 2], pdS[iDimS - 1], dDxLoc, dDeltax, dRho, 1, 
         pdGrille + iLastG, iUsed);
  iLastG += iUsed;

  iDimGrille = iLastG;
}

} // namespace mesh

} // namespace numeric

} // namespace ito33
