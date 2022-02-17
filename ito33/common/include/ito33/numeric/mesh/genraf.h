/**************************************************************************
 * file      : ito33/numeric/mesh/genraf.h
 * Purpose   : "genraf" mesh generator
 * Author    : ZHANG Yunzhi
 * RCS-ID    : $Id: genraf.h,v 1.5 2004/11/10 16:37:14 afrolov Exp $
 * Copyright : (c) 1999 - 2003 Trilemma LLP all rights reserved
 **************************************************************************/

/**
   @file ito33/numeric/mesh/genraf.h
   @brief functions for adapted mesh generartion
 */

#ifndef _ITO33_NUMERIC_MESH_GENRAF_H_
#define _ITO33_NUMERIC_MESH_GENRAF_H_

#include <cmath>
#include <stddef.h>

namespace ito33
{

namespace numeric
{

namespace mesh
{


inline size_t GenRafSize(double dWidth,
                         double dSmallStep, 
                         double dBigStep, 
                         double dStretch);

/**
   Generates a mesh for given interval.

   @param dA left point of the interval
   @param dB right point of the interval
   @param dDeltaMin the minimum delta
   @param dDeltaMax the maximum delta
   @param dRho multiplier used to grow from min delta to max delta
   @param iRafA refine at left if true, refine at right otherwise
   @param pdX output, array of the generated mesh
   @param iDimX output size of the generated mesh
 */
void GenRaf(double dA, double dB,
            double dDeltaMin, double dDeltaMax,
            double dRho,
            int iRafA,
            double *pdX, int &iDimX);

size_t GenMeshSize(int iNRaf,
                   double dDeltax1,
                   double dDeltax,
                   double dRho,
                   double *pdS,
                   int iDimS);


void GenMesh(int iNraf,
             double dDeltax1,
             double dDeltax,
             double dRho,
             double *pdS,
             int iDimS,
             double *pdGrille,
             int &iDimGrille);


// ----------------------------------------------------------------

size_t GenRafSize(double dWidth, 
                         double dSmallStep, 
                         double dBigStep, 
                         double dStretch)
{
  return 2 * (size_t) 
    (
      ceil( log(dBigStep / dSmallStep) / log(dStretch) + 1) + 
      ceil(dWidth / dBigStep + 1)
    );
}


} // namespace mesh

} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_MESH_GENRAF_H_
