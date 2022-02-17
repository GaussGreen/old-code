/////////////////////////////////////////////////////////////////////////////
// Name:        numeric/mesh/cbmesh.cpp
// Purpose:     Some useful function for the meshes in the CB case.
// Author:      Nabil
// Created:     2003/12/01
// RCS-ID:      $Id: cbmesh.cpp,v 1.3 2004/10/05 09:13:46 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////
#include "ito33/beforestd.h"
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/numeric/mesh/cbmesh.h"

namespace ito33
{

namespace numeric
{

namespace mesh
{

int WidenGrid(double *_pdOld, int _iNbOld,
               double _dRhoLeft, double _dDeltaMaxLeft, double _dMinLeft,
               double _dRhoRight, double _dDeltaMaxRight, double _dMaxRight,
               double **_ppdNew, int &_iNbNew, int &_iIndiceInitial)
{
  if(_iNbOld < 2)
    return 1;

  int
    piNbRef[2],
    I,
    iTmp;
  double 
    dTmp,
    dDelta,
    dDeltaLeft = _pdOld[1] - _pdOld[0],
    dDeltaRight = _pdOld[_iNbOld - 1] - _pdOld[_iNbOld - 2],
    dDeltaMaxLeft = (_dDeltaMaxLeft > dDeltaLeft) ? _dDeltaMaxLeft : dDeltaLeft,
    dDeltaMaxRight = (_dDeltaMaxRight > dDeltaRight) ? _dDeltaMaxRight : dDeltaRight,
    dRhoLeft = _dRhoLeft,
    dRhoRight = _dRhoRight,
    *pdX
    ;

  // calculate first the dimension of the new grid
  I = _iNbOld;

  piNbRef[0] = 0;
  if((dTmp = _pdOld[0] - _dMinLeft) > 0)
  {
    I += int(dTmp / dDeltaMaxLeft) + 1;
    if(_dDeltaMaxLeft > dDeltaLeft)
    {
      piNbRef[0] = int(log(_dDeltaMaxLeft / dDeltaLeft) / log(_dRhoLeft)) + 1;
      I += piNbRef[0];
      // we change dRhoLeft so that after dDelatSmall run  to dDeltaBig after iTmp times
      dRhoLeft = pow( (_dDeltaMaxLeft / dDeltaLeft), 1. / piNbRef[0]);
    }
  }

  piNbRef[1] = 0;
  if((dTmp = _dMaxRight - _pdOld[_iNbOld - 1]) > 0)
  {
    I += int(dTmp / dDeltaMaxRight) + 1;
    if(_dDeltaMaxRight > dDeltaRight)
    {
      piNbRef[1] = static_cast<int>(ceil(log(_dDeltaMaxRight / dDeltaRight) / log(_dRhoRight)));
      I += piNbRef[1];
      // we change dRhoRight so that after dDelatSmall run  to dDeltaBig after iTmp times
      dRhoRight = pow( (_dDeltaMaxRight / dDeltaRight), 1. / piNbRef[1]);
    }
  }

  pdX = new double [I];
  
  // add points on the left
  pdX[0] = _pdOld[0];

  dDelta = dDeltaLeft;
  for(I = 1; I <= piNbRef[0] && pdX[I - 1] > _dMinLeft; I++)
    {
    dDelta *= dRhoLeft;
    pdX[I] = pdX[I - 1] - dDelta;
    }
  for(; pdX[I - 1] > _dMinLeft; I++)
    pdX[I] = pdX[I - 1] - dDelta;

  _iNbNew = I;

  // exchange the points
  for(I = 0, iTmp = _iNbNew - 1; iTmp > I; I++, iTmp--)
    {
    dTmp       = pdX[I];
    pdX[I]     = pdX[iTmp];
    pdX[iTmp] = dTmp;
    }

  if(_iIndiceInitial >= 0)
    _iIndiceInitial += _iNbNew - 1;

  // add the original points
  for(I = 1; I < _iNbOld; )
    pdX[_iNbNew++] = _pdOld[I++];

  // add points on the right
  dDelta = dDeltaRight;
  for(I = 0; I < piNbRef[1] &&  pdX[_iNbNew - 1] < _dMaxRight; I++, _iNbNew++)
    {
    dDelta *= dRhoRight;
    if(dDelta > _dDeltaMaxRight)
      dDelta = _dDeltaMaxRight;
    pdX[_iNbNew] = pdX[_iNbNew - 1] + dDelta;
    }
  for(; pdX[_iNbNew - 1] < _dMaxRight; _iNbNew++)
    pdX[_iNbNew] = pdX[_iNbNew - 1] + dDelta;

  *_ppdNew = pdX;

  return 0;
}


int WidenGrid(double *_pdOld, int _iNbOld,
               double _dRhoLeft, double _dDeltaMaxLeft, double _dMinLeft,
               double _dRhoRight, double _dDeltaMaxRight, double _dMaxRight,
               double **_ppdNew, int &_iNbNew, int &_iIndiceInitial, double _dXRight)
{
  if(_iNbOld < 2)
    return 1;

  int
    piNbRef[2],
    I,
    iTmp;
  double 
    dTmp,
    dDelta,
    dDeltaLeft = _pdOld[1] - _pdOld[0],
    dDeltaRight = _pdOld[_iNbOld - 1] - _pdOld[_iNbOld - 2],
    dDeltaMaxLeft = (_dDeltaMaxLeft > dDeltaLeft) ? _dDeltaMaxLeft : dDeltaLeft,
    dDeltaMaxRight = (_dDeltaMaxRight > dDeltaRight) ? _dDeltaMaxRight : dDeltaRight,
    dRhoLeft = _dRhoLeft,
    dRhoRight = _dRhoRight,
    *pdX
    ;

  // calculate first the dimension of the new grid
  I = _iNbOld;

  piNbRef[0] = 0;
  if((dTmp = _pdOld[0] - _dMinLeft) > 0)
  {
    I += int(dTmp / dDeltaMaxLeft) + 1;
    if(_dDeltaMaxLeft > dDeltaLeft)
    {
      piNbRef[0] = int(log(_dDeltaMaxLeft / dDeltaLeft) / log(_dRhoLeft)) + 1;
      I += piNbRef[0];
      // we change dRhoLeft so that after dDelatSmall run  to dDeltaBig after iTmp times
      dRhoLeft = pow( (_dDeltaMaxLeft / dDeltaLeft), 1. / piNbRef[0]);
    }
  }

  piNbRef[1] = 0;
  if((dTmp = _dMaxRight - _pdOld[_iNbOld - 1]) > 0)
  {
    I += int(dTmp / dDeltaMaxRight) + 1;
    if(_dDeltaMaxRight > dDeltaRight)
    {
      piNbRef[1] = int(log(_dDeltaMaxRight / dDeltaRight) / log(_dRhoRight)) + 1;
      I += piNbRef[1];
      // we change dRhoRight so that after dDelatSmall run  to dDeltaBig after iTmp times
      dRhoRight = pow( (_dDeltaMaxRight / dDeltaRight), 1. / piNbRef[1]);
    }
  }

  pdX = new double [I];
  
  // add points on the left
  pdX[0] = _pdOld[0];

  dDelta = dDeltaLeft;
  for(I = 1; I <= piNbRef[0] && pdX[I - 1] > _dMinLeft; I++)
    {
    dDelta *= dRhoLeft;
    pdX[I] = pdX[I - 1] - dDelta;
    }
  for(; pdX[I - 1] > _dMinLeft; I++)
    pdX[I] = pdX[I - 1] - dDelta;

  _iNbNew = I;

  // exchange the points
  for(I = 0, iTmp = _iNbNew - 1; iTmp > I; I++, iTmp--)
    {
    dTmp       = pdX[I];
    pdX[I]     = pdX[iTmp];
    pdX[iTmp] = dTmp;
    }

  if(_iIndiceInitial >= 0)
    _iIndiceInitial += _iNbNew - 1;

  // add the original points
  for(I = 1; I < _iNbOld; )
    pdX[_iNbNew++] = _pdOld[I++];

  // add points on the right
  if(_dXRight < pdX[_iNbNew - 1] || _dXRight > _dMaxRight)
    {
    dDelta = dDeltaRight;
    for(I = 0; I < piNbRef[1] &&  pdX[_iNbNew - 1] < _dMaxRight; I++, _iNbNew++)
      {
      dDelta *= dRhoRight;
      if(dDelta > _dDeltaMaxRight - 1.e-8)
        dDelta = _dDeltaMaxRight;
      pdX[_iNbNew] = pdX[_iNbNew - 1] + dDelta;
      }
    for(; pdX[_iNbNew - 1] < _dMaxRight; _iNbNew++)
      pdX[_iNbNew] = pdX[_iNbNew - 1] + dDelta;
    }
  else
    {
    dDelta = dDeltaRight;

    for(; ;)
      {
      if(_dXRight < pdX[_iNbNew - 1] + dDelta)
        {
        pdX[_iNbNew - 1] = (_dXRight + pdX[_iNbNew - 2]) * 0.5;
        pdX[_iNbNew++] = _dXRight;
        break;
        }/*
      else if(_dXRight < pdX[_iNbNew - 1] + dDelta + dDelta)
        {
        pdX[_iNbNew++] = _dXRight;
        break;
        }*/
      else
        {
        dDelta *= dRhoRight;
        if(dDelta > _dDeltaMaxRight - 1.e-8)
          dDelta = _dDeltaMaxRight;
        pdX[_iNbNew] = pdX[_iNbNew - 1] + dDelta;
        _iNbNew++;
        }
      }

    dDelta = pdX[_iNbNew - 1] - pdX[_iNbNew - 2];
    for(; pdX[_iNbNew - 1] < _dMaxRight; _iNbNew++)
      {
      dDelta *= dRhoRight;
      if(dDelta > _dDeltaMaxRight - 1.e-8)
        dDelta = _dDeltaMaxRight;
      pdX[_iNbNew] = pdX[_iNbNew - 1] + dDelta;
      }
    for(; pdX[_iNbNew - 1] < _dMaxRight; _iNbNew++)
      pdX[_iNbNew] = pdX[_iNbNew - 1] + dDelta;

    }

  *_ppdNew = pdX;

  return 0;
}


} // namespace mesh

} // namespace numeric

} // namespace ito33
