//***********************************************************************************************
// file........: tridiagpre.h
// author......: WANG (ITO33)
// Objectif: Preconditionneur par la partie tridiagonale d'une mmatrice 
//***********************************************************************************************

#ifndef TRIDIAGPRE_H
#define TRIDIAGPRE_H

#include <math.h>
#include "mvvd.h"

#define VECTOR_double CVecteurDouble

class ILUPreconditioner_double {
private:
  const double 
    *mpdA, *mpdB, *mpdC; 
  double 
    *mpdZ, *mpdInvBMoinsAZ;
  int 
    miNbInconnu;

public:

  ILUPreconditioner_double() : miNbInconnu(0){};

  ILUPreconditioner_double(double *_pdA, double *_pdB, double *_pdC, int _iNbInconnu)
                         : mpdA(_pdA), mpdB(_pdB), mpdC(_pdC), miNbInconnu(_iNbInconnu)
  {
    int 
      iIdElem;

    mpdZ = new double [miNbInconnu];
    mpdInvBMoinsAZ = new double [miNbInconnu];

    mpdInvBMoinsAZ[0] = 1. / mpdB[0];
    mpdZ[0] = mpdC[0] * mpdInvBMoinsAZ[0];
    
    for (iIdElem = 1; iIdElem < miNbInconnu - 1; iIdElem++)
    {
      mpdInvBMoinsAZ[iIdElem] = 1. / (mpdB[iIdElem] - mpdA[iIdElem] * mpdZ[iIdElem - 1]);
      mpdZ[iIdElem] = mpdC[iIdElem] * mpdInvBMoinsAZ[iIdElem];
    }
    mpdInvBMoinsAZ[iIdElem] = 1. / (mpdB[iIdElem] - mpdA[iIdElem] * mpdZ[iIdElem - 1]);
  }


  ~ILUPreconditioner_double (void) 
  {
    delete [] mpdZ;
    delete [] mpdInvBMoinsAZ;
  };

  VECTOR_double solve (const VECTOR_double &b) const
  {
    VECTOR_double 
      x(b);  
    int 
      iIdElem;
  
    x[0] = x(0) * mpdInvBMoinsAZ[0];
    for (iIdElem = 1; iIdElem < miNbInconnu; iIdElem++)
      x[iIdElem] = (x[iIdElem] - mpdA[iIdElem] * x[iIdElem - 1]) * mpdInvBMoinsAZ[iIdElem];

    for (iIdElem = miNbInconnu - 2; iIdElem >= 0; iIdElem--)
      x[iIdElem] = x[iIdElem] - mpdZ[iIdElem] * x[iIdElem + 1];

    return x;
  }
};


#endif
