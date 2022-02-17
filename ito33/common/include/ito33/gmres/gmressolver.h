#ifndef GMRESSOLVER_H
#define GMRESSOLVER_H

/**
   @deprecated, use instead ito33/numeric/gmressolver.h
 */

#include "ito33/mv/mvmd.h"
#include "ito33/gmres/gmres.h"

class CGmresSolver
{ 
  int
    miNbUnknown;

  double
    mdTolerance;

  MATRIX_double
    *mpHessenberg;
        
  int
    miRestart, miMaxIteration;

public: 
  
  CGmresSolver(int _iNbUnknown, int _iRestart = 35, int _iMaxIteration = 1000, double _dTolerance = 1.e-15)
  {
    miNbUnknown = _iNbUnknown;
    miRestart = _iRestart;
    miMaxIteration = _iMaxIteration;
    mdTolerance = _dTolerance; 
    
    // memory allocation for the hessenberg matrix used in Gmres
    mpHessenberg = new MATRIX_double(miRestart + 1, miRestart);
  }

  ~CGmresSolver()
  {
    delete mpHessenberg;
  }
 
  template <class Matrix>

  void solve(Matrix *_pMatrix, double *_pdB, double *_pdX)
  {
    CVecteurDouble
      *pB, *pX;
    int
      iRestart = miRestart,
      iMaxIteration = miMaxIteration;
    double
      dTolerance = mdTolerance;

    pB = new CVecteurDouble(_pdB, miNbUnknown, MV_Vector_::ref);
    pX = new CVecteurDouble(_pdX, miNbUnknown, MV_Vector_::ref);

    GMRES(*_pMatrix, *pX, *pB, *(_pMatrix->getPreconditioner()), *mpHessenberg, iRestart, iMaxIteration, dTolerance);

    delete pB;
    delete pX;
  }
};

#endif
