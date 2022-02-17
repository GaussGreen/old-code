#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/arrayutils.h"

using namespace std;

void SetZero(double *_pdX, int _iDim)
{
  for (int iIdx = 0; iIdx < _iDim; iIdx++)
    _pdX[iIdx] = 0;
}

void CopyToArray(double *_pdA, double *_pdB, int _iDim)
{
  for (int iIdx = 0; iIdx < _iDim; iIdx++)
    _pdB[iIdx] = _pdA[iIdx];
}

void PrintArray(double *_pdA, int _iDim)
{
  for (int iIdx = 0; iIdx < _iDim; iIdx++)
    cout << _pdA[iIdx] << endl;
}

double InnerProduct(double *_pdA, double *_pdB, int _iDim)
{
  double
    dTmp = 0;

  for (int iIdx = 0; iIdx < _iDim; iIdx++)
    dTmp += _pdA[iIdx] * _pdB[iIdx];

  return dTmp;
}

double InnerProduct(double *_pdA, double *_pdB, double *_pdPoid, int _iDim)
{
  double
    dTmp = 0;

  for (int iIdx = 0; iIdx < _iDim; iIdx++)
    dTmp += _pdPoid[iIdx] * _pdA[iIdx] * _pdB[iIdx];

  return dTmp;
}

double InnerProduct(double **_ppdA, double **_ppdB, int _iDimA, int _iDimB)
{
  double
    dTmp = 0;

  for (int iIdx = 0; iIdx < _iDimA; iIdx++)
    dTmp += InnerProduct(_ppdA[iIdx], _ppdB[iIdx], _iDimB);

  return dTmp;
}

int AddToArray(int _iElem, int *_piArray, int &_iDim)
{
  int 
    iIdx;

  for (iIdx = 0; iIdx < _iDim; iIdx++)
    if (_iElem == _piArray[iIdx])
      return 1;

  _piArray[_iDim++] = _iElem;

  return 0;
}


