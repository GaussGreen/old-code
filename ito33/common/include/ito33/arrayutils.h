#ifndef ARRAYUTILS_H
#define ARRAYUTILS_H

void SetZero(double *_pdX, int _iDim);

void CopyToArray(double *_pdA, double *_pdB, int _iDim);

void PrintArray(double *_pdA, int _iDim);

double InnerProduct(double *_pdA, double *_pdB, int _iDim);

double InnerProduct(double *_pdA, double *_pdB, double *_pdPoid, int _iDim);

double InnerProduct(double **_ppdA, double **_ppdB, int _iDimA, int _iDimB);

int AddToArray(int _iElem, int *_piArray, int &_iDim);

#endif
