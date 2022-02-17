// file........: ilupre_double.hpp
// author......: ZHANG (ITO33)
//
#ifndef ILUPRE_H
#define ILUPRE_H

#include <math.h>
#include "ito33/mv/mvvd.h"
#include "ito33/mv/cmorse.hpp"
#define VECTOR_double CVecteurDouble

class ILUPreconditioner_double {
private:
  const CMorseStruct *mpMS;
  int miDim;
  int *mpiIndiceDiag;
  double **mppdElem;

public:
  ILUPreconditioner_double():miDim(0){};

   ILUPreconditioner_double(CMorseMatrix &_M)
  :  mpMS(_M.GetMorseStruct())
  {
    int
      i,
      j,
      ji,
      jp,
      ip,
      iNb = 0;
    double
      dMul;

    miDim = mpMS->Max() - mpMS->Min() + 1;
    mpiIndiceDiag = new int [miDim];

    for(i = 0; i < miDim; i++)
      {
      iNb += mpMS->mpiLongLignes[i];
      for(j = 0; mpMS->mppiCols[i][j] != i; j++)
        ;
      mpiIndiceDiag[i] = j;
      }

    mppdElem = new double* [miDim];
    mppdElem[0] = new double [iNb];
    for(i = 0; i < miDim - 1; i++)
      mppdElem[i + 1] = mppdElem[i] + mpMS->mpiLongLignes[i];

    for(i = 0; i < miDim; i++)
      for(j = 0; j < mpMS->mpiLongLignes[i]; j++)
        mppdElem[i][j] = _M.ElementDirect(i,j);

    for(i = 0; i < miDim - 1; i++)
    {
      dMul = 1. / mppdElem[i][mpiIndiceDiag[i]];
      for(j = mpiIndiceDiag[i] + 1; j < mpMS->mpiLongLignes[i]; j++)
        mppdElem[i][j] *= dMul;

      for(j = 0; j < mpiIndiceDiag[i + 1]; j++)
      {
        // le premier element de chaque ligne ne change pas, il est bien initie
        dMul = mppdElem[i + 1][j];  // l(I + 1,J)
        ji = mpMS->mppiCols[i + 1][j]; // J

        // nous calculons l'influence de cet element a ceux qui le suivent
        ip = j + 1;
        // calcule l(I+1,*). puisqu'on utlise u(J,*), il faut que *>J
        for(jp = mpiIndiceDiag[ji] + 1;
          jp < mpMS->mpiLongLignes[ji] && mpMS->mppiCols[ji][jp] <= i+1;
          jp++)
        {
          for(;
              ip <= mpiIndiceDiag[i+1] && mpMS->mppiCols[i+1][ip] < mpMS->mppiCols[ji][jp];
              ip++)
            ;
          if(mpMS->mppiCols[i+1][ip] == mpMS->mppiCols[ji][jp])
            mppdElem[i+1][ip] -= dMul * mppdElem[ji][jp];
        }

        // on utilise u(I+1,*)
        ip = mpiIndiceDiag[i+1] + 1;
        for(; jp < mpMS->mpiLongLignes[ji]; jp++)
        {
          for(;
            ip < mpMS->mpiLongLignes[i+1] && mpMS->mppiCols[i+1][ip] < mpMS->mppiCols[ji][jp];
            ip++)
            ;
          if(ip < mpMS->mpiLongLignes[i+1] && mpMS->mppiCols[i+1][ip] == mpMS->mppiCols[ji][jp])
            mppdElem[i + 1][ip] -= dMul * mppdElem[ji][jp];
        }
      }
    }
  }


  ~ILUPreconditioner_double (void) 
  {
    if(miDim)
    {
      delete [] mppdElem[0];
      delete [] mppdElem;
      delete [] mpiIndiceDiag;
    }
  };

  VECTOR_double solve (const VECTOR_double &b) const
  {
    VECTOR_double x(b);
    double dTmp;
    int i, k;

    i = 0;
    x(i) /= mppdElem[i][mpiIndiceDiag[i]];
    for(i++; i < miDim; i++)
      {
      dTmp = 0;
      for(k = 0; k < mpiIndiceDiag[i]; k++)
        dTmp += mppdElem[i][k] * x(mpMS->mppiCols[i][k]);
      x(i) = (x(i) - dTmp) / mppdElem[i][mpiIndiceDiag[i]];
      }
  /*
    for(i = 0; i < miDim; i++)
      {
      if(mpMS->mppiCols[i][mpiIndiceDiag[i]] != i)
        cout << i << endl;
      t[i] = x[i];
      dTmp = 0;
      for(k = 0; k <= mpiIndiceDiag[i]; k++)
        dTmp += mppdElem[i][k] * x(mpMS->mppiCols[i][k]);
      if((dTmp = fabs(dTmp - b(i))) > 1.e-7)
      {
        dTmp = 0;
        for(k = 0; k < mpiIndiceDiag[i]; k++)
        {
          cout << k << " " << x(mpMS->mppiCols[i][k]) << "\t" << mppdElem[i][k] << endl;
          dTmp += mppdElem[i][k] * x(mpMS->mppiCols[i][k]);
        }
        cout << "...\n";
        cout << (b(i) - dTmp) / mppdElem[i][mpiIndiceDiag[i]] << "\n" << x(i) << endl;
        exit(1);
      }
      }
  */
    for(i = miDim - 1; i >= 0; i--)
      {
      dTmp = 0;
      for(k = mpMS->mpiLongLignes[i] - 1; k > mpiIndiceDiag[i]; k--)
        dTmp += mppdElem[i][k] * x(mpMS->mppiCols[i][k]);
      x(i) = (x(i) - dTmp);
      }

    return x;
  }

  VECTOR_double trans_solve (const VECTOR_double &b) const
  {
    VECTOR_double x(b);
    int i, k;

    for(i = 0; i < miDim; i++)
      {
      for(k = mpiIndiceDiag[i] + 1; k < mpMS->mpiLongLignes[i]; k++)
        x(mpMS->mppiCols[i][k]) -= mppdElem[i][mpMS->mppiCols[i][k]] * x(i);
      }

    for(i = miDim - 1; i >= 0; i--)
      {
      x(i) /= mppdElem[i][mpiIndiceDiag[i]];
      for(k = 0; k < mpiIndiceDiag[i]; k++)
        x(mpMS->mppiCols[i][k]) -= mppdElem[i][mpMS->mppiCols[i][k]] * x(i);
      }

    return x;
  }

};


#endif
