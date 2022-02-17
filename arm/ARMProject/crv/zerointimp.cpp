/*
 * $Log: zerointimp.cpp,v $
 * Revision 1.13  2004/01/21 14:31:46  jpriaudel
 * modif dans le bump
 *
 * Revision 1.12  2003/08/20 18:32:27  jpriaudel
 * memory leaked dans shift curve
 *
 * Revision 1.11  2003/07/02 10:53:58  arm
 * Formatting
 *
 * Revision 1.10  2003/07/01 20:06:52  jpriaudel
 * modif for int i
 *
 * Revision 1.9  2003/07/01 18:56:28  arm
 * correction after for : int i : added
 *
 * Revision 1.8  2002/10/11 08:20:09  mab
 * Improvements
 *
 * Revision 1.7  2002/08/08 10:24:59  mab
 * formattage
 *
 * Revision 1.6  2002/05/30 13:36:30  mab
 * in char Terms[ARM_NB_TERMS][6] : 6 replaced by 12
 *
 * Revision 1.5  2001/06/21 12:12:25  smysona
 * Nouveau constructeur (inout ZCLI)
 *
 * Revision 1.4  2001/06/08 13:10:57  abizid
 * Nettoyage memory leaks
 *
 * Revision 1.3  2000/07/26 13:30:02  mab
 * debug construction de l'objet
 *
 * Revision 1.2  2000/06/26 16:57:02  mab
 * Modif constructeur "Swap"
 *
 * Revision 1.1  2000/06/22 15:41:51  mab
 * Initial revision
 *
 * Revision 1.3  1999/03/02 15:14:44  ypilchen
 * Rajout du comment. RCS
 *
 */


/*----------------------------------------------------------------------------*
 
    zerointerp.cpp
 
    This file implements the ARM_ZeroInterpolation class, a class for 
         computing a ARM_ZeroCurve implicit interpolation.

    Between data, we compute a regularized method in order 
    to have more precision on points. 
    And cubic splines are used between these points

*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "util.h"
#include "zerointimp.h"
#include "currency.h"
#include "expt.h"
#include "newton.h"




/*---------------------------------------------------------------------------*/


ARM_ZeroInterpolation::ARM_ZeroInterpolation(ARM_Date& asOf, 
                                             ARM_Vector* yearTerms,
                                             ARM_Vector* zeroYields,
                                             int compoundMeth, int lambda, 
                                             int nbPoints):ARM_ZeroCurve(asOf)
{
    // Check variables

    SetName(ARM_ZERO_INTERPOLATION);

    if ( yearTerms->GetSize() < K_MIN_NUM_YEAR_TERMS )
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
           "The input year terms and yields must have a minimum of 3 elements");
    }

    if ( yearTerms->GetSize() != zeroYields->GetSize() )
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
           "The input year terms and yields must have the same size");
    }

    // Set initial variables

    int size = yearTerms->GetSize();

    itsLambda = lambda;

    itsNbPoints = nbPoints;

    Init();

    // mise en place de la procedure d'interpolation

    ARM_Vector* newterms = new ARM_Vector(itsNbPoints);

    ARM_Vector* newyields = new ARM_Vector(itsNbPoints);

    /* construction de l'interpolation...*/

    NonParamCurve(yearTerms, zeroYields, itsLambda,
                  itsNbPoints, newterms, newyields);


    // setting des parametres pour les splines cubiques

    ARM_Vector* terms = new ARM_Vector(newterms, itsNbPoints+2, 0,
                                       itsNbPoints-1, 1);

    terms->Elt(0) = 0.0;
    terms->Elt(itsNbPoints+1) = 1000.0;

    ARM_Vector* yields = new ARM_Vector(newyields, itsNbPoints+2,
                                        0, itsNbPoints-1, 1);

    yields->Elt(0) = newyields->Elt(0)-(newyields->Elt(1)
                      -newyields->Elt(0))*newterms->Elt(0)
                      /(newterms->Elt(1)-newterms->Elt(0));

    yields->Elt(itsNbPoints+1) = newyields->Elt(itsNbPoints-1);

    SetYearTerms(terms);

    SetZeroRates(yields);

    delete newterms;
    delete newyields;

    int Size = GetYearTerms()->GetSize();

    ARM_Vector* BPShifts = new ARM_Vector(Size, 0.0);

    for (int i = 0; i < Size; i++)
    {
        if ( GetYearTerms()->Elt(i) >= GetBucketStartPeriod()
                && GetYearTerms()->Elt(i) < GetBucketEndPeriod())
        {
           BPShifts->Elt(i) = GetBPShift()*0.01; // using yield in percent
        }
    }

    SetBPShifts(BPShifts);

    itsCompoundMeth = compoundMeth;

    GenerateDateTerms();

    GenerateDiscountFactors(compoundMeth);

    ComputePolynomCoef();
}


ARM_ZeroInterpolation::ARM_ZeroInterpolation(ARM_ZeroLInterpol* inCurve,
                                             int lambda, 
                                             int nbPoints):
                                          ARM_ZeroCurve(inCurve->GetAsOfDate())
{
    ARM_Date asOf = inCurve->GetAsOfDate();

    ARM_Vector* yearTerms = inCurve->GetYearTerms();
    ARM_Vector* zeroYields = inCurve->GetZeroRates();

    int compoundMeth = inCurve->GetCompoundMeth();
    
    // Check variables

    SetName(ARM_ZERO_INTERPOLATION);

    if ( yearTerms->GetSize() < K_MIN_NUM_YEAR_TERMS )
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
        "The input year terms and yields must have a minimum of 3 elements");
    }

    if ( yearTerms->GetSize() != zeroYields->GetSize() )
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
        "The input year terms and yields must have the same size");
    }

    // Set initial variables

    int size = yearTerms->GetSize();

    itsLambda = lambda;

    itsNbPoints = nbPoints;

    Init();

    // mise en place de la procedure d'interpolation

    ARM_Vector* newterms = new ARM_Vector(itsNbPoints);

    ARM_Vector* newyields = new ARM_Vector(itsNbPoints);

    /* construction de l'interpolation...*/

    NonParamCurve(yearTerms, zeroYields, itsLambda,
                  itsNbPoints, newterms, newyields);


    // setting des parametres pour les splines cubiques

    ARM_Vector* terms = new ARM_Vector(newterms, itsNbPoints+2, 0,
                                       itsNbPoints - 1, 1);

    terms->Elt(0) = 0.0;
    terms->Elt(itsNbPoints+1) = 1000.0;

    ARM_Vector* yields = new ARM_Vector(newyields, itsNbPoints+2, 0,
                                        itsNbPoints - 1, 1);

    yields->Elt(0) = newyields->Elt(0) - (newyields->Elt(1) - newyields->Elt(0))
                     * newterms->Elt(0)/(newterms->Elt(1) - newterms->Elt(0));

    yields->Elt(itsNbPoints+1) = newyields->Elt(itsNbPoints-1);

    SetYearTerms(terms);

    SetZeroRates(yields);

    delete newterms;
    delete newyields;

    int Size = GetYearTerms()->GetSize();

    ARM_Vector* BPShifts = new ARM_Vector(Size, 0.0);

    for (int i = 0; i < Size; i++)
    {
        if ( GetYearTerms()->Elt(i) >= GetBucketStartPeriod()
                && GetYearTerms()->Elt(i) < GetBucketEndPeriod())
        {
           BPShifts->Elt(i) = GetBPShift()*0.01; // using yield in percent
        }
    }

    SetBPShifts(BPShifts);

    itsCompoundMeth = compoundMeth;

    GenerateDateTerms();

    GenerateDiscountFactors(compoundMeth);

    ComputePolynomCoef();
}



ARM_ZeroInterpolation::ARM_ZeroInterpolation(ARM_Date& asOf,
                                             ARM_CRV_TERMS& Terms,
                                             ARM_Vector* mktData, int MMVsFut,
                                             int SwapVsFut, int raw,
                                             int Cont_Lin, ARM_Currency* ccy,
                                             int lambda, int nbPoints) 
                     : ARM_ZeroCurve(asOf, Terms, mktData,
                                                     MMVsFut, SwapVsFut, raw,
                                                     Cont_Lin, ccy)
{
    // Check variables
    itsLambda = lambda;

    itsNbPoints = nbPoints;

    Init();

    SetName(ARM_ZERO_INTERPOLATION);

    // Set variables
    ARM_Vector *yearTerms = GetYearTerms();
    ARM_Vector *zeroYields = GetZeroRates();

    int size = GetYearTerms()->GetSize();

    // mise en place de la procedure d'interpolation

    ARM_Vector* newterms = new ARM_Vector(itsNbPoints);

    ARM_Vector* newyields = new ARM_Vector(itsNbPoints);

    /* construction de l'interpolation...*/

    NonParamCurve(yearTerms, zeroYields, itsLambda, itsNbPoints, newterms,
                  newyields);

    // setting des parametres pour les splines cubiques

    ARM_Vector* terms = new ARM_Vector(newterms, itsNbPoints+2, 0,
                                       itsNbPoints - 1, 1);

    terms->Elt(0) = 0.0;
    terms->Elt(itsNbPoints+1) = 1000.0;

    ARM_Vector* yields = new ARM_Vector(newyields, itsNbPoints+2, 0,
                                        itsNbPoints - 1, 1);

    yields->Elt(0) = newyields->Elt(0)-(newyields->Elt(1)-newyields->Elt(0))
                     *newterms->Elt(0)/(newterms->Elt(1)-newterms->Elt(0));

    yields->Elt(itsNbPoints+1) = newyields->Elt(itsNbPoints-1);

    SetYearTerms(terms);

    SetZeroRates(yields);

    delete newterms;
    delete newyields;

    int Size = GetYearTerms()->GetSize();
 
    ARM_Vector* BPShifts = new ARM_Vector(Size, 0.0);
 
    for (int i = 0; i < Size; i++)
    {
        if ( GetYearTerms()->Elt(i) >= GetBucketStartPeriod()
                && GetYearTerms()->Elt(i) < GetBucketEndPeriod())
        {
            BPShifts->Elt(i) = GetBPShift()*0.01; // using yield in percent
        }
    }

    SetBPShifts(BPShifts);

    itsCompoundMeth = K_COMP_CONT;

    GenerateDateTerms();

    GenerateDiscountFactors(itsCompoundMeth);

    ComputePolynomCoef();
}



ARM_ZeroInterpolation::ARM_ZeroInterpolation(ARM_CRV_TERMS& Terms,
                                             ARM_Date& asOf,
                                             ARM_Vector* mktData, int MMVsFut,
                                             int SwapVsFut, int raw,
                                             int Cont_Lin, ARM_Currency* ccy,
                                             int lambda, int nbPoints) 
                         : ARM_ZeroCurve(Terms, asOf, mktData,
                                         MMVsFut, SwapVsFut,
                                         raw, Cont_Lin, ccy)
{
    // Check variables
    itsLambda = lambda;

    itsNbPoints = nbPoints;

    Init();
 
    SetName(ARM_ZERO_INTERPOLATION);
 
    // Set variables
    ARM_Vector* yearTerms = GetYearTerms();
    ARM_Vector* zeroYields = GetZeroRates();
 
    int size = GetYearTerms()->GetSize();

    // mise en place de la procedure d'interpolation

    ARM_Vector* newterms = new ARM_Vector(itsNbPoints);

    ARM_Vector* newyields = new ARM_Vector(itsNbPoints);

    /* construction de l'interpolation...*/

    NonParamCurve(yearTerms, zeroYields, itsLambda, itsNbPoints, newterms,
                  newyields);

    // setting des parametres pour les splines cubiques

    ARM_Vector* terms = new ARM_Vector(newterms, itsNbPoints+2, 0,
                                       itsNbPoints-1, 1);

    terms->Elt(0) = 0.0;
    terms->Elt(itsNbPoints+1) = 1000.0;

    ARM_Vector* yields = new ARM_Vector(newyields, itsNbPoints+2, 0,
                                        itsNbPoints-1, 1);

    yields->Elt(0) = newyields->Elt(0)-(newyields->Elt(1)-newyields->Elt(0))
                     *newterms->Elt(0)/(newterms->Elt(1)-newterms->Elt(0));

    yields->Elt(itsNbPoints+1) = newyields->Elt(itsNbPoints-1);

    SetYearTerms(terms);

    SetZeroRates(yields);

    delete newterms;
    delete newyields;

    int Size = GetYearTerms()->GetSize();
 
    ARM_Vector* BPShifts = new ARM_Vector(Size, 0.0);
 
    for (int i = 0; i < Size; i++)
    {
        if ( GetYearTerms()->Elt(i) >= GetBucketStartPeriod()
                && GetYearTerms()->Elt(i) < GetBucketEndPeriod())
        {
            BPShifts->Elt(i) = GetBPShift()*0.01; // using yield in percent
        }
    }
 
    SetBPShifts(BPShifts);
 
    itsCompoundMeth = K_COMP_CONT;

    GenerateDateTerms();

    GenerateDiscountFactors(itsCompoundMeth);

    ComputePolynomCoef();
}



ARM_ZeroInterpolation::ARM_ZeroInterpolation(ARM_Date& asOf,
                                             ARM_CRV_TERMS& Terms, 
                                             ARM_Vector* mktData,
                                             ARM_Container* bonds, 
                                             ARM_Vector* Yields, int MMVsFut,
                                             ARM_Currency* ccy,
                                             int lambda, int nbPoints)
                    : ARM_ZeroCurve(asOf, Terms, mktData,
                                    bonds, Yields,
                                    MMVsFut, ccy)
{
    // Check variables
    itsLambda = lambda;

    itsNbPoints = nbPoints;

    Init();

    SetName(ARM_ZERO_INTERPOLATION);

    // Set variables
    ARM_Vector* yearTerms = GetYearTerms();
    ARM_Vector* zeroYields = GetZeroRates();

    int size = GetYearTerms()->GetSize();

    // mise en place de la procedure d'interpolation

    ARM_Vector* newterms = new ARM_Vector(itsNbPoints);

    ARM_Vector* newyields = new ARM_Vector(itsNbPoints);

    /* construction de l'interpolation...*/

    NonParamCurve(yearTerms, zeroYields, itsLambda, itsNbPoints, newterms,
                  newyields);

    // setting des parametres pour les splines cubiques

    ARM_Vector* terms = new ARM_Vector(newterms, itsNbPoints+2, 0,
                                       itsNbPoints-1, 1);

    terms->Elt(0) = 0.0;
    terms->Elt(itsNbPoints+1) = 1000.0;

    ARM_Vector* yields = new ARM_Vector(newyields, itsNbPoints+2, 0,
                                        itsNbPoints-1, 1);

    yields->Elt(0) = newyields->Elt(0)-(newyields->Elt(1)-newyields->Elt(0))
                     *newterms->Elt(0)/(newterms->Elt(1)-newterms->Elt(0));

    yields->Elt(itsNbPoints+1) = newyields->Elt(itsNbPoints-1);

    SetYearTerms(terms);

    SetZeroRates(yields);

    delete newterms;
    delete newyields;

    int Size = GetYearTerms()->GetSize();
 
    ARM_Vector* BPShifts = new ARM_Vector(Size, 0.0);
 
    for (int i = 0; i < Size; i++)
    {
        if ( GetYearTerms()->Elt(i) >= GetBucketStartPeriod()
             && GetYearTerms()->Elt(i) < GetBucketEndPeriod())
        {
           BPShifts->Elt(i) = GetBPShift()*0.01; // using yield in percent
        }
    }

    SetBPShifts(BPShifts);

    itsCompoundMeth = K_COMP_CONT;

    GenerateDateTerms();

    GenerateDiscountFactors(itsCompoundMeth);

    ComputePolynomCoef();
}



ARM_ZeroInterpolation::ARM_ZeroInterpolation(const ARM_ZeroInterpolation& zeroInterpolation)
                      : ARM_ZeroCurve(zeroInterpolation)
{
    Init();

    SetName(ARM_ZERO_INTERPOLATION);

    BitwiseCopy(&zeroInterpolation);
}



ARM_ZeroInterpolation& ARM_ZeroInterpolation::operator
                       = (const ARM_ZeroInterpolation&
                                   zeroInterpolation)
{
    (*this).ARM_ZeroCurve::operator = (zeroInterpolation);

    BitwiseCopy(&zeroInterpolation);

    return(*this);
}



/* fonction de lissage non parametrique                             */
/* en fait, on implemente (-grad_0) car la solution est plus simple */
/* la fonction de regularisation n'apparait que dans la hessienne   */
void ARM_ZeroInterpolation::NonParamCurve(ARM_Vector* terms, ARM_Vector* yields,
                                          int lambda, int nbPoints,
                                          ARM_Vector* newterms,
                                          ARM_Vector* newyields)
{
    int size = terms->GetSize();

    int k = (nbPoints-1)/(size-1); // points equi repartis
    int l = (nbPoints-1) - k*(size-1); // points restants

    int i,j;

    ARM_Vector* grad_0 = new ARM_Vector(nbPoints); // gradient en zero
    ARM_Vector* Sol = new ARM_Vector(nbPoints);    // solution

    ARM_Matrix* Hess = new ARM_Matrix(nbPoints,nbPoints); // hessienne

    ARM_Matrix* L = new ARM_Matrix(nbPoints,nbPoints); // matrice L
    ARM_Matrix* U = new ARM_Matrix(nbPoints,nbPoints); // matrice U

    if ( k < 0 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                       "Non Valid Number of Points");
    }

    /* construction des vecteurs de dates et de taux initiaux,
       calcul du gradient et du terme constant dans la hessienne */
    for (i = 0; i < size-1; i++)
    {
        grad_0->Elt(i*k) = 2.0*lambda*yields->Elt(i);
        Hess->Elt(i*k,i*k) = 2.0*lambda;

        for (j = 0; j < k; j++)
        {
            newterms->Elt(i*k+j) = (terms->Elt(i)*(k-j)+terms->Elt(i+1)*j)/k;
        }
    }

    newterms->Elt((size-1)*k) = terms->Elt(size-1);
    grad_0->Elt((size-1)*k) = 2.0*lambda*yields->Elt(size-1);
    Hess->Elt((size-1)*k,(size-1)*k) = 2.0*lambda;

    for (j = 1; j <= l; j++)
    {
        newterms->Elt(j+k*(size-1)) = terms->Elt(size-1) + j * 
                                      (newterms->Elt((size-1)*k) - 
                                      newterms->Elt((size-1)*k-1));
    }

    // normalisation du temps :
    double norm = newterms->Elt(nbPoints-1);

    for (i = 0 ; i < nbPoints ; i++)
    {
        newterms->Elt(i) = newterms->Elt(i) / norm;
    }

    // calcul du terme de regularisation dans la hessienne
    Calc_Hess(newterms,Hess);

    /* Factorisation de la matrice Hessienne en LU */
   
    factLU(Hess,L,U);

/*  Resolution du systeme. Calcul de la solution approchee u=-H^{-1}.grad(0) */

    resoud(L,U,grad_0,Sol);

    j = 0;
    for (i = 0; i < nbPoints; i++) 
    {
        newyields->Elt(i) = Sol->Elt(i);
        newterms->Elt(i) = newterms->Elt(i) * norm; // renormalisation
    }

    delete grad_0;
    delete Hess;
    delete Sol;
    delete L;
    delete U;
}


void ARM_ZeroInterpolation::resoud(ARM_Matrix* L, ARM_Matrix* U,
                                   ARM_Vector* Grad, ARM_Vector* Solution)
{
    int i,j;
    double som,*y;
    int size = Grad->GetSize();

    y = (double*) calloc(size, sizeof(double));

    y[0] = Grad->Elt(0);

    for (i = 1; i < size; i++)  /* descente */
    { 
        som = 0.0;

        for (j = 0; j < i; j++)
            som += L->Elt(i,j)*y[j];

        y[i] = (Grad->Elt(i)-som)/L->Elt(i,i); 
    }

    Solution->Elt(size-1)=y[size-1]/U->Elt(size-1,size-1);

    for (i = size-2; i >= 0; i--)
    {
        som=0.0;

        for(j = i+1; j < size; j++)
           som = som+U->Elt(i,j)*Solution->Elt(j);

        Solution->Elt(i)=(y[i]-som)/U->Elt(i,i);
    }

    free((double *)y);
}



void ARM_ZeroInterpolation::factLU(ARM_Matrix* Hess,
                                   ARM_Matrix* L, ARM_Matrix* U)
{
    int i,j,k,p,ll;
    double som;

    int size = Hess->GetNumLines();

    for (k = 0; k < size; k++)
    {
        L->Elt(k,k) = 1.0;

        ll = k+6;

        if ( ll > size ) 
           ll = size;

        for (j = k; j < ll; j++)
        {
            som=0.0;

            for (p = 0; p < k; p++)
                som += L->Elt(k,p)*U->Elt(p,j);

            U->Elt(k,j)=Hess->Elt(k,j)-som;  
        }

        for (i = k+1; i < ll; i++)
        {
            som = 0.0;

            for (p = 0; p < k; p++)
                som+=L->Elt(i,p)*U->Elt(p,k);

            if ( U->Elt(k,k) == 0 )
            {
               throw Exception(__LINE__, __FILE__, ERR_NUMERICAL_CALL,
                  "Division by Zero in ARM_ZeroInterpolation::factLU");
            }

            L->Elt(i,k)=(Hess->Elt(i,k)-som)/U->Elt(k,k); 
        }
    }
}


/* fonctions mathematiques */

double ARM_ZeroInterpolation::Regularisation(double x)
{
    double u1; 

    u1 = x*log(2.0+x); /* exp(-2*(T+xkm));  xkm*log(T+xkm);*/  

    return(u1);
}


void ARM_ZeroInterpolation::Remp_Theta(ARM_Vector* newterms,ARM_Matrix* theta)
{
    int size = newterms->GetSize();

    int i;

    for (i = 1; i < size-1 ; i++)
    {
        theta->Elt(i,0) = 2.0*Regularisation(newterms->Elt(i-1))
                          /((newterms->Elt(i-1)-newterms->Elt(i))
                          *(newterms->Elt(i-1)-newterms->Elt(i+1)));

        theta->Elt(i,1) = 2.0*Regularisation(newterms->Elt(i))/
                          ((newterms->Elt(i)-newterms->Elt(i-1))
                          *(newterms->Elt(i)-newterms->Elt(i+1)));

        theta->Elt(i,2) = 2.0*Regularisation(newterms->Elt(i+1))/
                          ((newterms->Elt(i+1)-newterms->Elt(i-1))
                          *(newterms->Elt(i+1)-newterms->Elt(i)));
    }

    i = size - 1;

    theta->Elt(i,0) = 2.0*Regularisation(newterms->Elt(i-1))
                      /((newterms->Elt(i-1)-newterms->Elt(i))*newterms->Elt(i-1));

    theta->Elt(i,1) = 2.0*Regularisation(newterms->Elt(i)) 
                      /((newterms->Elt(i-1)-newterms->Elt(i))*newterms->Elt(i));

    theta->Elt(i,2) = 0.0;
}



/* rem : le += dans Hess[i,i] apparait car il y a un terme constant,
         calcule auparavant */

void ARM_ZeroInterpolation::Calc_Hess(ARM_Vector* newterms, ARM_Matrix* Hess)
{
    int i,j;
    int size = newterms->GetSize();
    double u1,u2,u3;

    // remplissage de la matrice des theta_1,2,3
    ARM_Matrix *theta = new ARM_Matrix(size,3);

    Remp_Theta(newterms,theta);

    /* remplit hessienne j=i-2 */

    Hess->Elt(2,0) = (newterms->Elt(2) - newterms->Elt(1))
                     * theta->Elt(1,0) * theta->Elt(1,2);

    for (i=3; i<size-1; i++)
    {
        Hess->Elt(i,i-2) = (newterms->Elt(i) - newterms->Elt(i-2))
                           * theta->Elt(i-1,0) * theta->Elt(i-1,2);
    }

    Hess->Elt(size-1,size-3) = (newterms->Elt(size-1)-newterms->Elt(size-2))
                               * theta->Elt(size-2,0) * theta->Elt(size-2,2);

    /* remplit hessienne j=i-1 */

    Hess->Elt(1,0) = (newterms->Elt(2) - newterms->Elt(1)) * theta->Elt(1,0)
                     * theta->Elt(1,1);

    u1=(newterms->Elt(2)-newterms->Elt(1))*theta->Elt(1,1)*theta->Elt(1,2);
    u2=(newterms->Elt(3)-newterms->Elt(1))*theta->Elt(2,0)*theta->Elt(2,1);
    Hess->Elt(2,1)=u1+u2;

    for (i=3; i<size-2; i++)
    {
        u1 = (newterms->Elt(i+1)-newterms->Elt(i-1)) * theta->Elt(i,0)
             * theta->Elt(i,1) ;
        u2 = (newterms->Elt(i)-newterms->Elt(i-2)) * theta->Elt(i-1,1)
             * theta->Elt(i-1,2);
        Hess->Elt(i,i-1)=u1+u2;
    }

    u1 = (newterms->Elt(size-2)-newterms->Elt(size-3)) * theta->Elt(size-2,0)
         * theta->Elt(size-2,1);  
    u2 = (newterms->Elt(size-2)-newterms->Elt(size-4)) * theta->Elt(size-3,1)
         * theta->Elt(size-3,2);
    Hess->Elt(size-2,size-3)=u1+u2;

    u2 = (newterms->Elt(size-1)-newterms->Elt(size-2)) * theta->Elt(size-2,1)
         * theta->Elt(size-2,2);
    Hess->Elt(size-1,size-2)=u2;


    /* remplit hessienne i=j */

    Hess->Elt(0,0) += (newterms->Elt(2)-newterms->Elt(1)) * theta->Elt(1,0)
                      * theta->Elt(1,0);

    u1=(newterms->Elt(2)-newterms->Elt(1))*theta->Elt(1,1)*theta->Elt(1,1);
    u2=(newterms->Elt(3)-newterms->Elt(1))*theta->Elt(2,0)*theta->Elt(2,0);
    Hess->Elt(1,1)+=u1+u2;

    u1=(newterms->Elt(2)-newterms->Elt(1))*theta->Elt(1,2)*theta->Elt(1,2);
    u2=(newterms->Elt(3)-newterms->Elt(1))*theta->Elt(2,1)*theta->Elt(2,1);
    u3=(newterms->Elt(4)-newterms->Elt(2))*theta->Elt(3,0)*theta->Elt(3,0);   
    Hess->Elt(2,2)+=u1+u2+u3;

    for (i=3; i<size-3; i++)
    {
        u1 = (newterms->Elt(i+2)-newterms->Elt(i)) * theta->Elt(i+1,0)
             * theta->Elt(i+1,0) ;
        u2 = (newterms->Elt(i+1)-newterms->Elt(i-1)) * theta->Elt(i,1)
             * theta->Elt(i,1);
        u3 = (newterms->Elt(i)-newterms->Elt(i-2)) * theta->Elt(i-1,2)
             * theta->Elt(i-1,2);
        Hess->Elt(i,i)+=u1+u2+u3;
    }

    i=size-3;

    u1=(newterms->Elt(i+1)-newterms->Elt(i))*theta->Elt(i+1,0)
       *theta->Elt(i+1,0);
    u2=(newterms->Elt(i+1)-newterms->Elt(i-1))*theta->Elt(i,1)*theta->Elt(i,1);
    u3=(newterms->Elt(i)-newterms->Elt(i-2))*theta->Elt(i-1,2)
       *theta->Elt(i-1,2);
    Hess->Elt(i,i)+=u1+u2+u3;

    i=size-2;

    u1=(newterms->Elt(i)-newterms->Elt(i-1))*theta->Elt(i,1)*theta->Elt(i,1);
    u2=(newterms->Elt(i)-newterms->Elt(i-2))*theta->Elt(i-1,2)
       *theta->Elt(i-1,2);
    Hess->Elt(i,i)+=u1+u2;

    i=size-1;

    Hess->Elt(i,i) += (newterms->Elt(i)-newterms->Elt(i-1)) * theta->Elt(i-1,2)
                      * theta->Elt(i-1,2);

    /* matrice symetrique */
    for (i=0; i<size; i++)
    {
        for (j=0; j<i; j++)
        {
            Hess->Elt(j,i)=Hess->Elt(i,j);
        }
    }

    delete theta;
}


/*----------------------------------------------------------------------------*
   Returns the discount price with maturity yearTerm years from settlement, 
    computed from cubic splines.
*----------------------------------------------------------------------------*/

double ARM_ZeroInterpolation::DiscountFunction(double yearTerm)
{
    double z, zeroShift=0.0, intYield;
    

    if ( yearTerm < 0.0 )
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "YearTerm must be non negative");
    }

    if (GetBucketEndPeriod() > K_DOUBLE_TOL
        && yearTerm >= GetBucketStartPeriod() 
        && yearTerm < GetBucketEndPeriod())
    {
        zeroShift = 0.0001 * GetBPShift();
    }
    

    // cubic spline interpolation 
    
    intYield = ComputeInterpolValue(yearTerm)/100.0;

    // compute discount bond price

    if ( itsCompoundMeth == 0 )
    {
        z = exp(-yearTerm*(intYield + zeroShift));
    }

    if ( itsCompoundMeth == -1 )
    {
        z = 1.0/(1.0+yearTerm*(intYield + zeroShift));
    }
    
    if ( itsCompoundMeth > 0 )
    {
        z = pow(1.0 + (zeroShift + intYield) / itsCompoundMeth,
                 - yearTerm * itsCompoundMeth);
    }

    return(z);
}



/*----------------------------------------------------------------------------*
    Returns the d(discount price) / d(yearTerm) with maturity yearTerm 
    years from settlement, computed from linear interpolation.
*----------------------------------------------------------------------------*/

double ARM_ZeroInterpolation::D1DiscountFunction(double yearTerm)
{
    double    zp;

    if ( yearTerm <= 0.0 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                       "YearTerm must be non negative");
    }

    zp = (DiscountFunction(yearTerm+0.000001)
          -DiscountFunction(yearTerm-0.000001))*500000.0; 
    
    return(zp);
}


    
/*----------------------------------------------------------------------------*
    Returns the d2(discount price) / d(yearTerm)2 with maturity yearTerm 
    years from settlement, computed from splines.
*----------------------------------------------------------------------------*/
double ARM_ZeroInterpolation::D2DiscountFunction(double yearTerm)
{
    double zp;

    if ( yearTerm <= 0.0 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                       "YearTerm must be non negative");
    }

    zp = (D1DiscountFunction(yearTerm+0.000001)
           -D1DiscountFunction(yearTerm-0.000001))*500000.0; 

    return(zp);
}



/*------------------------ TAM/T4M Curve ------------------------------------*/


ARM_ZeroInterpolation::ARM_ZeroInterpolation(ARM_Date& asOf,
                                             char* Terms[ARM_NB_TERMS],
                                             ARM_Vector* mktData,
                                             double mean_rates,
                                             int raw, int Cont_Lin,
                                             ARM_Currency* ccy,
                                             int lambda, int nbPoints)
                     : ARM_ZeroCurve(asOf, Terms, mktData,
                                             mean_rates,
                                             raw, Cont_Lin, ccy)
 
{
    // Check variables
    itsLambda = lambda;

    itsNbPoints = nbPoints;

    Init();

    SetName(ARM_ZERO_INTERPOLATION);
 
    // Set variables
    ARM_Vector *yearTerms = GetYearTerms();
    ARM_Vector *zeroYields = GetZeroRates();
 
    int size = GetYearTerms()->GetSize();

    // mise en place de la procedure d'interpolation

    ARM_Vector* newterms = new ARM_Vector(itsNbPoints);

    ARM_Vector* newyields = new ARM_Vector(itsNbPoints);

    /* construction de l'interpolation...*/

    NonParamCurve(yearTerms, zeroYields, itsLambda, itsNbPoints, newterms,
                  newyields);

    // setting des parametres pour les splines cubiques

    ARM_Vector* terms = new ARM_Vector(newterms, itsNbPoints+2, 0,
                                       itsNbPoints-1, 1);

    terms->Elt(0) = 0.0;
    terms->Elt(itsNbPoints+1) = 1000.0;

    ARM_Vector* yields = new ARM_Vector(newyields, itsNbPoints+2, 0,
                                        itsNbPoints-1, 1);

    yields->Elt(0) = newyields->Elt(0) - (newyields->Elt(1)-newyields->Elt(0))
                     * newterms->Elt(0) / (newterms->Elt(1)-newterms->Elt(0));

    yields->Elt(itsNbPoints+1) = newyields->Elt(itsNbPoints-1);

    SetYearTerms(terms);

    SetZeroRates(yields);

    delete newterms;
    delete newyields;

    int Size = GetYearTerms()->GetSize();
 
    ARM_Vector* BPShifts = new ARM_Vector(Size, 0.0);
 
    for (int i = 0; i < Size; i++)
    {
        if ( GetYearTerms()->Elt(i) >= GetBucketStartPeriod()
             && GetYearTerms()->Elt(i) < GetBucketEndPeriod())
        {
            BPShifts->Elt(i) = GetBPShift()*0.01; // using yield in percent
        }
    }
 
    SetBPShifts(BPShifts);
 
    itsCompoundMeth = K_COMP_CONT;

    GenerateDateTerms();

    GenerateDiscountFactors(itsCompoundMeth);
}


/* interpolation splines cubiques */
double ARM_ZeroInterpolation::ComputeInterpolValue(double x)
{
    ARM_Vector *terms = GetYearTerms();
    ARM_Vector *yields = GetZeroRates();

    ARM_Matrix *polycoef = GetPolyCoef();

    int size = yields->GetSize();

    double result = 0.0;

    // Case out of array

    if ( x >= terms->Elt(size-2) )
    {
       result = yields->Elt(size-2);

       return(result);
    }

    if ( x <= terms->Elt(1) )
    {
       result = linInterpol(terms, x, yields);

       return(result);
    }

    // recherche de l'index

    int index = size-1;
    int trouve = 0;

    while (index > 1 && !trouve)
    {
        if (terms->Elt(index) >= x)
           index--;
        else
           trouve = 1;
    }

    result = polycoef->Elt(index, 0);
    result += polycoef->Elt(index, 1)*x;
    result += polycoef->Elt(index, 2)*x*x;
    result += polycoef->Elt(index, 3)*x*x*x;

    return(result);
}


void ARM_ZeroInterpolation::ComputePolynomCoef(void)
{
    int i, size;
    double t1, t2, t3, t4;
    double r1, r2, r3, r4;
    double det = 0;

    ARM_Vector *terms = GetYearTerms();
    ARM_Vector *yields = GetZeroRates();

    size = terms->GetSize();

    ARM_Matrix *coef = new ARM_Matrix(size, 4, 0.0);

    ARM_Matrix mat4(4, 4, 0.0);
    ARM_Vector *SndMb = new ARM_Vector(4);

    for (i=1; i < size-2; i++)
    {
       t1 = terms->Elt(i-1);
       t2 = terms->Elt(i);
       t3 = terms->Elt(i+1);
       t4 = terms->Elt(i+2);
     
       r1 = yields->Elt(i-1);
       r2 = yields->Elt(i);
       r3 = yields->Elt(i+1);
       r4 = yields->Elt(i+2);

       mat4.Elt(0,0) = t2*t2*t2;
       mat4.Elt(0,1) = t2*t2;
       mat4.Elt(0,2) = t2;
       mat4.Elt(0,3) = 1;
       mat4.Elt(1,0) = t3*t3*t3;
       mat4.Elt(1,1) = t3*t3;
       mat4.Elt(1,2) = t3;
       mat4.Elt(1,3) = 1;
       mat4.Elt(2,0) = 3*t2*t2;
       mat4.Elt(2,1) = 2*t2;
       mat4.Elt(2,2) = 1;
       mat4.Elt(2,3) = 0;
       mat4.Elt(3,0) = 3*t3*t3;
       mat4.Elt(3,1) = 2*t3;
       mat4.Elt(3,2) = 1;
       mat4.Elt(3,3) = 0;

       SndMb->Elt(0) = r2;
       SndMb->Elt(1) = r3;
       SndMb->Elt(2) = ((r2-r1)/(t2-t1) + (r3-r2)/(t3-t2))/2;
       SndMb->Elt(3) = ((r3-r2)/(t3-t2) + (r4-r3)/(t4-t3))/2;

       mat4.LinSolve(SndMb, det);

       if (det == 0.0)
       {
          throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Error invert matrix");
       }

       coef->Elt(i, 0) = SndMb->Elt(3);
       coef->Elt(i, 1) = SndMb->Elt(2);
       coef->Elt(i, 2) = SndMb->Elt(1);
       coef->Elt(i, 3) = SndMb->Elt(0);
    }


    SetPolyCoef(coef);

    delete SndMb;
}



ARM_ZeroCurve* ARM_ZeroInterpolation::GenerateShiftCurve(ARM_CRV_TERMS& Term, ARM_Vector* epsilon)
{
    if (GetMktData() == NULL)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
          "The input curve cant be shifted");
    }

    ARM_ZeroInterpolation* zc;
    ARM_Date tmpdat;
    ARM_MarketData* MktData = (ARM_MarketData*) GetMktData()->Clone();

    // verification de l'existence des datamarkets
    if (!MktData)
        return NULL;

    int szinp =0;
    szinp = epsilon->GetSize();

    int k=0;
    int l=0;

    for (l=0; l<szinp; l++)
    {
        for (k=0; k<ARM_NB_TERMS; k++)
        {
            if (!strcmp(Term[l],MktData->itsMktTerms[k]))
            {
               MktData->itsMktValue->Elt(k) = MktData->itsMktValue->Elt(k)
                                              +epsilon->Elt(l);

               k = ARM_NB_TERMS;
            }
        }
    }

    switch (MktData->itsConstructionMeth)
    {
        case 0 :
                zc = new ARM_ZeroInterpolation(GetAsOfDate(),
                                MktData->itsMktTerms,
                                MktData->itsMktValue,
                                MktData->itsMMVsFut,
                                MktData->itsSwapVsFut,
                                MktData->itsraw,
                                MktData->itsCont_Lin,
                                GetCurrencyUnit(),
                                this->itsLambda,
                                this->itsNbPoints);
        break;

        case 1 :
                zc = new ARM_ZeroInterpolation(MktData->itsMktTerms,
                                GetAsOfDate(),
                                MktData->itsMktValue,
                                MktData->itsMMVsFut,
                                MktData->itsSwapVsFut,
                                MktData->itsraw,
                                MktData->itsCont_Lin,
                                GetCurrencyUnit(),
                                this->itsLambda,
                                this->itsNbPoints);

        default:
        { 
               if (zc)
                  delete zc;

               zc = NULL;

               delete MktData;

               return NULL;
        }
    }

    delete MktData;

    return(zc);
}


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
