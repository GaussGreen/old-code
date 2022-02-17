/*
 * $Log: markovtreepricer.h,v $
 * Revision 1.4  2003/07/21 16:00:55  jmprie
 * modifs pour calcul des probas d'exercice
 *
 * Revision 1.3  2003/06/24 14:46:00  jmprie
 * ajout des fcts pour le calcul analytique des esperances via regression
 * polynomiale
 *
 * Revision 1.2  2003/05/15 09:15:24  ykhlif
 * Add computePrices for pricing option and underlyings
 *
 * Revision 1.1  2003/03/07 17:40:17  mab
 * Initial revision
 *
 * Revision 1.1  2003/02/10 08:51:18  mab
 * Initial revision
 *
 */ 


/*----------------------------------------------------------------------------*
 
  markovtreepricer.h

  Header of the Markov Tree Pricer 

*----------------------------------------------------------------------------*/



#ifndef _MARKOVTREEPRICER_H
#define _MARKOVTREEPRICER_H

#include "treepricer.h"




class ARM_MarkovTreePricer : public ARM_TreePricer
{
    public :


    ARM_MarkovTreePricer(void)
    {
    }

    ARM_MarkovTreePricer(ARM_Security* sec, ARM_Model* mod);

    ARM_MarkovTreePricer(const ARM_MarkovTreePricer& pricer)
    {
        this->BitwiseCopy(&pricer);
    }



    double Price(void); 

    ARM_Vector* ComputePrices(void);
    ARM_Vector* ComputePrices(bool isWithExerProba);

    ARM_Vector* NodePolynomialInterpolation(int nodeIdx,bool isTrunc,
                    int nbPt,ARM_Vector* state,
                    int payOffIdx,ARM_Matrix* payOff);

    ARM_Vector* NodePolynomialInterpolation(int nodeIdx,
                    double stateMin,double stateMax,
                    int nbPt,ARM_Vector* state,
                    int payOffIdx,ARM_Matrix* payOff);

    double ComputeNodeAnalyticalExpectation(double mean,double stdDev,
                   ARM_Vector* polyInf,double xstar,ARM_Vector* polySup);

    double ComputeNodeAnalyticalExpectation(double mean,double stdDev,
                    ARM_Vector* polyCoef);

    double ComputeNodeNumericalExpectation(int nodeIdx,bool isTrunc,
                ARM_Vector* proba,int payOffIdx,ARM_Matrix* payOff);

    void ComputeSliceAnalyticalExpectation(int treeIdx,ARM_Object* treeState,
                    int prevTreeIdx,ARM_GenMatrix** treePrice);


    // Services

    ARM_MarkovTreePricer &operator = (const ARM_MarkovTreePricer& pricer)
    {
        (*this).ARM_TreePricer::operator = (pricer);

        BitwiseCopy(&pricer);

        return (*this);
    }

    
    void BitwiseCopy(const ARM_Object* opricer)
    {
    }

    void Copy(const ARM_Object* src)
    {    
        ARM_TreePricer::Copy(src);

        BitwiseCopy(src);
    }

};

#endif
/*------------------------------------------------------------------------------*/
/*---- End Of File ----*/
