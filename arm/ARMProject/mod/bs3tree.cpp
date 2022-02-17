/*----------------------------------------------------------------------------*
    bs3tree.cc
 
    This file implements the ARM_BSTrinoTree class, an abstract class 
    for solving parabolic
    PDEs using a lattice method. 
    Always define subclasses of this class. Subclasses should 
    implement their own
    propagation methods using a relevant finite difference scheme. Check for
    methods that must be overriden.

*----------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "util.h"
#include "armglob.h"
#include "crrtree.h"



/*---------------------------------------------------------------------------*
    Set
*----------------------------------------------------------------------------*/

void ARM_BSTrinoTree::Set(ARM_Date& startDate, ARM_Date& horizon, 
                          int numSteps, 
                          double spot, double dividend, 
                          double discountRate, double volatility)
{
    ARM_PDESolver::Set(startDate, horizon, numSteps);

    itsSpot = spot;
    itsDividendYield = dividend;
    itsDiscountRate = discountRate;
    itsVolatility = volatility;
}



/*----------------------------------------------------------------------------*
    Constructor        
*----------------------------------------------------------------------------*/

ARM_BSTrinoTree::ARM_BSTrinoTree(ARM_Date &startDate, 
                         ARM_Date &horizon, 
                         int numSteps, 
                         double spot,
                         double dividend, 
                         double discountRate, 
                         double volatility)
            :ARM_PDESolver(startDate, horizon, numSteps)
{
    Init();

    SetName(ARM_BSTrinoTree);

    itsSpot = spot;
    itsDividendYield = dividend;
    itsDiscountRate = discountRate;
    itsVolatility = volatility;
}



/*----------------------------------------------------------------------------*
    Constructor    (copy).
*----------------------------------------------------------------------------*/

ARM_BSTrinoTree::ARM_BSTrinoTree(ARM_BSTrinoTree &tree) : ARM_PDESolver(tree)
{
    Init();

    SetName(ARM_BSTrinoTree);

    itsSpot = tree.itsSpot;
    itsDividendYield = tree.itsDividendYield;
    itsDiscountRate = tree.itsDiscountRate;
    itsVolatility = tree.itsVolatility;
}



/*----------------------------------------------------------------------------*
    Assignment operator
*----------------------------------------------------------------------------*/

ARM_BSTrinoTree& ARM_BSTrinoTree::operator = (ARM_BSTrinoTree& tree)
{
    (*this).ARM_PDESolver::operator = (tree);

    itsSpot = tree.itsSpot;
    itsDividendYield = tree.itsDividendYield;
    itsDiscountRate = tree.itsDiscountRate;
    itsVolatility = tree.itsVolatility;

    return(*this);
}



/*----------------------------------------------------------------------------*
    Compute Prices in the binomial tree at settlement 
*----------------------------------------------------------------------------*/

ARM_Vector* ARM_BSTrinoTree::CFPrices(ARM_Date& settlement, 
                                 ARM_Vector* cfDates, 
                                 ARM_Vector* cfValues, 
                                 int DomOrFrg) 
{
    int i, t;

    
    if ( settlement > GetHorizon() )
    {
       throw Exception(__LINE__, __FILE__, ERR_YEAR_TERMS_DISCRETE_TIME,
            "settlement does not fit in discrete time");
    }
    
    try
    {
        //    get last time index before settlement

        t = GetDateIndex(settlement);
    }

    catch(Exception& m)
    {
        m.DebugPrint();
        throw Exception(__LINE__, __FILE__, ERR_YEAR_TERMS_DISCRETE_TIME,
                             "Year term does not fit in discrete time");
    }

    ARM_Vector *prices;
    
    double dt, upStepLength;

    if ( t == 0 )
    {
        prices = new ARM_Vector(1, itsSpot);

        return(prices);
    }
    else
    {
        dt = (*GetDiscreteTime())[t] - (*GetDiscreteTime())[t-1];
    
        upStepLength = exp(0.01*itsVolatility*sqrt(dt));
    
        prices = new ARM_Vector(t+1, itsSpot*pow(upStepLength, t));

        for (i = 1; i < t+1; i++) 
        {
            (*prices)[i] = (*prices)[i-1]/SQR(upStepLength);
        }
    }
    
    return(prices);
}



/*----------------------------------------------------------------------------*
    Compute MaxPrices at stateIndex-th node of timeIndex 
*----------------------------------------------------------------------------*/

ARM_Vector* ARM_BSTrinoTree::GetMaxPrices(int timeIndex, int stateIndex, 
                                          double upStepLength)
{
    ARM_Vector* MaxPrices;



    if ( timeIndex == 0 )
    {
        MaxPrices = new ARM_Vector(1, itsSpot);
    }
    else
    {
        int nMax, j;

        nMax=1 + MIN(stateIndex, timeIndex-stateIndex);

        MaxPrices = new ARM_Vector(nMax, 0.0);

        for (j=0; j<nMax; j++)
        {
           //LES MAXIMUMS SONT DANS L'ORDRE DECROISSANTS
           //L'INDICE 0 CORRESPOND AU PLUS GRAND MAXIMUM POSSIBLE
           //L'INDICE NumberMax CORRESPOND AU PLUS PETIT MAXIMUM POSSIBLE

           (*MaxPrices)[j]=itsSpot*pow(upStepLength,timeIndex-stateIndex-j);
        }
    }

    return(MaxPrices);
}



/*----------------------------------------------------------------------------*
    Compute next nodes Indexes (at timeIndex+1)
    from stateIndex-th node max prices at timeIndex
    
    The 1st elt corresponds to up mvt, the 2nd to down mvt
/*---------------------------------------------------------------------------*/

ARM_Vector* ARM_BSTrinoTree::GetNextNodes(int tIndex, int stateId)
{
    ARM_Vector *nextNodes = new ARM_Vector(2, (double) stateId);

    nextNodes->Elt(1) = stateId+1.0;
    
    return(nextNodes);
}



ARM_Vector* ARM_BSTrinoTree::GetProbaVector(int stateIndex, int timeIndex)
{
    double dt, upStepLength, rnDrift, pUp;
    
    ARM_Vector* prob = new ARM_Vector(2, 0.0);
    dt = (*GetDiscreteTime())[timeIndex+1] - (*GetDiscreteTime())[timeIndex];
    
    upStepLength = exp(0.01*itsVolatility*sqrt(dt));
    
    rnDrift = exp(0.01*dt*(itsDiscountRate-itsDividendYield));

    pUp = (rnDrift-1.0/upStepLength)/(upStepLength-1.0/upStepLength); 

    prob->Elt(0) = pUp;
    prob->Elt(1) = 1-pUp;
    
    return(prob);
}



/*----------------------------------------------------------------------------*
    Compute next nodes MaxPrices (at timeIndex+1)
    from stateIndex-th node max prices at timeIndex
    
    The 1st elt corresponds to up mvt, the 2nd to down mvt
*----------------------------------------------------------------------------*/

ARM_Vector* ARM_BSTrinoTree::NextMaxPrices(double maxPrice, 
                                           double upStepLength)
{
    ARM_Vector *MaxPrices;
    
    MaxPrices = new ARM_Vector(2, maxPrice);

    MaxPrices->Elt(0) = maxPrice*upStepLength;
    
    return(MaxPrices);
}



/*----------------------------------------------------------------------------*
    Compute next nodes MinPrices (at timeIndex+1)
    from stateIndex-th node max prices at timeIndex 
    
    The 1st elt corresponds to up mvt, the 2nd to down mvt
/*---------------------------------------------------------------------------*/

ARM_Vector* ARM_BSTrinoTree::NextMinPrices(double minPrice, 
                                           double upStepLength)
{
    ARM_Vector *MinPrices;
    
    MinPrices = new ARM_Vector(2, minPrice);

    MinPrices->Elt(0) = minPrice/upStepLength;
    
    return(MinPrices);
}



/*----------------------------------------------------------------------------*
    Compute MinPrices at stateIndex-th node of timeIndex 
/*---------------------------------------------------------------------------*/

ARM_Vector* ARM_BSTrinoTree::GetMinPrices(int timeIndex, int stateIndex, 
                                          double upStepLength)
{
    ARM_Vector *MinPrices;

    if (timeIndex==0)
    {
        MinPrices = new ARM_Vector(1, itsSpot);
    }
    else
    {
        int nMin, j;

        nMin=1+MIN(stateIndex, timeIndex-stateIndex);

        MinPrices = new ARM_Vector(nMin, 0.0);

        for (j=0; j<nMin; j++)
        {
            //LES MINIMUMS SONT DANS L'ORDRE CROISSANTS
            //L'INDICE 0 CORRESPOND AU PLUS PETIT MINIMUM POSSIBLE
            //L'INDICE NumberMax CORRESPOND AU PLUS GRAND MAXIMUM POSSIBLE

            (*MinPrices)[j]=itsSpot*pow(upStepLength, j-stateIndex);
        }
    }

    return(MinPrices);
}



void ARM_BSTrinoTree::DiscountValues(int timeIndex, double dt, 
                                     ARM_GenMatrix** value, int DomOrFrg)
{
    int l, k, nL, nC;

    nL = (*value)->GetNumLines();
    nC = (*value)->GetNumCols();
    
    for (k=0; k<nC; k++)
    {
        for (l=0; l<nL; l++)
        {
            (*value)->Elt(l, k) *= exp(-0.01*itsDiscountRate*dt);
        }
    }
}



/*----------------------------------------------------------------------------*
    Do backward propagation of value, e.g. solve Kolmogorov
    Backward equation from time2 back to time1. 
    
      Note that contrarily to yield curve tree :
        1) back propagation will begin at tIndex2+1 if time2 does not fit 
           in a discrete date of the tree : obsolete rmk?
        2) time2>itsHorizon will provoque an error  
*----------------------------------------------------------------------------*/
        
void ARM_BSTrinoTree::PropagateBackward(ARM_Date time1,
                ARM_Date time2, 
                ARM_GenMatrix** value,
                int kDiscount, int DomOrFrg)
{
    int    tIndex1, tIndex2, t;
    double    dt;
    
    //  check time2 is less than itsHorizon
    if (time2 > GetHorizon())
    {
        throw Exception(__LINE__, __FILE__, ERR_YEAR_TERMS_DISCRETE_TIME,
            "time1 and/or time2 do not fit in discrete time");
    }

     try
    {
        //    get last time index before time2
        tIndex2 = GetDateIndex(time2);

        //    get first time index before time1
        tIndex1 = GetDateIndex(time1);
    }

    catch(Exception& m)
    {
        m.DebugPrint();
        throw Exception(__LINE__, __FILE__, ERR_YEAR_TERMS_DISCRETE_TIME,
             "time1 and/or time2 do not fit in discrete time");
    }

    //    check value has correct size
    if ((*value)->GetNumLines() != GetTotalStateSize(tIndex2))
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Invalid size of value vector at the specified 2nd time index");
    }

    //    case time1 and time2 within same time step

    if (tIndex1==tIndex2)
    {
        if (kDiscount)
        {
            dt = CountYears(KACTUAL_365, time1, time2);

            DiscountValues(tIndex1, dt, value);
        }
    }

    //    general case    
    
    else
    {
        //    discount back to tIndex2 if discountRate != NULL
    
        if (kDiscount) 
        {
            dt = CountYears(KACTUAL_365, time2, 
                      ARM_Date((*GetDiscreteDayTerms())[tIndex2]));
            
            if ( dt > 0.0 )
            {
                DiscountValues(tIndex2, dt, value);
            }
        }

        // propagate from tIndex2 back to tIndex1+1
    

        for (t=tIndex2-1; t>=tIndex1+1; t--)
        {

            try
            {
                PropagateOneStepBackward(t, value, kDiscount);
            }

            catch(Exception& m)
            {
                m.DebugPrint();

                throw Exception(__LINE__, __FILE__, ERR_PROPAGATION_PB,
                             "Backward propagation problem");
            }
        }
    
        // propagate back to time1

        try
        {
            PropagateOneStepBackward(tIndex1, value, 0);
        }

        catch(Exception& m)
        {
            m.DebugPrint();
            throw Exception(__LINE__, __FILE__, ERR_PROPAGATION_PB,
                             "Backward propagation problem");
        }

        if (kDiscount) 
        {
            dt = CountYears(KACTUAL_365, time1, 
                            ARM_Date((*GetDiscreteDayTerms())[tIndex1+1]));

            if (dt>0.0)
            {
                DiscountValues(tIndex1, dt, value);
            }
        }
    }
}


 
/*----------------------------------------------------------------------------*
    NAME    ARM_BSTrinoTree::PropagateOneStepBackward

    Do one step backward propagation of value, e.g. solve Kolmgorov
    Backward equation from timeIndex+1 discrete time back to timeIndex 
    discrete time. 
        
*----------------------------------------------------------------------------*/

void ARM_BSTrinoTree::PropagateOneStepBackward(int timeIndex, 
                                            ARM_GenMatrix** value,
                                            int kDiscount, int DomOrFrg)
{
    int i, k, nCols;
    double dt, upStepLength, pUp, rnDrift;


    // Copy value into propagation array
    ARM_GenMatrix* propArray = NULL;

    if ((*value)->GetNumLines() != timeIndex + 2)
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "Invalid size of value vector at the specified time index");
    }

    // Resize value

    nCols = (*value)->GetNumCols();

    if ( (*value)->GetName() == ARM_MATRIX )
    {
       ARM_Matrix* newArray = (ARM_Matrix *) (*value)->Clone();
 
       delete *value;

       *value = new ARM_Matrix(timeIndex+1, nCols);

       propArray = newArray;
    }
    else
    {
       ARM_Vector* newArray = (ARM_Vector *) (*value)->Clone();
 
       delete *value;

       *value = new ARM_Vector(timeIndex+1, 0.0);

       propArray = newArray;
    }

    dt = (*GetDiscreteTime())[timeIndex+1]-(*GetDiscreteTime())[timeIndex];
    
    upStepLength = exp(0.01*itsVolatility*sqrt(dt));
    
    rnDrift = exp(0.01*dt*(itsDiscountRate-itsDividendYield));

    pUp = (rnDrift-1.0/upStepLength)/(upStepLength-1.0/upStepLength); 

    // Take expectation
    
    for (k = 0; k < nCols; k++)
    {
        for (i = 0; i < timeIndex+1; i++)
        {
            (*value)->Elt(i, k) = pUp*propArray->Elt(i, k)+ 
                                   (1.0-pUp)*propArray->Elt(i+1, k);
        }
    }
    
    delete propArray;
    propArray = NULL;
    
    // Discount if discountRate != NULL
    
    if (kDiscount) 
    {
       DiscountValues(timeIndex, dt, value);
    }
}



/*----------------------------------------------------------------------------*
    NAME    ARM_BSTrinoTree::HistPropagateOneStepBackward

    Do one step backward propagation of value, e.g. solve Kolmgorov
    Backward equation from timeIndex+1 discrete time back to timeIndex 
    discrete time. 
    Returns 0 if failed.

*----------------------------------------------------------------------------*/

void ARM_BSTrinoTree::HistPropagateOneStepBackward(int timeIndex, 
                                            ARM_GenMatrix** value,
                                            int kDiscount, int DomOrFrg)
{
    int i, k;
    double dt, upStepLength, pUp, histDrift;

    //    copy value into propagation array

    ARM_GenMatrix *propArray = *value;


    if ((*value)->GetNumLines() != timeIndex + 2)
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
           "Invalid size of value vector at the specified time index");
    }

    
    // resize value
    
    if ((*value)->GetName() == ARM_MATRIX)
    {
        *value = new ARM_Matrix(timeIndex+1, (*value)->GetNumCols());
    }
    else
    {
        *value = new ARM_Vector(timeIndex+1, 0.0);
    }
    
    dt = (*GetDiscreteTime())[timeIndex+1] - (*GetDiscreteTime())[timeIndex];
    
    upStepLength = exp(0.01*itsVolatility*sqrt(dt));
    
    histDrift = exp(0.01*dt*(itsHistDrift-itsDividendYield));

    pUp = (histDrift-1.0/upStepLength)/(upStepLength-1.0/upStepLength); 

    //    take expectation
    
    for (i=0; i<timeIndex+1; i++)
    {
        for (k=0; k<(*value)->GetNumCols(); k++)
        {
            (*value)->Elt(i, k) = pUp*propArray->Elt(i, k) + 
                (1.0-pUp)*propArray->Elt(i+1, k);
        }
    }
    
    delete propArray;
    propArray = NULL;
}


 
/*----------------------------------------------------------------------------*
    NAME    ARM_BSTrinoTree::PropagateOneStepForward

    Do one step forward propagation of value, e.g. solve Kolmogorov Forward
    equation from timeIndex discrete time to timeIndex+1 discrete time.
    Returns 0 if failed.

    INPUT
        int     timeIndex    - propagation is from timeIndex+1 back totimeIndex
    INPUT / OUTPUT
        ARM_Vector     *value    - ARM_Vector of values to propagate
        
*----------------------------------------------------------------------------*/

void ARM_BSTrinoTree::PropagateOneStepForward(int timeIndex,
                         ARM_Vector *&value,
                         int kDiscount)
{
    int i;
    double dt, upStepLength;

    ARM_Vector * propArray;


    if ( value->GetSize() != timeIndex + 19 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
           "Invalid size of value vector at the specified time index");
    }

    
    // Copy value into propagation array
        
    propArray = value;
    
    value = new ARM_Vector(timeIndex+2, 0.0);

    dt = (*GetDiscreteTime())[timeIndex+1] - (*GetDiscreteTime())[timeIndex];
    
    upStepLength = exp(0.01*itsVolatility*sqrt(dt));
    
    //    propagate forward
    
    for (i=0; i<timeIndex+1; i++) 
    {
        (*value)[i] = (*propArray)[i] * upStepLength;
    }
    
    (*value)[timeIndex+1] = (*value)[timeIndex+1]/SQR(upStepLength);
    
    delete propArray;
        
}


 
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
