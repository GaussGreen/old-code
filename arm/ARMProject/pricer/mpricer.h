/*
 * $Log: mpricer.h,v $
 * Revision 1.5  2003/07/01 19:37:59  jpriaudel
 * suppression of include vcmarray
 *
 * Revision 1.4  2003/06/18 15:13:17  ebenhamou
 * remove include
 *
 * Revision 1.3  2003/06/18 11:31:44  ebenhamou
 * remove STL warning
 *
 *
 */


#ifndef _MPRICER_H
#define _MPRICER_H

#include "treepricer.h"




/**********************************************/
/*         Version 1D                         */
/**********************************************/

class ARM_MTreePricer1 : public ARM_TreePricer
{
private :

    bool itsSimplePayoff;

    int itsNbHoldValue;
    int itsNbExtraInfo;     
    

public :

	ARM_MTreePricer1(void);

	ARM_MTreePricer1(ARM_Security *sec, ARM_Model *mod);

    
    double Price(void) ;

	// Services

	ARM_MTreePricer1(const ARM_MTreePricer1 &pricer);

    ARM_MTreePricer1 &operator = (const ARM_MTreePricer1 &pricer);

    
    void BitwiseCopy(const ARM_Object* opricer) {};

    void Copy(const ARM_Object* src)
    {
        ARM_TreePricer::Copy(src);

        BitwiseCopy(src);
    }

    void Init(void)
    {
        itsSimplePayoff    = true;
        itsNbHoldValue     = 1;
        itsNbExtraInfo = 0;
    }

    void BeFittedTo(ARM_Security *sec);

};



/**********************************************/
/*         Version Template                   */
/**********************************************/
/*

// P is a payoff type :
// double for simple payoff, ARM_Vector otherwise
// N is the dimension of the pricer (tree)

//template <int N>
class ARM_MTreePricer : public ARM_TreePricer
{
private :

    bool itsSimplePayoff;

    int itsNbHoldValue;
    int itsNbExtraInfo;     
    

public :

	ARM_MTreePricer(void);

	ARM_MTreePricer(ARM_Security *sec, ARM_Model *mod);

    
    // Price  descendra des slices de type 
    // MArray <P,N>
    double Price(void) ;

	// Services

	ARM_MTreePricer(const ARM_MTreePricer &pricer);

    ARM_MTreePricer &operator = (const ARM_MTreePricer &pricer);

    
    void BitwiseCopy(const ARM_Object* opricer) {};

    void Copy(const ARM_Object* src)
    {
        ARM_TreePricer::Copy(src);

        BitwiseCopy(src);
    }

    void Init(void)
    {
        itsSimplePayoff    = true;
        itsNbHoldValue     = 1;
        itsNbExtraInfo = 0;
    }

    void BeFittedTo(ARM_Security *sec);

};



*/

#endif
