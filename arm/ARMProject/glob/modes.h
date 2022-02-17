/*
 * $Log: modes.h,v $
 * Revision 1.10  2002/05/30 13:46:15  mab
 * Improvement YANLI
 *
 * Revision 1.9  2001/08/17 15:42:15  smysona
 * modifs accesseurs.
 *
 * Revision 1.8  2001/08/14 13:11:03  smysona
 * Ajout DIAGCALIB
 *
 * Revision 1.7  2001/07/30 09:15:15  smysona
 * modif des enums
 *
 * Revision 1.6  2001/05/23 17:11:30  smysona
 * Ajout de MC modes supplementaires
 *
 * Revision 1.5  2001/04/27 09:25:36  smysona
 * IsCallOrPut
 *
 * Revision 1.4  2001/04/23 09:30:46  smysona
 * Supression des modes a 0
 *
 * Revision 1.3  2001/04/03 12:08:15  nicolasm
 * masques de CALIB_PORTFOLIO_MODE
 *
 * Revision 1.2  2001/03/06 17:09:03  smysona
 * Ajout de BS dans les expectation modes
 *
 * Revision 1.1  2001/01/30 09:52:14  smysona
 * Initial revision
 *
 */

#ifndef _CLASS_MODES_H
#define _CLASS_MODES_H




#include "vector"
#include "armdef.h"





/*-------- FRM Tree management Mode --------*/
/* UTILS defined in frmutils.               */


class TREE_MANAGEMENT_MODE
{
public :
    enum modes
    {
        precompute = 1,             // data (probas, forwards) are precomputed
        prebuild = precompute << 1,               // the tree is prebuilt
        memory = prebuild << 1      // do not destroy the tree when backward propagating
    };

    static inline bool IsPrecomputeMode(int mode)
    {
        return ((mode & TREE_MANAGEMENT_MODE::precompute) 
                 == TREE_MANAGEMENT_MODE::precompute);
    }

    static inline bool IsPrebuildMode(int mode)
    {
        return ((mode & TREE_MANAGEMENT_MODE::prebuild) == TREE_MANAGEMENT_MODE::prebuild);
    }

    static inline bool IsMemoryMode(int mode)
    {
        return ((mode & TREE_MANAGEMENT_MODE::memory) 
                 == TREE_MANAGEMENT_MODE::memory);
    }
};



class BRANCH_MODE
{
    public :

        enum modes
        {
           Uninomial = 1,
           Binomial = 2,
           Trinomial = 3
        };

    static inline bool IsUninomialMode(int mode)
    {
        return((mode & BRANCH_MODE::Uninomial) == BRANCH_MODE::Uninomial);
    }

    static inline bool IsBinomialMode(int mode)
    {
        return((mode & BRANCH_MODE::Binomial) == BRANCH_MODE::Binomial);
    }

    static inline bool IsTrinomialMode(int mode)
    {
        return((mode & BRANCH_MODE::Trinomial) == BRANCH_MODE::Trinomial);
    }

};



class VALUATION_CONTEXT
{
    public :

        enum contexts 
        {
            Price = 1,
            Hedge = Price << 1
        }; 


    static inline bool IsContextPrice(int Type)
    {
        return ((Type & VALUATION_CONTEXT::Price) == VALUATION_CONTEXT::Price);
    }

    static inline bool IsContextHedge(int Type)
    {
        return ((Type & VALUATION_CONTEXT::Hedge) == VALUATION_CONTEXT::Hedge);
    }
};



class CORRECTION_TYPE
{
    public :

        enum types 
        {
            None    = 1 ,
            PathDep = None << 1,    // for path dependancy
            Barrier = PathDep << 1, // for barriers
            FRM     = Barrier << 1, // for FRM model
            BS      = FRM << 1      // Available correction for euro option (last but one layer)
        };

    inline static bool IsCorrPathDep(int Type)
    {
        return ((Type & CORRECTION_TYPE::PathDep) == CORRECTION_TYPE::PathDep);
    }

    inline static bool IsCorrNone(int Type)
    {
        return ((Type & CORRECTION_TYPE::None) == CORRECTION_TYPE::None);
    }

    inline static bool IsCorrBarrier(int Type)
    {
        return ((Type & CORRECTION_TYPE::Barrier) == CORRECTION_TYPE::Barrier);
    }

    inline static bool IsCorrFRM(int Type)
    {
        return ((Type & CORRECTION_TYPE::FRM) == CORRECTION_TYPE::FRM);
    }

    inline static bool IsCorrBS(int Type)
    {
        return ((Type & CORRECTION_TYPE::BS) == CORRECTION_TYPE::BS);
    }
};



class CALIB_PORTFOLIO_MODE
{
    public :

         enum types 
         {
             None  = 1,
             Strip = None << 1, // strip bermuda product into europen
        Early   = Strip   << 1, // adds early calibration products (before lockout)
        Control = Early   << 1, // builds for control portfolio
        Level1  = Control << 1, // strip sur un seul niveau de profondeur, sinon strip en profondeur
        Under   = Level1  << 1, // adds the underlying of the option in the portfolio
        Call    = Under   << 1, // we build PF for a call on the underlying
        Put     = Call    << 1, // we build PF for a put on the underlying
        Hedge   = Put     << 1,
        Sort    = Hedge   << 1  // sort by price
    };    
};



class PAYOFF_MODE
{
    public :

        enum types
        {
            Payoff = 1,
            Probas = Payoff << 1,
            Correction = Probas << 1
        };
};



class AUTO_MODE
{
     public :
    
         enum _type
         {
             NONE = 1,
             IRG  = NONE << 1, // calibrate on swaptions
             IDX  = IRG << 1,  // calibrate on caps
             SWOPT_IRG = IDX <<1,   // calibrate on swaptions, optimize on caps
        IRG_SWOPT = SWOPT_IRG << 1, // calibrate on caps, optimize on swaptions
        DUAL      = IRG_SWOPT << 1, // calibrate swaptions and captions, 2 vols curves
        AUTOMEAN  = DUAL << 1,      // optimize the mean reversion rate
        BASISDISC = AUTOMEAN << 1,  // basis discounting
        CAPDIAG   = BASISDISC << 1,
        DIAG      = CAPDIAG << 1
    };

    static inline bool IsAutoModeNone(int mode)
    {
        return ((mode & AUTO_MODE::NONE) == AUTO_MODE::NONE);
    }

    static inline bool IsAutoModeIdx(int mode)
    {
        return ((mode & AUTO_MODE::IDX) == AUTO_MODE::IDX);
    }

    static inline bool IsAutoModeIrg(int mode)
    {
        return ((mode & AUTO_MODE::IRG) == AUTO_MODE::IRG);
    }

    static inline bool IsAutoModeSwoptIrg(int mode)
    {
        return ((mode & AUTO_MODE::SWOPT_IRG) == AUTO_MODE::SWOPT_IRG);
    }

    static inline bool IsAutoModeIrgSwopt(int mode)
    {
        return ((mode & AUTO_MODE::IRG_SWOPT) == AUTO_MODE::IRG_SWOPT);
    }

    static inline bool IsAutoModeDual(int mode)
    {
        return ((mode & AUTO_MODE::DUAL) == AUTO_MODE::DUAL);
    }

    static inline bool IsAutoModeAutoMean(int mode)
    {
        return ((mode & AUTO_MODE::AUTOMEAN) == AUTO_MODE::AUTOMEAN);
    }

    static inline bool IsAutoModeBasisDisc(int mode)
    {
        return ((mode & AUTO_MODE::BASISDISC) == AUTO_MODE::BASISDISC);
    }

    static inline bool IsAutoModeCapDiag(int mode)
    {
        return ((mode & AUTO_MODE::CAPDIAG) == AUTO_MODE::CAPDIAG);
    }

    static inline bool IsAutoModeDiag(int mode)
    {
        return ((mode & AUTO_MODE::DIAG) == AUTO_MODE::DIAG);
    }

    static inline int ValAutoShapeMode(int mode)
    {
        return (IsAutoModeDiag(mode)?K_DIAG:K_ROW);
    }



};



class MAV_MODE
{
    public :

        enum mode
        {
            Fill = 1,
            Resize = Fill << 1,
            Alloc = Resize << 1
        };
};



class EXPECTATION_MODE
{
    public :

        enum mode
        {
            envelopeSup = 1,
            forceExe    = envelopeSup << 1,
            forceHold   = forceExe    << 1,
            pathDepCorr = forceHold   << 1,
            BS          = pathDepCorr << 1
        };
};



class MC_MODE
{
    public :

        enum modes
        {
            control = 1,             
            OLS     = control << 1,  
            entropy = OLS << 1   
        };

    static inline bool IsModeContol(int mode)
    {
        return((mode & MC_MODE::control) == MC_MODE::control);
    }

    static inline bool IsModeOLS(int mode)
    {
        return ((mode & MC_MODE::OLS) == MC_MODE::OLS);
    }

    static inline bool IsModeEntropy(int mode)
    {
        return ((mode & MC_MODE::entropy) == MC_MODE::entropy);
    }

};


class XCCY_MODE
{
     public :

         enum modes
         {
             none            = 0,
             basisDiscount   = 1,
             translate       = basisDiscount << 1
         };
};



class VOL_TYPE
{
    public :

        enum modes
        {
            IRG = 1,             
            IDX     = IRG << 1  
        };

    static inline bool IsVolIDX(int mode)
    {
        return ((mode & VOL_TYPE::IDX) == VOL_TYPE::IDX);
    }

    static inline bool IsVolIRG(int mode)
    {
        return ((mode & VOL_TYPE::IRG) == VOL_TYPE::IRG);
    }
};





/***************************************************************/
/*                                                             */
/*            MASQUES DE LECTURE                               */
/*                                                             */
/***************************************************************/








/*****************************************************************/
/*                 Masques de Lecture des branch modes           */
/*****************************************************************/



inline bool IsPriceMode(int mode)
{
    return((mode & VALUATION_CONTEXT::Price) == VALUATION_CONTEXT::Price);
}



inline bool IsHedgeMode(int mode)
{
    return((mode & VALUATION_CONTEXT::Hedge) == VALUATION_CONTEXT::Hedge);
}


/*****************************************************************/
/*                 Masques de Lecture des correction modes       */
/*****************************************************************/





/*****************************************************************/
/*          Masques de Lecture portefeuilles de calibration      */
/*****************************************************************/


inline bool IsNonePortfolio(int Type)
{
    return((Type & CALIB_PORTFOLIO_MODE::None) == CALIB_PORTFOLIO_MODE::None);
}



inline bool IsStripPortfolio(int Type)
{
    return((Type & CALIB_PORTFOLIO_MODE::Strip) == CALIB_PORTFOLIO_MODE::Strip);
}


inline bool IsLevel1(int Type)
{
    return((Type & CALIB_PORTFOLIO_MODE::Level1) == CALIB_PORTFOLIO_MODE::Level1);
}

inline bool IsControl(int Type)
{
    return((Type & CALIB_PORTFOLIO_MODE::Control) == CALIB_PORTFOLIO_MODE::Control);
}

inline bool IsEarly(int Type)
{
    return((Type & CALIB_PORTFOLIO_MODE::Early) == CALIB_PORTFOLIO_MODE::Early);
}

inline bool IsUnderlying(int Type)
{
    return((Type & CALIB_PORTFOLIO_MODE::Under) == CALIB_PORTFOLIO_MODE::Under);
}


inline bool IsHedge(int Type)
{
    return((Type & CALIB_PORTFOLIO_MODE::Hedge) == CALIB_PORTFOLIO_MODE::Hedge);
}


inline bool IsSort(int Type)
{
    return((Type & CALIB_PORTFOLIO_MODE::Sort) == CALIB_PORTFOLIO_MODE::Sort);
}


int IsCallOrPut(int Type);









/*****************************************************************/
/*          Masques de Lecture des payoffs modes                 */
/*****************************************************************/


inline bool IsPayoffPayoff(int Type)
{
    return ((Type & PAYOFF_MODE::Payoff) == PAYOFF_MODE::Payoff);
}


inline bool IsPayoffProbas(int Type)
{
    return ((Type & PAYOFF_MODE::Probas) == PAYOFF_MODE::Probas);
}

inline bool IsPayoffCorrection(int Type)
{
    return ((Type & PAYOFF_MODE::Correction) == PAYOFF_MODE::Correction);
}






/***************************************************************/
/*                                                             */
/*                   AUTO_MODE for AUTO_CALIB                  */
/*                                                             */
/***************************************************************/




/***************************************************************/
/*                                                             */
/*                   MAV MODE                                  */
/*                                                             */
/***************************************************************/



inline bool IsMAVFill(int mode)
{
    return ((mode & MAV_MODE::Fill) == MAV_MODE::Fill);
}

inline bool IsMAVResize(int mode)
{
    return ((mode & MAV_MODE::Resize) == MAV_MODE::Resize);
}

inline bool IsMAVAlloc(int mode)
{
    return ((mode & MAV_MODE::Alloc) == MAV_MODE::Alloc);
}



/***************************************************************/
/*                                                             */
/*                   EXPECTATION MODE                          */
/*                                                             */
/***************************************************************/


inline bool IsExpectationSupEnv(int mode)
{
    return((mode & EXPECTATION_MODE::envelopeSup) == EXPECTATION_MODE::envelopeSup);
}

inline bool IsExpectationForceExe(int mode)
{
    return((mode & EXPECTATION_MODE::forceExe) == EXPECTATION_MODE::forceExe);
}

inline bool IsExpectationForceHold(int mode)
{
    return((mode & EXPECTATION_MODE::forceHold) == EXPECTATION_MODE::forceHold);
}

inline bool IsExpectationPathDepCorr(int mode)
{
    return((mode & EXPECTATION_MODE::pathDepCorr ) == EXPECTATION_MODE::pathDepCorr);
}

inline bool IsExpectationBS(int mode)
{
    return((mode & EXPECTATION_MODE::BS ) == EXPECTATION_MODE::BS);
}



/***************************************************************/
/*                                                             */
/*                   XCCY MODE                                 */
/*                                                             */
/***************************************************************/

inline bool IsModeBasisDiscount(int mode)
{
    return ((mode & XCCY_MODE::basisDiscount) == XCCY_MODE::basisDiscount);
}

inline bool IsModeTranslate(int mode)
{
    return ((mode & XCCY_MODE::translate) == XCCY_MODE::translate);
}




#endif
/*---------------------------------------------------------------*/
/*--- End Of File ---*/
