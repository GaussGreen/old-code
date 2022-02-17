/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file modelfitter.h
 *
 *  \brief base class for all model fitter
 *
 *	\author  E.M Ezzine
 *	\version 1.0
 *	\date November 2003
 */


#ifndef _INGPCALIB_MODELFITTER_H
#define _INGPCALIB_MODELFITTER_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix

/// gpbase
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/timer.h"
#include "gpbase/numericconstant.h"
#include "gpbase/rootobject.h"
#include "gpbase/warning.h"
#include "gpbase/warningkeeper.h"

//gpinfra
#include "gpinfra/calibdirection.h"

/// gpcalib
#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_PricingModel;


class ARM_ModelFitter : public ARM_RootObject, public ARM_Timer, public ARM_WarningKeeper
{    
private:
	bool itsGetDetails;				/// whether we give detail about the method!
	string* itsFileName;		/// corresponding fileName
    ARM_PricingModel* itsPricingModel;		/// not cloned because shared object
    ARM_StdPortfolioPtr itsPortfolio;
    ARM_ModelParamVector itsCalibParam;
	ARM_ParamTypeVector  itsCalibParamsType;
    ARM_VanillaArgVector itsArgsVector;
    ARM_CalibDirection itsCalibDirection;
    ARM_MktTargetType itsTargetType;

    /// to combine  the limped or successes calibration
    ARM_ModelFitterPtr itsLinkedModelFitter;
    ARM_ModelFitterPtr itsPreviousModelFitter;

    /// for model fitter, there is a vector of single errors
    CC_IS_MUTABLE ARM_GP_MatrixPtr itsError;

	/// for multi-factor model
	size_t itsFactorNb; 

	/// for iterative Calibration
	size_t itsNbCalibIter;

    /// max iter for optimiser routine
    size_t itsMax_iter;

    void CleanUp();
    void CopyNoCleanUp(const ARM_ModelFitter& rhs );

	void Validate();
public:
    ARM_ModelFitter( ARM_PricingModel* model,
        const ARM_StdPortfolioPtr portfolio,
        const ARM_ModelParamVector& calibParam,
        ARM_ModelFitterPtr& linkedModelFitter	= ARM_ModelFitterPtr(NULL ),
        ARM_ModelFitterPtr& previousModelFitter = ARM_ModelFitterPtr(NULL ),
        const size_t max_iter					= ARM_NumericConstants::ARM_GP_MAX_ITER ,
		bool getDetails							= false,
		size_t factorNb							= 0, 
		size_t NbCalibIter						= 1,
        ARM_MktTargetType						= ARM_CalibrationTarget::PriceTarget);

	ARM_ModelFitter();
    ARM_ModelFitter& operator = (const ARM_ModelFitter& rhs);
    ARM_ModelFitter(const ARM_ModelFitter& rhs);
    virtual ~ARM_ModelFitter();   

    /// the core of the calibration
	virtual void Calibrate() = 0;

	/// std method (non virtual so not to be redefined!)
	void SetCalibParamToInitialValue();
    void PostProcessing();
	void DoCalibrationProcess();
	void SetUpError() const;
    
    /// accessors
    inline ARM_PricingModel* GetPricingModel() const            { return itsPricingModel;           } 
	inline void SetPricingModel(ARM_PricingModel* model)        {  itsPricingModel = model;         }
    inline const ARM_StdPortfolioPtr GetPortfolio() const       { return itsPortfolio;              }
    inline void SetPortfolio(ARM_StdPortfolioPtr Pf)            { itsPortfolio = Pf;                }                                                                                                                        
    inline ARM_ModelParamVector GetCalibParams() const          { return itsCalibParam;             }
    inline void SetCalibParams(ARM_ModelParamVector CalibParams){ itsCalibParam = CalibParams;      }
                                                                                        
    ARM_ModelParamVector::const_iterator FindCalibParamWType( int type ) const;                                         
    inline ARM_ModelParamVector::const_iterator UnknownCalibParamIterator() const   { return itsCalibParam.end();       }
                                                                                                            
    inline const ARM_ModelParam* GetCalibParam(int index = 0) const		{ return itsCalibParam[index];          }
    inline ARM_ModelParam* GetCalibParam(int index = 0)					{ return itsCalibParam[index];          }
    inline void SetCalibParam(ARM_ModelParam* calibparam ,int index)    { itsCalibParam[index] = calibparam;    }
    inline ARM_ModelFitterPtr GetLinkedModelFitter() const              { return itsLinkedModelFitter;          }
    inline ARM_ModelFitterPtr GetPreviousModelFitter() const            { return itsPreviousModelFitter;        }
    inline ARM_VanillaArgVector GetVanillaArgVector() const             { return itsArgsVector;                 }
    void SetVanillaArgVector(ARM_VanillaArgVector argsVector); 
	inline ARM_ParamTypeVector GetCalibParamsType() const               { return itsCalibParamsType;		    }
    inline ARM_MktTargetType GetTargetType() const                      { return itsTargetType;                 }
    
    inline ARM_VanillaArg* GetVanillaArgVector(int index) const         { return itsArgsVector[index];          }
    inline void SetVanillaArgVector(ARM_VanillaArg* arg,int index)      { itsArgsVector[index] = arg;           }
    inline const size_t GetMax_Iter() const                             { return itsMax_iter;                   }
    inline void SetMax_Iter(size_t max_iter)                            { itsMax_iter = max_iter;               }
    inline const ARM_GP_MatrixPtr GetError() const                      { return itsError;                      }
    inline bool GetDetails() const                                      { return itsGetDetails;                 }
    inline ARM_CalibDirection GetCalibDirection() const                 { return itsCalibDirection;             }
    inline void SetCalibDirection(ARM_CalibDirection calibDirection)    { itsCalibDirection = calibDirection;   }
	inline size_t GetFactorNb() const									{ return itsFactorNb;                   }
	inline size_t GetNbCalibIter() const								{ return itsNbCalibIter;                }
	inline void SetNbCalibIter(int nbCalibIter)							{ itsNbCalibIter = nbCalibIter;			}

	virtual void setStartAndEndCalibrationTimes( double startTime, double endTime){ ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": setStartAndEndCalibrationTimes not implemented!"); }


	virtual void Initialise(){}
	inline string* GetFileName() const									{ return itsFileName;	                }
    /// to force redefinition
	virtual void StoreDetailsInFileIfRequired(); 
	virtual string detailedString( const string& indent = "" ) const ;

	void AddModelFitterWarning( const ARM_ModelFitterPtr& modelFitter );
};



CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
