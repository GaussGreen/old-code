/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \file hwfactory.cpp
 *  \brief
 *	\author  E Ezzine
 *	\version 1.0
 *	\date July 2005
 *
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/hwfactory.h"

///gpbase
#include "gpbase/singleton.h"
#include "gpbase/ostringstream.h"
#include "gpbase/numericconstant.h"

///gpmodels
#include "gpmodels/HW1F.h"
#include "gpmodels/HW2F.h"
#include "gpmodels/ModelParamsHW1F.h"
#include "gpmodels/ModelParamsHW2F.h"

///gpinfra
#include "gpinfra/modelparamutil.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_HWModelFactoryImp
///	Routine: CreateHWModel
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_HullWhite* ARM_HWModelFactoryImp::CreateHWModel(const ARM_ZeroCurvePtr& zc,
        const ARM_ModelParamVector& params)
{

    size_t sizeparam = params.size();
    if(sizeparam>2)
        return  new ARM_HullWhite2F(zc, &(ARM_ModelParamsHW2FStd(params)));
	else if(sizeparam == 2)
	{
        ARM_ModelParamVector::const_iterator foundMRS =  (CC_NS( std, find_if) ( params.begin(), params.end(), 
                                FindModelParamWEnumUnaryVersion( ARM_ModelParamType::MeanReversion )));
        if(foundMRS != params.end())
        {
            if((*foundMRS)->size()>1)
                return  new ARM_HullWhite1F(zc, &ARM_ModelParamsHW1FExt(params));
            else if ((*foundMRS)->size() == 1)
                return  new ARM_HullWhite1F(zc, &ARM_ModelParamsHW1FStd(params));
            else
                ARM_THROW( ERR_INVALID_ARGUMENT, "Mean Reversion model param is empty" );
        }
        else
            ARM_THROW( ERR_INVALID_ARGUMENT, "Building HW Model needs a Mean reversion model param, please advice" );
    }
	else
	{
		CC_Ostringstream os;
		os << " HW Model consctuction needs a vector of modelparams(at least 2), please advice" ;
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
}

ARM_SingletonHolder<ARM_HWModelFactoryImp> ARM_HWModelFactory;


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
