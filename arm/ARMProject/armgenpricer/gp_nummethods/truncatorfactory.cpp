/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file truncatorfactory.cpp
 *
 *  \brief
 *
 *	\author  JM Prie E Benhamou R. Guillemot
 *	\version 1.0
 *	\date October 2005
 */


#include "gpnummethods/truncatorfactory.h"

/// gpbase
#include "gpbase/singleton.h"
#include "gpbase/autocleaner.h"


/// gpinfra
#include "gpinfra/pricingstates.h"

/// gpnummethods
#include "gpnummethods/truncator.h"
#include "gpnummethods/argconvdefault.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_TruncatorFactoryData
///	Routine: toString
///	Returns: 
///	Action : Object dump
////////////////////////////////////////////////////
string ARM_TruncatorFactoryData::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;

	os << indent << "Truncator Factory Data\n";
	os << indent << "--------------------\n";
	os << indent << "\nTruncator=" << ARM_ArgConvReverse_TruncatorType.GetString(itsTruncatorType) << "\n";
    if(itsTruncatorDatas.size()>1)
        os << indent << "   NbStdDev="<<itsTruncatorDatas[0] << "\n";

	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_TruncatorFactoryImp
///	Routine: CreateTruncator
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_TruncatorBase* ARM_TruncatorFactoryImp::CreateTruncator(const ARM_TruncatorFactoryData& truncatorFactoryData)
{
    return CreateTruncator(
		truncatorFactoryData.GetNbDims(),
		truncatorFactoryData.GetTruncatorType(), 
		truncatorFactoryData.GetTruncatorDatas());
}

////////////////////////////////////////////////////
///	Class  : ARM_TruncatorFactoryImp
///	Routine: CreateScheduler
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_TruncatorBase* ARM_TruncatorFactoryImp::CreateTruncator(
		int nbDims,
        int truncatorType, 
		const std::vector<double>& truncatorDatas)
{
    ARM_TruncatorBase* truncator        = NULL;

    /// Truncator creation
    size_t nbDatas = truncatorDatas.size();
    switch(truncatorType)
    {
    case ARM_TruncatorBase::StandardDeviation:
        if(truncatorDatas.size()<1)
        {
			if(nbDims==1)
                truncator = new ARM_Truncator1D;
            else
                truncator = new ARM_TruncatorND(nbDims);
		}
        else if(nbDatas==1 || nbDims==1)
        {
            if(nbDims > 1)
                /// MaxStdDev are all equal to the input parameter
			    truncator = new ARM_TruncatorND(nbDims,truncatorDatas[0]);
            else
			    truncator = new ARM_Truncator1D(truncatorDatas[0]);
		}
        else if(nbDatas== nbDims)
        {
            /// MaxStdDev are set to the input parameter vector
			truncator = new ARM_TruncatorND(truncatorDatas);
        }
		else 
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : cannot interpret truncator data");
		}

        break;

    case ARM_TruncatorBase::ArrowDebreu:
        if(nbDatas<2)
        {
			if(nbDims==1)
                truncator = new ARM_Truncator1DArrowDebreu;
            else
                truncator = new ARM_TruncatorNDArrowDebreu(nbDims);
		}
        else if(nbDatas==4)
        {
            if(nbDims>1)
			    truncator = new ARM_TruncatorNDArrowDebreu(nbDims,
                                truncatorDatas[0],static_cast<size_t>(truncatorDatas[1]),
                                truncatorDatas[2],truncatorDatas[3] );
            else
			    truncator = new ARM_Truncator1DArrowDebreu(
                                truncatorDatas[0],static_cast<size_t>(truncatorDatas[1]),
                                truncatorDatas[2],truncatorDatas[3] );
		}
        else if(nbDatas== 2*nbDims+2)
        {
			std::vector<double> stdDevs(nbDims);
			ARM_IntVector maxMaxIndex(nbDims);
			for(size_t i=0;i<nbDims;++i)
			{
				stdDevs[i]=truncatorDatas[i];
				maxMaxIndex[i]=static_cast<int>(truncatorDatas[nbDims+i]);
			}

            truncator = new ARM_TruncatorNDArrowDebreu(
				stdDevs,maxMaxIndex,truncatorDatas[2*nbDims],
                truncatorDatas[2*nbDims+1]);
		}
        break;

    default:
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : unknown truncator type");
    }

	/// return the result
    return truncator;
}


ARM_SingletonHolder<ARM_TruncatorFactoryImp> ARM_TruncatorFactory;


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/