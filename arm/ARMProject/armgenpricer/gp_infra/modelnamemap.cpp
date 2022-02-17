/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file modelnamemap.cpp
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

#include "gpinfra/modelnamemap.h"

/// gpbase
#include "gpbase/ostringstream.h"

/// gpinfra
#include "gpinfra/pricingmodelir.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingmodeltype.h"

/// STL
#include <iomanip>		/// necessary for setw
#include <vector>		/// necessary for vector
#include <algorithm>	/// necessary to find quickly things
CC_USING( std::find );
CC_USING( std::sort );

/// ARM Kernel
#include "glob/expt.h"		/// necessary for exception throwing


CC_BEGIN_NAMESPACE( ARM )



////////////////////////////////////////////////////
/// Class  : ARM_ModelNameMap
///	Routine: constructor
///	Returns: built object
///	Action : 
////////////////////////////////////////////////////
ARM_ModelNameMap::ARM_ModelNameMap(
	const ARM_StringVector& Names, 					/// Names of the models.
	const vector< ARM_PricingModelPtr >& Models,	/// The models.
	const ARM_StringVectorVector& OtherModelNames	/// The names of other models for diffusion
)
:
	ARM_RootObject(),
	itsData(),
	itsSortedData()
{
    SetName(ARM_MODELNAMEMAP);

	/// some validation on the size!
	if( Names.size() != Models.size() )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
		"ARM_ModelNameMap(): name and model arrays are different sizes... Please advise" );
	
	if( Names.size() != OtherModelNames.size() )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
		"ARM_ModelNameMap(): name and other model names arrays are different sizes... Please advise" );
	
	
	ARM_StringVector::const_iterator 
		Name    = Names.begin(),
		NameEnd = Names.end();
	
    CC_NS( std, vector )< ARM_PricingModelPtr >::const_iterator	
		Model = Models.begin();
	
	ARM_StringVectorVector::const_iterator
		OtherModelNameRow = OtherModelNames.begin();
	
	/// size with the appropriate number the vector of DataElem
	itsData = DataT( Names.size() );
	
	DataT::iterator
		Data = itsData.begin();
	

	/// construct the single DataElem
	while( Name < NameEnd )
	{
		ARM_IntVector
			OtherModelRowIndices( OtherModelNameRow->size() );
		
		ARM_StringVector::const_iterator
			OtherModelNameCol    = OtherModelNameRow->begin(),
			OtherModelNameColEnd = OtherModelNameRow->end();
		
		ARM_IntVector::iterator
			OtherModelRowIndex = OtherModelRowIndices.begin();
	
		/// find the corresponding index for each otherModel
		/// from the name
		while( OtherModelNameCol < OtherModelNameColEnd )
		{
			/// Could be a bit more efficient if we build a map from
			/// names to rows, but these tables are so small it is not
			/// worth it.
			ARM_StringVector::const_iterator
				FoundName = find( Names.begin(), Names.end(), *OtherModelNameCol );
			
			/// validation that we do not include another model
			/// that does not exist!
			if( FoundName == Names.end() )
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				string( "Diffusion other model named \"" ) + *OtherModelNameCol 
				+ string( "\" not listed in model name table." ) );

			*OtherModelRowIndex = FoundName - Names.begin();
			
			/// increment everything
			++OtherModelNameCol;
			++OtherModelRowIndex;
		}
		
		*Data = DataElem( *Name, *Model, OtherModelRowIndices );
		Data->Model()->SetModelName( *Name );

		++Model;
		++Data;
		++Name;
		++OtherModelNameRow;
	}
	
	/// Mark those that are used as other models.
	/// and set the correct name!
	DataT::iterator
		DataEnd = itsData.end();
	
	Data = itsData.begin();
	
	while( Data < DataEnd )
	{
		ARM_IntVector::const_iterator
			OtherModelElemNumber    = Data->OtherModelRefNb().begin(),
			OtherModelElemNumberEnd = Data->OtherModelRefNb().end();
		
		while( OtherModelElemNumber != OtherModelElemNumberEnd )
			itsData[ *OtherModelElemNumber++ ].SetUsedAsOtherModel( true );

		++Data;
	}

	//
	// Set if model is used according to model type
	//
	Data = itsData.begin();
	
	while( Data < DataEnd )
	{
		if ( (Data->Model()->GetType() & MT_NON_STOCHASTIC_MODEL) == MT_NON_STOCHASTIC_MODEL )
			Data->SetUsedInPricing(false);
		++Data;
	}


	itsSortedData = itsData;
	sort( itsSortedData.begin(), itsSortedData.end() );
}


////////////////////////////////////////////////////
/// Class  : ARM_ModelNameMap
///	Routine: Copy Constructor, Desctructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_ModelNameMap::ARM_ModelNameMap(	const ARM_ModelNameMap &ModelNameMap)
:	ARM_RootObject(ModelNameMap), itsData(ModelNameMap.itsData), itsSortedData(ModelNameMap.itsSortedData)
{
    size_t nbData = itsData.size();
    for(size_t i=0;i<nbData;++i)
    {
        /// Clone each model
		itsData[i].Model() = ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(ModelNameMap.itsData[i].Model()->Clone()));
    }
	itsSortedData = itsData;
	sort( itsSortedData.begin(), itsSortedData.end() );
}


ARM_ModelNameMap::~ARM_ModelNameMap()
{}


////////////////////////////////////////////////////
/// Class  : ARM_ModelNameMap
///	Routine: Assignment operator
///	Returns:
///	Action : 
////////////////////////////////////////////////////
ARM_ModelNameMap & ARM_ModelNameMap::operator=(	const ARM_ModelNameMap& rhs)
{
	if( this != &rhs )
	{
		/// delegate to the internal object!
		ARM_RootObject::operator=(rhs);
		itsData = rhs.itsData;
		itsSortedData = rhs.itsSortedData;
	}
	return *this;
}

////////////////////////////////////////////////////
/// Class  : ARM_ModelNameMap
///	Routine: operator[]
///	Returns: The mapped element or an error if not found
///	Action : 
////////////////////////////////////////////////////
ARM_ModelNameMap::const_iterator ARM_ModelNameMap::operator[] ( const string& Name ) const
{	
	DataT::const_iterator
			Found = itsData.begin(),
			End   = itsData.end();

	while( Found != End )
	{
		if( Found->ModelName() == Name )
			return Found;
			
		++Found;
	}

    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
		"The model name is not found in the map" );
}

ARM_ModelNameMap::iterator ARM_ModelNameMap::operator[] ( const string& Name ) 
{	
	DataT::iterator
			Found = itsData.begin(),
			End   = itsData.end();

	while( Found != End )
	{
		if( Found->ModelName() == Name )
			return Found;
			
		++Found;
	}

    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
		ARM_USERNAME + ": model with name " + Name + " not found in the map" );
}

////////////////////////////////////////////////////
/// Class  : ARM_ModelNameMap
///	Routine: toString
///	Returns:
///	Action : dump in a string for easy debbuging
////////////////////////////////////////////////////

bool ARM_ModelNameMap::TestIfModelExisting ( const string& Name ) const
{	
	DataT::const_iterator
			Found = itsData.begin(),
			End   = itsData.end();

	while( Found != End )
	{
		if( Found->ModelName() == Name )
			return true;
			
		++Found;
	}

	return false;   
}


////////////////////////////////////////////////////
/// Class  : ARM_ModelNameMap
///	Routine: toString
///	Returns:
///	Action : dump in a string for easy debbuging
////////////////////////////////////////////////////
string ARM_ModelNameMap::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	/*const_iterator
		Iter,
		Begin= begin(),
		End  = end();

	os << "\n\n";
	os << indent << "ARM_ModelNameMap\n";
	os << indent << "----------------\n\n\n";

	os << indent << "---------------------------------------------\n";
	os << indent << "Quick View\n";
	os << indent << "---------------------------------------------\n\n";

	size_t i=0;
	os << "-" << size() << " Model(s)\n";
	for( Iter=Begin; Iter!=End; ++Iter, ++i )
	{
		CC_Ostringstream ost;
		os << indent << i+1 << ") Name: " <<  CC_NS(std,left) << CC_NS(std,setw)(25) << Iter->ModelName()  << " Ccy :";
		if( (*Iter).Model()->GetZeroCurve() != ARM_ZeroCurvePtr(NULL) )
			os << CC_NS(std,left) << CC_NS(std,setw)(5) << (*Iter).Model()->GetZeroCurve()->GetCurrencyUnit()->GetCcyName();
		else
			os << CC_NS(std,left) << CC_NS(std,setw)(5) << "NA";

		os << " Other Models[" << Iter->OtherModelRefNb().size() << "] = ";
		ost << "{";
		
		if( Iter->OtherModelRefNb().size() )
		{
			
			ARM_IntVector::const_iterator 
				OtherModel	  = Iter->OtherModelRefNb().begin(), 
				OtherModelEnd = Iter->OtherModelRefNb().end();

			while( OtherModel<OtherModelEnd )
			{
				ost << itsData[ *OtherModel ].ModelName();
				if( ++OtherModel != OtherModelEnd )
					ost << ",";
			}
		}
		ost <<"}";
		os << CC_NS(std,left) << CC_NS(std,setw)(35) << ost.str();
		os <<"-Used in pricing? " << ((*Iter).UsedInPricing()? "true  " : "false ")
			<< " -Used As Other Model? " << ( (*Iter).UsedAsOtherModel()? "true  " : "false ") << "\n";
	}

	os << indent << "\n\n\nDetailled View\n";
	os << indent << "----------------\n";

	i=0;
	for( Iter=Begin; Iter!=End; ++Iter, ++i )
	{
		os << i+1 << ") Model\n";
		os << indent << "Name: " <<  Iter->ModelName() << "\n" << indent << "Ccy ";
        if( (*Iter).Model()->GetZeroCurve() != ARM_ZeroCurvePtr(NULL) )
            os << (*Iter).Model()->GetZeroCurve()->GetCurrencyUnit()->GetCcyName() << "\n";
		else
			os << " Not Specified\n";
		os << indent << "Model Details : \n" << (*Iter).Model()->toString(indent+nextIndent,nextIndent) << "\n\n";
	}*/

	/// return the final result
	return os.str();
}

////////////////////////////////////////////////////
/// Class  : ARM_ModelNameMap
///	Routine: getNextUsedIter
///	Returns: const_iterator&
///	Action : skip for a given iterator unused model
////////////////////////////////////////////////////

ARM_ModelNameMap::const_iterator& ARM_ModelNameMap::getNextUsedIter( ARM_ModelNameMap::const_iterator& iter ) const
{
	/// increment at least by one.
	++iter;
	/// check if this is used!
	while( iter != itsData.end() && !(*iter).UsedInPricing() )
		++iter;
	return iter;
}


////////////////////////////////////////////////////
/// Class  : ARM_ModelNameMap
///	Routine: getNextUsedIter
///	Returns: iterator&
///	Action : skip for a given iterator unused model
////////////////////////////////////////////////////

ARM_ModelNameMap::iterator& ARM_ModelNameMap::getNextUsedIter( ARM_ModelNameMap::iterator& iter )
{
	/// increment at least by one.
	++iter;
	/// check if this is used!
	while( iter != itsData.end() && !(*iter).UsedInPricing() )
		++iter;
	return iter;
}

////////////////////////////////////////////////////
/// Class  : ARM_ModelNameMap
///	Routine: getSortedNextUsedIter
///	Returns: const_iterator&
///	Action : skip for a given iterator unused model
////////////////////////////////////////////////////

ARM_ModelNameMap::const_iterator& ARM_ModelNameMap::getSortedNextUsedIter( ARM_ModelNameMap::const_iterator& iter ) const
{
	/// increment at least by one.
	++iter;
	/// check if this is used!
	while( iter != itsSortedData.end() && !(*iter).UsedInPricing() )
		++iter;
	return iter;
}


////////////////////////////////////////////////////
/// Class  : ARM_ModelNameMap
///	Routine: getSortedNextUsedIter
///	Returns: iterator&
///	Action : skip for a given iterator unused model
////////////////////////////////////////////////////

ARM_ModelNameMap::iterator& ARM_ModelNameMap::getSortedNextUsedIter( ARM_ModelNameMap::iterator& iter )
{
	/// increment at least by one.
	++iter;
	/// check if this is used!
	while( iter != itsSortedData.end() && !(*iter).UsedInPricing() )
		++iter;
	return iter;
}



////////////////////////////////////////////////////
///	Class   : ARM_ModelNameMap
///	Routine : UsedModelsSize
///	Returns : size_t
///	Action  : computes the number of used models
////////////////////////////////////////////////////
size_t ARM_ModelNameMap::UsedModelsSize() const
{
	size_t i=0;
	for( const_iterator iter=begin(); iter!=end(); getNextUsedIter(iter), ++i )
		;
	return i;
}

///////////////////////////////////////////////////////
///	Class   : ARM_ModelNameMap::DataElem
///	Routine : operator<
///	Returns : bool
///	Action  : compares two DataElem (by model.GetType)
///////////////////////////////////////////////////////
bool ARM_ModelNameMap::DataElem::operator< ( const DataElem &rhs ) 
{ 
	return ( itsModel->GetType() > rhs.itsModel->GetType() );
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

