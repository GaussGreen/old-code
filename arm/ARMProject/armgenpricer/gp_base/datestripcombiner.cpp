/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file datestripcombiner.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#include "gpbase/datestripcombiner.h"
#include "gpbase/ostringstream.h"
#include "gpbase/stringmanip.h"
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/gpvector.h"

#include <iomanip>
#include <algorithm>

CC_BEGIN_NAMESPACE( ARM )

/// linkt this static constant to the one of the datestrip
const double ARM_DateStripCombiner::DateStripCombiner_BlankData = ARM_DateStrip::DateStripCombiner_BlankData;

/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStripCombiner
///	Routine: Constructor
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner::~ARM_DateStripCombiner( )
{}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStripCombiner
///	Routine: Constructor
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner::ARM_DateStripCombiner( const ARM_DateStripVector& dateStrips, 
	const string& mergeDateFuncString )
:	ARM_RootObject(), itsContents(dateStrips.size()), itsMergeData()
{
	/// first copies all the dateStrip cloning them!
	for(size_t i=0; i<dateStrips.size(); ++i )
		itsContents[i] = ARM_DateStripPtr( static_cast<ARM_DateStrip*>( dateStrips[i]->Clone() ) );
	
	/// mergeData
	MergeData( GetMergeFuncFromString(mergeDateFuncString ) );
	
	CC_ARM_SETNAME(ARM_DATESTRIP_COMBINER);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStripCombiner
///	Routine: MergeData
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////
void ARM_DateStripCombiner::MergeData( const MergeFunc& func ) 
{
	/// first merge the data and insert blank afterwards in the other dateStrips
	if( itsContents.size() == 1 )
	{
		std::vector<double>* result = ((*itsContents[0]).*func)();
		itsMergeData = ARM_VectorPtr( static_cast<std::vector<double>*>( result ) );
	}
	else
	{
		/// 1) create the mergeData
		size_t totalSize = 0;
		size_t i;
		for(i=0;i<itsContents.size(); ++i )
			totalSize += ((*itsContents[i]).*func)()->size();
		
		vector<double> tmpMergeData(totalSize);
		CC_NS(std,vector)<double>::iterator mergeDataIter = tmpMergeData.begin();
		CC_NS(std,vector)<double>::iterator iter,begin,end;

		for(i=0;i<itsContents.size(); ++i )
		{
			begin	= ((*itsContents[i]).*func)()->begin();
			end		= ((*itsContents[i]).*func)()->end();
			for(iter=begin; iter<end; ++iter, ++mergeDataIter )
				*mergeDataIter = *iter;
		}

#if defined( __GP_STRICT_VALIDATION )
		if( mergeDataIter-tmpMergeData.begin() != totalSize )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "totalSize different from real size" );
#endif

		CC_NS(std,sort)(tmpMergeData.begin(),tmpMergeData.end());
		CC_NS(std,vector)<double>::iterator pos = CC_NS(std,unique)(tmpMergeData.begin(),tmpMergeData.end());
		tmpMergeData.resize( pos-tmpMergeData.begin() );
		size_t mergeDataSize = tmpMergeData.size();
		itsMergeData = ARM_VectorPtr( new std::vector<double>( tmpMergeData.begin(), tmpMergeData.end() ) );
	
		/// 2) add DateStripCombiner_BlankData to all dateStrip if necessary!
		size_t j,k;
		for(i=0;i<itsContents.size(); ++i )
		{
			/// create data
			std::vector<double>* newFlowStartDates	= new std::vector<double>(mergeDataSize);
			std::vector<double>* newFlowEndDates		= new std::vector<double>(mergeDataSize);
			std::vector<double>* newFwdStartDates	= new std::vector<double>(mergeDataSize);
			std::vector<double>* newFwdEndDates		= new std::vector<double>(mergeDataSize);
			std::vector<double>* newResetDates		= new std::vector<double>(mergeDataSize);
			std::vector<double>* newPaymentDates		= new std::vector<double>(mergeDataSize);
			std::vector<double>* newInterestDays		= new std::vector<double>(mergeDataSize);
			std::vector<double>* newInterestTerms	= new std::vector<double>(mergeDataSize);
	
			std::vector<double>* initialVector	= ((*itsContents[i]).*func)();
			size_t 	initialVectorSize	= initialVector->size();
			

			for(j=0,k=0; j<mergeDataSize; ++j )
			{	
				/*if(k<initialVectorSize && initialVector->Elt(k) == itsMergeData->Elt(j) )
				{
					newFlowStartDates->[j]	= (*itsContents[i]).GetFlowStartDates()->Elt(k);
					newFlowEndDates->Elt(j)		= (*itsContents[i]).GetFlowEndDates()->Elt(k);
					newFwdStartDates->Elt(j)	= (*itsContents[i]).GetFwdRateStartDates()->Elt(k);
					newFwdEndDates->Elt(j)		= (*itsContents[i]).GetFwdRateEndDates()->Elt(k);
					newResetDates->Elt(j)		= (*itsContents[i]).GetResetDates()->Elt(k);
					newPaymentDates->Elt(j)		= (*itsContents[i]).GetPaymentDates()->Elt(k);
					newInterestDays->Elt(j)		= (*itsContents[i]).GetInterestDays()->Elt(k);
					newInterestTerms->Elt(j)	= (*itsContents[i]).GetInterestTerms()->Elt(k);
					++k;
				}
				else
				{
					newFlowStartDates->Elt(j)	= ARM_DateStripCombiner::DateStripCombiner_BlankData;
					newFlowEndDates->Elt(j)		= ARM_DateStripCombiner::DateStripCombiner_BlankData;
					newFwdStartDates->Elt(j)	= ARM_DateStripCombiner::DateStripCombiner_BlankData;
					newFwdEndDates->Elt(j)		= ARM_DateStripCombiner::DateStripCombiner_BlankData;
					newResetDates->Elt(j)		= ARM_DateStripCombiner::DateStripCombiner_BlankData;
					newPaymentDates->Elt(j)		= ARM_DateStripCombiner::DateStripCombiner_BlankData;
					newInterestDays->Elt(j)		= ARM_DateStripCombiner::DateStripCombiner_BlankData;
					newInterestTerms->Elt(j)	= ARM_DateStripCombiner::DateStripCombiner_BlankData;
				}*/
			}

			/// 3) update the corresponding datestrip!
			(*itsContents[i]).SetFlowStartDatesNoClone( newFlowStartDates );
			(*itsContents[i]).SetFlowEndDatesNoClone( newFlowEndDates );
			(*itsContents[i]).SetFwdRateStartDatesNoClone( newFwdStartDates );
			(*itsContents[i]).SetFwdRateEndDatesNoClone( newFwdEndDates );
			(*itsContents[i]).SetResetDatesNoClone( newResetDates );
			(*itsContents[i]).SetPaymentDatesNoClone( newPaymentDates );
			(*itsContents[i]).SetInterestDaysNoClone( newInterestDays );
			(*itsContents[i]).SetInterestTermsNoClone( newInterestTerms );
		}
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStripCombiner
///	Routine: Copy constructor
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////

ARM_DateStripCombiner::ARM_DateStripCombiner( const ARM_DateStripCombiner& rhs )
:	ARM_RootObject(rhs), 
	itsContents( rhs.itsContents ),
	itsMergeData( rhs.itsMergeData )
{}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStripCombiner
///	Routine: Operator=
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////

ARM_DateStripCombiner& ARM_DateStripCombiner::operator=( const ARM_DateStripCombiner& rhs )
{

	if( this != &rhs )
	{
		ARM_RootObject::operator =(rhs);
		itsContents		= rhs.itsContents;
		itsMergeData	= rhs.itsMergeData;
	}
	return *this;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStripCombiner
///	Routine: GetMergeFuncFromString
///	Returns: 
///	Action : function to return the appropriate function to use for the merge
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner::MergeFunc ARM_DateStripCombiner::GetMergeFuncFromString( 
	const string& mergeDateFuncString ) const
{
	/// convert to upper case
	string mergeDateFuncStringUpper = stringGetUpper( mergeDateFuncString );

// FIXMEFRED: mig.vc8 (22/05/2007 15:50:36): function pointer explicit cast

	if( mergeDateFuncStringUpper == "STARTDATE" )
		return &ARM_DateStrip::GetFlowStartDates;
	if( mergeDateFuncStringUpper == "ENDDATE" )
		return &ARM_DateStrip::GetFlowEndDates;
	if( mergeDateFuncStringUpper == "RESETDATE" )
		return &ARM_DateStrip::GetResetDates;
	if( mergeDateFuncStringUpper == "PAYMENTDATE" )
		return &ARM_DateStrip::GetPaymentDates;
	if( mergeDateFuncStringUpper == "FWDSTARTDATE" )
		return &ARM_DateStrip::GetFwdRateStartDates;
	if( mergeDateFuncStringUpper == "FWDENDDATE" )
		return &ARM_DateStrip::GetFwdRateEndDates;
	/// other cases
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Unknown type for merge" );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStripCombiner
///	Routine: PrintOneDateStripDoubleElem
///	Returns: 
///	Action : print a double elem that is a date
/////////////////////////////////////////////////////////////////
string ARM_DateStripCombiner::PrintOneDateStripDateElem( double elem )
{
	if(elem != ARM_DateStripCombiner::DateStripCombiner_BlankData )
		return ARM_Date(elem).toString();
	else
		return "    ";
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStripCombiner
///	Routine: PrintOneDateStripDoubleElem
///	Returns: 
///	Action : print a double elem that is a double
/////////////////////////////////////////////////////////////////
string ARM_DateStripCombiner::PrintOneDateStripDoubleElem( double elem )
{
	char msg[10];
	
	if(elem != ARM_DateStripCombiner::DateStripCombiner_BlankData )
		sprintf( msg, "%f", elem );
	else
		sprintf( msg, "     ");
	return string(msg);
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStripCombiner
///	Routine: toString
///	Returns: 
///	Action : stringify the object
/////////////////////////////////////////////////////////////////
string ARM_DateStripCombiner::toString(const string& indent, const string& nextIndent) const
{
	size_t i,j;

	CC_Ostringstream os;
	os << "DateStrip Combiner\n";
	
	os << " Row Nb\t Merge Data\t ";

	for( i=0; i<itsContents.size(); ++i )
		os << "Start Dates\t End Dates\t Fixing Dates\t Fwd Start Dates Fwd End Dates\tPayment Dates\t Interest Days\t Interest Terms\t";
	os << "\n";

	/*for( i=0; i<itsMergeData->size(); ++i )
	{
		os << " " << CC_NS(std,left) << i+1 << "\t " << ((ARM_Date) (*itsMergeData)[i]).toString() << "\t";

		for( j=0; j<itsContents.size(); ++j )
		{
			os << CC_NS(std,setw)(14) << ARM_DateStripCombiner::PrintOneDateStripDateElem( itsContents[j]->GetFlowStartDates()->Elt(i) ) << "\t ";
			os << CC_NS(std,setw)(14) << ARM_DateStripCombiner::PrintOneDateStripDateElem( itsContents[j]->GetFlowEndDates()->Elt(i) ) << "\t ";
			os << CC_NS(std,setw)(14) << ARM_DateStripCombiner::PrintOneDateStripDateElem( itsContents[j]->GetResetDates()->Elt(i) ) << "\t ";
			os << CC_NS(std,setw)(14) << ARM_DateStripCombiner::PrintOneDateStripDateElem( itsContents[j]->GetFwdRateStartDates()->Elt(i) ) << "\t ";
			os << CC_NS(std,setw)(14) << ARM_DateStripCombiner::PrintOneDateStripDateElem( itsContents[j]->GetFwdRateEndDates()->Elt(i) ) << "\t ";
			os << CC_NS(std,setw)(14) << ARM_DateStripCombiner::PrintOneDateStripDateElem( itsContents[j]->GetPaymentDates()->Elt(i) ) << "\t ";

			/// double part!
			os << CC_NS(std,setw)(14) << ARM_DateStripCombiner::PrintOneDateStripDoubleElem( itsContents[j]->GetInterestDays()->Elt(i) ) << "\t ";
			os << CC_NS(std,setw)(14) << ARM_DateStripCombiner::PrintOneDateStripDoubleElem( itsContents[j]->GetInterestTerms()->Elt(i) )<< "\t ";
		}

		/// end of line mark
		os << CC_NS(std,endl);
	}

	for( i=0; i<itsContents.size(); ++i )
    {
		os << "\n\n========== > Date Strip Nb: " << i+1 << "\n";
		os << itsContents[i]->toString(indent,nextIndent);
	}*/

	return os.str();
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStripCombiner
///	Routine: Clone, View
///	Returns: ARM_Object*
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_DateStripCombiner::Clone() const
{
	return new ARM_DateStripCombiner( *this );
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/



