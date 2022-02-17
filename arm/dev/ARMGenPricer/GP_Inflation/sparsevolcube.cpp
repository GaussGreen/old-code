/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file sparsevolcube.h
 *
 *  \brief object to store the inflation vol data
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date August 2003
 */

#include "gpinflation/sparsevolcube.h"

/// gpbase
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/stringconvert.h"
#include "gpbase/gpvector.h"
#include "gpbase/gplinalgconvert.h"
#include "gpbase/comparisonfunctor.h"

//// kernel 
#include <crv/volcube.h>
#include <glob/expt.h>

CC_USING_NS( std, pair )
CC_USING_NS( std, vector )
//CC_USING_NS_T( std, vector, ARM_VolCurve* )
//CC_USING_NS_T( std, vector, CC_USING_BI_T( pair, int, int ) )


CC_BEGIN_NAMESPACE( ARM )

///////////////////////////////////////////
/// constant to do the conversion between
/// Zero coupon and expired Year to Year
///////////////////////////////////////////
const double ARM_SparseVolCube::spotTime = StringMaturityToYearTerm( "1d" );

///////////////////////////////////////////////
//// Sparse VolCube
/// the validation function make sure that input are making
/// sense and are correctly sorted
///////////////////////////////////////////////
void ARM_SparseVolCube::Validation(
	ARM_GP_Vector* dim1,
	ARM_GP_Vector* strikes,
	int dim1Type,
	int otherDimType )
{
	/// first do various validation tests
	if( dim1->sort() != (*dim1) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
		"Dim1 has to be sorted");
	
	/// first do various validation tests
	if( strikes->sort() != (*strikes) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
		"strikes has to be sorted");

	if( otherDimType == K_TENOR && 
		( dim1Type != K_TIMETOEXPIRY && dim1Type != K_TIMETOSTART ) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
		"otherDimType used Tenor while dim1Type is neither time to expiry or time to start");
	
	if( dim1Type == K_TENOR && otherDimType != K_TIMETOSTART )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
		"dim1Type used Tenor while otherDimType is not time to start");

	if( itsAsOf < itsLastKnownDate )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
		string( "Asof Date: " ) + itsAsOf.toString() + string( " smaller than lastKnownDate: " ) 
		+ itsLastKnownDate.toString() );
}


////////////////////////////////////////////////////
///	Class  : ARM_SparseVolCube
///	Routine: ARM_SparseVolCube
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_SparseVolCube::ARM_SparseVolCube(
	const ARM_Date& asOf,
	const ARM_Date& lastKnownDate,
	const string& indexName,
	int dim1Type,
	ARM_GP_Vector* dim1,
	int otherDimType,
	double otherDimValue,
	ARM_GP_Vector* strikes,
	ARM_Matrix* volatilities,
	int strikeType,
	int volType )
: itsAsOf( asOf ), itsLastKnownDate( lastKnownDate )
{
	Validation( dim1, strikes, dim1Type, otherDimType );

	switch( otherDimType )
	{
		
	// standard case
	/// easy case is that the sparse vol cube is correctly set up
	case K_TENOR:
		{
			int nbOfWrongTimeToStart = 0;
			/// no size as we do not know it
			/// will use thererefore push_back functions
			itsVols			= new vector < ARM_VolCurve* >;
			vector<double> underlyingsTmp;

			/// if we specify the vol in term of the time to expiry
			/// we need to convert into time To Start
			/// as the vol cube dimension is time To Start!
			if( dim1Type == K_TIMETOEXPIRY )
			{
				int size = dim1->size();
				int colSize = volatilities->GetNumCols();
				underlyingsTmp.reserve(size);
				
				for( int i=0; i<size; ++i )
				{
					if( dim1->Elt(i) - otherDimValue < ARM_SparseVolCube::spotTime )
					{
						/// expired vol are similar to the zero coupon vols
						/// in this case we take 1d maturities
						ARM_GP_Vector* maturities = new ARM_GP_Vector( 1, ARM_SparseVolCube::spotTime );
						
						/// and a tenor corresponding to the time to expiry
						/// forces it to be more than one Day!
						double tenor = dim1->Elt(i);
						if( tenor < ARM_SparseVolCube::spotTime )
							tenor = ARM_SparseVolCube::spotTime;

						ARM_Vector* valuesVect  = volatilities->GetRow(i);
						ARM_Matrix* data = new ARM_Matrix( 1.0, colSize, &(*valuesVect)[0] );
						delete valuesVect;

						ARM_Vector* maturitiesVec = To_pARM_Vector(maturities);
						delete maturities;
						ARM_Vector* strikesVec = To_pARM_Vector(strikes);

						ARM_VolCurve* curv = new ARM_VolLInterpol( 
							const_cast<ARM_Date&>(asOf), 
							maturitiesVec,
							strikesVec,
							data, 
							strikeType, 
							volType );
						
						itsVols->push_back( curv );
						underlyingsTmp.push_back( tenor );
						nbOfWrongTimeToStart++;

						delete maturitiesVec;
						delete strikesVec;
					}
					else
						dim1->Elt(i) -= otherDimValue;
				}
			}
			
			if( nbOfWrongTimeToStart )
			{
				int dim1Size = dim1->size();
				ARM_GP_Vector* newDim1 = new ARM_GP_Vector( dim1Size - nbOfWrongTimeToStart );
				
				/// first get the correct dim1 number
				int i;
				for( i=nbOfWrongTimeToStart; i<dim1Size; ++i )
					newDim1->Elt( i- nbOfWrongTimeToStart ) = dim1->Elt(i);
				delete dim1;
				dim1 = newDim1;
				
				/// second gets correct volatilities
				int colSize = volatilities->GetNumCols();
				ARM_Matrix* newVolatilities = new ARM_Matrix( dim1Size - nbOfWrongTimeToStart, colSize );
				
				int j;
				for( i=nbOfWrongTimeToStart; i<dim1Size; ++i )
					for( j=0; j<colSize; ++j )
						newVolatilities->Elt(i- nbOfWrongTimeToStart ,j) = volatilities->Elt(i,j);

				/// in the case of already defined volcurves
				/// we need to add these!!! so use the interface addVol
				ARM_Vector* tmpDim1 = To_pARM_Vector(dim1);
				delete dim1;

				ARM_VolCurve* curv = new ARM_VolLInterpol( 
					const_cast<ARM_Date&>(asOf), 
					tmpDim1,
					To_pARM_Vector(strikes), 
					newVolatilities, 
					strikeType, 
					volType );

				double tenor = otherDimValue;
				/// assign all the result of underlyings
				itsUnderlyings = new ARM_GP_Vector(underlyingsTmp.size(), &underlyingsTmp[0] );

				AddOneVolCurve( curv, tenor );
			}
			else
			{
				ARM_Vector* tmpDim1 = To_pARM_Vector(dim1);
				delete dim1;

				/// otherwise simply just add these
				ARM_VolCurve* curv = new ARM_VolLInterpol( 
					const_cast<ARM_Date&>(asOf), 
					tmpDim1,
					To_pARM_Vector(strikes), 
					static_cast<ARM_Matrix*>(volatilities->Clone()), 
					strikeType, 
					volType );
				itsVols->push_back( curv );
				underlyingsTmp.push_back( otherDimValue );

				/// assign all the result of underlyings
				itsUnderlyings = new ARM_GP_Vector(underlyingsTmp.size(), &underlyingsTmp[0] );
			}
			break;
		}
		
		
	/// in this case we specify in the other dim
	/// hence we need to create a series of vol curves
	case K_TIMETOSTART:
		{
			int size	= dim1->size();
			int colSize = volatilities->GetNumCols();
			
			ARM_GP_Vector* maturities = new ARM_GP_Vector(1.0, otherDimValue );
			itsVols			= new vector< ARM_VolCurve* >( size );
			itsUnderlyings = dim1;
			
			int i = 0;

			ARM_Vector* valuesVect = volatilities->GetRow(i);
			ARM_Matrix* data = new ARM_Matrix( 1.0, colSize, &(*valuesVect)[0] );
			delete valuesVect;

			(*itsVols)[i] = new ARM_VolLInterpol( 
				const_cast<ARM_Date&>(asOf), 
				To_pARM_Vector(maturities), 
				To_pARM_Vector(strikes), 
				data, 
				strikeType, 
				volType );
			
			/// start after one for the cloning of necessary data
			for( i=1; i<size; ++i )
			{
				/// the cloning is necessary to give
				/// to each object its own version of maturities and strikes
				/// so that they can delete objects properly
				ARM_Vector* valuesVect = volatilities->GetRow(i);
				ARM_Matrix* data = new ARM_Matrix( 1.0, colSize, &(*valuesVect)[0]);
				delete valuesVect;

				(*itsVols)[i] = new ARM_VolLInterpol( 
					const_cast<ARM_Date&>(asOf), 
					To_pARM_Vector(maturities), 
					To_pARM_Vector(strikes),
					data, 
					strikeType, 
					volType );
			}
			delete maturities;
			break;
		}

	/// time to expiry is not permitted for the other Dim type
	default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, 
			"Unknown other Dim type... be aware that time to expiry is not permitted for the other dim type");
	}
	
	SetName(ARM_SPARSE_VOL_CUBE);
}


////////////////////////////////////////////////////
///	Class  : ARM_SparseVolCube
///	Routine: intersectPos
///	Returns: vector< pairIntInt >
///	Action : routine to get the common indexes of the 
///				two ARM_GP_Vector
/// 
///				the result is a vector of pair of positions
///				for each common number 
///				hence first int is position in v1
///				second int is position in v2
///
///				ONE OF THE BIG ASSUMTION IS THAT V1 AND V2 ARE
///				SORTED IN INCREASING ORDER
//////////////////////////////////////////////////
vector< pairIntInt > ARM_SparseVolCube::intersectPos( const ARM_GP_Vector& v1, const ARM_GP_Vector& v2 )
{
	vector< pairIntInt > result;
	int i=0, j=0;
	int v1Size = v1.size();
	int v2Size = v2.size();

	for( ; i<v1Size; ++i )
	{
		// increase while v2 is smaller
		while(  j < v2Size  && v2.Elt(j) < v1.Elt(i) )
			++j;
		// are we at the end
		if( j == v2Size )
			break;
		else
			// otherwise are we equal
			if( fabs( v2.Elt(j) - v1.Elt(i) ) < PRECISION )
				result.push_back( pairIntInt( i, j ) );
	}
	return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_SparseVolCube
///	Routine: mergePos
///	Returns: vector< pairIntInt >
///	Action : routine to get the index and the corresponding 
///				vector for the merged elements
///
///				the result is a vector of pair of integer where
///					- the first integer is the position of the first or second vector
///					- the second integer is 
///						- 1 if only found in V1
///						- 2 if only found in V2
///						- 3+index of the second vector if found in both
///
///				ONE OF THE BIG ASSUMTION IS THAT V1 AND V2 ARE
///				SORTED IN INCREASING ORDER
//////////////////////////////////////////////////

vector< pairIntInt >  ARM_SparseVolCube::mergePos( const ARM_GP_Vector& v1, const ARM_GP_Vector& v2, ARM_GP_Vector& writeResult )
{
	/// use a temp vector because we want to use functionalities 
	/// of STL vector
	vector< pairIntInt > result;
	vector< double > tempWriteResult;
	tempWriteResult.reserve( v1.size() + v2.size() );

	int i = 0,j = 0;
	int v1Size = v1.size();
	int v2Size = v2.size();

	for( ; i< v1Size; ++i )
	{
		/// add all the elements from v2 lower than v1.Elt(i)
		while( j < v2Size && v2.Elt(j) < v1.Elt(i) )
		{
			/// copy from v2 
			/// put pos in V2 in the first int
			/// and 2 in the second int
			result.push_back( pairIntInt( j, 2 ) );
			tempWriteResult.push_back( v2.Elt(j) );
			++j;
		}

		if( j == v2Size )
		{
			/// copy from v1 
			/// put pos in V1 in the first int
			/// with 1 in the second int
			result.push_back( pairIntInt( i, 1 ) );
			tempWriteResult.push_back( v1.Elt(i) );
		}
		else
		{
			if( fabs( v2.Elt(j) - v1.Elt(i) ) < PRECISION )
			{
				/// v1 and v2 contains both the element

				/// copy using value in V1
				/// put pos in V1 in the first int
				/// and 3+pos in V2 in the second int
				
				result.push_back( pairIntInt( i, 3 + j ) );
				tempWriteResult.push_back( v1.Elt(i) );
			++j;
			}
			else
			{
				/// copy from v2 
				/// put pos in V2 in the first int
				/// and 2 in the second int
				result.push_back( pairIntInt( j, 2 ) );
				tempWriteResult.push_back( v2.Elt(j) );
			}
		}
	}


	/// we now need to check the elements from v2
	/// that could potentially be left
	for( ; j< v2Size; ++j )
	{
		if( fabs( v2.Elt(j) - v1.Elt(v1Size-1) ) < PRECISION )
		{
			/// v1 and v2 contains both the element

			/// copy using value in V1
			/// put pos in V1 in the first int
			/// and 3+pos in V2 in the second int
			
			result.push_back( pairIntInt(v1Size-1 , 3 + j ) );
			tempWriteResult.push_back( v1.Elt(v1Size-1) );
		}
		else
		{
			/// copy from v2 
			/// put pos in V2 in the first int
			/// and 2 in the second int
			result.push_back( pairIntInt( j, 2 ) );
			tempWriteResult.push_back( v2.Elt(j) );
		}
	}

	writeResult = ARM_GP_Vector( tempWriteResult );
	return result;
}



////////////////////////////////////////////////////
///	Class  : ARM_SparseVolCube
///	Routine: compare
///	Returns: void
///	Action : helper function
///				compares two vol curve to see that overlapping
///				data are similar
////////////////////////////////////////////////
void ARM_SparseVolCube::compare( ARM_VolCurve* nVolCrv, ARM_VolCurve* oVolCrv, double tenor )
{
	// this works only for ARM_VolLInterpol
	ARM_VolLInterpol* newVolCrv = dynamic_cast<ARM_VolLInterpol*>( nVolCrv );
	ARM_VolLInterpol* oldVolCrv = dynamic_cast<ARM_VolLInterpol*>( oVolCrv );

	if( !newVolCrv || !oldVolCrv )
		return ;

	/// to compare  we check that 
	///		the overlap between the old vol curve and the new one
	///		is the same
	
	/// if this is true, return true
	/// else return false
	
	/// a vol curve is composed of strike and maturities
	/// hence we need to see for each strike and maturities if this exists in the oldCurv
	ARM_Vector* newExpiries = newVolCrv->GetExpiryTerms();
	ARM_Vector* oldExpiries = oldVolCrv->GetExpiryTerms();
	
	/// get strikes
	ARM_Vector* newStrikes = newVolCrv->GetStrikes();
	ARM_Vector* oldStrikes = oldVolCrv->GetStrikes();

	/// get the intersectExpiries position
	vector< pairIntInt > intersectExpiries = intersectPos( To_ARM_GP_Vector(*newExpiries), To_ARM_GP_Vector(*oldExpiries) );
	vector< pairIntInt > intersectStrikes  = intersectPos( To_ARM_GP_Vector(*newStrikes ), To_ARM_GP_Vector(*oldStrikes ) );
	
	vector< pairIntInt >::const_iterator intersectExpiriesIter;
	vector< pairIntInt >::const_iterator intersectStrikesIter; 
	
	vector< pairIntInt >::const_iterator intersectExpiriesEnd= intersectExpiries.end();
	vector< pairIntInt >::const_iterator intersectStrikesEnd = intersectStrikes.end();
	
	ARM_Matrix* newVol	   = newVolCrv->GetVolatilities();
	ARM_Matrix* oldVol	   = oldVolCrv->GetVolatilities();

	/// once we have found the intersect position
	/// we need to compare the corresponding vol
	for( intersectExpiriesIter = intersectExpiries.begin();
	intersectExpiriesIter != intersectExpiriesEnd;
	++intersectExpiriesIter )
	{
		for( intersectStrikesIter = intersectStrikes.begin();
		intersectStrikesIter!= intersectStrikesEnd; 
		++intersectStrikesIter )
		{

			/// we do not allow to have different data
			/// vol matrices should be consistent
			if(	fabs(
				newVol->Elt( (*intersectExpiriesIter).second, (*intersectStrikesIter).second ) -
				oldVol->Elt( (*intersectExpiriesIter).first, (*intersectStrikesIter).first ) ) > PRECISION )
			{
				char msg[255];
				sprintf( msg, "the existing data are not the same as the new ones\n" );
				sprintf( msg, "inconsistency for tenor %f maturity %f and strike %f, old value %f, new value %f",
					tenor,
					oldExpiries->Elt( (*intersectExpiriesIter).first ), 
					oldStrikes->Elt( (*intersectStrikesIter).first ), 
					oldVol->Elt( (*intersectExpiriesIter).first, (*intersectStrikesIter).first ),
					newVol->Elt( (*intersectExpiriesIter).second, (*intersectStrikesIter).second ) );
				throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, msg );
			}
		}
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_SparseVolCube
///	Routine: compare
///	Returns: void
///	Action : helper function. This merge two vol curve 
////////////////////////////////////////////////////

ARM_VolCurve* ARM_SparseVolCube::merge( ARM_VolCurve* nVolCrv, ARM_VolCurve* oVolCrv, double tenor )
{
	// this works only for ARM_VolLInterpol
	ARM_VolLInterpol* newVolCrv = dynamic_cast<ARM_VolLInterpol* >( nVolCrv );
	ARM_VolLInterpol* oldVolCrv = dynamic_cast<ARM_VolLInterpol* >( oVolCrv  );

	// if it is not an ARM_VolLInterpol
	// don't do anything
	if( !newVolCrv || !oldVolCrv )
		return NULL;

	/// a vol curve is composed of maturities and strikes
	/// hence we need to see for each strike and maturities if this exists in the oldCurv
	ARM_Vector* newExpiries = newVolCrv->GetExpiryTerms();
	ARM_Vector* oldExpiries = oldVolCrv->GetExpiryTerms();

	ARM_Vector* newStrikes = newVolCrv->GetStrikes();
	ARM_Vector* oldStrikes = oldVolCrv->GetStrikes();
	
	ARM_GP_Vector* mergedExpiriesVec = new ARM_GP_Vector;
	ARM_GP_Vector* mergedStrikesVec  = new ARM_GP_Vector;
	
	/// we find the common expiries and strikes
	/// the routine mergePos fill as well the ARM_GP_Vector merged... blahblah
	vector< pairIntInt > mergedExpiries = mergePos( To_ARM_GP_Vector(*newExpiries), To_ARM_GP_Vector(*oldExpiries), *mergedExpiriesVec );
	vector< pairIntInt > mergedStrikes  = mergePos( To_ARM_GP_Vector(*newStrikes ), To_ARM_GP_Vector(*oldStrikes ), *mergedStrikesVec );

	// we create an empty vol matrix with the correct dimension
	ARM_Matrix* mergedVol = new ARM_Matrix( mergedExpiriesVec->size(), mergedStrikesVec->size() );
	
	vector< pairIntInt >::const_iterator mergedExpiriesIter;
	vector< pairIntInt >::const_iterator mergedStrikesIter; 
	
	vector< pairIntInt >::const_iterator mergedExpiriesEnd= mergedExpiries.end();
	vector< pairIntInt >::const_iterator mergedStrikesEnd = mergedStrikes.end();

	int i1 = 0, j1 = 0;
	ARM_Matrix* newVol	   = newVolCrv->GetVolatilities();
	ARM_Matrix* oldVol	   = oldVolCrv->GetVolatilities();

	/// do a loop to get data from submatrices
	/// if missing data, throw an exception
	for( mergedExpiriesIter = mergedExpiries.begin();
	mergedExpiriesIter != mergedExpiriesEnd;
	++mergedExpiriesIter )
	{
		j1 = 0;
		for( mergedStrikesIter = mergedStrikes.begin();
		mergedStrikesIter!= mergedStrikesEnd; 
		++mergedStrikesIter )
		{
			// test that the data are missing
			bool missing =  ((*mergedExpiriesIter).second == 1 && (*mergedStrikesIter).second == 2)
			||  ((*mergedExpiriesIter).second == 2 && (*mergedStrikesIter).second == 1);
			

			///if missing data, throw an exception
			if( missing )
			{
				char msg[255];
				double expiry = (*mergedExpiriesIter).second == 1 ?
					newExpiries->Elt( (*mergedExpiriesIter).first )
					: oldExpiries->Elt( (*mergedExpiriesIter).first );

				double strike = (*mergedStrikesIter).second == 1 ?
					newStrikes->Elt( (*mergedStrikesIter).first )
					: oldStrikes->Elt( (*mergedStrikesIter).first );

				sprintf( msg, "data are missing in the new vol curve\n" );
				sprintf( msg, "misssing for tnore %f maturity %f and strike %f", tenor, expiry, strike );
					
				throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, msg );
			
			}

			/// in this case, it does exist in the submatrix and we need to find them
			else
			{
				/// is it in the first sub matrix or second?

				/// since the convention is 
				/// when there is a common element it is taken from the first vector
				/// we need to test when there was a second element at least somewhere

				/// if there is a second element
				if( ( (*mergedExpiriesIter).second == 2 )
					|| (  (*mergedStrikesIter).second == 2 ) )
				{
					/// in this case, whether the index pos given by the first int of the 
					/// pair is in V2 and we are fine
					/// otherwise we substract 3 to the second int
					/// since the convention is for a common number
					/// to put the V1 pos in the first int
					/// and V2 pos + 3 in the second int
					int pos1 = (*mergedExpiriesIter).second == 2 ? 
						(*mergedExpiriesIter).first
						: (*mergedExpiriesIter).second - 3;

					int pos2 = (*mergedStrikesIter).second == 2 ? 
						(*mergedStrikesIter).first
						: (*mergedStrikesIter).second - 3;

					/// just need to read now
					mergedVol->Elt(i1,j1) = oldVol->Elt( pos1, pos2 );
				}
				else
				{
					/// in this case, the data are in the first matrix
					mergedVol->Elt(i1,j1) = newVol->Elt( (*mergedExpiriesIter).first,
						(*mergedStrikesIter).first );
				}
			}
			++j1;
		}
		++i1;
	}

	/// once we are done
	/// we create the ARM_VolCurve
	ARM_Date asOf	= oldVolCrv->GetAsOfDate();
	int strikeType	= oldVolCrv->GetStrikeType();
	int volType		= oldVolCrv->GetOptionType();

	ARM_Vector* tmpMergedExpiriesVec = To_pARM_Vector(mergedExpiriesVec);
	delete mergedExpiriesVec;
	ARM_Vector* tmpMergedStrikesVec = To_pARM_Vector(mergedStrikesVec);
	delete mergedStrikesVec;
	ARM_VolCurve* result = new ARM_VolLInterpol( 
		asOf, 
		tmpMergedExpiriesVec,
		tmpMergedStrikesVec,
		mergedVol, 
		strikeType, 
		volType );

	return result;
}




////////////////////////////////////////////////////
///	Class  : ARM_SparseVolCube
///	Routine: AddOneVolCurve
///	Returns: void
///	Action : AddOneVolCurve
////////////////////////////////////////////////////

void ARM_SparseVolCube::AddOneVolCurve( ARM_VolCurve* volCurv, double underlying )
{
	/// first check whether there is already an underlying corresponding to it
	int index = lower_boundPosWithPrecision(*itsUnderlyings, underlying)-1;

	// equal
	if( index > -1 && itsUnderlyings->Elt(index) == underlying )
	{
		/// check that we are consistent
		ARM_VolCurve* oldVolCurv = (*itsVols)[index];
		
		/// the first check is to see that we have consitent data
		compare( volCurv, oldVolCurv, underlying );
		
		/// merge the volCurve
		ARM_VolCurve* mergedVolCurv = merge( volCurv, oldVolCurv, underlying );

		/// delete old curves
		delete oldVolCurv;
		delete volCurv;

		/// insert the new data
		(*itsVols)[index]=mergedVolCurv;
	}
	else
	{
		/// insert the element at the correct place
		++index;

		/// resize the vols and underlyings
		int size = itsVols->size() + 1;
		itsVols->resize( size );

		/// basically resize the ARM_GP_Vector by its previous size + 1
		/// putting a zero at the end ... for STL developpers,
		/// we cannot use Resize on ARM_GP_Vector because it will
		/// reinitialise everything... 
		ARM_GP_Vector* underlyingsTmp = new ARM_GP_Vector(itsUnderlyings->size()+1);
		int i;
		for(i=0;i<itsUnderlyings->size();++i)
			underlyingsTmp->Elt(i)=itsUnderlyings->Elt(i);
		delete itsUnderlyings;
		itsUnderlyings=underlyingsTmp;

		// and shift the data by one
		for(i=size-1; i>index; --i)
		{
			(*itsVols)[i] = (*itsVols)[i-1];
			itsUnderlyings->Elt(i) = itsUnderlyings->Elt(i-1);
		}
		
		/// insert the new data
		(*itsVols)[index]=volCurv;
		itsUnderlyings->Elt(index) = underlying;
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_SparseVolCube
///	Routine: Add Vol Curves
///	Returns: void
///	Action : Add Vol Curves
////////////////////////////////////////////////////

void ARM_SparseVolCube::AddVolCurves(
	int dim1Type,
	ARM_GP_Vector* dim1,
	int otherDimType,
	double otherDimValue,
	ARM_GP_Vector* strikes,
	ARM_Matrix* volatilities )
{
	Validation( dim1, strikes, dim1Type, otherDimType );

	ARM_Date asOf	= (*itsVols)[0]->GetAsOfDate();
	int strikeType	= (*itsVols)[0]->GetStrikeType();
	int volType		= (*itsVols)[0]->GetOptionType();

	switch( otherDimType )
	{
	// standard case
	/// easy case is that the sparse vol cube is correctly set up
	case K_TENOR:
		{
			int nbOfWrongTimeToStart = 0;
			
			/// if we specify the vol in term of the time to expiry
			/// we need to convert into time To Start
			/// as the vol cube dimension is time To Start!
			if( dim1Type == K_TIMETOEXPIRY )
			{
				int size = dim1->size();
				int colSize = volatilities->GetNumCols();
				
				for( int i=0; i<size; ++i )
				{
					if( dim1->Elt(i) - otherDimValue < ARM_SparseVolCube::spotTime )
					{
						/// expired vol are similar to the zero coupon vols
						/// in this case we take 1d maturities
						ARM_GP_Vector* maturities = new ARM_GP_Vector( 1, ARM_SparseVolCube::spotTime  );
						
						/// and a tenor corresponding to the time to expiry
						/// forces it to be more than one Day!
						double tenor = dim1->Elt(i);
						if( tenor < ARM_SparseVolCube::spotTime )
							tenor = ARM_SparseVolCube::spotTime;

						ARM_Vector* valuesVect = volatilities->GetRow(i);
						ARM_Matrix* data = new ARM_Matrix( 1.0, colSize, &(*valuesVect)[0]);
						delete valuesVect;

						ARM_VolCurve* curv = new ARM_VolLInterpol( const_cast<ARM_Date&>(asOf), 
							To_pARM_Vector(maturities),
							To_pARM_Vector(strikes),
							data, strikeType, volType );

						delete maturities;
						
						AddOneVolCurve( curv, tenor );
						nbOfWrongTimeToStart++;
					}
					else
						dim1->Elt(i) -= otherDimValue;
				}
			}
			
			ARM_Matrix* newVolatilities = NULL;
			if( nbOfWrongTimeToStart )
			{
				int dim1Size  = dim1->size();
				ARM_GP_Vector* newDim1 = new ARM_GP_Vector( dim1Size  - nbOfWrongTimeToStart );
				
				/// first get the correct dim1 number
				int i;
				for( i=nbOfWrongTimeToStart; i<dim1Size; ++i )
					newDim1->Elt( i- nbOfWrongTimeToStart ) = dim1->Elt(i);
				delete dim1;
				dim1 = newDim1;
				
				/// second gets correct volatilities
				int colSize = volatilities->GetNumCols();
				newVolatilities = new ARM_Matrix( dim1Size - nbOfWrongTimeToStart, colSize );
				
				int j;
				for( i=nbOfWrongTimeToStart; i<dim1Size; ++i )
					for( j=0; j<colSize; ++j )
						newVolatilities->Elt(i- nbOfWrongTimeToStart ,j) = volatilities->Elt(i,j);
			}
			else
			{
				newVolatilities = (ARM_Matrix*)volatilities->Clone();
			}
			
			ARM_Vector* tmpDim1 = To_pARM_Vector(dim1);
			delete dim1;

			/// UGLY but forced to const cast
			ARM_VolCurve* curv = new ARM_VolLInterpol( const_cast<ARM_Date&>(asOf), 
				tmpDim1, 
				To_pARM_Vector(strikes),
				newVolatilities, 
				strikeType, volType );
			AddOneVolCurve( curv, otherDimValue );
			break;
		}

	case K_TIMETOSTART:
		{
			int size	= dim1->size();
			int colSize = volatilities->GetNumCols();
			ARM_GP_Vector* maturities = new ARM_GP_Vector(1.0, otherDimValue );

			ARM_Vector* valuesVect = (const_cast<ARM_Matrix*>(volatilities))->GetRow(0);
			ARM_Matrix* data	= new ARM_Matrix( 1.0, colSize, &(*valuesVect)[0] );
			delete valuesVect;

			ARM_VolCurve* curv	= new ARM_VolLInterpol(
				asOf, 
				To_pARM_Vector(maturities),
				To_pARM_Vector(strikes),
				data,
				strikeType,
				volType );
			AddOneVolCurve( curv, dim1->Elt(0) );

			int i;
			/// start after one for the cloning of necessary data
			for( i=1; i<size; ++i )
			{
				/// the cloning is necessary to give
				/// to each object its own version of maturities and strikes
				/// so that they can delete objects properly
				ARM_Vector* valuesVect = (const_cast<ARM_Matrix*>(volatilities))->GetRow(i);
				ARM_Matrix* data	= new ARM_Matrix( 1.0, colSize, &(*valuesVect )[0]);
				delete valuesVect;

				ARM_VolCurve* curv	= new ARM_VolLInterpol( asOf, 
					To_pARM_Vector(maturities),
					To_pARM_Vector(strikes),
					data, 
					strikeType, 
					volType );
				AddOneVolCurve( curv, dim1->Elt(i) );
			}

			delete maturities;

			break;
		}
	/// time to expiry is not permitted for the other Dim type
	default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, 
			"Unknown other Dim type... be aware that time to expiry is not permitted for the other dim type");
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_SparseVolCube
///	Routine: Copy constructor
///	Returns: void
///	Action : Copy constructor
////////////////////////////////////////////////////

ARM_SparseVolCube::ARM_SparseVolCube( const ARM_SparseVolCube& rhs )
:	ARM_RootObject( rhs )
{
	CopyNoCleanUp( rhs );
}



////////////////////////////////////////////////////
///	Class  : ARM_SparseVolCube
///	Routine: CleanUp
///	Returns: void
///	Action : Copy constructor
///				CleanUp frees all the allocated memory
///				because itsVols is a vector of pointor to vol
///				curve, we have first to browse to delete vols
///				and then delete the vector
////////////////////////////////////////////////////

void ARM_SparseVolCube::CleanUp()
{
	/// first delete the vols
	int size = itsVols->size();
	int i;

	for( i=0; i<size; ++i )
		delete (*itsVols)[i];
	delete itsVols;

	/// delete the corresponding underlyings
	delete itsUnderlyings;
}


////////////////////////////////////////////////////
///	Class  : ARM_SparseVolCube
///	Routine: Destructor
///	Returns: void
///	Action : Destructor to free the vector of VolCurves
////////////////////////////////////////////////////

ARM_SparseVolCube::~ARM_SparseVolCube()
{
	CleanUp();
}



////////////////////////////////////////////////////
///	Class  : ARM_SparseVolCube
///	Routine: Assignment operator
///	Returns: ARM_SparseVolCube&
///	Action : Assignment operator
////////////////////////////////////////////////////

ARM_SparseVolCube& ARM_SparseVolCube::operator=( const ARM_SparseVolCube& rhs )
{
	if( this != &rhs )
	{
	    ARM_RootObject::operator = ( rhs );
		CleanUp();
		CopyNoCleanUp( rhs );
	}

	return *this;
}


/*!
 * Copy no clean up: deep copy of the extra layer
 */
void ARM_SparseVolCube::CopyNoCleanUp( const ARM_SparseVolCube& rhs )
{
	/// Copies all the volatilities by cloning them
	itsVols = new vector<ARM_VolCurve*>( rhs.itsVols->size() );
	vector< ARM_VolCurve* >::const_iterator itInput;
	vector< ARM_VolCurve* >::iterator itOutput;
	
	for( itInput = rhs.itsVols->begin(), itOutput = itsVols->begin();
			itInput < rhs.itsVols->end();
			++itInput, ++itOutput )
		*itOutput = (ARM_VolCurve*) (*itInput)->Clone();

	itsUnderlyings = (ARM_GP_Vector*) rhs.itsUnderlyings->Clone();

	/// copy the asOf and lastKnownDate
	itsLastKnownDate = rhs.itsLastKnownDate;
	itsAsOf = rhs.itsAsOf;
}


////////////////////////////////////////////////////
///	Class  : ARM_SparseVolCube
///	Routine: Clone
///	Returns: void
///	Action : Clones the object
////////////////////////////////////////////////////

ARM_Object* ARM_SparseVolCube::Clone() const
{
	return new ARM_SparseVolCube(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_SparseVolCube
///	Routine: GetVols
///	Returns: vector < ARM_VolCurve* >
///	Action : Get volatility accessor method
////////////////////////////////////////////////////
vector < ARM_VolCurve* > * ARM_SparseVolCube::GetVols() const
{
	return itsVols;
}


////////////////////////////////////////////////////
///	Class  : ARM_SparseVolCube
///	Routine: GetUnderlyings
///	Returns: ARM_GP_Vector*
///	Action : Get underlying accessor method
////////////////////////////////////////////////////
ARM_GP_Vector*	ARM_SparseVolCube::GetUnderlyings() const
{
	return itsUnderlyings;
}



////////////////////////////////////////////////////
///	Class  : ARM_SparseVolCube
///	Routine: ConvertToVolCube
///	Returns: ARM_VolCube*
///	Action : Conversion to a volcube
////////////////////////////////////////////////////

ARM_VolCube* ARM_SparseVolCube::ConvertToVolCube() const
{
	ARM_Vector* tmpUnderlyings = To_pARM_Vector(itsUnderlyings);
	ARM_VolCube* result = new ARM_VolCube( itsVols, tmpUnderlyings, itsLastKnownDate );
	delete tmpUnderlyings;
	return result;
}



////////////////////////////////////////////////////
///	Class  : ARM_SparseVolCube
///	Routine: View
///	Returns: void
///	Action : View the object
////////////////////////////////////////////////////

void ARM_SparseVolCube::View(char* id, FILE* ficOut) const
{
    FILE* fOut;
    char fOutName[200];
	
	/// first determine that the file is not already opened
    if ( NULL == ficOut )
    {
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w");
    }
    else
		fOut = ficOut;

    fprintf(fOut, "\n\n >>>>>>> Sparse volatility Cube <<<<<<<\n\n");

	/// AsOfDate
    char strDate[20];
	ARM_Date tempDate = itsAsOf;
    tempDate.JulianToStrDateDay(strDate);
	fprintf(fOut, "\t AsOfDate  : %s \n\n", strDate);

	/// last known date
	tempDate = itsLastKnownDate;
    tempDate.JulianToStrDateDay(strDate);
	fprintf(fOut, "\t LastKnownDate  : %s \n\n", strDate);

	/// loop over the vol curve and underlying
	if( itsVols )
	{
		size_t sz = itsUnderlyings->size();
		for (size_t i = 0; i < sz; i++)
		{
			fprintf(fOut, "\n\n\n------------------> %o) vol Curve\n", i+1 );
			fprintf(fOut, "\t Tenor : %lf\n", (*itsUnderlyings)[i]);
			(*itsVols)[i]->View(id, fOut);
		}
	}

    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}



CC_END_NAMESPACE()


/*---------------------------------------------------------------*/
/*---- End Of File ----*/

