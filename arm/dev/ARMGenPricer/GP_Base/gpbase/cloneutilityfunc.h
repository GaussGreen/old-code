/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *
 *	\file cloneutilityfunc.h
 *
 *  \brief files to holds any function that is very general
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date December 2003
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPBASE_CLONEUTILITYFUNC_H
#define _INGPBASE_CLONEUTILITYFUNC_H

/// use our macro for namespace
#include "port.h"
#include "env.h"
#include "ostringstream.h"
#include "countedptr.h"
#include "typedef.h"
#include <glob/expt.h>

#include <vector>
CC_USING_NS( std, vector )

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Routine: DuplicateCloneablePointorVectorInPlace
///	Returns: void
///	Action : for a vector of cloneable objects (basically can call
///             Clone on the pointor
///         creates an output element that is similar 
///         to the input ....
///
///         basically to use for any vector< ARM_ObjectLike* > object
////////////////////////////////////////////////////
template <typename T>
     void DuplicateCloneablePointorVectorInPlace( const vector< T* >& input, vector< T *>& output )
{
    size_t inputSize = input.size();
    output.resize( inputSize );
    size_t i;
    for( i=0; i<inputSize; ++i )
    {
        output[i] = input[i] ? (T*) input[i]->Clone(): NULL;
    }
}




////////////////////////////////////////////////////
///	Routine: DuplicateCloneablePointorVectorVectorInPlace
///	Returns: void
///	Action : for a vector of vector of cloneable objects (basically can call
///             Clone on the pointor
///         creates an output element that is similar 
///         to the input ....
///
///         basically to use for any vector< ARM_ObjectLike* > object
////////////////////////////////////////////////////
template <typename T>
     void DuplicateCloneablePointorVectorVectorInPlace( const vector< vector< T* > >& input, vector< vector< T* > >& output )
{
    size_t inputSize = input.size();
    output.resize( inputSize );
    size_t i,j;
    for( i=0; i<inputSize; ++i )
    {
		size_t inputSize2 = input[i].size();
		output[i].resize(inputSize2);

		for( j=0; j<inputSize2; ++j )
		{
	#if defined(__GP_STRICT_VALIDATION)
			/// for efficiency reason assumed that the input elements are not NULL
			if( !input[i][j] )
			{
				CC_Ostringstream os;
				os << "argument " << i << << "," << j << " is NULL! cannot be cloned!";
				ARM_THROW( ERR_INVALID_ARGUMENT, os.str() );
			}
	#endif
			output[i][j] = (T*) input[i][j]->Clone();
		}
    }
}


////////////////////////////////////////////////////
///	Routine: DuplicateCloneablePointorVector
///	Returns: void
///	Action : same as DuplicateCloneablePointorVectorInPlace
///				except that it returns the result!
////////////////////////////////////////////////////

template <typename T>
     vector<T*> DuplicateCloneablePointorVector( const vector< T* >& input )
{
    vector<T*> output;
    DuplicateCloneablePointorVectorInPlace<T>(input, output );
    return output;
}


////////////////////////////////////////////////////
///	Routine: DuplicateCloneablePointorAndNullVectorInPlace
///	Returns: void
///	Action : same as DuplicateCloneablePointorVectorInPlace
///				except that it accepts null pointor
////////////////////////////////////////////////////
template <typename T>
     void DuplicateCloneablePointorAndNullVectorInPlace( const vector< T* >& input, vector< T *>& output )
{
    size_t inputSize = input.size();
    output.resize( inputSize );
    size_t i;
    for( i=0; i<inputSize; ++i )
        output[i] = input[i] ? (T*) input[i]->Clone() : NULL;
}


////////////////////////////////////////////////////
///	Routine: DuplicateCloneablePointorAndNullVector
///	Returns: void
///	Action : same as DuplicateCloneablePointorAndNullVectorInPlace
///				except that it accepts null pointor
////////////////////////////////////////////////////

template <typename T>
     vector<T*> DuplicateCloneablePointorAndNullVector( const vector< T* >& input )
{
    vector<T*> output;
    DuplicateCloneablePointorAndNullVectorInPlace<T>(input, output );
    return output;
}


////////////////////////////////////////////////////
///	Routine: DuplicateCloneablePtrVectorInPlace
///	Returns: void
///	Action : for a vector of cloneable pointor objects 
///         creates an output element that is similar 
///         to the input ....
///
///         basically to use for any vector< ARM_ObjectPtrLike > object
////////////////////////////////////////////////////
template <typename T>
     void DuplicateCloneablePtrVectorInPlace( const vector< ARM_CountedPtr< T > >& input, vector<  ARM_CountedPtr< T >  >& output )
{
    size_t inputSize = input.size();
    output.resize( inputSize );
    size_t i;
    for( i=0; i<inputSize; ++i )
    {
#if defined(__GP_STRICT_VALIDATION)
        /// for efficiency reason assumed that the input elements are not NULL
        if( input[i] == ARM_CountedPtr<T>(NULL) )
        {
            CC_Ostringstream os;
            os << "argument " << i << " is NULL! cannot be cloned!";
            ARM_THROW( ERR_INVALID_ARGUMENT, os.str() );
        }
#endif
        output[i] = ARM_CountedPtr<T>( (T*) input[i]->Clone() );
    }
}



////////////////////////////////////////////////////
///	Routine: DuplicateCloneablePtrVector
///	Returns: void
///	Action : same as DuplicateCloneablePtrVectorInPlace
///				except that it accepts null pointor
////////////////////////////////////////////////////

template <typename T>
     vector<T*> DuplicateCloneablePtrVector(  const vector< ARM_CountedPtr< T > >& input )
{
    vector< ARM_CountedPtr< T > > output;
    DuplicateCloneablePointorAndNullVectorInPlace<T>(input, output );
    return output;
}



////////////////////////////////////////////////////
///	Routine: DuplicateCloneablePtrAndNullVectorInPlace
///	Returns: void
///	Action : similar to DuplicateCloneablePtrVectorInPlace
///				except that it handles ARM_CountedPtr<T>(NULL) object! 
////////////////////////////////////////////////////
template <typename T>
     void DuplicateCloneablePtrAndNullVectorInPlace( const vector< ARM_CountedPtr< T > >& input, vector<  ARM_CountedPtr< T >  >& output )
{
    size_t inputSize = input.size();
    output.resize( inputSize );
    size_t i;
    for( i=0; i<inputSize; ++i )
        output[i] = input[i] != ARM_CountedPtr<T>(NULL) ? ARM_CountedPtr<T>( (T*) input[i]->Clone() ): ARM_CountedPtr<T>(NULL) ;
}



////////////////////////////////////////////////////
///	Routine: DuplicateCloneablePtrAndNullVector
///	Returns: void
///	Action : same as DuplicateCloneablePtrAndNullVectorInPlace
///				except that it accepts null pointor
////////////////////////////////////////////////////

template <typename T>
     vector<T*> DuplicateCloneablePtrAndNullVector(  const vector< ARM_CountedPtr< T > >& input )
{
    vector< ARM_CountedPtr< T > > output;
    DuplicateCloneablePtrAndNullVectorInPlace<T>(input, output );
    return output;
}




////////////////////////////////////////////////////
///	Routine: DeletePointorVector
///	Returns: void
///	Action : for an input that is a vector< ARM_ObjectLike* >
///             deletes the corresponding pointor and sets
///             them to NULL!
////////////////////////////////////////////////////

template <typename T>
     void DeletePointorVector( vector< T* >& vec )
{
    for( size_t i=0; i<vec.size(); ++i )
    {
        delete vec[i];
#if defined( __ARM_VECTOR_NO_RANGE_CHECK )
        vec[i] = NULL;
#endif
    }
}


////////////////////////////////////////////////////
///	Routine: DeletePointorVectorVector
///	Returns: void
///	Action : for an input that is a vector< ARM_ObjectLike* >
///             deletes the corresponding pointor and sets
///             them to NULL!
////////////////////////////////////////////////////

template <typename T>
     void DeletePointorVectorVector( vector< vector< T* > > & vec )
{

    for( size_t i=0; i<vec.size(); ++i )
	{
		for( size_t j=0; j<vec[i].size(); ++j )
		{
			delete vec[i][j];
#if defined( __ARM_VECTOR_NO_RANGE_CHECK )
			vec[i][j] = NULL;
#endif	
		}
    }
}


template <typename T>
    ARM_CountedPtr<T> CreateClonedPtr( T* p)
{
	return ARM_CountedPtr<T>( p? (T*) p->Clone() : NULL );
}

/// create clone and handle NULL case
template <typename T>
	T* CreateClone( T* p)
{
	return p ? (T*) p->Clone() : NULL ;
}





CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/


