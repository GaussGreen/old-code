/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: gramfunctorargcheck.cpp,v $
 * Revision 1.1  2003/10/20 19:43:31  ebenhamou
 * Initial revision
 *
 */


/*! \file gramfunctorargcheck.cpp
 *
 *  \brief checker class that are only used at debug time
 *		solves to nothing in release mode
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#include "gpinfra/gramfunctorargcheck.h"
#include "gpbase/ostringstream.h"
#include "gpbase/env.h"
#include "gpinfra/gramfunctor.h"
#include "gpinfra/gramnode.h"


/// for easy debugging of shared nodes!
#if defined(__GP_SHOW_SHARED_NODE_COORDINATES)
	#include "gpbase/pair.h"
#endif

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : GPAF_Pb
///	Routine: UnexpectecType
///	Returns: void
///	Action : throw an exception if type is not correct!
////////////////////////////////////////////////////

void GPAF_Pb::UnexpectedType( ARM_GramFctorArg::Type FoundType, 
	ARM_GramFctorArg::Type RequiredType, const string& funcName )
{
	CC_Ostringstream os;
	os << "In " << funcName << " function: unexpected type! expected " 
		<< ARM_GramFctorArg::TypeToString(RequiredType) << " found "
		<< ARM_GramFctorArg::TypeToString(FoundType)<< " " << ARM_USERNAME
		<< " please advise";
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
}


#ifdef __GP_STRICT_VALIDATION
	/// Checker class for Generic Pricer Argument Functor 
	/// in DEBUG mode, check the size of the vector of arguments
	/// is equal to nb
	GPAF_CheckArgSize::GPAF_CheckArgSize( const ARM_GramFctorArgVector& arg, int nb, const string& funcName )
	{
		/// checking of the size
		if( arg.size() != nb )
		{
			CC_Ostringstream os;
			os << "In " << funcName << " function: expected " << nb << " arguments  " 
				<< " found " << arg.size() << " " << ARM_USERNAME 
				<< " please advise";

			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}
	}
	/// Checker class for Generic Pricer Argument Functor 
	/// Check the arguments type
	GPAF_CheckArgType::GPAF_CheckArgType( const ARM_GramFctorArg& arg, ARM_GramFctorArg::Type t, const string& funcName )
	{
		switch(t)
		{
		case GFAT_DOUBLE_TYPE:
		case GFAT_STRING_TYPE:
		case GFAT_DATE_TYPE:
		case GFAT_CURVE_TYPE:
			{
				if( t!= arg.GetType() )
					GPAF_Pb::UnexpectedType( arg.GetType(), t, funcName );
				break;
			}
		case GFAT_MULTITOKENSTRING_TYPE:
		case GFAT_MATURITY_TYPE:
			{
				if( GFAT_STRING_TYPE != arg.GetType() )
					GPAF_Pb::UnexpectedType( arg.GetType(), t, funcName );
				break;
			}
		case GFAT_VECTOR_TYPE:
			{
				/// a little special as vector can be both double
				/// or double!
				if(		GFAT_DOUBLE_TYPE != arg.GetType() 
					&&	GFAT_VECTOR_TYPE != arg.GetType() )
					GPAF_Pb::UnexpectedType( arg.GetType(), t, funcName );
				break;
			}
		case GFAT_DATEORMATU_TYPE:
			{
				if(		GFAT_STRING_TYPE != arg.GetType() 
					&&	GFAT_DATE_TYPE   != arg.GetType() )
					GPAF_Pb::UnexpectedType( arg.GetType(), t, funcName );
				break;
			}
		case GFAT_DATE_VECTOR_OR_DOUBLE_TYPE: case GFAT_DATE_OR_VECTOR_TYPE:
        case GFAT_DATE_OR_DOUBLE_TYPE :
			{
				/// a little special as vector can be both double
				/// or double!
				if(		GFAT_DATE_TYPE   != arg.GetType() 
					&&	GFAT_VECTOR_TYPE != arg.GetType() 
					&&  GFAT_DOUBLE_TYPE != arg.GetType() )
					GPAF_Pb::UnexpectedType( arg.GetType(), t, funcName );
				break;
			}
		case GFAT_VECTOR_OR_CURVE_TYPE :
			{
				/// a little special as vector can be both double
				/// or double!
				if(		GFAT_DOUBLE_TYPE   != arg.GetType() 
					&&	GFAT_VECTOR_TYPE != arg.GetType() 
					&&  GFAT_CURVE_TYPE != arg.GetType() )
					GPAF_Pb::UnexpectedType( arg.GetType(), t, funcName );
				break;
			}
		case GFAT_STRING_OR_CURVE_TYPE :
			{
				if(		GFAT_STRING_TYPE   != arg.GetType() 
					&&	GFAT_CURVE_TYPE != arg.GetType() )
					GPAF_Pb::UnexpectedType( arg.GetType(), t, funcName );
				break;
			}
		case GFAT_DATESTRIP_TYPE:
			{
                /// Double case to handle arbitrary double default value
				if( GFAT_DOUBLE_TYPE != arg.GetType() && GFAT_DATESTRIP_TYPE != arg.GetType())
					GPAF_Pb::UnexpectedType( arg.GetType(), t, funcName );
				break;
			}

		default:
			{
				CC_Ostringstream os;
				os << "In " << funcName << " unknown type :"
					<< ARM_GramFctorArg::TypeToString(t)
					<< " " << ARM_USERNAME << " please advise";
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
			}
		}
	}

	/// Checker class for Generic Pricer Argument Functor 
	/// Check that the various arguments are vector of same size
	GPAF_CheckArgVecSameSize::GPAF_CheckArgVecSameSize( const ARM_GramFctorArgVector& arg, const string& funcName )
	{
		int i = 0, sizeMax = arg.size();
		
		/// found first vector
		while( i<sizeMax && GFAT_VECTOR_TYPE!= arg[i].GetType() )
			++i;

		/// not the end
		if( i != sizeMax )
		{
			size_t sizeStdt = arg[i].GetVector()->size();
			int iStdt = i;

			while(i<sizeMax)
			{
				if(		GFAT_VECTOR_TYPE == arg[i].GetType()
					&&	arg[i].GetVector()->size() != sizeStdt
					)
				{
					CC_Ostringstream os;
					os << "In " << funcName << " function: expected same size's vectors! Found vector "
						<< i << " of size " << arg[i].GetVector()->size()
						<< " and vector " << iStdt << "th of size " << sizeStdt << " "
						<< ARM_USERNAME << " please advise";
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
				}
				++i;
			}
		}
	}


	GPAF_CheckArgReturnType::GPAF_CheckArgReturnType( const ARM_GramFctorArg& arg, ARM_GramFctorArg::Type t, const ARM_ExpNodePtr& node, const string& funcName )
	{
		// check first special type like GFAT_FUTUREREF_TYPE
		if(GFAT_FUTUREREF_TYPE == t)
		{
			/// node is a smartPointor to ARM_ExpNode not ARM_ExpNodeRef
			/// hence the casting!
			if( ! ((ARM_ExpNodeRef&) *node).IsPreComputed() )
			{
				CC_Ostringstream os;
				os << "In " << funcName << "function: expected precomputed node "
					<< node->toString() << " " << ARM_USERNAME
					<< " please advise";
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
			}
		}
        else if( GFAT_VECTOR_OR_FUTUREREF_TYPE == t)
        {
            const ARM_ExpNodeRef* nodeRef = dynamic_cast<ARM_ExpNodeRef*>(&*node);
            if( nodeRef && !nodeRef->IsPreComputed() )
			    GPAF_CheckArgType argCheck( arg, ARM_GramFctorArg::GetResultType(t), funcName );
        }
		else
			/// other cases!
			GPAF_CheckArgType argCheck( arg, ARM_GramFctorArg::GetResultType(t), funcName );
	}
#else
	/// All check functions resolve to nothing

	GPAF_CheckArgSize::GPAF_CheckArgSize( const ARM_GramFctorArgVector& arg, int nb, const string& funcName ){}

	GPAF_CheckArgVecSameSize::GPAF_CheckArgVecSameSize( const ARM_GramFctorArgVector& arg, const string& funcName ){}

	GPAF_CheckArgType::GPAF_CheckArgType( const ARM_GramFctorArg& arg, ARM_GramFctorArg::Type, const string& funcName ){}

	GPAF_CheckArgReturnType::GPAF_CheckArgReturnType( const ARM_GramFctorArg& arg, ARM_GramFctorArg::Type t, const ARM_ExpNodePtr& node, const string& funcName ){}
#endif




CC_END_NAMESPACE()



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

