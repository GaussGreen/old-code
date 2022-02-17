/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *
 *	\file modelnamemap.h
 *
 *  \brief model map is a dictionary containing all the models
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPINFRA_MODELNAMEMAP_H
#define _INGPINFRA_MODELNAMEMAP_H

/// this headers has to come first
/// as it contains pragma for stupid warning
#include "gpbase/removeidentifiedwarning.h"

#include "gpbase/port.h"		/// to make code platform independent
#include "gpbase/rootobject.h"

#include "typedef.h"

#include "gpbase/gpvector.h"

#include <string>
CC_USING_NS( std, string )

/// use of a namespace
/// macro for cross platforms code
CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_PricingModel;
class ARM_PricingModelIR;



/// \class   ARM_ModelNameMap 
///
///	\brief  a class to map string to model. It is more than a simple map
/// between models and their corresponding names as some models
/// may need to consult one or more other models in order to diffuse
/// themselves properly. For example, the inflation model needs
/// to consult another model to get its drift. In order to make 
/// the difference between real model and model used only for other model
/// there is a boolean called itsUsedAsOtherModel
/// if a model needs another or some model(s), the corresponding number(s)
/// as stored in itsOtherModelRefNb

/// the design is also to use a nested class to do all the map work
/// and have a simple vector of this nested class
/// iterator const and non const version provides easy access to the
/// model name map!
///
///	\author  Eric Benhamou
///	\version 1.0
///


class ARM_ModelNameMap : public ARM_RootObject
{
/// forward declaration
/// the conventional ARM way of writting is to specify 
/// private variable before public
/// hence the need to do this forward declaration
public:
	/// \class DataElem
	/// nested class to do all the work
	/// this class is the equivalent of a map
	/// but we do not really need the STL here
	/// so not used!
    struct DataElem 
	{ 
		/// default constructor
		DataElem() : itsModelName(), itsModel( NULL ), itsOtherModelRefNb() {}

		/// real constructor
		DataElem( const string& Name, const ARM_PricingModelPtr& Model, 
			const  ARM_IntVector &OtherElems, bool UsedInPricing = true, bool UsedAsOtherModel = false ) 
		: 	itsModelName( Name ), 
			itsModel( Model ), 
			itsOtherModelRefNb( OtherElems ),
			itsUsedInPricing( UsedInPricing ),
			itsUsedAsOtherModel( UsedAsOtherModel )
		{}

		/// assignment operator
		DataElem &operator=( const DataElem &rhs )
		{
			if( this != &rhs )
			{
				itsModelName		= rhs.itsModelName;
				itsModel			= rhs.itsModel;  
				itsOtherModelRefNb	= rhs.itsOtherModelRefNb;
				itsUsedInPricing	= rhs.itsUsedInPricing;
				itsUsedAsOtherModel	= rhs.itsUsedAsOtherModel;
			}
			return *this;
		}

		bool operator< ( const DataElem &rhs );
		
		/// accessors
		const string&	ModelName() const { return itsModelName;  }
		const ARM_PricingModelPtr& Model() const { return itsModel; }
		ARM_PricingModelPtr& Model() { return itsModel; }
		const ARM_IntVector& OtherModelRefNb() const { return itsOtherModelRefNb; }
		bool UsedInPricing() const { return itsUsedInPricing; }
		bool UsedAsOtherModel()	const { return itsUsedAsOtherModel; }
		void ResetOtherModels() {itsOtherModelRefNb.clear();}

		/// and Set methods
		void SetUsedInPricing( bool UsedInPricing ) { itsUsedInPricing = UsedInPricing; }
		void SetUsedAsOtherModel( bool UsedAsOtherModel ) { itsUsedAsOtherModel = UsedAsOtherModel; }

	private:
		/// the model name
		string itsModelName; 
	
		/// the pointor to the real model (shared model then use ARM_CountedPtr)
		ARM_PricingModelPtr itsModel;

		/// if a model does need some other models and their corresponding number!
		ARM_IntVector itsOtherModelRefNb;

		/// it used in the pricing?
		bool itsUsedInPricing;

		/// is it use only for another model?
		bool itsUsedAsOtherModel;
	};

	/// for the cpp file
	/// for easyness of the code
	typedef CC_STL_VECTOR( DataElem ) DataT;	

private:
	/// the vector of elements
	/// we use a vector as this should not 
	/// grow to much.. hence even a linear
	/// search will do the job properly!
	CC_STL_VECTOR( DataElem ) itsData;

	/// vector of elements sorted by order
	/// of priority of the model
	CC_STL_VECTOR( DataElem ) itsSortedData;

public:
	
	/// constructor
	ARM_ModelNameMap( const ARM_StringVector& Names		= ARM_StringVector(),
		const vector< ARM_PricingModelPtr >& Models		= vector< ARM_PricingModelPtr >(), 
		const ARM_StringVectorVector& OtherModelNames	= ARM_StringVectorVector() );

	/// copy constructor
	ARM_ModelNameMap( const ARM_ModelNameMap& ModelNameMap );

	/// assignment operator
	ARM_ModelNameMap& operator=( const ARM_ModelNameMap& rhs );

	/// destructor
	virtual ~ARM_ModelNameMap();
	
	/// iterators are the ones of the STL
	/// this class wraps up the STL iterator corresponding to the
	/// vector

	/// Iterators:
	typedef CC_STL_VECTOR( DataElem )::const_iterator const_iterator;
	typedef CC_STL_VECTOR( DataElem )::iterator iterator;

	/// iterator support
	/// const version
	inline const_iterator begin() const { return itsData.begin(); }
	inline const_iterator end() const { return itsData.end(); }

	/// non const version
	inline iterator begin() { return itsData.begin(); }
	inline iterator end() { return itsData.end(); }

	/// sorted iterator support
	/// const version
	inline const_iterator sortedBegin() const { return itsSortedData.begin(); }
	inline const_iterator sortedEnd() const { return itsSortedData.end(); }

	/// non const version
	inline iterator sortedBegin() { return itsSortedData.begin(); }
	inline iterator sortedEnd() { return itsSortedData.end(); }

    /// accessor operator with operator []
	iterator operator[] ( const string& Name );
	const_iterator operator[] ( const string& Name ) const;
	inline iterator operator[] ( size_t Index ) { return itsData.begin()+ Index; }
	inline const_iterator operator[] ( size_t Index ) const { return itsData.begin()+ Index; }

	/// Test if the model exist
	bool TestIfModelExisting ( const string& Name ) const; 

	const_iterator& getNextUsedIter( const_iterator& iter ) const;
	iterator& getNextUsedIter( iterator& iter );

	const_iterator& getSortedNextUsedIter( const_iterator& iter ) const;
	iterator& getSortedNextUsedIter( iterator& iter );

	/// to access the size of the data
	size_t size() const { return itsData.size(); }
	size_t UsedModelsSize() const;

	/// ARM_Object stuff
	virtual ARM_Object* Clone() const { return new ARM_ModelNameMap(*this); }
	/// for easy debugging in main
	virtual string toString(const string& indent="", const string& nextIndent="") const;
};

CC_END_NAMESPACE()

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/


