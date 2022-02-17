/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pdenumericalschemes.h
 *	\author  A Schauly
 *	\version 1.0
 *	\date July 2005
 */

#ifndef _INGPINFRA_PDENUMSCHEME_H
#define _INGPINFRA_PDENUMSCHEME_H

#include "gpbase/port.h"
#include "gpbase/gpvector.h"
#include "gpbase/gpmatrix.h"
#include "gpinfra/typedef.h"
#include "gpbase/typedef.h"
#include "gpbase/mempoolmatrix.h"
#include "scheduler.h"//for the timesteps in the init of pde3Fnumericalscheme


CC_BEGIN_NAMESPACE( ARM )

class PDE_Truncator;

class ARM_PDENumericalScheme : public ARM_RootObject
{
public: 
	enum NumericalSchemeTypes
	{
		None = 0,
		Explicit1F,
		CN1F,
		Explicit2F,
		ADI2F,
		CS3F
	};

protected:
	ARM_MatrixVector itsGlobalVariances;
private: 
	/// Vector that changes toTimeIndex into matrixIndexes
	ARM_GP_VectorPtr itsTimeIndexesToMatrixIndexes;

	size_t itsSpaceDiscretizationPointsNb;  // Last Date SpaceSteps
	double itsSpaceDiscretizationStep;	  // How many StdDev LastEvent
	double itsStdDevNb;
	int LastTimeIdx;
	size_t itsPriceIndex;				/// Index of the 0 NumMethodSates

	PDE_Truncator * itsTruncator;

public:
	/// Functions Used for Pricing. To be reimplented for each PDENumerical Scheme
	virtual ARM_PricingStatesPtr Init( const ARM_PricingModel& model ) = 0;
	virtual void Induct( const ARM_PricingModel& model, ARM_PricingStatesPtr& states, int toTimeIdx) = 0;
	/// With truncation
	/// virtual void BuildTransitionMatrixes( const ARM_GP_VectorPtr& timeSteps,const ARM_GP_MatrixPtr& relativeDrifts,const ARM_GP_MatrixPtr& absoluteDrifts,const ARM_GP_MatrixPtr& vols,const ARM_GP_MatrixPtr& d1Vols,const ARM_GP_MatrixPtr& correls, const ARM_GP_T_Vector<size_t>& TruncationSizeVector, const ARM_PricingStatesPtr& states) {}
	/// Without truncation
	virtual void BuildTransitionMatrixes( const ARM_PricingModel& model,  const ARM_GP_VectorPtr& timeSteps, const ARM_PricingStatesPtr& states, double* timeStepStart = NULL, double* timeStepEnd = NULL) = 0;

	virtual bool IsOneDim() {return false;};
	/// OtherStuff
	virtual bool CheckCompatibilityWithModel( const ARM_PricingModel& model ) = 0;
	virtual NumericalSchemeTypes getNumericalSchemeType() const = 0;
	static ARM_PDENumericalScheme* getNumericalSchemeInstanceById( int NumSchemeId, size_t NX=-1, size_t NY = -1, size_t NZ = -1, double Theta1 = -1.0, double Theta2 = -1.0, double Theta3 = -1.0, int BoundConditionId = -1, double lambda = 0.0, int gridType = -1, ARM_GP_Matrix gridData=ARM_GP_Matrix(0), ARM_GP_Matrix schedulerData=ARM_GP_Matrix(0) );


	/// Accessors
	inline size_t getSpaceDiscretizationPointsNb() const { return itsSpaceDiscretizationPointsNb; }
	inline void setSpaceDiscretizationPointsNb( size_t discretizationPointsNb ) { itsSpaceDiscretizationPointsNb=discretizationPointsNb; }
	inline double getSpaceDiscretizationStep() const { return itsSpaceDiscretizationStep; }
	inline void setSpaceDiscretizationStep( double SpaceDiscretizationStep ) { itsSpaceDiscretizationStep = SpaceDiscretizationStep; }
	inline double getStdDevNb() const { return itsStdDevNb; }
	inline void setStdDevNb( double StdDevNb ) { itsStdDevNb = StdDevNb; }
	inline int getLastTimeIdx() const { return LastTimeIdx; }
	inline void setLastTimeIdx( int LTI ) { LastTimeIdx = LTI; }
	inline size_t getPriceIndex() const { return itsPriceIndex; }
	inline void setPriceIndex( size_t PriceIndex) { itsPriceIndex = PriceIndex; }
	inline void setTimeToMatrixIndex( size_t TimeIdx, size_t MatrixIdx ) { itsTimeIndexesToMatrixIndexes->Elt(TimeIdx)=MatrixIdx;}
	inline size_t getMatrixIdxFromTimeIdx( size_t TimeIdx ) { return itsTimeIndexesToMatrixIndexes->Elt(TimeIdx); }
	inline void setTimeIndexesToMatrixIndexes( const ARM_GP_VectorPtr& TimeIndexesToMatrixIndexes ) { itsTimeIndexesToMatrixIndexes = TimeIndexesToMatrixIndexes;}
	inline ARM_GP_VectorPtr getTimeIndexesToMatrixIndexes() const { return itsTimeIndexesToMatrixIndexes; }
	inline PDE_Truncator * getTruncator() const { return itsTruncator; }
	inline void setTruncator( PDE_Truncator * Truncator ) {  itsTruncator=Truncator; }

	const ARM_MatrixVector& GetGlobalVars() const { return itsGlobalVariances; };

	/// Standard ARM Support
	virtual string toString(const string& indent="", const string& nextIndent="") const = 0;
	virtual ARM_Object* Clone() const = 0;

	/// Constructor/Destructor
	ARM_PDENumericalScheme( size_t DiscretizationPointsNb = 50, double StdDevNb = 6., PDE_Truncator* Truncator =NULL);
	ARM_PDENumericalScheme( const ARM_PDENumericalScheme& rhs );
	virtual ~ARM_PDENumericalScheme();
};


class ARM_PDE1FNumericalScheme : public ARM_PDENumericalScheme
{
protected: 
	ARM_GP_VectorPtr itsVolChangeTimeSteps;
	ARM_GP_VectorPtr itsVolChangeIndexes;

#if defined(__GP_STRICT_VALIDATION)
	ARM_GP_MatrixPtr StoredVols;
#endif

public:
	virtual bool IsOneDim() {return true;};
	/// Init is the Same for Every 1F PDEs
	virtual ARM_PricingStatesPtr Init( const ARM_PricingModel& model );

	/// Basically checks that model is 1F
	virtual bool CheckCompatibilityWithModel( const ARM_PricingModel& model );

	/// Updates a vector of payoffs with Tridiagonal Matrixes (used both by Explicit and CN)
	void UpdateVectorWithMatrixByProduct( const ARM_GP_VectorPtr& DiagTerm, const ARM_GP_VectorPtr& UpperTerm, const ARM_GP_VectorPtr& LowerTerm, ARM_GP_Matrix& TransposedVector );
	void UpdateVectorWithMatrixByProduct( const ARM_GP_VectorPtr& DiagTerm, const ARM_GP_VectorPtr& UpperTerm, const ARM_GP_VectorPtr& LowerTerm, std::vector<double>& Vector );
	void UpdateVectorWithMatrixByProduct( const ARM_GP_VectorPtr& DiagTerm, const ARM_GP_VectorPtr& UpperTerm, const ARM_GP_VectorPtr& LowerTerm, ARM_MemPool_Matrix::iterator vecBegin );

	/// Standard ARM Support
	virtual string toString(const string& indent="", const string& nextIndent="") const = 0;
	virtual ARM_Object* Clone() const = 0;

	/// Constructors
#if defined(__GP_STRICT_VALIDATION)
	ARM_PDE1FNumericalScheme( size_t DiscretizationPointsNb, double StdDevNb, PDE_Truncator* Truncator ) : ARM_PDENumericalScheme( DiscretizationPointsNb, StdDevNb, Truncator ), StoredVols(NULL) {}
#else
	ARM_PDE1FNumericalScheme( size_t DiscretizationPointsNb, double StdDevNb, PDE_Truncator* Truncator ) : ARM_PDENumericalScheme( DiscretizationPointsNb, StdDevNb, Truncator ) {}
#endif

	ARM_PDE1FNumericalScheme( ) : ARM_PDENumericalScheme() {}
	ARM_PDE1FNumericalScheme( const ARM_PDE1FNumericalScheme& rhs ) : ARM_PDENumericalScheme(rhs) {}
	virtual ~ARM_PDE1FNumericalScheme(){};
};

class ARM_PDE2FNumericalScheme : public ARM_PDENumericalScheme
{
protected: 
	ARM_GP_VectorPtr itsVolChangeTimeSteps;
	ARM_GP_VectorPtr itsVolChangeIndexes;

	double itsSpaceDiscretizationStepX;
	double itsSpaceDiscretizationStepY;

public:
	/// Init is the Same for Every 2F PDEs
	virtual ARM_PricingStatesPtr Init( const ARM_PricingModel& model );

	/// Basically checks that model is 2F
	virtual bool CheckCompatibilityWithModel( const ARM_PricingModel& model );

	/// Accessors for 2F Models
	inline double getSpaceDiscretizationStepX() const { return itsSpaceDiscretizationStepX; }
	inline void setSpaceDiscretizationStepX( double SpaceDiscretizationStepX ) { itsSpaceDiscretizationStepX = SpaceDiscretizationStepX; }
	inline double getSpaceDiscretizationStepY() const { return itsSpaceDiscretizationStepY; }
	inline void setSpaceDiscretizationStepY( double SpaceDiscretizationStepY ) { itsSpaceDiscretizationStepY = SpaceDiscretizationStepY; }

	/// Standard ARM Support
	virtual string toString(const string& indent="", const string& nextIndent="") const = 0;
	virtual ARM_Object* Clone() const = 0;

	/// Constructors
	ARM_PDE2FNumericalScheme( size_t DiscretizationPointsNb, double StdDevNb, PDE_Truncator* Truncator  ) : ARM_PDENumericalScheme( DiscretizationPointsNb, StdDevNb, Truncator  ) {}
	ARM_PDE2FNumericalScheme( ) : ARM_PDENumericalScheme() {}
	ARM_PDE2FNumericalScheme( const ARM_PDE2FNumericalScheme& rhs ) : ARM_PDENumericalScheme(rhs) {}
	virtual ~ARM_PDE2FNumericalScheme(){};
};

class ARM_PDE3FNumericalScheme : public ARM_PDENumericalScheme
{
public:
	enum GridType
	{
		StdDev,
		Fixed,
		BiReg
	};

	/// Init is the Same for Every 3F PDEs
	virtual ARM_PricingStatesPtr Init( const ARM_PricingModel& model, double lambda= 0);

	/// Basically checks that model is 3F
	virtual bool CheckCompatibilityWithModel( const ARM_PricingModel& model );

	/// Standard ARM Support
	virtual string toString(const string& indent="", const string& nextIndent="") const = 0;
	virtual ARM_Object* Clone() const = 0;

	/// Constructors
	ARM_PDE3FNumericalScheme(size_t NX, size_t NY, size_t NZ, GridType gridType, const ARM_GP_Matrix& gridData, const ARM_GP_Matrix& schedulerData);
	ARM_PDE3FNumericalScheme( const ARM_PDE3FNumericalScheme& rhs );
	virtual ~ARM_PDE3FNumericalScheme(){};

	// accessors
	ARM_GP_VectorPtr GetTimeSteps() const { return itsTimeSteps; }
	const ARM_PricingModel* GetModel() const { return itsModel; }

	const size_t GetNX() const { return itsNX; }
	const size_t GetNY() const { return itsNY; }
	const size_t GetNZ() const { return itsNZ; }

	inline void SetNX( size_t NX ) { itsNX=NX; }
	inline void SetNY( size_t NY ) { itsNY=NY; }
	inline void SetNZ( size_t NZ ) { itsNZ=NZ; }

	const ARM_GP_VectorPtr& GetDeltaX() const { return itsDeltaX; }
	const ARM_GP_VectorPtr& GetDeltaY() const { return itsDeltaY; }
	const ARM_GP_VectorPtr& GetDeltaZ() const { return itsDeltaZ; }

	inline void SetDeltaX( const ARM_GP_VectorPtr& deltaX) { itsDeltaX = deltaX; }
	inline void SetDeltaY( const ARM_GP_VectorPtr& deltaY) { itsDeltaY = deltaY; }
	inline void SetDeltaZ( const ARM_GP_VectorPtr& deltaZ) { itsDeltaZ = deltaZ; }

	const ARM_GP_VectorPtr& GetXGrid() const { return itsXGrid; }
	const ARM_GP_VectorPtr& GetYGrid() const { return itsYGrid; }
	const ARM_GP_VectorPtr& GetZGrid() const { return itsZGrid; }

	inline void SetXGrid( const ARM_GP_VectorPtr& XGrid) { itsXGrid = XGrid; }
	inline void SetYGrid( const ARM_GP_VectorPtr& YGrid) { itsYGrid = YGrid; }
	inline void SetZGrid( const ARM_GP_VectorPtr& ZGrid) { itsZGrid = ZGrid; }


private:

	GridType itsGridType;
	ARM_GP_Matrix itsGridData;
	ARM_GP_Matrix itsSchedulerData;

	size_t itsNX;
	size_t itsNY;
	size_t itsNZ;

	ARM_GP_VectorPtr itsDeltaX;
	ARM_GP_VectorPtr itsDeltaY;
	ARM_GP_VectorPtr itsDeltaZ;

	ARM_GP_VectorPtr itsXGrid;
	ARM_GP_VectorPtr itsYGrid;
	ARM_GP_VectorPtr itsZGrid;

	ARM_GP_VectorPtr itsTimeSteps;

	// Bidouille
	const ARM_PricingModel* itsModel;
};

CC_END_NAMESPACE()

#endif