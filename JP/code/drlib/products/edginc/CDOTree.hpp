//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : CDOTree.hpp
//
//   Description : SpreadLossTree product for GeneralisedCDO
//
//   Author      : Matthias Arnsdorf
//
//   Date        : August, 2006
//
//----------------------------------------------------------------------------
#ifndef QR_CDO_TREE_HPP
#define QR_CDO_TREE_HPP

#include "edginc/SpreadLossTree.hpp"
#include "edginc/GeneralisedCDO.hpp"
DRLIB_BEGIN_NAMESPACE


/////////////////////////////////////////////////////////////////////////////////////
// CDOTree class
////////////////////////////////////////////////////////////////////////////////////
class CDOTree : public FDProduct {
public:

	/******************** variables ********************/
	GeneralisedCDOConstSP inst;

	SpreadLossTree*  crTree;

	/******************** methods ********************/
	// constructor
	CDOTree(GeneralisedCDOConstSP inst, FDModel* model);

	/** initialisation, called ONCE only after InitState() for each new model instance */
	virtual void init(Control*) const;

	/** initialising and setting product variables */
	// this is called per pricing call before each pricing
	virtual void initProd(void);

	/** update payoff slice at each tree time step */
	virtual void update(int & step, FDProduct::UpdateType);

	/** get value of product at each tree node at model timeline step */
	virtual const TreeSlice & getValue(int step) const;

	virtual const TreeSlice & getValue( int step, DateTime eventDate ) const
	{
		return getValue(step);
	}

   
	virtual DateTime getStartDate(void) const { return model->getDate(0);}

	virtual string getOutputName() const;

	void printInfo(ostream& outputStream) const;

	virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results); 
	
private:

	// FIELDS ///////////////////

	/** slice that holds output obtained with getValue(step) */
	TreeSliceSP getValueSlice;

    /** slice to hold conditional contLegValues */
    TreeSliceSP condContLegSlice;

    /** slice to hold conditional fee leg values */
    TreeSliceSP condFeeLegSlice;

	// to help debugging and know which product instance this is, take copy of the instrument's
	// outputName field as the debuggers are unable to get from instrumentSP in the watch windows
	string instOutputName;  

	/** protection leg product */
	FDProductSP protLegProd;

	/** fee leg product */
	FDProductSP feeLegProd;

    /** flag to indicate whether conditional fwds are to be calculated */
    // mutable since has to be set in init()
    mutable bool calcCondFwds;

    /** date at which to calc fwds. default is set to valueDate */
    // mutable since has to be set in init()
    mutable DateTime condFwdDate;

    /** step in timeline corresponding to condFwdDate */
    int condFwdStep;

    /** value date */
    DateTime valueDate;

};

DRLIB_END_NAMESPACE

#endif


