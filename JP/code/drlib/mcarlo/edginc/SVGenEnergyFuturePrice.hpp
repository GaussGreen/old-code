// ---------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : SVGenEnergyFuturePrice.hpp
//
//   Description : State Variable for Energy Futures prices
//
//   Author      : Spyridon Schismenos
//
//   Date        : June 22, 2005
//
// ---------------------------------------------------------


#ifndef EDR_MCENERGYFUTUREPRICE_HPP
#define EDR_MCENERGYFUTUREPRICE_HPP

#include "edginc/ElemStateVariableGen.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/EnergyContractLabel.hpp"
#include "edginc/DoubleMatrix.hpp"

DRLIB_BEGIN_NAMESPACE

class MCARLO_DLL SVGenEnergyFuturePrice:   virtual public IElemStateVariableGen,
                                        public virtual VirtualDestructorBase
{
public:

	// State Variable interface
	class MCARLO_DLL IStateVar: public virtual IStateVariable
	{
	public:
		virtual ~IStateVar() {}
        virtual double getSpot() = 0;
		virtual DoubleMatrix getFullPath() = 0;
	};
    typedef smartPtr<IStateVar> IStateVarSP;

	SVGenEnergyFuturePrice(const DateTimeArray& dates,const EnergyContractLabelArray& labels);

    class MCARLO_DLL PastSV : public IStateVar {
    public:
        double getSpot() {  return 0; }
        bool doingPast() const { return true; }
		DoubleMatrix getFullPath() {
			DoubleMatrix result;
			result=DoubleMatrix(1,1);
			result[0][0]=1.0;
			return result;            //no use, just implemented
		}
    };

    virtual IStateVariableSP create(IStateVariableSP oldStateVar,
                                    IStateVariableGen::IStateGen* pathGen) const;
    IStateVarSP getEnergyFuturePriceSV(IStateVariableGen::IStateGen* pathGen) const;
	DateTimeArray getPriceDates() const;
    void attachSVGen(IElemStateVariableGenVisitor*) const;


    // used with past path generator
    IStateVarSP pastSV() const;
private:
	DateTimeArray                 priceDates;
	EnergyContractLabelArray      contractLabels;
};

typedef smartPtr<SVGenEnergyFuturePrice> SVGenEnergyFuturePriceSP;

DRLIB_END_NAMESPACE

#endif
