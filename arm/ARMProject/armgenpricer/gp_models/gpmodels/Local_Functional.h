
#ifndef _INGPMODELS_LOCAL_FUNCTIONAL_H
#define _INGPMODELS_LOCAL_FUNCTIONAL_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "gpbase/rootobject.h"
#include "gpbase/gpvector.h"
#include "gpbase/typedef.h"


CC_BEGIN_NAMESPACE( ARM )

class ARM_DensityFunctor;

class ARM_HestonModel_Fx;

class ARM_LocalFunctional : public ARM_RootObject
{

private:
	ARM_VectorPtrVector					itsFunctionals;
	ARM_VectorPtrVector					itsGrids;
	std::vector<double>						itsFwds;
	std::vector<double>						itsVols;
	std::vector<double>						itsResetTimes;
	int									itsGridSize;
	int									itsNbStdDev;

public:
	///constructors/destructors
	ARM_LocalFunctional (int sizeGrid = 501,double nbStdDev = 6.);
	ARM_LocalFunctional (const ARM_LocalFunctional& rhs);
	ASSIGN_OPERATOR(ARM_LocalFunctional)
	virtual ~ARM_LocalFunctional() {};
	ARM_Object* Clone() const { return new ARM_LocalFunctional(*this); };

	/// utilities
	virtual string toString(const string& indent="",const string& nextIndent="") const;

	virtual void Calibrate (
		double resetTime,
		double fwd,
		double vol,
		ARM_DensityFunctor& density,
		bool rescalling);

	virtual void CalibrateHestonFx (
		double resetTime,
		double fwd,
		const ARM_HestonModel_Fx* hestonModelFx,
		ARM_DensityFunctor& density);

	inline void setFunc(ARM_GP_VectorPtr func)			{ itsFunctionals.push_back(func); }
	inline void setTime(double time)					{ itsResetTimes.push_back(time); }

	inline const ARM_GP_VectorPtr& Grid(size_t idx) const			{ return itsGrids[idx]; }
	inline const std::vector<double>& ResetTimes() const		{ return itsResetTimes; }
	inline size_t GetGridSize()	const					{ return itsGridSize; }
	inline double GetGrid(size_t idx1, size_t idx2)	const				{ return (*itsGrids[idx1])[idx2]; }

	inline void setVol ( double vol )					{ itsVols.push_back(vol); }
	inline void setFwd ( double fwd )					{ itsFwds.push_back(fwd); }
	inline double GetFwd(size_t i) const				{ return itsFwds[i]; }
	inline double GetVol(size_t i) const				{ return itsVols[i]; }

	inline double FuncValue(size_t idx,size_t j) const
											{ return (*itsFunctionals[idx])[j]; }

	virtual ARM_VectorPtr Func(double evalTime,const ARM_GP_VectorPtr& values) const;
};



CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
