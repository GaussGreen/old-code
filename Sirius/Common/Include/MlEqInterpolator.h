#ifndef INTERPOLATORH
#define INTERPOLATORH

#include "cMatrix.h"
#include "edcubicspline.h"
#include "smart.h"
#include "MLEqMaths.h"


using namespace std;



/***************************************************************
**	Class   : MlEqInterpolator 
**	Function: 
**	Returns : 
**	Comment : 
****************************************************************/

class MlEqInterpolator :public RCObject
{

protected:

	GVector < CVector >  m_xData;//[whichdata][idata]
	GVector < CVector >  m_yData;//[whichdata][idata]

	int		m_lagPoints;////how many points before/after needed for this interpolation to work
	
	
public:
	
	MlEqInterpolator(int lagPoints=1);
	MlEqInterpolator(const GVector <CVector>& xData,const GVector<CVector>& yData,int lagPoints=1);
	MlEqInterpolator(const CVector& xData,const GVector<CVector>& yData,int lagPoints=1);
	MlEqInterpolator(const CVector& xData,const CVector& yData);
	MlEqInterpolator(const CVector& xData,const CVector& yData,int n);


	MlEqInterpolator(const MlEqInterpolator& rhs);
	
	virtual void initialize(const GVector<CVector>& xData,const GVector<CVector>& yData);
	virtual void initialize(const CVector& xData,const CVector& yData);
	void		 initialize(const CVector& xData,const CVector& yData,int n);
	
	virtual double getValue(double xVal,int whichData=0) ;
	virtual void   getValues(CVector& yVals,double xVal) ;

	virtual const CVector& getXData(int whichData=0)const;
	virtual const CVector& getYData(int whichData=0)const;
	
	virtual getNumberOfSlices()const{return m_yData.getsize();};
	virtual int lagPoints()const{return m_lagPoints;};

	InterpolatorTypeEnum m_type;

	virtual MlEqInterpolatorHandle copy() const{ return new MlEqInterpolator( *this ); };
	virtual void reinitialize(const CVector& xData,const CVector& yData,int whichData);

	virtual getNumberOfSlices(){return m_yData.getsize();};
	virtual getNumberOfXDataSlices(){return m_xData.getsize();};
	
};	


/***************************************************************
**	Class   : MlEqCubicSplineInterpolator 
**	Function: 
**	Returns : 
**	Comment : 
****************************************************************/

class MlEqCubicSplineInterpolator :public MlEqInterpolator
{

protected:

	vector<EDCubicSpline> m_cubicSpline;

public:

	MlEqCubicSplineInterpolator(){};

	MlEqCubicSplineInterpolator( 
					const GVector<CVector>& 	xData,			// Input Data, on the x dimension 
					const GVector<CVector>& 	yData,			// Input Data, on the y dimension
					int					addTanhWings,	// addTanhWings
					double				cL,				// cL
					double				cR,				// cR
					int					yPower			// yPower
					);

	MlEqCubicSplineInterpolator(const MlEqCubicSplineInterpolator& rhs);

	virtual MlEqInterpolatorHandle copy() const{ return new MlEqCubicSplineInterpolator( *this ); };

	virtual double getValue(double xVal,int whichData=0);
	virtual getNumberOfSlices(){return m_cubicSpline.size();};

	virtual void initialize(vector<EDCubicSpline>&m);
	virtual void initialize(GVector<CVector>& xData,GVector<CVector>& yData);
	virtual void initialize(const CVector& xData,const CVector& yData);

	virtual void reinitialize(const CVector& xData,const CVector& yData,int whichData);

	virtual const CVector& getXData(int whichData=0)const;
	virtual const CVector& getYData(int whichData=0)const;

};



/***************************************************************
**	Class   : MlEqMonotonicSplineInterpolator 
**	Function: 
**	Returns : 
**	Comment : 
****************************************************************/

class MlEqMonotonicSplineInterpolator :public MlEqInterpolator
{

protected:

	vector<MlEqMonoSpline> m_monoSpline;

public:

	MlEqMonotonicSplineInterpolator(){};
	MlEqMonotonicSplineInterpolator( 
					const GVector<CVector>& 	xData,			// Input Data, on the x dimension 
					const GVector<CVector>& 	yData,			// Input Data, on the y dimension
					int					addTanhWings,	// addTanhWings
					double				cL,				// cL
					double				cR,				// cR
					int					yPower			// yPower
					);


	virtual void initialize(const GVector<CVector>& 	xData,			// Input Data, on the x dimension 
					const GVector<CVector>& 	yData,			// Input Data, on the y dimension
					int					addTanhWings,	// addTanhWings
					double				cL,				// cL
					double				cR,				// cR
					int					yPower			// yPower
					);

	MlEqMonotonicSplineInterpolator(const MlEqMonotonicSplineInterpolator& rhs);

	virtual MlEqInterpolatorHandle copy() const{ return new MlEqMonotonicSplineInterpolator( *this ); };

	virtual double getValue(double xVal,int whichData=0) ;
	virtual getNumberOfSlices(){return m_monoSpline.size();};

	virtual void initialize(vector<MlEqMonoSpline>&m);
	virtual void initialize(const GVector<CVector>& xData,const GVector<CVector>& yData);
	virtual void initialize(const CVector& xData,const CVector& yData);

	virtual void reinitialize(const CVector& xData,const CVector& yData,int whichData);


	virtual const CVector& getXData(int whichData=0)const;
	virtual const CVector& getYData(int whichData=0)const;
};

/***************************************************************
**	Class   : MlEqMonotonicSplineInterpolator 
**	Function: 
**	Returns : 
**	Comment : 
****************************************************************/

class MlEqConstantInterpolator :public MlEqMonotonicSplineInterpolator
{


public:

	MlEqConstantInterpolator(){};

	MlEqConstantInterpolator(
					const GVector<CVector>& 	xData,			// Input Data, on the x dimension 
					const GVector<CVector>& 	yData,			// Input Data, on the y dimension
					int					addTanhWings,	// addTanhWings
					double				cL,				// cL
					double				cR,				// cR
					int					yPower			// yPower
					);

	MlEqConstantInterpolator(
					const GVector<CVector>& 	xData,			// Input Data, on the x dimension 
					const GVector<CVector>& 	yData			// Input Data, on the y dimension
					);

	MlEqConstantInterpolator(const MlEqConstantInterpolator& rhs);

	virtual MlEqInterpolatorHandle copy() const{ return new MlEqConstantInterpolator( *this ); };

	virtual double getValue(double xVal,int whichData=0) ;


};


/***************************************************************
**	Class   : MlEqLinearInterpolator 
**	Function: 
**	Returns : 
**	Comment : 
****************************************************************/

class MlEqLinearInterpolator :public MlEqMonotonicSplineInterpolator
{

public:

	MlEqLinearInterpolator(){};

	MlEqLinearInterpolator(
					const GVector<CVector>& 	xData,			// Input Data, on the x dimension 
					const GVector<CVector>& 	yData,			// Input Data, on the y dimension
					int					addTanhWings,	// addTanhWings
					double				cL,				// cL
					double				cR,				// cR
					int					yPower			// yPower
					);

	MlEqLinearInterpolator(
					const GVector<CVector>& 	xData,			// Input Data, on the x dimension 
					const GVector<CVector>& 	yData			// Input Data, on the y dimension
					);

	MlEqLinearInterpolator(const MlEqLinearInterpolator& rhs);

	virtual MlEqInterpolatorHandle copy() const{ return new MlEqLinearInterpolator( *this ); };
	virtual double getValue(double xVal,int whichData=0) ;
};



class MlEq2DInterpolator : public MlEqInterpolator
{

	MlEqInterpolatorHandle m_interpolateAcrossY;
	MlEqInterpolatorHandle m_interpolateAcrossYTemp;

	GVector< MlEqInterpolatorHandle > m_interpolateAcrossX;

	void Locate(const CVector& xx, const double x, int& j);

public:

	void	initialize(	MlEqInterpolatorHandle interpolateAcrossY,
						MlEqInterpolatorHandle interpolateAcrossX);

	void	initialize(	MlEqInterpolatorHandle interpolateAcrossY,
						vector < MlEqInterpolatorHandle > & interpolateAcrossX);

	MlEq2DInterpolator(	MlEqInterpolatorHandle interpolateAcrossY,
						MlEqInterpolatorHandle interpolateAcrossX);

	MlEq2DInterpolator(	MlEqInterpolatorHandle interpolateAcrossY,
						vector < MlEqInterpolatorHandle > & interpolateAcrossX);


	MlEq2DInterpolator(const MlEq2DInterpolator& rhs);

	MlEqConstInterpolatorHandle getXInterpolator(int whichslice);

	virtual MlEqInterpolatorHandle copy() const{ return new MlEq2DInterpolator( *this ); };
	virtual void reinitialize(const CVector& xData,const CVector& yData,int whichData);

	double	getValue(double xVal,double yVal);
};


/***************************************************************
**	Class   : MlEqAnalyticCurveWithTanhWingEdge 
**	Function: 
**	Returns : 
**	Comment : 
****************************************************************/

class MlEqAnalyticCurveWithTanhWingEdge : public MlEqInterpolator
{

//  base class implementation corresponds to a polynomial curve

private:

	CVector m_coeff;

protected:


	MlEqInterpolatorHandle m_interp;

	double	m_lowerEdge;
	double	m_upperEdge;
	bool	m_addTanhWing; 
	CVector m_initialGuess;
	CMatrix m_xValBounds;
	double  m_cL;
	double  m_cR;
	double  m_finalTol;
	double  m_stopTol;
	double  m_yPower;
	int     m_npoints;
	CVector m_xfitData;


public:


	const CVector& getInitialGuess() const{return m_initialGuess;};
	const CVector& getCoeff() const{return m_coeff;};
	void setCoeff(CVector& newCoeff);


	MlEqAnalyticCurveWithTanhWingEdge(CVector& xfitData,CVector& yfitData,CVector& initialGuess,CMatrix& xValBounds,int addTanhWing,double cL,double cR,double yPower,double finalTol,double stopTol);
	MlEqAnalyticCurveWithTanhWingEdge(double lowerEdge,double upperEdge,int addTanhWing,double cL,double cR,double yPower,int npoints,CVector& coeff);
	MlEqAnalyticCurveWithTanhWingEdge(){};

	void initialize(MlEqInterpolatorHandle& interp);
	void initialize(double lowerEdge,double upperEdge,int addTanhWing,double cL,double cR,double yPower,int npoints,CVector& coeff);
	void initialize(CMatrix& fitData,CVector& initialGuess,CMatrix& xValBounds,double lowerEdge,double upperEdge,int addTanhWing,double cL,double cR,double yPower,double finalTol=1e-12,double stopTol=1e-12);
	void initialize(const CVector& xfitData,const CVector& yfitData,CVector& initialGuess,CMatrix& xValBounds,int addTanhWing,double cL,double cR,double yPower,double finalTol,double stopTol);
	void initialize(const CVector& xfitData,CVector& initialGuess,CMatrix& xValBounds,int addTanhWing,double cL,double cR,double yPower,double finalTol,double stopTol);// fit is done later


	void reinitialize(const CVector& coeff);
	virtual void reinitialize(const CVector& xData,const CVector& yData,int whichData);

	void FittoData(const CVector& xData,const CVector& yData,CVector& initial_guess,CMatrix& limits,double fitTol,double stopTol);
	void FittoData( CMatrix& fitData,CVector& initialGuess,CMatrix& xValBounds,double finalTol,double stopTol);

	virtual double getAnalyticValueInside(double xVal);
	double  getValue(double xVal,int whichData=0) ;
	void    getValues(CVector& yVals,double xVal); 

};


#endif

