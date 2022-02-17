
#include "stdafx.h"
#include "MlEqInterpolator.h"
#include "MlEqMaths.h"
#include "utility.h"
#include "smart.h"


/***************************************************************
**	Class   : MlEqInterpolator 
**	Function: 
**	Returns : nothing
**	Comment : constructor
****************************************************************/

MlEqInterpolator::MlEqInterpolator(int lagPoints)
{
	m_lagPoints = lagPoints;
	m_type = Linear;
}


/***************************************************************
**	Class   : MlEqInterpolator 
**	Function: 
**	Returns : nothing
**	Comment : constructor
****************************************************************/


MlEqInterpolator::MlEqInterpolator(const GVector<CVector>& xData,const GVector<CVector>& yData,int lagPoints)
{
	m_lagPoints = lagPoints;
	m_type		= Linear;

	initialize(xData,yData);
}

/***************************************************************
**	Class   : MlEqInterpolator 
**	Function: initialize
**	Returns : nothing
**	Comment : initialization routine
****************************************************************/


void MlEqInterpolator::initialize(const CVector& xData,const CVector& yData)
{
// constructor for one timeslice only


	m_xData.resize(1);;
	m_yData.resize(1);

	m_xData[0] = xData;
	m_yData[0] = yData;

	m_lagPoints = 1;
}

/***************************************************************
**	Class   : MlEqInterpolator 
**	Function: initialize
**	Returns : nothing
**	Comment : initialization routine
****************************************************************/

void MlEqInterpolator::reinitialize(const CVector& xData,const CVector& yData,int whichData)
{


	if ( m_xData.getsize() > 1 ){

		if ( m_xData.getsize()-1 < whichData ){
			throw("incorrect dimensioning encountered in re-initializing interpolaotr data");
		}

		m_xData[whichData]	= xData;
	}
	else{
		m_xData[0]	= xData;
	}


	if ( m_yData.getsize() > 1 ){

		if ( m_yData.getsize()-1 < whichData ){
			throw("incorrect dimensioning encountered in re-initializing interpolaotr data");
		}

		m_yData[whichData]	= yData;
	}
	else{
		m_yData[0]	= yData;
	}


}

/***************************************************************
**	Class   : MlEqInterpolator 
**	Function: initialize
**	Returns : nothing
**	Comment : initialization routine
****************************************************************/


void MlEqInterpolator::initialize(const CVector& xData,const CVector& yData,int n)
{
// constructor for one timeslice only
// take only first n out of array


	if ( n+1 >= yData.getsize() ){
		initialize(xData,yData);
		return;
	}

	m_xData.resize(1);;
	m_yData.resize(1);

	m_xData[0].resize(n+1);
	m_yData[0].resize(n+1);

	for ( int i = 0 ; i <= n ; i++ )
	{
		m_xData[0][i] = xData[i];
		m_yData[0][i] = yData[i];
	}

	m_lagPoints = 1;
}


/***************************************************************
**	Class   : MlEqInterpolator 
**	Function: MlEqInterpolator
**	Returns : nothing
**	Comment : initialization routine
****************************************************************/


MlEqInterpolator::MlEqInterpolator(const CVector& xData,const CVector& yData)
{
	initialize(xData,yData);
}

/***************************************************************
**	Class   : MlEqInterpolator 
**	Function: MlEqInterpolator
**	Returns : nothing
**	Comment : initialization routine
****************************************************************/


MlEqInterpolator::MlEqInterpolator(const CVector& xData,const CVector& yData,int n)
{
	initialize(xData,yData,n);
}

/***************************************************************
**	Class   : MlEqInterpolator 
**	Function: initialize
**	Returns : nothing
**	Comment : initialization routine
****************************************************************/


void MlEqInterpolator::initialize(const GVector<CVector>& xData,const GVector<CVector>& yData)
{

	if ( !( xData.getsize() == yData.getsize()  || xData.getsize() == 1 ))
	{
		throw("inconsistent data encountered in interpolator object");
	}

	int whichXData=0;
	for ( int i = 0; i < xData.getsize(); i++ )
	{
		if ( xData.getsize() > 1 )
		{
			whichXData = i;
		}

		if ( xData[whichXData].getsize() != yData[i].getsize() )
		{
			throw(" inconsistent data encountered in interpolator");
		}
	}

// check if inputs are simply 2 column vectors and transpose in this case

	int columnFlag = 1;
	for ( int i = 0; i < xData.getsize(); i++ )
	{
		if ( xData[i].getsize() != 1 )
		{
			columnFlag = 0;
			break;
		}
	
	}

	if ( columnFlag ) 
	{
		m_xData.resize(1);
		m_yData.resize(1);

		m_xData[0].resize(xData.getsize());
		m_yData[0].resize(xData.getsize());

		for ( int i = 0 ; i < xData.getsize(); i++ )
		{
			m_xData[0][i] =	xData[i][0]; 
			m_yData[0][i] = yData[i][0];
		}

	}
	else
	{
		m_xData = xData;
		m_yData = yData;
	}
}


/***************************************************************
**	Class   : MlEqInterpolator
**	Function: getXData
**	Returns : vector of x data
**	Comment : 
****************************************************************/


const CVector& MlEqInterpolator::getXData(int whichData)const
{
	const CVector& xData =  getObjectFromVector(m_xData,whichData);//sos

	return xData;

}


/***************************************************************
**	Class   : MlEqInterpolator 
**	Function: getYData
**	Returns : vector of y data
**	Comment : 
****************************************************************/


const CVector& MlEqInterpolator::getYData(int whichData)const
{
	return m_yData[whichData];
}


/***************************************************************
**	Class   : MlEqInterpolator  
**	Function: getValue
**	Returns : interpolated value
**	Comment : 
****************************************************************/


double MlEqInterpolator::getValue(double xVal,int whichData) 
{
	int method = 0;
	int whichXData = whichData;
	if ( m_xData.getsize() == 1 )
	{
		whichXData = 0;
	}
	double val = MlEqMaths::linearInterp(m_xData[whichXData], m_yData[whichData],xVal,method);
	return val;
}

/***************************************************************
**	Class   : MlEqInterpolator  
**	Function: getValues
**	Returns : interpolated value
**	Comment : 
****************************************************************/


void MlEqInterpolator::getValues(CVector& yVals,double xVal) 
{
	int N = getNumberOfSlices();
	yVals.resize(N);

	for ( int i = 0 ; i < N; i++ )
	{
		yVals[i] = getValue(xVal,i);
	}
}



/***************************************************************
**	Class   : MlEqCubicSplineInterpolator 
**	Function: MlEqCubicSplineInterpolator
**	Returns : 
**	Comment : 
****************************************************************/


MlEqCubicSplineInterpolator::MlEqCubicSplineInterpolator(
					const GVector<CVector>& 	xData,			// Input Data, on the x dimension 
					const GVector<CVector>& 	yData,			// Input Data, on the y dimension
					int					addTanhWings,	// addTanhWings
					double				cL,				// cL
					double				cR,				// cR
					int					yPower			// yPower
					)
{

	m_lagPoints = 2;
	MlEqInterpolator::initialize(xData,yData);

	int nsize = yData.getsize();

	if ( xData.getsize() != 1 )
	{
		if ( nsize != xData.getsize() )
		{
			throw ( " inconsistent xData,yData entered in cubicspline interpolator");
		}
	}

	int whichXdata=0;
	m_cubicSpline.resize(nsize);
	for ( int i = 0 ; i < nsize; i++ )
	{
		if ( xData.getsize() != 1 )
		{
			whichXdata = i;
		}

		m_cubicSpline[i].initialize(xData[whichXdata],yData[i],addTanhWings,cL,cR,yPower);
	}

}


/***************************************************************
**	Class   : MlEqCubicSplineInterpolator 
**	Function: initialize
**	Returns : 
**	Comment : 
****************************************************************/

void MlEqCubicSplineInterpolator::initialize(GVector<CVector>& xData,GVector<CVector>& yData)
{

	m_type		= CubicSpline;

	if ( xData.getsize() == 0 || m_cubicSpline.size() == 0 )
	{
		throw("at least one cubicspline interpolator must be defined");
		return;
	}

	if ( yData.getsize() != xData.getsize() )
	{
		if ( xData.getsize() != 1 )
		{
			throw ( " inconsistent xData,yData entered in cubicspline interpolator");
		}
	}


	int addTanhWing =	m_cubicSpline[0].addTanhWing();
	double cL		=	m_cubicSpline[0].cL();
	double cR		=	m_cubicSpline[0].cR();
	double yPower	=	m_cubicSpline[0].yPower();

	MlEqInterpolator::initialize(xData,yData);

	m_cubicSpline.resize(xData.getsize());

	for ( int i = 0 ; i < xData.getsize(); i++ )
	{
		m_cubicSpline[i].initialize(xData[i],yData[i],addTanhWing,cL,cR,yPower);
	}

}

/***************************************************************
**	Class   : MlEqCubicSplineInterpolator 
**	Function: initialize
**	Returns : 
**	Comment : 
****************************************************************/

void MlEqCubicSplineInterpolator::initialize(const CVector& xData,const CVector& yData)
{

	m_cubicSpline.resize(1);
	GVector<CVector> vxData(1),vyData(1);

	vxData[0] = xData;
	vyData[0] = yData;

	initialize(vxData,vyData);
}


/***************************************************************
**	Class   : MlEqCubicSplineInterpolator 
**	Function: initialize
**	Returns : 
**	Comment : 
****************************************************************/

void MlEqCubicSplineInterpolator::reinitialize(const CVector& xData,const CVector& yData,int whichData)
{

	if ( m_xData.getsize() > 1 ){

		if ( m_xData.getsize()-1 < whichData ){
			throw("incorrect dimensioning encountered in re-initializing interpolaotr data");
		}

		m_xData[whichData]	= xData;
	}
	else{
		m_xData[0]	= xData;
	}


	if ( m_yData.getsize() > 1 ){

		if ( m_yData.getsize()-1 < whichData ){
			throw("incorrect dimensioning encountered in re-initializing interpolaotr data");
		}

		m_yData[whichData]	= yData;
	}
	else{
		m_yData[0]	= yData;
	}



	MlEqInterpolator::reinitialize(xData,yData,whichData);

	int n;
	if ( m_cubicSpline.size() > 1 )
	{
		if ( m_cubicSpline.size() < whichData){
			throw("indexing problem encountered in reinitialising interpolator");
		}

		n = whichData;
	}
	else
	{
		n = 0;
	}


	m_cubicSpline[n].initialize(xData,yData,m_cubicSpline[n].addTanhWing (),m_cubicSpline[n].cL(),m_cubicSpline[n].cR(),m_cubicSpline[n].yPower());


}

/***************************************************************
**	Class   : MlEqCubicSplineInterpolator 
**	Function: MlEqCubicSplineInterpolator
**	Returns : 
**	Comment : 
****************************************************************/


MlEqCubicSplineInterpolator::MlEqCubicSplineInterpolator(const MlEqCubicSplineInterpolator& rhs)
:
MlEqInterpolator( rhs ),
m_cubicSpline(rhs.m_cubicSpline )
{
	m_type = rhs.m_type;
}


/***************************************************************
**	Class   : MlEqCubicSplineInterpolator 
**	Function: getXData
**	Returns : vector of xData
**	Comment : 
****************************************************************/

const CVector& MlEqCubicSplineInterpolator::getXData(int whichData)const
{
	return m_cubicSpline[whichData].getXdata();
}


/***************************************************************
**	Class   : MlEqCubicSplineInterpolator 
**	Function: getYData
**	Returns : Vector of yData
**	Comment : 
****************************************************************/


const CVector& MlEqCubicSplineInterpolator::getYData(int whichData)const
{
	return m_cubicSpline[whichData].getYdata();
}



/***************************************************************
**	Class   : MlEqCubicSplineInterpolator 
**	Function: initialize
**	Returns : nothing
**	Comment : 
****************************************************************/


void MlEqCubicSplineInterpolator::initialize(vector<EDCubicSpline>&m)
{
//  fill in data into base class as well sos
	m_cubicSpline	=	m;
}


/***************************************************************
**	Class   : MlEqCubicSplineInterpolator  
**	Function: getValue
**	Returns : interpolated value
**	Comment : 
****************************************************************/


double MlEqCubicSplineInterpolator::getValue(double xVal,int whichData) 
{

	double val = getObjectFromVector(m_cubicSpline,whichData).getValue(xVal); 

//ouble val = m_cubicSpline[whichData].getValue(xVal);

//	double val = cubicSpline.getValue(xVal);
	return val;	

}






/***************************************************************
**	Class   : MlEqMonotonicSplineInterpolator 
**	Function: MlEqMonotonicSplineInterpolator
**	Returns : 
**	Comment : 
****************************************************************/

MlEqMonotonicSplineInterpolator::MlEqMonotonicSplineInterpolator(
					const GVector<CVector>& 	xData,			// Input Data, on the x dimension 
					const GVector<CVector>& 	yData,			// Input Data, on the y dimension
					int					addTanhWings,	// addTanhWings
					double				cL,				// cL
					double				cR,				// cR
					int					yPower			// yPower
					)
{
	initialize(	xData,yData,addTanhWings,cL,cR,yPower);
}


/***************************************************************
**	Class   : MlEqMonotonicSplineInterpolator 
**	Function: initialize
**	Returns : 
**	Comment : 
****************************************************************/


void MlEqMonotonicSplineInterpolator::initialize(	const GVector<CVector>& 	xData,			// Input Data, on the x dimension 
					const GVector<CVector>& 	yData,			// Input Data, on the y dimension
					int					addTanhWings,	// addTanhWings
					double				cL,				// cL
					double				cR,				// cR
					int					yPower			// yPower
					)

{

	m_lagPoints = 2;
	MlEqInterpolator::initialize(xData,yData);

	int nsize = yData.getsize();

	if ( xData.getsize() != 1 )
	{
		if ( nsize != xData.getsize() )
		{
			throw ( " inconsistent xData,yData entered in monotonic spline interpolator");
		}
	}

	int whichXdata=0;
	m_monoSpline.resize(nsize);
	for ( int i = 0 ; i < nsize; i++ )
	{
		if ( xData.getsize() != 1 )
		{
			whichXdata = i;
		}

		m_monoSpline[i].initialize(xData[whichXdata],yData[i],addTanhWings,cL,cR,yPower);
	}

}


/***************************************************************
**	Class   : MlEqMonotonicSplineInterpolator 
**	Function: initialize
**	Returns : 
**	Comment : 
****************************************************************/

void MlEqMonotonicSplineInterpolator::initialize(const GVector<CVector>& xData,const GVector<CVector>& yData)
{

	m_type		= MonotonicSpline;

	if ( xData.getsize() == 0 || m_monoSpline.size() == 0 )
	{
		throw("at least one monotonic spline interpolator must be defined");
		return;
	}

	if ( yData.getsize() != xData.getsize() )
	{
		if ( xData.getsize() != 1 )
		{
			throw ( " inconsistent xData,yData entered in monotonic spline interpolator");
		}
	}


	int addTanhWing =	m_monoSpline[0].addTanhWing();
	double cL		=	m_monoSpline[0].cL();
	double cR		=	m_monoSpline[0].cR();
	double yPower	=	m_monoSpline[0].yPower();

	MlEqInterpolator::initialize(xData,yData);

	m_monoSpline.resize(xData.getsize());

	for ( int i = 0 ; i < xData.getsize(); i++ ){
		m_monoSpline[i].initialize(xData[i],yData[i],addTanhWing,cL,cR,yPower);
	}

}

/***************************************************************
**	Class   : MlEqMonotonicSplineInterpolator 
**	Function: initialize
**	Returns : 
**	Comment : 
****************************************************************/

void MlEqMonotonicSplineInterpolator::initialize(const CVector& xData,const CVector& yData)
{

	m_monoSpline.resize(1);
	GVector<CVector> vxData(1),vyData(1);

	vxData[0] = xData;
	vyData[0] = yData;

	initialize(vxData,vyData);
}

/***************************************************************
**	Class   : MlEqCubicSplineInterpolator 
**	Function: initialize
**	Returns : 
**	Comment : 
****************************************************************/

void MlEqMonotonicSplineInterpolator::reinitialize(const CVector& xData,const CVector& yData,int whichData)
{
	if ( m_xData.getsize() > 1 ){

		if ( m_xData.getsize()-1 < whichData ){
			throw("incorrect dimensioning encountered in re-initializing interpolaotr data");
		}

		m_xData[whichData]	= xData;
	}
	else{
		m_xData[0]	= xData;
	}


	if ( m_yData.getsize() > 1 ){

		if ( m_yData.getsize()-1 < whichData ){
			throw("incorrect dimensioning encountered in re-initializing interpolaotr data");
		}

		m_yData[whichData]	= yData;
	}
	else{
		m_yData[0]	= yData;
	}



	MlEqInterpolator::reinitialize(xData,yData,whichData);

	int n;
	if ( m_monoSpline.size() > 1 )
	{
		if ( m_monoSpline.size() < whichData){
			throw("indexing problem encountered in reinitialising interpolator");
		}

		n = whichData;
	}
	else
	{
		n = 0;
	}

	m_monoSpline[n].initialize(xData,yData,m_monoSpline[n].addTanhWing (),m_monoSpline[n].cL(),m_monoSpline[n].cR(),m_monoSpline[n].yPower());

}


/***************************************************************
**	Class   : MlEqMonotonicSplineInterpolator 
**	Function: MlEqMonotonicSplineInterpolator
**	Returns : 
**	Comment : 
****************************************************************/


MlEqMonotonicSplineInterpolator::MlEqMonotonicSplineInterpolator(const MlEqMonotonicSplineInterpolator& rhs)
:
MlEqInterpolator( rhs ),
m_monoSpline(rhs.m_monoSpline )
{
	m_type = rhs.m_type;
}


/***************************************************************
**	Class   : MlEqMonotonicSplineInterpolator 
**	Function: getXData
**	Returns : vector of xData
**	Comment : 
****************************************************************/

const CVector& MlEqMonotonicSplineInterpolator::getXData(int whichData)const
{
	return m_monoSpline[whichData].getXdata();
}


/***************************************************************
**	Class   : MlEqMonotonicSplineInterpolator 
**	Function: getYData
**	Returns : Vector of yData
**	Comment : 
****************************************************************/


const CVector& MlEqMonotonicSplineInterpolator::getYData(int whichData)const
{
	return m_monoSpline[whichData].getYdata();
}



/***************************************************************
**	Class   : MlEqMonotonicSplineInterpolator 
**	Function: initialize
**	Returns : nothing
**	Comment : 
****************************************************************/


void MlEqMonotonicSplineInterpolator::initialize(vector<MlEqMonoSpline>&m)
{
//  fill in data into base class as well sos
	m_monoSpline	=	m;
}


/***************************************************************
**	Class   : MlEqMonotonicSplineInterpolator  
**	Function: getValue
**	Returns : interpolated value
**	Comment : 
****************************************************************/

/*
double MlEqMonotonicSplineInterpolator::getValue(double xVal,int whichData) 
{
	MlEqMonoSpline monoSpline  = getObjectFromVector(m_monoSpline,whichData); 

	double val = monoSpline.getValue(xVal);
	return val;
}
*/


double MlEqMonotonicSplineInterpolator::getValue(double xVal,int whichData) 
{
	double val = getObjectFromVector(m_monoSpline,whichData).getValue(xVal);
	return val;
}






/***************************************************************
**	Class   : MlEqConstantInterpolator  
**	Function: getValue
**	Returns : interpolated value
**	Comment : 
****************************************************************/


double MlEqConstantInterpolator::getValue(double xVal,int whichData) 
{

	int whichXData = whichData;
	if ( m_xData.getsize() == 1 ){
		whichXData = 0;
	}

	bool isEdge = false;
	double val =  MlEqMaths::constantInterp(m_xData[whichXData], m_yData[whichData],xVal,isEdge);

	if ( m_monoSpline.size() == 0 ){
		return val;
	}

	if ( isEdge )
	{
		// this is an edge point
		val  = getObjectFromVector(m_monoSpline,whichData).getValue(xVal);
		return val;
	}
	else
	{
		return val;
	}
	
}


/***************************************************************
**	Class   : MlEqConstantInterpolator  
**	Function: MlEqConstantInterpolator
**	Returns : 
**	Comment : 
****************************************************************/


MlEqConstantInterpolator::MlEqConstantInterpolator( 
					const GVector<CVector>& 	xData,			// Input Data, on the x dimension 
					const GVector<CVector>& 	yData,			// Input Data, on the y dimension
					int					addTanhWings,	// addTanhWings
					double				cL,				// cL
					double				cR,				// cR
					int					yPower			// yPower
					)
:MlEqMonotonicSplineInterpolator(xData,yData,addTanhWings,cL,cR,yPower)
{}


/***************************************************************
**	Class   : MlEqConstantInterpolator  
**	Function: MlEqConstantInterpolator
**	Returns : 
**	Comment : 
****************************************************************/


MlEqConstantInterpolator::MlEqConstantInterpolator(const MlEqConstantInterpolator& rhs)
:
MlEqMonotonicSplineInterpolator( rhs )
{
	m_type = rhs.m_type;
}


/***************************************************************
**	Class   : MlEqConstantInterpolator  
**	Function: MlEqConstantInterpolator
**	Returns : 
**	Comment : 
****************************************************************/
MlEqConstantInterpolator::MlEqConstantInterpolator(
					const GVector<CVector>& 	xData,			// Input Data, on the x dimension 
					const GVector<CVector>& 	yData			// Input Data, on the y dimension
					)
{
	MlEqInterpolator::initialize(xData,yData);
}




/***************************************************************
**	Class   : MlEqLinearInterpolator  
**	Function: getValue
**	Returns : interpolated value
**	Comment : 
****************************************************************/


double MlEqLinearInterpolator::getValue(double xVal,int whichData) 
{

	int whichXData = whichData;
	if ( m_xData.getsize() == 1 ){
		whichXData = 0;
	}

	bool isEdge = false;
	int method = 0;

	double val = MlEqMaths::linearInterp(m_xData[whichXData], m_yData[whichData],xVal,method,isEdge);

	if ( m_monoSpline.size() == 0 ){
		return val;
	}

	if ( isEdge )
	{
		// this is an edge point
		val  = getObjectFromVector(m_monoSpline,whichData).getValue(xVal);
		return val;
	}
	else
	{
		return val;
	}
	
}


/***************************************************************
**	Class   : MlEqLinearInterpolator  
**	Function: MlEqLinearInterpolator
**	Returns : 
**	Comment : 
****************************************************************/


MlEqLinearInterpolator::MlEqLinearInterpolator( 
					const GVector<CVector>& 	xData,			// Input Data, on the x dimension 
					const GVector<CVector>& 	yData,			// Input Data, on the y dimension
					int					addTanhWings,	// addTanhWings
					double				cL,				// cL
					double				cR,				// cR
					int					yPower			// yPower
					)
:MlEqMonotonicSplineInterpolator(xData,yData,addTanhWings,cL,cR,yPower)
{}


/***************************************************************
**	Class   : MlEqLinearInterpolator  
**	Function: MlEqLinearInterpolator
**	Returns : 
**	Comment : 
****************************************************************/


MlEqLinearInterpolator::MlEqLinearInterpolator(const MlEqLinearInterpolator& rhs)
:
MlEqMonotonicSplineInterpolator( rhs )
{
	m_type = rhs.m_type;
}


/***************************************************************
**	Class   : MlEqLinearInterpolator  
**	Function: MlEqLinearInterpolator
**	Returns : 
**	Comment : 
****************************************************************/
MlEqLinearInterpolator::MlEqLinearInterpolator(
					const GVector<CVector>& 	xData,			// Input Data, on the x dimension 
					const GVector<CVector>& 	yData			// Input Data, on the y dimension
					)
{
	MlEqInterpolator::initialize(xData,yData);
}


/***************************************************************
**	Class   : MlEq2DInterpolator  
**	Function: MlEq2DInterpolator
**	Returns : 
**	Comment : 
****************************************************************/


double	MlEq2DInterpolator::getValue(double xVal,double yVal)
{

	// find where this date is relative to m_date
	// get the m_dates that is of interest based on the lagPoints

	int whichPosition = -1;
	const CVector& yDat = m_interpolateAcrossY->getXData();

	Locate( yDat, yVal, whichPosition );

	int lapPoints = m_interpolateAcrossY->lagPoints();

	// do the m_InterpAcrossY
	int startPoint = 0;
	int endPoint = yDat.getsize() - 1;

	startPoint	= MlEqMaths::Max(whichPosition - lapPoints, startPoint );
	endPoint	= MlEqMaths::Min(whichPosition + lapPoints, endPoint );	
	
	int numPoints = endPoint - startPoint + 1;

	CVector yData(numPoints);
	CVector zData(numPoints);

	for(int i = 0; i < numPoints; i++ ){
		yData[i] = yDat[ startPoint + i ] ;
	}

	for(int i = 0; i < numPoints; i++ ){

		if ( m_interpolateAcrossX.getsize() == 1 ){
			zData[i] =  m_interpolateAcrossX[0]->getValue(xVal,startPoint + i);
		}
		else{
			zData[i] =  m_interpolateAcrossX[startPoint + i]->getValue(xVal);
		}
	}

	m_interpolateAcrossYTemp->initialize( yData, zData );	
	double res = m_interpolateAcrossYTemp->getValue(yVal);

	return res;

}



/***************************************************************
**	Class   : MlEq2DInterpolator  
**	Function: MlEq2DInterpolator
**	Returns : 
**	Comment : 
****************************************************************/


void	MlEq2DInterpolator::initialize(	MlEqInterpolatorHandle interpolateAcrossY,
									    MlEqInterpolatorHandle interpolateAcrossX)
{
	m_interpolateAcrossY = interpolateAcrossY;
	m_interpolateAcrossYTemp=interpolateAcrossY->copy();
	m_interpolateAcrossX.resize(1); 
	m_interpolateAcrossX[0] = interpolateAcrossX;

	int n = interpolateAcrossY->getNumberOfSlices();
	const CVector& yDat = m_interpolateAcrossY->getXData();
}

/***************************************************************
**	Class   : MlEq2DInterpolator  
**	Function: MlEq2DInterpolator
**	Returns : 
**	Comment : 
****************************************************************/


void	MlEq2DInterpolator::initialize(	MlEqInterpolatorHandle interpolateAcrossY,
									    vector < MlEqInterpolatorHandle > & interpolateAcrossX)
{

	m_interpolateAcrossY = interpolateAcrossY;
	m_interpolateAcrossYTemp	= interpolateAcrossY->copy();
	m_interpolateAcrossX.resize(interpolateAcrossX.size());
	
	for ( int i = 0; i < interpolateAcrossX.size(); i++ ){
		m_interpolateAcrossX[i] = interpolateAcrossX[i];
	}

}

/***************************************************************
**	Class   : MlEq2DInterpolator  
**	Function: MlEq2DInterpolator
**	Returns : 
**	Comment : 
****************************************************************/


// searches array xx to find i such that xx[ii] <= x <= xx[i+1]	
void MlEq2DInterpolator::Locate(const CVector& xx, const double x, int& j)
{
	int ju,jm,jl;
	bool ascnd;

    int n=xx.getsize();
	jl=-1;						// initialise lower
	ju=n;						// and upper limits
	ascnd=(xx[n-1]>xx[0]);		// true if ascending order of table, false otherwise
	while (ju-jl > 1)			// if we are not yet done
	{
		jm=(ju+jl)>>1;			// compute a midpoint
		if(x>=xx[jm]==ascnd)	
			jl=jm;				// and replace either the lower limit
		else
			ju=jm;				// or upper limit, as appropriate
	}							// repeat until test condition satisfied
	
	if( x == xx[0] ) j = 0;		// set the output
	else if( x == xx[n-1] ) j = n-2;
	else j = jl;
}

/***************************************************************
**	Class   : MlEq2DInterpolator  
**	Function: MlEq2DInterpolator
**	Returns : 
**	Comment : 
****************************************************************/

MlEq2DInterpolator::MlEq2DInterpolator(	MlEqInterpolatorHandle interpolateAcrossY,
						MlEqInterpolatorHandle interpolateAcrossX)
{
	initialize(	interpolateAcrossY,interpolateAcrossX);
}


/***************************************************************
**	Class   : MlEq2DInterpolator  
**	Function: MlEq2DInterpolator
**	Returns : 
**	Comment : 
****************************************************************/


MlEq2DInterpolator::MlEq2DInterpolator(	MlEqInterpolatorHandle interpolateAcrossY,
					vector < MlEqInterpolatorHandle > & interpolateAcrossX)
{
	initialize(	interpolateAcrossY,
				interpolateAcrossX);
}


/***************************************************************
**	Class   : MlEqInterpolator  
**	Function: MlEqInterpolator
**	Returns : 
**	Comment : 
****************************************************************/
		
		
MlEqInterpolator::MlEqInterpolator(const MlEqInterpolator& rhs)
{			
	m_xData		= rhs.m_xData;
	m_yData		= rhs.m_yData;
	m_lagPoints	= rhs.m_lagPoints;
	m_type		= rhs.m_type;	
}		



/***************************************************************
**	Class   : MlEq2DInterpolator  
**	Function: MlEq2DInterpolator
**	Returns : 
**	Comment : 
****************************************************************/


MlEq2DInterpolator::MlEq2DInterpolator(const MlEq2DInterpolator& rhs)
{	
//  or should I make a deep copy here ???
	
	m_interpolateAcrossY		=	rhs.m_interpolateAcrossY->copy();
	m_interpolateAcrossYTemp	=   rhs.m_interpolateAcrossYTemp->copy();
	
	
	m_interpolateAcrossX.resize(m_interpolateAcrossX.getsize());
	
	for ( int i = 0 ; i < m_interpolateAcrossX.getsize(); i++ ){
		m_interpolateAcrossX[i]		=	rhs.m_interpolateAcrossX[i]->copy();
	}
	
}	
	
	
/***************************************************************
**	Class   : MlEq2DInterpolator  
**	Function: MlEq2DInterpolator
**	Returns : 
**	Comment : 
****************************************************************/
	
	
MlEqConstInterpolatorHandle MlEq2DInterpolator::getXInterpolator(int whichslice)
{
	MlEqInterpolatorHandle pInterp;
	pInterp = getObjectFromVector(m_interpolateAcrossX,whichslice);
	return pInterp;
}	


/***************************************************************
**	Class   : MlEq2DInterpolator  
**	Function: MlEq2DInterpolator
**	Returns : 
**	Comment : 
****************************************************************/

void MlEq2DInterpolator::reinitialize(const CVector& xData,const CVector& yData,int whichData)
{

	MlEqInterpolatorHandle pInterp;
	pInterp = getObjectFromVector(m_interpolateAcrossX,whichData);



	if ( pInterp->getNumberOfSlices() > 1 ){
		pInterp->reinitialize(xData,yData,whichData);
	}
	else{
		pInterp->reinitialize(xData,yData,0);
	}


}


void MlEqAnalyticCurveWithTanhWingEdge::initialize(MlEqInterpolatorHandle& interp)
{

	int n = interp->getNumberOfSlices();

	if ( n > 1 ){
		throw("interpolator should contain one data set only");
	}

	m_interp	= interp;

	const CVector& xData = interp->getXData();

	m_lowerEdge = xData[0];
	m_upperEdge = xData[xData.getsize()-1];

}


void MlEqAnalyticCurveWithTanhWingEdge::initialize(double lowerEdge,double upperEdge,int addTanhWing,double cL,double cR,double yPower,int npoints,CVector & coeff)
{

	m_coeff			= coeff;
	m_npoints		= npoints;

	if ( npoints < 2 ){
		throw("increase number of points in MlEqAnalyticCurveWithTanhWingEdge");
	}

	m_lowerEdge		= lowerEdge;
	m_upperEdge		= upperEdge;

	if ( addTanhWing )
		m_addTanhWing	= true;
	else
		m_addTanhWing = false;


	m_cL		= cL;
	m_cR		= cR;
	m_finalTol	= 1e-10;
	m_stopTol	= 1e-10;
	m_yPower    = yPower;

	double x		= lowerEdge;
	double delta	= (upperEdge-lowerEdge)/(double)(npoints-1);
	double y;
	
	if ( m_addTanhWing )
	{
		GVector<CVector> xData(1);
		GVector<CVector> yData(1);
		xData[0].resize(npoints);
		yData[0].resize(npoints);

		for ( int i = 0 ; i < npoints; i++ )
		{
			xData[0][i]	= x;
			y			= getAnalyticValueInside(x);
			yData[0][i]	= y;
			x			+= delta;
		}

		m_interp = new MlEqMonotonicSplineInterpolator(
						xData,yData,addTanhWing,cL,cR,yPower);
	}

}


double MlEqAnalyticCurveWithTanhWingEdge::getValue(double xVal,int whichData)
{

	if ( whichData != 0 ){
		throw("only one data set is implemented");
	}

	double res;

	if ( (xVal < m_upperEdge && xVal > m_lowerEdge) || !m_addTanhWing){
		res = getAnalyticValueInside(xVal);
		return res;
	}

	res = m_interp->getValue(xVal);
	return res;
}


void MlEqAnalyticCurveWithTanhWingEdge::getValues(CVector& yVals,double xVal)
{
	for ( int i = 0 ; i < yVals.getsize(); i++ )
	{
		yVals[i] = getValue(xVal);
	}
}


double MlEqAnalyticCurveWithTanhWingEdge::getAnalyticValueInside(double xVal)
{
	int n = m_coeff.getsize();

	if ( !n ){
		return 1e99;
	}

	double res = 0.0;
	for ( int i = 0 ; i < n; i++ ){
		res += m_coeff[i]*pow(xVal,(double)i);
	}

	return res;
}


void MlEqAnalyticCurveWithTanhWingEdge::FittoData(const CVector& xData,const CVector& yData,CVector& initialGuess,CMatrix& xValBounds,double finalTol=1e-12,double stopTol=1e-12)
{
	CMatrix fitData(xData.getsize(),2);

	for ( int i = 0 ; i < xData.getsize(); i++ )
	{
		fitData[i][0] = xData[i];
		fitData[i][1] = yData[i];
	}

	FittoData(fitData,initialGuess,xValBounds,finalTol,stopTol);
}

void MlEqAnalyticCurveWithTanhWingEdge::FittoData( CMatrix& fitData,CVector& initialGuess,CMatrix& xValBounds,double finalTol=1e-12,double stopTol=1e-12)
{

	fitPolynomial fittingCurve;
	
	CMatrix NonZeroPartialDerivativesSpecification;
	double FinalTolerance		= finalTol;
	double StoppingTolerance	= 1e-12;
	double InitialTolerance		= 2*FinalTolerance;
	int outputFlag				= 0;

	m_finalTol		= FinalTolerance;
	m_stopTol		= StoppingTolerance;
	m_initialGuess	= initialGuess;
	m_xValBounds	= xValBounds;

	CMatrix ObjectiveBounds;
	fittingCurve.initialize(initialGuess,fitData,
					    xValBounds,ObjectiveBounds,
					    InitialTolerance,FinalTolerance,StoppingTolerance,
					    NonZeroPartialDerivativesSpecification,outputFlag);
	
	int maximizeFlag = 0;
	int returnCode;
	
	m_coeff.resize(3);
	fittingCurve.solve( m_coeff, maximizeFlag, returnCode);

	m_initialGuess = m_coeff;
}

void MlEqAnalyticCurveWithTanhWingEdge::reinitialize(const CVector& coeff)
{
	CVector c;
	c = coeff;

	initialize(m_lowerEdge,m_upperEdge,m_addTanhWing,m_cL,m_cR,m_yPower,m_npoints,c);
//	initialize(m_xfitData,vols,m_initialGuess,m_xValBounds,m_addTanhWing,m_cL,m_cR,m_yPower,m_finalTol,m_stopTol);
}


MlEqAnalyticCurveWithTanhWingEdge::MlEqAnalyticCurveWithTanhWingEdge(double lowerEdge,double upperEdge,int addTanhWing,double cL,double cR,double yPower,int npoints,CVector& coeff)
{
	initialize(lowerEdge,upperEdge,addTanhWing,cL,cR,yPower,npoints,coeff);
}

void MlEqAnalyticCurveWithTanhWingEdge::initialize(CMatrix& fitData,CVector& initialGuess,CMatrix& xValBounds,double lowerEdge,double upperEdge,int addTanhWing,double cL,double cR,double yPower,double finalTol,double stopTol)
{

	CVector xData(fitData.rows());
	CVector yData(fitData.rows());

	for ( int i = 0 ; i < fitData.rows(); i++ )
	{
		xData[i] = fitData[i][0];
		yData[i] = fitData[i][1];		
	}

	m_xfitData = xData;

	initialize(xData,yData,initialGuess,xValBounds,addTanhWing,cL,cR,yPower,finalTol,stopTol);
}

MlEqAnalyticCurveWithTanhWingEdge::MlEqAnalyticCurveWithTanhWingEdge(CVector& xfitData,CVector& yfitData,CVector& initialGuess,CMatrix& xValBounds,int addTanhWing,double cL,double cR,double yPower,double finalTol,double stopTol)
{
	initialize(xfitData,yfitData,initialGuess,xValBounds,addTanhWing,cL,cR,yPower,finalTol,stopTol);
}



void MlEqAnalyticCurveWithTanhWingEdge::initialize(const CVector& xfitData,const CVector& yfitData,CVector& initialGuess,CMatrix& xValBounds,int addTanhWing,double cL,double cR,double yPower,double finalTol,double stopTol)
{

	MlEqInterpolator::initialize(xfitData,yfitData);

	int npoints;
	int ndim = 10;

	if ( xfitData.getsize() < ndim ){
		npoints = xfitData.getsize();
	}
	else{
		npoints = ndim;
	}

	m_npoints	= npoints;
	m_xfitData	= xfitData;

	FittoData(xfitData,yfitData,initialGuess,xValBounds,finalTol,stopTol);

	m_lowerEdge = xfitData[0];
	m_upperEdge	= xfitData[xfitData.getsize()-1];

	m_initialGuess = initialGuess;
	initialize(m_lowerEdge,m_upperEdge,addTanhWing,cL,cR,yPower,npoints,m_coeff);
}


void MlEqAnalyticCurveWithTanhWingEdge::setCoeff(CVector& newCoeff)
{
	if ( m_coeff.getsize() != newCoeff.getsize() ){
		throw("incorrect dimesion encountered in resetting cooeff");
	}

	m_coeff = newCoeff;
}

void MlEqAnalyticCurveWithTanhWingEdge::initialize(const CVector& xfitData,CVector& initialGuess,CMatrix& xValBounds,int addTanhWing,double cL,double cR,double yPower,double finalTol,double stopTol)
{

	MlEqInterpolator::initialize(xfitData,xfitData);

	int npoints;
	int ndim = 10;

	if ( xfitData.getsize() < ndim ){
		npoints = xfitData.getsize();
	}
	else{
		npoints = ndim;
	}

	m_npoints	= npoints;
	m_xfitData	= xfitData;


	m_lowerEdge = xfitData[0];
	m_upperEdge	= xfitData[xfitData.getsize()-1];

	m_initialGuess = initialGuess;
	initialize(m_lowerEdge,m_upperEdge,addTanhWing,cL,cR,yPower,npoints,m_coeff);
}

void MlEqAnalyticCurveWithTanhWingEdge::reinitialize(const CVector& xData,const CVector& yData,int whichData)
{
	if ( xData.getsize() != yData.getsize() ){
		throw("error in setting up MlEqAnalyticCurveWithTanhWingEdge");
	}

	initialize(xData,yData,m_initialGuess,m_xValBounds,m_addTanhWing,m_cL,m_cR,m_yPower,m_finalTol,m_stopTol);
}






