#include<iostream>
#include "stdafx.h"
#include "edcubicspline.h"
#include "nrutil.h"

/////////////////////////////////////////////////////////////////////////////////////////


#undef  SQR
#define SQR(x)         ( (x)*(x) )
#undef  MAX
#define MAX(a, b)      ( (a) > (b) ? a : b )
#undef  MIN
#define MIN(a, b)      ( (a) < (b) ? a : b )


#undef   MYPRINT
#define  MYPRINT(x) 	std::cout << " " << #x << "		= " << x << endl << std::cout.flush()

/////////////////////////////////////////////////////////////////////////////////////////


/****************************************************************
**	Class  : EDCubicSpline
**	Routine: EDCubicSpline
**	Returns: nothing
**	Action : constructor
****************************************************************/

EDCubicSpline::EDCubicSpline( const EDCubicSpline& rhs )
:
pars_		( rhs.pars_ ),
xData_		( rhs.xData_ ),
yData_		( rhs.yData_ ),
hasData_	( rhs.hasData_ ),
hasCurve_	( rhs.hasCurve_ ),
Npars_		( rhs.Npars_ ),
Ndata_		( rhs.Ndata_ ),
Prt_		( rhs.Prt_ ),
PrtErr_		( rhs.PrtErr_ ),
yD1_		( rhs.yD1_ ),
yDN_		( rhs.yDN_ ),
addTanhWings_( rhs.addTanhWings_ ),
bL_			( rhs.bL_ ),
cL_			( rhs.cL_ ),	
bR_			( rhs.bR_ ),
cR_			( rhs.cR_ ),
xL_			( rhs.xL_ ),
yL_			( rhs.yL_ ),
xR_			( rhs.xR_ ),
yR_			( rhs.yR_ ),
yPower_		( rhs.yPower_)
{
}


EDCubicSpline& EDCubicSpline::operator = ( const EDCubicSpline& rhs )
{

	if( this == &rhs )
		return * this;

	pars_		= rhs.pars_ ;
	xData_		= rhs.xData_ ;
	yData_		= rhs.yData_ ;
	hasData_	= rhs.hasData_ ;
	hasCurve_	= rhs.hasCurve_ ;
	Npars_		= rhs.Npars_ ;
	Ndata_		= rhs.Ndata_ ;
	Prt_		= rhs.Prt_ ;
	PrtErr_		= rhs.PrtErr_ ;
	yD1_		= rhs.yD1_ ;
	yDN_		= rhs.yDN_ ;
	addTanhWings_= rhs.addTanhWings_ ;
	bL_			= rhs.bL_ ;
	cL_			= rhs.cL_ ;	
	bR_			= rhs.bR_ ;
	cR_			= rhs.cR_ ;
	xL_			= rhs.xL_ ;
	yL_			= rhs.yL_ ;
	xR_			= rhs.xR_ ;
	yR_			= rhs.yR_ ;
	yPower_		= rhs.yPower_;
	return *this;
}


////////// 	Set Input pars etc:

status	EDCubicSpline::setInputData( const CVector& 	xData,
									 const CVector& 	yData,
									 int				yPower)
{
  int i, ok = 1;

  yPower_	= yPower;

  hasData_	= 0;
  Ndata_ 	= xData.getsize();

  if( Ndata_ != yData.getsize() )
  {
	  if( PrtErr_ ) std::cout << " EDCubicSpline::setInputData: data size mismatch!"<< std::endl;
    return 0;
  }

  if( Ndata_ < 1 )
  {
    if( PrtErr_ ) std::cout << " EDCubicSpline::setInputData: need at least one data point!"<< std::endl;
    return 0;
  }

  int incr = xData[Ndata_-1] - xData[0] > 0 ? 1 : -1;

  // TK: because of tanhWings I will for now insist on increasing X after all:

  if( incr == -1 ) 
  {
    if( PrtErr_ ) std::cout << " EDCubicSpline::setInputData: xData have to be strictly increasing!"<<endl;
    return 0;      

  }


  for ( i = 0; i < xData.getsize(); i++ )
  if( i > 0 )
  {
    if( incr * ( xData[i] - xData[i-1] ) <= 0.0 )
    {
      if( PrtErr_ ) std::cout << " EDCubicSpline::setInputData: xData have to be strictly monotonic!"<<endl;
      return 0;      
    }
  }


  xData_ 	= xData;
  yData_ 	= yData;

  // change yData to yData to the power of yPower
  for(int yIndex = 0; yIndex < yData_.getsize(); yIndex++ )
	  yData_[ yIndex ] = pow( yData_[ yIndex ], (double) yPower_ );

  hasData_	= 1;

  Npars_	= Ndata_;
  pars_.resize( Npars_ );

  setAlgPars();

  return ok;
}


int	EDCubicSpline::buildCurve( const double& yD1, const double& yDN )
{
  int ok = 1;

  hasCurve_ = 0;

  if( !hasData_ )
  { 
    if( PrtErr_ ) std::cout << " EDCubicSpline::buildCurve: set proper input data first!"<<endl;

    throw ( "EDCubicSpline::buildCurve: set proper input data first!",
			  "EDCubicSpline::buildCurve()" );

    return 0;
  }

  yD1_ = yD1;
  yDN_ = yDN;

  if( addTanhWings_ ) // use natural BC for TanhWings: !!
  {
    yD1_ = 1.e30;
    yDN_ = 1.e30;
  }

  // write results for 2nd derivs into pars_:

  ok = spline( xData_.getPtr()-1, yData_.getPtr()-1, Ndata_, yD1_, yDN_, pars_.getPtr()-1 );

  if( !ok ) return 0;


  if( addTanhWings_ )
  {

 


    // Need derivatives at x1, xN to proceed: have to hack it here: clean up a bit later
    double x, y, slopeL, slopeR;

    // Is the sufficiently safe choice (or should we take curvature into account?):
    double eps = 0.001 * MAX( fabs(xData_[0]), fabs(xData_[Ndata_-1]) )/Ndata_;		// 0.0001;

    x = xData_[0] + eps;
    splint( xData_.getPtr()-1, yData_.getPtr()-1, pars_.getPtr()-1, Ndata_, x, &y );
    slopeL = ( y - yData_[0] ) / eps;

    if( Prt_ > 3 )
    {
      MYPRINT( xData_[0] );
      MYPRINT( yData_[0] );
      MYPRINT( x );
      MYPRINT( y );
      MYPRINT( slopeL );
      MYPRINT( xL_ );
      MYPRINT( yL_ );
      std::cout << endl;
    }

    x = xData_[Ndata_-1] - eps;
    splint( xData_.getPtr()-1, yData_.getPtr()-1, pars_.getPtr()-1, Ndata_, x, &y );
    slopeR = ( yData_[Ndata_-1] - y ) / eps;

    if( Prt_ > 3 )
    {
      MYPRINT( xData_[Ndata_-1] );
      MYPRINT( yData_[Ndata_-1] );
      MYPRINT( x );
      MYPRINT( y );
      MYPRINT( slopeR );
      MYPRINT( xR_ );
      MYPRINT( yR_ );
      std::cout << endl;
    }



//    if( slopeL * (yData_[0] - yL_) <= 0 )
//    {
//      if( PrtErr_ ) std::cout << " EDCubicSpline::buildCurve: bad choice of yL !"<<endl;
//      return 0;      
//    }

//    if( slopeR * (yR_ - yData_[Ndata_-1]) <= 0 )
//    {
//      if( PrtErr_ ) std::cout << " EDCubicSpline::buildCurve: bad choice of yR !"<<endl;
//      return 0;      
//    }

    // TK: what if slopes are = 0 !!!???

//    cL_ = solveForC( xData_[0], yData_[0], xL_, yL_, slopeL );
//    bL_ = slopeL / cL_;

//    cR_ = solveForC( xData_[Ndata_-1], yData_[Ndata_-1], xR_, yR_, slopeR );
//    bR_ = slopeR / cR_;
 
    // Now that cL, cR are given, replace above by:

    if( cL_ == 0.0 ) 	bL_ = 0.0;
    else		bL_ = slopeL / cL_;

    if( cR_ == 0.0 ) 	bR_ = 0.0;
    else		bR_ = slopeR / cR_;

    ok = 1;

    if( Prt_ > 3 )
    {
      MYPRINT( bL_ );
      MYPRINT( cL_ );
      MYPRINT( bR_ );
      MYPRINT( cR_ );
      std::cout << endl;
    }

  } // end of TanhWings


  if( ok ) 	hasCurve_ = 1;

  return ok;
}


double	EDCubicSpline::operator()( const double& x )	
{
  int 	 ok;
  double y;

  if( !hasCurve_ ) return -999.0;

  if( !addTanhWings_ )
  {
    ok = splint( xData_.getConstPtr()-1, yData_.getConstPtr()-1, pars_.getConstPtr()-1, Ndata_, x, &y );
    // can anything really go wrong here !!?? Not with strictly monotonic xData, I think....
    
    if( !ok ) return -999.0;
    
    return y;
  }

  else	// addTanhWings:
  {

    if( x < xData_[0] ) // remember we assume x are increasing now...
    {
      return yData_[0] + bL_ * tanh( cL_*(x-xData_[0]) );
    }

    if( x > xData_[Ndata_-1] )
    {
      return yData_[Ndata_-1] + bR_ * tanh( cR_*(x-xData_[Ndata_-1]) );
    }

    else
    {
      ok = splint( xData_.getConstPtr()-1, yData_.getConstPtr()-1, pars_.getConstPtr()-1, Ndata_, x, &y );
      //  can anything really go wrong here !!?? Not with strictly monotonic xData, I think....
    
      if( !ok ) return -999.0;
    
      return y;    
    }
  }
}





status	EDCubicSpline::setAlgPars( int 	  addTanhWings	,
				   double cL		,
				   double cR		)
{
  int ok = 1;

  // checks...

  addTanhWings_	= addTanhWings;
  cL_		= cL;
  cR_		= cR;

  hasCurve_	= 0;

  return ok;
}




/////////////////////////////////////////////////////////////////////////////////////////////

// Call this only if  abs( slope ) > abs( (y1-y0)/(x1-x0) ):

double	 EDCubicSpline::solveForC( double x0, double y0, double x1, double y1, double slope ) 
{
  // check above condition....

  int i=0, maxIts = 1000;

  double c = 1.0;	// TK: always iterate for positive soln ( -c is also a soln...)
  double cNew;

  while( i < maxIts )
  {
    cNew = slope * tanh( c*(x1-x0) ) / ( y1 - y0 );

    if( fabs( cNew/c - 1.0 ) < 1.e-8 ) break;

    c = cNew;
    i++;
  }

  return cNew;
}


////////////////////////////////////////////////////////////////////////////////////////////////
/////////   Numerical Recipes code:

int EDCubicSpline::spline( double* x, double* y, int n, double yp1, double ypn, double* y2 )


{
  int ok = 1;

  int i,k;
  double p,qn,sig,un,*u;

  u = dvector(1,n-1);

  if (yp1 > 0.99e+30)
    y2[1]=u[1]=0.0;
  else {
    y2[1] = -0.5;
    u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }

  for (i=2;i<=n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }

  if (ypn > 0.99e+30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }

  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];

  free_dvector(u,1,n-1);
	
  return ok;
}

int	EDCubicSpline::splint( const double* xa, const double* ya, const double* y2a, int n, double x, double *y ) 
{
  int ok = 1;

  //        void nrerror(char error_text[]);
  int klo,khi,k;
  double h,b,a;

  klo = 1;
  khi = n;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }

  h = xa[khi] - xa[klo];
  // TK: this should never happen anymore:
  if (h == 0.0) 
  { 
	  throw("Bad xa input to routine splint"); 
	  return ok;
  } // nrerror always rets 0 !?

  a = (xa[khi]-x)/h;
  b = (x-xa[klo])/h;
  *y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;

  return ok;
}




void EDCubicSpline::initialize( 
	const CVector& 	xData,			// Input Data, on the x dimension 
	const CVector& 	yData,			// Input Data, on the y dimension
	int					addTanhWings,	// addTanhWings
	double				cL,				// cL
	double				cR,				// cR
	int					yPower			// yPower
)
{
	setInputData( xData, yData, yPower );

	setAlgPars(addTanhWings, cL, cR );
	
	buildCurve( );
}


EDCubicSpline:: EDCubicSpline( 
	const CVector& 	xData,			// Input Data, on the x dimension 
	const CVector& 	yData,			// Input Data, on the y dimension
	int					addTanhWings,	// addTanhWings
	double				cL,				// cL
	double				cR,				// cR
	int					yPower			// yPower
)
: 
Prt_( 0 ), 
PrtErr_( 0 )
{
	initialize( xData,yData,addTanhWings,cL,cR,yPower);			
}


double EDCubicSpline::getValue(const double xVal)
{
	double value = (*this)(xVal);
	
	if ( yPower_ == 1 )
	{
		return value;
	}

	if(value < 0 && (yPower_%2 == yPower_/2) ){
		throw(" value should be positive in cubic spline get value");
	}

	value = pow(value,1.0/(double)yPower_);
	return value;

}

/****************************************************************
**	Name			MlEqCubicSplineCreate
**	Summary:		create a cubic spline interpolator
**	References:		
**	Description:	
**			
****************************************************************/


 EDCubicSpline* MlEqCubicSplineCreate( 
	const CVector& 	xData,			// Input Data, on the x dimension 
	const CVector& 	yData,			// Input Data, on the y dimension
	int					addTanhWings,	// addTanhWings
	double				cL,				// cL
	double				cR,				// cR
	int					yPower			// yPower
)
{
	
	EDCubicSpline* retVal = new EDCubicSpline();

	MlEqCubicSplineInitialize(*retVal,xData,yData,addTanhWings,cL,cR,yPower);

	return retVal;
}


/****************************************************************
**	Adlib Name:		MlEqCubicSplineCreate
**	Summary:		initializes a cubic spline interpolator
**	References:		
**	Description:	
**	Example:		
****************************************************************/


 void MlEqCubicSplineInitialize( 

	EDCubicSpline& cubicSp,				// cubic spline object to be initialized
	const CVector& 	xData,			// Input Data, on the x dimension 
	const CVector& 	yData,			// Input Data, on the y dimension
	int					addTanhWings,	// addTanhWings
	double				cL,				// cL
	double				cR,				// cR
	int					yPower			// yPower
)
{



	cubicSp.setInputData( xData, yData, yPower );

	cubicSp.setAlgPars(addTanhWings, cL, cR );
	
	cubicSp.buildCurve( );


}


/****************************************************************
**					MlEqCubicSplineGetValue
**	Summary:		get F(x)
**	References:		
**	Description:	
**	Example:		
****************************************************************/

double MlEqCubicSplineGetValue( 
	EDCubicSpline& cSpline,		// the spline 
	double				x				// the x value of which you want F(x)
)
{
	double retValue = cSpline(x);

	if( cSpline.yPower() % 2 == 0 && retValue < 0.0)
	{
		char buff[1024];
		sprintf( buff, "%.3lf%s%i",retValue, " returned is negative while power passed in is ", cSpline.yPower());
		throw(  buff);
	}
	
	double exponent_ =  1.0 / (( double )cSpline.yPower() );

	retValue = pow( retValue, exponent_ ); 
	return retValue;
}







MlEqMonoSpline::MlEqMonoSpline	( const MlEqMonoSpline& ms )
{
  copyDeep( ms );
}

MlEqMonoSpline&	MlEqMonoSpline::operator=( const MlEqMonoSpline& ms )
{
  return copyDeep( ms );
}


/****************************************************************
**      Class  : MlEqMonoSpline
**      Routine: copy
**      Returns: copy of self
**      Action : 
****************************************************************/
 
MlEqMonoSpline *MlEqMonoSpline::copy() const
{
        return new MlEqMonoSpline( *this );
}



MlEqMonoSpline&	MlEqMonoSpline::copyDeep( const MlEqMonoSpline& ms )
{
  if( this == &ms ) return *this;
 
  Ndata_	= ms.Ndata_;
  xData_	= ms.xData_;
  yData_	= ms.yData_;
  yPower_   = ms.yPower_;
  y1_		= ms.y1_;	
  y2_ 	       	= ms.y2_;
  y3_	        = ms.y3_;
  addTanhWings_	= ms.addTanhWings_;
  cL_		= ms.cL_; 
  cR_		= ms.cR_;
  bL_		= ms.bL_;
  bR_		= ms.bR_;
  Prt_		= ms.Prt_;
  PrtErr_	= ms.PrtErr_;

  return *this;
}

/****************************************************************
**      Class  : MlEqMonoSpline
**      Routine: MlEqMonoSpline
**      Returns: nothing
**      Action : Constructor
****************************************************************/

MlEqMonoSpline::MlEqMonoSpline	( const CVector&  	xData,
				  const CVector&  	yData,			  
				  int	 		addTanhWings,
				  const double&		cL,
				  const double&		cR,
				  const double&     yPower)

{
	initialize(xData,yData,addTanhWings,cL,cR,yPower );
}

/****************************************************************
**      Class  : MlEqMonoSpline
**      Routine: initialize
**      Returns: nothing
**      Action : Constructor
****************************************************************/

void MlEqMonoSpline::initialize	( const CVector&  	xData,
				  const CVector&  	yData,			  
				  int	 		addTanhWings,
				  const double&		cL,
				  const double&		cR,
				  const double&     yPower)

{  
	Prt_ = 0;
	PrtErr_ = 0;

  short ok = build( xData, yData, addTanhWings, cL, cR,yPower );

  if( !ok )
  {
    throw( "build failed" );
  }
}




short	MlEqMonoSpline::build 	( const CVector&  	xData,
				  const CVector&  	yData,			  
				  short 		addTanhWings,
				  const double&		cL,
				  const double&		cR ,
				  const double&     yPower)

{
  int i;
  short ok = 1;

  if( xData.getsize() != yData.getsize() ) return 0;

  ok = increasing( xData );

  if( !ok ) 
  {
    if( PrtErr_ ) std::cout << " MlEqMonoSpline::build(): xData not properly ordered!" << endl;

    throw ( "MlEqMonoSpline::build(): xData not properly ordered!",
			  "MlEqMonoSpline::build()" );
    return 0;
  }

  xData_	= xData;
  yData_	= yData;
  yPower_	= yPower;

  // change yData to yData to the power of yPower
  for(int yIndex = 0; yIndex < yData_.getsize(); yIndex++ )
	  yData_[ yIndex ] = pow( yData_[ yIndex ], (double) yPower_ );

  Ndata_	= xData.getsize();

  if( Ndata_ == 0 )	return 0;	// require at least 1 data point:

  y1_.resize(Ndata_);
  y2_.resize(Ndata_);
  y3_.resize(Ndata_);

  cL_		= cL;
  cR_		= cR;
  yPower_	= yPower;

  // NOTE: I checked in some examples that also for Mono splines do we need tanh wings!!

  addTanhWings_	= addTanhWings;

  // Now calculate y1_, y2_, y3_, bL_, bR_:
  // For tanh wings can now use y1_ for first-order derivatives:

  int lo = 0;
  int hi = Ndata_-1;


  if( Ndata_ == 1 )
  {
    y1_[lo]	= 0.0;
    y2_[lo]	= 0.0;
    y3_[lo]	= 0.0;
  }

  else if( Ndata_ == 2 )
  {
    y1_[lo]    	= ( yData_[hi] - yData_[lo] ) / ( xData_[hi] - xData_[lo] );
    y2_[lo]    	= 0.0;
    y3_[lo]	= 0.0;

    y1_[hi]    	= y1_[lo];
    y2_[hi]    	= y2_[lo];
    y3_[hi]	= y3_[lo];

  }

  else		// 3 or more points:
  {
    y1_[lo]	= deriv1( xData_, yData_, lo+1, -1 );
    y1_[hi]	= deriv1( xData_, yData_, hi-1,  1 );

    for(i=lo+1; i<hi; ++i)
    {
      y1_[i]	= deriv2( xData_, yData_, i );
    }

    if( Ndata_ > 3 )	// 4 or more points:
    {
      // adjust derivatives at bdy:
      if( y1_[lo] * y1_[lo+1] < 0.0 )	y1_[lo] = 0.0;
      if( y1_[hi] * y1_[hi-1] < 0.0 )	y1_[hi] = 0.0;
     
      // adjust y1_ according to Fritsch-Carlson:
      adjustDerivs( y1_, xData_, yData_ );
    }

    // calculate remaining y2_, y3_ coefficients:
    calcCubicSplineCoeffs( xData_, yData_, y1_, y2_, y3_ );
  }


  double slopeL	= y1_[lo];
  double slopeR	= y1_[hi];

  if( cL_ == 0.0 ) 	bL_ = 0.0;
  else			bL_ = slopeL / cL_;

  if( cR_ == 0.0 ) 	bR_ = 0.0;
  else			bR_ = slopeR / cR_;


  if( Prt_ > 3 )
  {
      MYPRINT( bL_ );
      MYPRINT( cL_ );
      MYPRINT( bR_ );
      MYPRINT( cR_ );
      std::cout << endl;
  }

  
  return ok;
}




short	MlEqMonoSpline::increasing( const CVector& x ) const
{
  int i;

  for( i=1; i<x.getsize(); ++i)
  {
      if( x[i] <= x[i-1] )
	  {
		  return 0;
	  }
  }

  return 1;
}






double	MlEqMonoSpline::operator()( const double& x ) const
{ 
  if( !addTanhWings_ )	
  {
	  double retval = calcCubicSplineHorner( x, xData_, yData_, y1_, y2_, y3_ );  
    return retval;
  }


  else	// addTanhWings:
  {

    if( x < xData_[0] ) // remember we assume x are increasing now...
    {
      return yData_[0] + bL_ * tanh( cL_*(x-xData_[0]) );
    }

    if( x > xData_[Ndata_-1] )
    {
      return yData_[Ndata_-1] + bR_ * tanh( cR_*(x-xData_[Ndata_-1]) );
    }

    else
    {
		double retval = calcCubicSplineHorner( x, xData_, yData_, y1_, y2_, y3_ ); 
      return retval;
    }
  }
}


double	MlEqMonoSpline::operator()( const double& x ,int index) const
{ 
  if( !addTanhWings_ )	
  {
	  double retval = calcCubicSplineHorner( x, xData_, yData_, y1_, y2_, y3_,index ); 
    return retval;
  }


  else	// addTanhWings:
  {

    if( x < xData_[0] ) // remember we assume x are increasing now...
    {
      return yData_[0] + bL_ * tanh( cL_*(x-xData_[0]) );
    }

    if( x > xData_[Ndata_-1] )
    {
      return yData_[Ndata_-1] + bR_ * tanh( cR_*(x-xData_[Ndata_-1]) );
    }

    else
    {
      return calcCubicSplineHorner( x, xData_, yData_, y1_, y2_, y3_ ,index);
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////
//
// Stuff specific to Monotonic Splines:


//----------------------------------------------------------------------------//
//
// Calculate the spline coefficients y2[i] and y3[i] for a natural cubic
// spline, given the abscissa x[i], the ordinate y[i], and the 1st 
// derivative y1[i].
//
// The spline interpolation can then be evaluated using Horner's rule
//
//      P(u) = y[i] + dx * (y1[i] + dx * (y2[i] + dx * y3[i]))
//
// where  x[i] <= u < x(i+1) and dx = u - x[i].
//
//----------------------------------------------------------------------------//

void MlEqMonoSpline::calcCubicSplineCoeffs ( const CVector& 	x, 
					     const CVector& 	y, 
					     const CVector& 	y1, 
					     CVector& 		y2, 
					     CVector& 		y3) const
{
  int lo = 0;
  int hi = x.getsize()-1;

  for (int i = lo; i < hi; i++) {
    double h  = x[i+1] - x[i],
           mi = (y[i+1] - y[i]) / h;
    y2[i] = (3 * mi - 2 * y1[i] - y1[i+1]) / h;
    y3[i] = (y1[i] + y1[i+1] - 2 * mi) / (h * h);
  }

  y2[hi] = y3[hi] = 0.0;
}



//-----------------------------------------------------------------------------//
// 
// Initial derivatives at boundaries of data set using
// quadratic Newton interpolation
// 
//-----------------------------------------------------------------------------//


double MlEqMonoSpline::deriv1 ( const CVector& x, 
				const CVector& y, 
				int i, int sgn ) 	const 
{
  double di, dis, di2, his;
  int 	 i1,i2;

  i1  = i+1;
  i2  = i-1;
  his = x[i1]-x[i2];
  dis = (y[i1]-y[i2]) / his;
  di  = (y[i1]-y[i]) / (x[i1]-x[i]);
  di2 = (di-dis) / (x[i]-x[i2]);

  return dis + sgn * di2 * his;
}


//-----------------------------------------------------------------------------//
// 
// Initial derivatives within data set using
// quadratic Newton interpolation
// 
//-----------------------------------------------------------------------------//


double MlEqMonoSpline::deriv2 ( const CVector& x, 
				const CVector& y, 
				int i )	const
{
  double di0, di1, di2, hi0;
  int 	 i1, i2;

  i1  = i+1;
  i2  = i-1;
  hi0 = x[i]-x[i2];
  di0 = (y[i]-y[i2]) / hi0;
  di1 = (y[i1]-y[i]) / (x[i1]-x[i]);
  di2 = (di1-di0) / (x[i1]-x[i2]);

  return di0 + di2 * hi0;
}


//-----------------------------------------------------------------------------//
// 
// Fritsch-Carlson iteration to adjust the monotone
// cubic interpolant. The iteration converges with cubic order.
// 
//-----------------------------------------------------------------------------//

void MlEqMonoSpline::adjustDerivs (CVector& d, const CVector& x, const CVector& y) const
{
  int i, i1;

  int stop;

  double d1,r2,t;

  const int max_loop = 20; // should never happen! Note, that currently it
			   // can happen when the curve is not strictly 
  			   // monotone. In future this case should be handled
			   // more gracefully without wasting CPU time.
  int loop = 0;

  do 
  {
    stop = 1;
 
    for(i=0; i<x.getsize()-1; ++i)
    {
      i1 = i + 1;
      d1 = (y[i1]-y[i]) / (x[i1]-x[i]);

      if (d1 == 0.0)
	d[i] = d[i1] = 0.0;

      else 
      {
	t 	= d[i]/d1;
	r2 	= t*t;
	t 	= d[i1]/d1;
	r2     += t*t;

	if (r2 > 9.0) 
	{
	  t 		= 3.0 / sqrt(r2);
	  d[i] 	       *= t;
	  d[i1]        *= t;
	  stop 	       	= 0;
	}
      }
    }

  } while (!stop && ++loop < max_loop);

}


//----------------------------------------------------------------------------//
//
// Return the interpolation value P(u) for a piecewise cubic curve determined
// by the abscissa vector x, the ordinate vector y, the 1st derivative
// vector y1, the 2nd derivative vector y2, and the 3rd derivative vector y3,
// using the Horner scheme. 
// All vectors must have conformant dimenions.
// The abscissa x[i] values must be strictly increasing.
// In the special case of empty data vectors (x,y) a value of 0.0 is returned.
//
// This subroutine evaluates the function
//
//    P(u) = y[i] + dx * (y1[i] + dx * (y2[i] + dx * y3[i]))
//
// where  x[i] <= u < x(i+1) and dx = u - x[i], using Horner's rule
//
//    lo <= i <= hi is the index range of the vectors.
//    if  u <  x[lo] then  i = lo  is used.
//    if  u >= x[hi] then  i = hi  is used.
//
// Input arguments:
//
//    double u       = the abscissa at which the interpolation is to be evaluated
//
//    vector& x
//    vector& y      = the vectors (lo,hi) of data abscissas and ordinates
//
//    vector& y1
//    vector& y2
//    vector& y3     = vectors (lo,hi) of coefficients
//                     y1[i] contains the 1st derivative y'(x[i])
//                     y2[i] contains the 2nd derivative y''(x[i])
//
//  Note: A binary search is performed to determine the proper interval.
//
//----------------------------------------------------------------------------//

double MlEqMonoSpline::calcCubicSplineHorner ( double u, 
					       const CVector& x, 
					       const CVector& y, 
					       const CVector& y1, 
					       const CVector& y2, 
					       const CVector& y3 ) const
{
  // if no data, return 0:

  if( x.getsize() <= 0 ) return 0.0;

  int i = locate( u, x );  

  if( i  < 0 	      ) i = 0;  	// extrapolate to the left
  if( i == x.getsize()-1 ) i--; 	 	// extrapolate to the right

  double dx = u - x[i];
  return ( y[i] + dx * ( y1[i] + dx * ( y2[i] + dx * y3[i] ) ) );
}

double MlEqMonoSpline::calcCubicSplineHorner ( double u, 
					       const CVector& x, 
					       const CVector& y, 
					       const CVector& y1, 
					       const CVector& y2, 
					       const CVector& y3 ,int  i) const
{
  // if no data, return 0:

  if( x.getsize() <= 0 ) return 0.0;

//  int i = locate( u, x );  

  if( i  < 0 	      ) i = 0;  	// extrapolate to the left
  if( i == x.getsize()-1 ) i--; 	 	// extrapolate to the right

  double dx = u - x[i];
  return ( y[i] + dx * ( y1[i] + dx * ( y2[i] + dx * y3[i] ) ) );
}


//----------------------------------------------------------------------------//
//
// Find index of largest element in the increasingly ordered vector x, 
// which is smaller than u. If u is smaller than the smallest value in 
// the vector then the lowest index minus one is returned. 
// A fast binary search is performed.
//
//----------------------------------------------------------------------------//

int MlEqMonoSpline::locate (double u, const CVector& x) const
{
  int i, j;
  int lo = 0;
  int hi = x.getsize()-1;

  if( u < x[lo] ) 
  {
    i = lo - 1;	// attention: return index below smallest index
  } 
  else if( u >= x[hi] ) 
  {
    i = hi;	// attention: return highest index
  } 
  else 
  {
    i = lo;
    j = hi;
    do 
    {
      int k = (i + j) / 2;
      if (u <  x[k])  j = k;
      if (u >= x[k])  i = k;
    } while (j > i + 1);
  }
  return i;
}


double MlEqMonoSpline::getValue(const double xVal)
{
	double value = (*this)(xVal);
	
	if ( yPower_ == 1 )
	{
		return value;
	}

	if( value < 1e-12 )	{
		value = 1e-12;	// added 22/2/5 Romain...local vol may crash otherwise
	}
//	assert(value > 0.0 );

	value = pow(value,1.0/(double)yPower_);
	return value;

}


/////////////////////////////////////////////////////////////////////////////////////////











#undef SQR
#undef MAX
#undef MIN

#undef MYPRINT


