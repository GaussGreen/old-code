#ifndef __EDCubicSpline_H__
#define __EDCubicSpline_H__

#include "cMatrix.h"

/*************************************************************************************/


typedef	 int    status;

/*************************************************************************************


  If we write for the tanh wing  y(x) = A + B * tanh( C*(x-x0) )

  we can match an inflection point at y0 = y(x0) with A = y0 and B*C determined by
  the slope. In the setAlgPars fct you can set a value for C in the left and right
  wing: cL, cR. The bL, bR are then calculated internally using the (currently
  numerically calculated...) slope on the left and right boundary.


*************************************************************************************/


class   EDCubicSpline
{

public:
	//////////   constructor:
	EDCubicSpline(int 	Prt		= 0, 
				  int 	PrtErr	= 0 )
				: Prt_( Prt ), PrtErr_( PrtErr ), pars_(0), 
				  xData_(0), yData_(0), Npars_(0), Ndata_(0), yPower_(1) {}


	EDCubicSpline( 
					const CVector& 	xData,			// Input Data, on the x dimension 
					const CVector& 	yData,			// Input Data, on the y dimension
					int					addTanhWings,	// addTanhWings
					double				cL,				// cL
					double				cR,				// cR
					int					yPower			// yPower
					);


	void initialize( 
					const CVector& 	xData,			// Input Data, on the x dimension 
					const CVector& 	yData,			// Input Data, on the y dimension
					int					addTanhWings,	// addTanhWings
					double				cL,				// cL
					double				cR,				// cR
					int					yPower			// yPower
					);


	double getValue(const double xVal);

  //////////   destructor:
  virtual ~EDCubicSpline() {}
  
  // 	make private, do not allow copying for Calculator:
  EDCubicSpline	( const EDCubicSpline& cs );
  EDCubicSpline& operator = 	( const EDCubicSpline& cs );

  ////////// 	Set Inputs:
  status	setInputData( const CVector& 	xData, 
						  const CVector& 	yData,
						  int				yPower = 1);

  status	setAlgPars	( int addTanhWings = 0, 
						  double cL	= 1,
						  double cR	= 1 );

  /*****************************************************************
  status	setAlgPars	( int 	 addTanhWings 	= 0, 
				  double xL	      	= 0,
				  double yL		= 0,
				  double xR		= 0,
				  double yR		= 0 );
  ******************************************************************/


  status	buildCurve	(	const double& yD1 	= 1.e+30, 
							const double& yDN 	= 1.e+30 );

  status 	setPrt  	( int Prt ) { Prt_ = Prt; return 1; }  // 1,2,3,4,5,6
  status 	setPrtError	( int Prt ) { PrtErr_ = Prt; return 1; } 

  //////////	Getting of results/info:

  double	operator()	( const double& x )	;
  //  int		getValues	( CVector& yy, const CVector& xx )	const;

  //  kill this: only allow operator() access to results !??
  const	CVector& pars	() const { return pars_; }		// what if !hasCurve_ ??

  int		Npars		() const { return Npars_; }
  int		Ndata		() const { return Ndata_; }
  int		yPower		() const { return yPower_; }
  double    cL			() const  { return cL_;};
  double    cR			() const { return cR_;};
  int       addTanhWing () const { return addTanhWings_;};
  
  const CVector&		getXdata()const{return xData_;};
  const CVector&		getYdata()const{return yData_;};


  protected:


  int 		spline( double* x, double* y, int n, double yp1, double ypn, double* y2 );
  int		splint( const double* xa, const double* ya, const double* y2a, int n, double x, double *y ) ;

protected:

  CVector	pars_;		// or should we call this yDDs or so, since it is 2nd derivatives... !?
  CVector	xData_;
  CVector	yData_;	

  status	hasData_;
  status	hasCurve_;

  int		Npars_;
  int		Ndata_;
  int		Prt_;		// whether/how much to print to stdout
  int		PrtErr_;	// whether to print error messages

  double	yD1_, yDN_;	// derivatives at boundaries

  int		addTanhWings_;
  double	bL_, cL_;	// b,c are coeffts of tanh in the wings...
  double	bR_, cR_;


  double	xL_, yL_;	// xL, yL, xR, yR are points needed to determine tanh coeffts...
  double	xR_, yR_;

  int		yPower_;	// before doing anything, make yData_ be yData_^yPower_



private: 

  double     	solveForC( double x0, double y0, double x1, double y1, double slope );


};

//////////////////////////////////////////////////////////////////////////////////


 EDCubicSpline* MlEqCubicSplineCreate( 
	const CVector& 	xData,			// Input Data, on the x dimension 
	const CVector& 	yData,			// Input Data, on the y dimension
	int					addTanhWings,	// addTanhWings
	double				cL,				// cL
	double				cR,				// cR
	int					yPower			// yPower
);


 double MlEqCubicSplineGetValue( 
	EDCubicSpline& cSpline,		// the spline 
	double				x				// the x value of which you want F(x)
);

 void MlEqCubicSplineInitialize( 
	EDCubicSpline& cubicSp,				// cubic spline object to be initialized
	const CVector& 	xData,			// Input Data, on the x dimension 
	const CVector& 	yData,			// Input Data, on the y dimension
	int					addTanhWings,	// addTanhWings
	double				cL,				// cL
	double				cR,				// cR
	int					yPower			// yPower
);






class   MlEqMonoSpline
{

public:

  MlEqMonoSpline	( )
    :	Ndata_(0), addTanhWings_(1), cL_(1), cR_(1), Prt_(0), PrtErr_(0),yPower_(1) {};


  MlEqMonoSpline( const CVector&  	xData,
			  const CVector&  	yData,			  
			  int 		addTanhWings,
			  const double&		cL=1,	// how fast tanh-wing asymptotes on left
			  const double&		cR=1, 	// how fast tanh-wing asymptotes on right
			  const double&     yPower=1

			  );


  void MlEqMonoSpline::initialize	( const CVector&  	xData,
				  const CVector&  	yData,			  
				  int	 		addTanhWings,
				  const double&		cL,
				  const double&		cR,
				  const double&     yPower);

  virtual ~MlEqMonoSpline() {};

  MlEqMonoSpline	( const MlEqMonoSpline& ms );

  MlEqMonoSpline&	operator = ( const MlEqMonoSpline& ms );


  MlEqMonoSpline*	copy() const;



  short			build 	( const CVector&  	xData,
				  const CVector&  	yData,			  
				  short 		addTanhWings,
				  const double&		cL=1,	// how fast tanh-wing asymptotes on left
				  const double&		cR=1, 	// how fast tanh-wing asymptotes on right
				  const double&     yPower=1
				  );

  double		operator()	( const double& x )	const;
  double		operator()	( const double& x,int index )	const;


  double		getValue(double xVal) ;

  short 		setPrt  	( int Prt ) { Prt_ = Prt; return 1; } 
  short 		setPrtError	( int Prt ) { PrtErr_ = Prt; return 1; } 

  // useful, at least for debugging:
  double		bL() const { return bL_; }
  double		bR() const { return bR_; }
  double		cL() const { return cL_; }
  double		cR() const { return cR_; }
  int			addTanhWing () const { return addTanhWings_;};
  int			yPower		() const { return yPower_; }
 

  const CVector&		getXdata()const{return xData_;};
  const CVector&		getYdata()const{return yData_;};

private:

  int			Ndata_;
  CVector		xData_;
  CVector		yData_;
  CVector		y1_, y2_, y3_;		// 1st, 2nd, 3rd derivatives @ xData

  short			addTanhWings_;
  double		cL_, cR_;
  double		bL_, bR_;
  double		yPower_;

  short			Prt_, PrtErr_;

  MlEqMonoSpline&	copyDeep( const MlEqMonoSpline& ms);

  // all fcts below could be static:

  short			increasing	( const CVector& x ) const;

  double 		calcCubicSplineHorner	( double x, 
						  const CVector& xData, 
						  const CVector& yData, 
						  const CVector& y1, 
						  const CVector& y2, 
						  const CVector& y3 ) const;

  double 		calcCubicSplineHorner	( double x, 
						  const CVector& xData, 
						  const CVector& yData, 
						  const CVector& y1, 
						  const CVector& y2, 
						  const CVector& y3 ,int index) const;

  void 			calcCubicSplineCoeffs	( const CVector& 	x, 
						  const CVector& 	y, 
						  const CVector& 	y1, 
						  CVector& 	y2,	// y2,y3=output
						  CVector& 	y3 ) const;

  double 		deriv1 	( const CVector& 	x, 	
				  const CVector& 	y, 
				  int i, 
				  int sgn ) const;

  double 		deriv2 	( const CVector& 	x, 
				  const CVector& 	y, 
				  int i ) const;

  void 			adjustDerivs	( CVector& 		d, 	// output
					  const CVector& 	x, 	
					  const CVector& 	y) const;

  int 			locate	(double x, const CVector& xData) const;

};



#endif

