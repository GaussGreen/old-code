/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		BRENTSOLVER.H
	PROJECT:	UTIL
	
	DESCRIPTION:	this class provides methods to use a BRENT Solver
					implemented following the numerical recepies.

   -----------------------------------------------------------------

	CREATION:	October 8, 2004

	LAST MODIF:	October 8, 2004
   
	-----------------------------------------------------------------
   
	ICM Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

# ifndef _BRENTSOLVER_H
# define _BRENTSOLVER_H

# include <math.h> 

# include "ICMKernel\glob\icm_constants.h"
# include "ICMKernel\glob\icm_types.h"

class BrentSolver
{
	// ***********************************************************
	// THE DATA
	// ***********************************************************

protected:

	int		NbMaxIter;	// Max iteration Number
	double	Target;	 	// Goal in One Dimension

	// ***********************************************************
	// THE METHODS
	// ***********************************************************

public :

	// Constructors
	BrentSolver()	{Init();}
	BrentSolver(const BrentSolver* Data) {Init();Copy(Data);}
	
	// Destructor
	~BrentSolver()	{Reset();}

	// ------------------------------------------------------------------------------
	// Constructors & alike
	// ------------------------------------------------------------------------------
public:

	void	BitwiseCopy(const BrentSolver* src);
	void	Copy(const BrentSolver* src);
	
	BrentSolver* Clone() const;

	void	Init();
	void	Reset(); 

	// ------------------------------------------------------------------------------
	// METHODS
	// ------------------------------------------------------------------------------

public:
	
	void	SetNbMaxIter(int Data) {NbMaxIter = Data;}
	void	GetNbMaxIter(int& Data) const {Data = NbMaxIter;}

	void	SetTarget(double Data) {Target = Data;}
	void	GetTarget(double& Data) const {Data = Target;}

	ReturnCode	Solve(	double x1,
					double x2,
					double tol,
					int&	nb_iter,
					ReturnCode (*fct)(void*, double, double&),
					void* param,
					double& root,
					bool	NeedZeroBracket = false,
					int		ZeroBracketChoice = 0,
					int		NbIntervals	=	10) const;	// 0 - outward // else inward -> number of roots

	ReturnCode	Solve(	double x1,
					double x2,
					double fx1,
					double fx2,
					double tol,
					int&	nb_iter,
					ReturnCode (*fct)(void*, double, double&),
					void* param,
					double& root) const;

protected:

	double ZBrent(double Src1,double Dest1,double Src2,double Dest2,double ,double fc) const;

	ReturnCode ZeroBracket(	double x1,
							double x2,
							double&	fx1,
							double&	fx2,
							unsigned int NbIterMax,
							double factor,
							ReturnCode (*fct)(void*, double, double&),
							void* param,
							bool& Result) const;

	ReturnCode ZeroBracket_In(	double x1,
										double x2,
										ReturnCode (*fct)(void*, double, double&),
										void* param,
										int	n,
										double	xb1[],
										double	xb2[],
										int*		nb,
										double	fx1[],
										double	fx2[]
										) const;

};

# endif

