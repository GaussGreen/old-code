/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file basic_distribution.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpbase/numericconstant.h"

#include <algorithm>
#include <cmath>	/// for fabs

#include "gpclosedforms/basic_distributions.h"

using namespace std;
CC_BEGIN_NAMESPACE(ARM)
///////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Argument: an int that specifies the number of arguments of the 
/// associated expression.
///	Returns: nothing
///	Action : Constructor
///////////////////////////////////////////////////////////////////////////////////////

ArgumentList::ArgumentList(int n) : itsContainer(n)
{}


ArgumentList::ArgumentList(int n,int p) : itsContainer(n),itsVectorContainer(p)
{}
////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 1 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,double x0) : itsContainer(1,x0)
,itsVectorContainer(2)
{	itsContainer[0]=x0;
	itsVectorContainer[0]=x;
	itsVectorContainer[1]=y;}


////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 2 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,double x0,double x1):itsContainer(2)
,itsVectorContainer(2)
{
	itsContainer[0]=x0;
	itsContainer[1]=x1;
	itsVectorContainer[0]=x;
	itsVectorContainer[1]=y;
}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 3 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////

ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,double x0,double x1,double x2) : itsContainer(3)
,itsVectorContainer(2)
{
	itsContainer[0]=x0;
	itsContainer[1]=x1;
	itsContainer[2]=x2;
	itsVectorContainer[0]=x;
	itsVectorContainer[1]=y;
		
}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 4 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,double x0,double x1,double x2,double x3) : itsContainer(4)
,itsVectorContainer(2)
{
		itsContainer[0]=x0;
		itsContainer[1]=x1;
		itsContainer[2]=x2;
		itsContainer[3]=x3;
		itsVectorContainer[0]=x;
		itsVectorContainer[1]=y;
		
}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 5 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,double x0,double x1,double x2,double x3,double x4) : itsContainer(5)
,itsVectorContainer(2)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsVectorContainer[0]=x;
			itsVectorContainer[1]=y;
			
	}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 6 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,double x0,double x1,double x2,double x3,double x4,double x5)
 : itsContainer(6),itsVectorContainer(2)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsVectorContainer[0]=x;
			itsVectorContainer[1]=y;
		
	}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 7 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6)
						    : itsContainer(7),itsVectorContainer(2)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsVectorContainer[0]=x;
			itsVectorContainer[1]=y;
			
	}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 1 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,ARM_GP_Vector*  xv,double x0) : itsContainer(1,x0)
,itsVectorContainer(3)
{	itsContainer[0]=x0;
	itsVectorContainer[0]=x;
	itsVectorContainer[1]=y;
	itsVectorContainer[2]=xv;
}


////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 2 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,ARM_GP_Vector*  xv,double x0,double x1):itsContainer(2)
,itsVectorContainer(3)
{
	itsContainer[0]=x0;
	itsContainer[1]=x1;
	itsVectorContainer[0]=x;
	itsVectorContainer[1]=y;
	itsVectorContainer[2]=xv;
}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 3 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////

ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,ARM_GP_Vector*  xv,double x0,double x1,double x2) : itsContainer(3)
,itsVectorContainer(3)
{
	itsContainer[0]=x0;
	itsContainer[1]=x1;
	itsContainer[2]=x2;
	itsVectorContainer[0]=x;
	itsVectorContainer[1]=y;
	itsVectorContainer[2]=xv;
		
}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 4 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,ARM_GP_Vector*  xv,double x0,double x1,double x2,double x3) : itsContainer(4)
,itsVectorContainer(3)
{
		itsContainer[0]=x0;
		itsContainer[1]=x1;
		itsContainer[2]=x2;
		itsContainer[3]=x3;
		itsVectorContainer[0]=x;
		itsVectorContainer[1]=y;
		itsVectorContainer[2]=xv;
		
}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 5 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,ARM_GP_Vector*  xv,double x0,double x1,double x2,double x3,double x4) : itsContainer(5)
,itsVectorContainer(3)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsVectorContainer[0]=x;
			itsVectorContainer[1]=y;
			itsVectorContainer[2]=xv;
			
			
	}


////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 7 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,ARM_GP_Vector*  xv,double x0,double x1,double x2,double x3,double x4,double x5
						  )
						    : itsContainer(6),itsVectorContainer(3)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsVectorContainer[0]=x;
			itsVectorContainer[1]=y;
			itsVectorContainer[2]=xv;
			
	}
////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 7 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,ARM_GP_Vector*  xv,double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6)
						    : itsContainer(7),itsVectorContainer(3)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsVectorContainer[0]=x;
			itsVectorContainer[1]=y;
			itsVectorContainer[2]=xv;
			
	}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 1 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,ARM_GP_Vector*  xv,ARM_GP_Vector*  y2,double x0) : itsContainer(1,x0)
,itsVectorContainer(4)
{	itsContainer[0]=x0;
	itsVectorContainer[0]=x;
	itsVectorContainer[1]=y;
	itsVectorContainer[2]=xv;
	itsVectorContainer[3]=y2;}


////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 2 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,ARM_GP_Vector*  xv,ARM_GP_Vector*  y2,double x0,double x1):itsContainer(2)
,itsVectorContainer(4)
{
	itsContainer[0]=x0;
	itsContainer[1]=x1;
	itsVectorContainer[0]=x;
	itsVectorContainer[1]=y;
	itsVectorContainer[2]=xv;
	itsVectorContainer[3]=y2;
}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 3 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////

ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,ARM_GP_Vector*  xv,ARM_GP_Vector*  y2,double x0,double x1,double x2) : itsContainer(3)
,itsVectorContainer(4)
{
	itsContainer[0]=x0;
	itsContainer[1]=x1;
	itsContainer[2]=x2;
	itsVectorContainer[0]=x;
	itsVectorContainer[1]=y;
	itsVectorContainer[2]=xv;
	itsVectorContainer[3]=y2;
		
}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 4 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,ARM_GP_Vector*  xv,ARM_GP_Vector*  y2,double x0,double x1,double x2,double x3) : itsContainer(4)
,itsVectorContainer(4)
{
		itsContainer[0]=x0;
		itsContainer[1]=x1;
		itsContainer[2]=x2;
		itsContainer[3]=x3;
		itsVectorContainer[0]=x;
		itsVectorContainer[1]=y;
		itsVectorContainer[2]=xv;
		itsVectorContainer[3]=y2;
		
}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 5 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,ARM_GP_Vector*  xv,ARM_GP_Vector*  y2,double x0,double x1,double x2,double x3,double x4) : itsContainer(5)
,itsVectorContainer(4)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsVectorContainer[0]=x;
			itsVectorContainer[1]=y;
			itsVectorContainer[2]=xv;
			itsVectorContainer[3]=y2;
			
	}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 6 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,ARM_GP_Vector*  xv,ARM_GP_Vector*  y2,double x0,double x1,double x2,double x3,double x4,double x5)
 : itsContainer(6),itsVectorContainer(4)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsVectorContainer[0]=x;
			itsVectorContainer[1]=y;
			itsVectorContainer[2]=xv;
			itsVectorContainer[3]=y2;
		
	}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 7 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,ARM_GP_Vector*  xv,ARM_GP_Vector*  y2,double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6)
						    : itsContainer(7),itsVectorContainer(4)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsVectorContainer[0]=x;
			itsVectorContainer[1]=y;
			itsVectorContainer[2]=xv;
			itsVectorContainer[3]=y2;
			
	}
////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 7 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,ARM_GP_Vector*  xv,ARM_GP_Vector*  y2,double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14,double x15,double x16,double x17)
						    : itsContainer(18),itsVectorContainer(4)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			itsContainer[15]=x15;
			itsContainer[16]=x16;
			itsContainer[17]=x17;
			itsVectorContainer[0]=x;
			itsVectorContainer[1]=y;
			itsVectorContainer[2]=xv;
			itsVectorContainer[3]=y2;
			
	}
////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 7 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,ARM_GP_Vector*  xv,ARM_GP_Vector*  y2,double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14,double x15,double x16,double x17,double x18)
						    : itsContainer(19),itsVectorContainer(4)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			itsContainer[15]=x15;
			itsContainer[16]=x16;
			itsContainer[17]=x17;
			itsContainer[18]=x18;
			itsVectorContainer[0]=x;
			itsVectorContainer[1]=y;
			itsVectorContainer[2]=xv;
			itsVectorContainer[3]=y2;
			
	}
////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 7 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(ARM_GP_Vector*  x,ARM_GP_Vector*  y,ARM_GP_Vector*  xv,ARM_GP_Vector*  y2,double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19)
						    : itsContainer(20),itsVectorContainer(4)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			itsContainer[15]=x15;
			itsContainer[16]=x16;
			itsContainer[17]=x17;
			itsContainer[18]=x18;
			itsContainer[19]=x19;
			itsVectorContainer[0]=x;
			itsVectorContainer[1]=y;
			itsVectorContainer[2]=xv;
			itsVectorContainer[3]=y2;
			
	}


////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 1 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0) : itsContainer(1,x0)
{}


////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 2 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1):itsContainer(2)
{
	itsContainer[0]=x0;
	itsContainer[1]=x1;
}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 3 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////

ArgumentList::ArgumentList(double x0,double x1,double x2) : itsContainer(3)
{
	itsContainer[0]=x0;
	itsContainer[1]=x1;
	itsContainer[2]=x2;
		
}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 4 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3) : itsContainer(4)
{
		itsContainer[0]=x0;
		itsContainer[1]=x1;
		itsContainer[2]=x2;
		itsContainer[3]=x3;
		
}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 5 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4) : itsContainer(5)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			
	}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 6 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5)
 : itsContainer(6)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
		
	}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 7 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6)
						    : itsContainer(7)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			
	}



////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 8 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7)
						    : itsContainer(8)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			
	}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 9 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8)
						    : itsContainer(9)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			
	}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 10 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9)
						    : itsContainer(10)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			
	}


////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 11 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10)
						    : itsContainer(11)
	{
			itsContainer=*new vector<double>(11);
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			
	}


////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 12 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11)
						    : itsContainer(12)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			
	}


////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 13 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12)
						    : itsContainer(13)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			
	}


////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 14 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13)
						    : itsContainer(14)
	{
			itsContainer=*new vector<double>(14);
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			
	}


////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 15 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14)
						    : itsContainer(15)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			
	}
////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 16 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14,double x15)
						    : itsContainer(16)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			itsContainer[15]=x15;
			
	}
////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 17 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14,double x15,double x16)
						    : itsContainer(17)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			itsContainer[15]=x15;
			itsContainer[16]=x16;
			
}
////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 18 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14,double x15,double x16,double x17)
						    : itsContainer(18)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			itsContainer[15]=x15;
			itsContainer[16]=x16;
			itsContainer[17]=x17;
			
	}
////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 19 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14,double x15,double x16,double x17,double x18)
						    : itsContainer(19)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			itsContainer[15]=x15;
			itsContainer[16]=x16;
			itsContainer[17]=x17;
			itsContainer[18]=x18;
			
	}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 20 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19)
						    : itsContainer(20)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			itsContainer[15]=x15;
			itsContainer[16]=x16;
			itsContainer[17]=x17;
			itsContainer[18]=x18;
			itsContainer[19]=x19;
			
	}
////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 21 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,double x20)
						    : itsContainer(21)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			itsContainer[15]=x15;
			itsContainer[16]=x16;
			itsContainer[17]=x17;
			itsContainer[18]=x18;
			itsContainer[19]=x19;
			itsContainer[20]=x20;
			
	}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 22 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,double x20,double x21)
						    : itsContainer(22)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			itsContainer[15]=x15;
			itsContainer[16]=x16;
			itsContainer[17]=x17;
			itsContainer[18]=x18;
			itsContainer[19]=x19;
			itsContainer[20]=x20;
			itsContainer[21]=x21;
			
	}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 23 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,
						   double x20,double x21,double x22)
						    : itsContainer(23)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			itsContainer[15]=x15;
			itsContainer[16]=x16;
			itsContainer[17]=x17;
			itsContainer[18]=x18;
			itsContainer[19]=x19;
			itsContainer[20]=x20;
			itsContainer[21]=x21;
			itsContainer[22]=x22;
			
	}
////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 24 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,
						   double x20,double x21,double x22,double x23)
						    : itsContainer(24)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			itsContainer[15]=x15;
			itsContainer[16]=x16;
			itsContainer[17]=x17;
			itsContainer[18]=x18;
			itsContainer[19]=x19;
			itsContainer[20]=x20;
			itsContainer[21]=x21;
			itsContainer[22]=x22;
			itsContainer[23]=x23;
			
	}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 25 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,
						   double x20,double x21,double x22,double x23,double x24)
						    : itsContainer(25)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			itsContainer[15]=x15;
			itsContainer[16]=x16;
			itsContainer[17]=x17;
			itsContainer[18]=x18;
			itsContainer[19]=x19;
			itsContainer[20]=x20;
			itsContainer[21]=x21;
			itsContainer[22]=x22;
			itsContainer[23]=x23;
			itsContainer[24]=x24;
			
	}
////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 26 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,
						   double x20,double x21,double x22,double x23,double x24,double x25)
						    : itsContainer(26)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			itsContainer[15]=x15;
			itsContainer[16]=x16;
			itsContainer[17]=x17;
			itsContainer[18]=x18;
			itsContainer[19]=x19;
			itsContainer[20]=x20;
			itsContainer[21]=x21;
			itsContainer[22]=x22;
			itsContainer[23]=x23;
			itsContainer[24]=x24;
			itsContainer[25]=x25;
			
	}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 27 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,
						   double x20,double x21,double x22,double x23,double x24,double x25,double x26)
						    : itsContainer(27)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			itsContainer[15]=x15;
			itsContainer[16]=x16;
			itsContainer[17]=x17;
			itsContainer[18]=x18;
			itsContainer[19]=x19;
			itsContainer[20]=x20;
			itsContainer[21]=x21;
			itsContainer[22]=x22;
			itsContainer[23]=x23;
			itsContainer[24]=x24;
			itsContainer[25]=x25;
			itsContainer[26]=x26;
			
	}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 28 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,
						   double x20,double x21,double x22,double x23,double x24,double x25,double x26,double x27)
						    : itsContainer(28)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			itsContainer[15]=x15;
			itsContainer[16]=x16;
			itsContainer[17]=x17;
			itsContainer[18]=x18;
			itsContainer[19]=x19;
			itsContainer[20]=x20;
			itsContainer[21]=x21;
			itsContainer[22]=x22;
			itsContainer[23]=x23;
			itsContainer[24]=x24;
			itsContainer[25]=x25;
			itsContainer[26]=x26;
			itsContainer[27]=x27;
			
	}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 29 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,
						   double x20,double x21,double x22,double x23,double x24,double x25,double x26,double x27,
						   double x28)
						    : itsContainer(29)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			itsContainer[15]=x15;
			itsContainer[16]=x16;
			itsContainer[17]=x17;
			itsContainer[18]=x18;
			itsContainer[19]=x19;
			itsContainer[20]=x20;
			itsContainer[21]=x21;
			itsContainer[22]=x22;
			itsContainer[23]=x23;
			itsContainer[24]=x24;
			itsContainer[25]=x25;
			itsContainer[26]=x26;
			itsContainer[27]=x27;
			itsContainer[28]=x28;
			
	}


////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 30 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,
						   double x20,double x21,double x22,double x23,double x24,double x25,double x26,double x27,
						   double x28,double x29)
						    : itsContainer(30)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			itsContainer[15]=x15;
			itsContainer[16]=x16;
			itsContainer[17]=x17;
			itsContainer[18]=x18;
			itsContainer[19]=x19;
			itsContainer[20]=x20;
			itsContainer[21]=x21;
			itsContainer[22]=x22;
			itsContainer[23]=x23;
			itsContainer[24]=x24;
			itsContainer[25]=x25;
			itsContainer[26]=x26;
			itsContainer[27]=x27;
			itsContainer[28]=x28;
			itsContainer[29]=x29;
			
	}
////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 31 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,
						   double x20,double x21,double x22,double x23,double x24,double x25,double x26,double x27,
						   double x28,double x29,double x30)
						    : itsContainer(31)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			itsContainer[15]=x15;
			itsContainer[16]=x16;
			itsContainer[17]=x17;
			itsContainer[18]=x18;
			itsContainer[19]=x19;
			itsContainer[20]=x20;
			itsContainer[21]=x21;
			itsContainer[22]=x22;
			itsContainer[23]=x23;
			itsContainer[24]=x24;
			itsContainer[25]=x25;
			itsContainer[26]=x26;
			itsContainer[27]=x27;
			itsContainer[28]=x28;
			itsContainer[29]=x29;
			itsContainer[30]=x30;
			
	}
////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ArgumentList 
///	Routine: ArgumentList
/// Arguments: 32 double that need to be provided as argument of an associated expression 
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////////////////////////////////////////
ArgumentList::ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,
						   double x6,double x7,double x8,double x9,double x10,double x11,
						   double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,
						   double x20,double x21,double x22,double x23,double x24,double x25,double x26,double x27,
						   double x28,double x29,double x30,double x31)
						    : itsContainer(32)
	{
			itsContainer[0]=x0;
			itsContainer[1]=x1;
			itsContainer[2]=x2;
			itsContainer[3]=x3;
			itsContainer[4]=x4;
			itsContainer[5]=x5;
			itsContainer[6]=x6;
			itsContainer[7]=x7;
			itsContainer[8]=x8;
			itsContainer[9]=x9;
			itsContainer[10]=x10;
			itsContainer[11]=x11;
			itsContainer[12]=x12;
			itsContainer[13]=x13;
			itsContainer[14]=x14;
			itsContainer[15]=x15;
			itsContainer[16]=x16;
			itsContainer[17]=x17;
			itsContainer[18]=x18;
			itsContainer[19]=x19;
			itsContainer[20]=x20;
			itsContainer[21]=x21;
			itsContainer[22]=x22;
			itsContainer[23]=x23;
			itsContainer[24]=x24;
			itsContainer[25]=x25;
			itsContainer[26]=x26;
			itsContainer[27]=x27;
			itsContainer[28]=x28;
			itsContainer[29]=x29;
			itsContainer[30]=x30;
			itsContainer[31]=x31;
			
	}

ArgumentList* ArgumentList::sublist(int n) const
{ 
	int p=size2();
	ArgumentList* new_arg;
	if(p==0)
	{
		new_arg=new ArgumentList(n);
		copy(itsContainer.begin(), itsContainer.begin()+n,
			new_arg->get_itsContainer().begin()	);
		return new_arg;
	}
	else
	{
		new_arg=new ArgumentList(n,p);
		copy(itsContainer.begin(), itsContainer.begin()+n,
			new_arg->get_itsContainer().begin()	);
		copy(itsVectorContainer.begin(), itsVectorContainer.begin()+n,
			new_arg->get_itsVectorContainer().begin()	);
		return new_arg;
	}

}
////////////////////////////////////////////////////////////////////////////////////////
///	Class  : Expression 
///	Routine: standard_first_derivative
/// Arguments:  double(* f)(ArgumentList * a), int i,ArgumentList* a,double s
///	Returns: double
///	Action :  computes the first derivative given by the integer i of the associated expression
/// and the given shift s
////////////////////////////////////////////////////////////////////////////////////////

double standard_first_derivative(double(* f)(const ArgumentList& a),int i,const ArgumentList& a,double s)
{
	ArgumentList* shifted_args=new ArgumentList(a);
	double v1=(*f)(a);
	double shift= fabs((a.get_itsContainer())[i]) < ARM_NumericConstants::ARM_TOLERENCE ? s : (a.get_itsContainer())[i]*s;
	double beforeshift = (shifted_args->get_itsContainer())[i];
	shifted_args->set_nth(i,shift+a.get_itsContainer()[i]);
	double aftershift = (shifted_args->get_itsContainer())[i];
	double valf0=a.get_itsContainer()[i];
	double valf=(*shifted_args)[i];
	double v2=(*f)(*shifted_args);
	delete shifted_args;
	return (v2-v1)/shift;
}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : Expression 
///	Routine: centred_standard_first_derivative
/// Arguments:  double(* f)(ArgumentList * a), int i,ArgumentList* a,double s
///	Returns: double
///	Action :  computes the first derivative given by the integer i of the associated expression
/// and the given shift s
////////////////////////////////////////////////////////////////////////////////////////

double centred_standard_first_derivative(double(* f)(const ArgumentList& a),int i,const ArgumentList& a,double s)
{
	ArgumentList* shifted_args1=new ArgumentList(a);
	ArgumentList* shifted_args2=new ArgumentList(a);

	double shift= fabs((a.get_itsContainer())[i]) < ARM_NumericConstants::ARM_TOLERENCE ? s : (a.get_itsContainer())[i]*s;
	double value=a.get_itsContainer()[i];
	shifted_args1->set_nth(i,-shift+value);
	shifted_args2->set_nth(i,shift+value);

	double v1=(*f)(*shifted_args1);
	double v2=(*f)(*shifted_args2);
	delete shifted_args1;
	delete shifted_args2;
	return (v2-v1)/(2.*shift);
}

////////////////////////////////////////////////////////////////////////////////////////
///	Class  : Expression 
///	Routine: standard_second_derivative
/// Arguments:  double(* f)(ArgumentList * a),int i,int j,ArgumentList* a,double s1, double s2
///	Returns: double
///	Action :  computes the second derivative given by the integer i and j of the associated expression
/// and the given shift s
////////////////////////////////////////////////////////////////////////////////////////

double standard_second_derivative(double(* f)(const ArgumentList& a),int i,int j, const ArgumentList& a,double s1,double s2)
{
	ArgumentList shifted_args1(a);
	ArgumentList shifted_args2(a);
	double shift1=(a.get_itsContainer())[i] * s1;
	double shift2=(a.get_itsContainer())[j] * s2;
	shifted_args1.set_nth(i,shifted_args1.get_itsContainer()[i] + shift1);
	shifted_args2.set_nth(j,shifted_args2.get_itsContainer()[j] + shift2);
	ArgumentList shifted_args21(shifted_args2);
	shifted_args21.set_nth(i,shifted_args21.get_itsContainer()[i] + shift1);
	double diff1=(*f)(shifted_args1)-(*f)(a);
	double diff2=(*f)(shifted_args21)-(*f)(shifted_args2);

	return (diff2-diff1)/(shift1*shift2);
}

ArgumentList_Checking_Result::ArgumentList_Checking_Result(bool res,string trb): checkresult(res) 
{
	trouble=trb;
}

ArgumentList_Checking_Result operator&&(ArgumentList_Checking_Result arg1,ArgumentList_Checking_Result arg2)
{
	if((arg1.checkresult)&&(arg2.checkresult)) return ArgumentList_Checking_Result(true,"");
	else
	return ArgumentList_Checking_Result(false,arg1.trouble+string(" and ")+arg2.trouble);
}

ArgumentList_Checking_Result operator||(ArgumentList_Checking_Result arg1,ArgumentList_Checking_Result arg2)
{
	if((arg1.checkresult)||(arg2.checkresult)) return ArgumentList_Checking_Result(true,"");
	else
	return ArgumentList_Checking_Result(false,arg1.trouble+string(" and ")+arg2.trouble);
}



CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/