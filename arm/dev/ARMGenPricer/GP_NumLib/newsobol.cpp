/*!
 *
 * Copyright (c) CDC IXIS CM March 2005 Paris
 *
 *	\file newsobol.cpp
 *  \brief 
 *	\author  A. Triki
 *	\version 1.0
 */


#include "gpbase/removeidentifiedwarning.h"
#include "gpnumlib/newsobol.h"
//#include "gpnumlib/sobolinit.h"
#include "gpbase/ostringstream.h"
#include <cstdlib>
#include <ctime>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////////////////
///             ARM_Sobol multidimentionnal sequence		 ///
////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_NewSobol
///	Routine: constructor
///	Returns: 
///	Action : build the object
////////////////////////////////////////////////////
ARM_NewSobol::ARM_NewSobol(long Seed, int firstSimulations)
:	ARM_QuasiRandom(firstSimulations),
	itsIndex(0),
	itsValues(),
	itsSeed(Seed)
{
}
////////////////////////////////////////////////////
///	Class  : ARM_NewSobol
///	Routine: D_UNIFORM_01
///	Returns: 
///	Action : returns a unit double precision pseudorandom number (to be changed by better generator)
////////////////////////////////////////////////////


double ARM_NewSobol::d_uniform_01 ( int *seed )

{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}

////////////////////////////////////////////////////
///	Class  : ARM_NewSobol
///	Routine: i4_bit_hi1
///	Returns: 
///	Action :  I4_BIT_HI1 returns the position of the high 1 bit base 2 in an integer
////////////////////////////////////////////////////
//  Example:
//
//       N    Binary    Hi 1
//    ----    --------  ----
//       0           0     0
//       1           1     1
//       2          10     2
//       3          11     2 
//       4         100     3
//       5         101     3
//       6         110     3
//       7         111     3
//       8        1000     4
//       9        1001     4
//      10        1010     4
//      11        1011     4
//      12        1100     4
//      13        1101     4
//      14        1110     4
//      15        1111     4
//      16       10000     5
//      17       10001     5
//    1023  1111111111    10
//    1024 10000000000    11
//    1025 10000000001    11

int ARM_NewSobol::i4_bit_hi1 ( int n )

{
  int bit;

  bit = 0;

  while ( 0 < n )
  {
    bit = bit + 1;
    n = n / 2;
  }

  return bit;
}

////////////////////////////////////////////////////
///	Class  : ARM_NewSobol
///	Routine: i4_bit_lo0
///	Returns: 
///	Action :  returns the position of the low 0 bit base 2 in an integer
////////////////////////////////////////////////////
int ARM_NewSobol::i4_bit_lo0 ( int n )
{
  int bit;
  int n2;

  bit = 0;

  while ( true )
  {
    bit = bit + 1;
    n2 = n / 2;

    if ( n == 2 * n2 )
    {
      break;
    }

    n = n2;

  }

  return bit;
}
////////////////////////////////////////////////////
///	Class  : ARM_NewSobol
///	Routine: i4_sobol
///	Returns: 
///	Action : generates a new quasirandom Sobol vector with each call
////////////////////////////////////////////////////
//  Parameters:
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//    DIM_NUM must satisfy 2 <= DIM_NUM <= 40.
//
//    Input/output, int *SEED, the "seed" for the sequence.
//    This is essentially the index in the sequence of the quasirandom
//    value to be generated.  On output, SEED has been set to the
//    appropriate next value, usually simply SEED+1.
//    If SEED is less than 0 on input, it is treated as though it were 0.
//    An input value of 0 requests the first (0-th) element of the sequence.
//
//    Output, float QUASI(DIM_NUM), the next quasirandom vector.

void ARM_NewSobol::i4_sobol ( int dim_num, int *seed, float quasi[ ] )
{
# define DIM_MAX 60

  static int atmost = 1073741823;
  static int dim_num_save = 0;
  int i;
  bool includ[8];
  static bool initialized = false;
  int j;
  int j2;
  int k;
  int l;
  static int lastq[DIM_MAX];
  int m;
  static int maxcol;
  int newv;
  static int poly[DIM_MAX] =
  {
        1,   3,   7,  11,  13,  19,  25,  37,  59,  47,
       61,  55,  41,  67,  97,  91, 109, 103, 115, 131,
      193, 137, 145, 143, 241, 157, 185, 167, 229, 171,
      213, 191, 253, 203, 211, 239, 247, 285, 369, 299 
  };
  static float recipd;
  static int seed_save = 0;
  int seed_temp;
  static int v[DIM_MAX][30];
//
  if ( !initialized || dim_num != dim_num_save )
  {
    initialized = true;
//
//  Initialize (part of) V.
//
    v[ 0][0] = 1;
    v[ 1][0] = 1;
    v[ 2][0] = 1;
    v[ 3][0] = 1;
    v[ 4][0] = 1;
    v[ 5][0] = 1;
    v[ 6][0] = 1;
    v[ 7][0] = 1;
    v[ 8][0] = 1;
    v[ 9][0] = 1;
    v[10][0] = 1;
    v[11][0] = 1;
    v[12][0] = 1;
    v[13][0] = 1;
    v[14][0] = 1;
    v[15][0] = 1;
    v[16][0] = 1;
    v[17][0] = 1;
    v[18][0] = 1;
    v[19][0] = 1;
    v[20][0] = 1;
    v[21][0] = 1;
    v[22][0] = 1;
    v[23][0] = 1;
    v[24][0] = 1;
    v[25][0] = 1;
    v[26][0] = 1;
    v[27][0] = 1;
    v[28][0] = 1;
    v[29][0] = 1;
    v[30][0] = 1;
    v[31][0] = 1;
    v[32][0] = 1;
    v[33][0] = 1;
    v[34][0] = 1;
    v[35][0] = 1;
    v[36][0] = 1;
    v[37][0] = 1;
    v[38][0] = 1;
    v[39][0] = 1;

    v[ 2][1] = 1;
    v[ 3][1] = 3;
    v[ 4][1] = 1;
    v[ 5][1] = 3;
    v[ 6][1] = 1;
    v[ 7][1] = 3;
    v[ 8][1] = 3;
    v[ 9][1] = 1;
    v[10][1] = 3;
    v[11][1] = 1;
    v[12][1] = 3;
    v[13][1] = 1;
    v[14][1] = 3;
    v[15][1] = 1;
    v[16][1] = 1;
    v[17][1] = 3;
    v[18][1] = 1;
    v[19][1] = 3;
    v[20][1] = 1;
    v[21][1] = 3;
    v[22][1] = 1;
    v[23][1] = 3;
    v[24][1] = 3;
    v[25][1] = 1;
    v[26][1] = 3;
    v[27][1] = 1;
    v[28][1] = 3;
    v[29][1] = 1;
    v[30][1] = 3;
    v[31][1] = 1;
    v[32][1] = 1;
    v[33][1] = 3;
    v[34][1] = 1;
    v[35][1] = 3;
    v[36][1] = 1;
    v[37][1] = 3;
    v[38][1] = 1;
    v[39][1] = 3;

    v[ 3][2] = 7;
    v[ 4][2] = 5;
    v[ 5][2] = 1;
    v[ 6][2] = 3;
    v[ 7][2] = 3;
    v[ 8][2] = 7;
    v[ 9][2] = 5;
    v[10][2] = 5;
    v[11][2] = 7;
    v[12][2] = 7;
    v[13][2] = 1;
    v[14][2] = 3;
    v[15][2] = 3;
    v[16][2] = 7;
    v[17][2] = 5;
    v[18][2] = 1;
    v[19][2] = 1;
    v[20][2] = 5;
    v[21][2] = 3;
    v[22][2] = 3;
    v[23][2] = 1;
    v[24][2] = 7;
    v[25][2] = 5;
    v[26][2] = 1;
    v[27][2] = 3;
    v[28][2] = 3;
    v[29][2] = 7;
    v[30][2] = 5;
    v[31][2] = 1;
    v[32][2] = 1;
    v[33][2] = 5;
    v[34][2] = 7;
    v[35][2] = 7;
    v[36][2] = 5;
    v[37][2] = 1;
    v[38][2] = 3;
    v[39][2] = 3;

    v[ 5][3] =  1;
    v[ 6][3] =  7;
    v[ 7][3] =  9;
    v[ 8][3] = 13;
    v[ 9][3] = 11;
    v[10][3] =  1;
    v[11][3] =  3;
    v[12][3] =  7;
    v[13][3] =  9;
    v[14][3] =  5;
    v[15][3] = 13;
    v[16][3] = 13;
    v[17][3] = 11;
    v[18][3] =  3;
    v[19][3] = 15;
    v[20][3] =  5;
    v[21][3] =  3;
    v[22][3] = 15;
    v[23][3] =  7;
    v[24][3] =  9;
    v[25][3] = 13;
    v[26][3] =  9;
    v[27][3] =  1;
    v[28][3] = 11;
    v[29][3] =  7;
    v[30][3] =  5;
    v[31][3] = 15;
    v[32][3] =  1;
    v[33][3] = 15;
    v[34][3] = 11;
    v[35][3] =  5;
    v[36][3] =  3;
    v[37][3] =  1;
    v[38][3] =  7;
    v[39][3] =  9;
  
    v[ 7][4] =  9;
    v[ 8][4] =  3;
    v[ 9][4] = 27;
    v[10][4] = 15;
    v[11][4] = 29;
    v[12][4] = 21;
    v[13][4] = 23;
    v[14][4] = 19;
    v[15][4] = 11;
    v[16][4] = 25;
    v[17][4] =  7;
    v[18][4] = 13;
    v[19][4] = 17;
    v[20][4] =  1;
    v[21][4] = 25;
    v[22][4] = 29;
    v[23][4] =  3;
    v[24][4] = 31;
    v[25][4] = 11;
    v[26][4] =  5;
    v[27][4] = 23;
    v[28][4] = 27;
    v[29][4] = 19;
    v[30][4] = 21;
    v[31][4] =  5;
    v[32][4] =  1;
    v[33][4] = 17;
    v[34][4] = 13;
    v[35][4] =  7;
    v[36][4] = 15;
    v[37][4] =  9;
    v[38][4] = 31;
    v[39][4] =  9;

    v[13][5] = 37;
    v[14][5] = 33;
    v[15][5] =  7;
    v[16][5] =  5;
    v[17][5] = 11;
    v[18][5] = 39;
    v[19][5] = 63;
    v[20][5] = 27;
    v[21][5] = 17;
    v[22][5] = 15;
    v[23][5] = 23;
    v[24][5] = 29;
    v[25][5] =  3;
    v[26][5] = 21;
    v[27][5] = 13;
    v[28][5] = 31;
    v[29][5] = 25;
    v[30][5] =  9;
    v[31][5] = 49;
    v[32][5] = 33;
    v[33][5] = 19;
    v[34][5] = 29;
    v[35][5] = 11;
    v[36][5] = 19;
    v[37][5] = 27;
    v[38][5] = 15;
    v[39][5] = 25;

    v[19][6] =  13;
    v[20][6] =  35;
    v[21][6] = 115;
    v[22][6] =  41;
    v[23][6] =  79;
    v[24][6] =  17;
    v[25][6] =  29;
    v[26][6] = 119;
    v[27][6] =  75;
    v[28][6] =  73;
    v[29][6] = 105;
    v[30][6] =   7;
    v[31][6] =  59;
    v[32][6] =  65;
    v[33][6] =  21;
    v[34][6] =   3;
    v[35][6] = 113;
    v[36][6] =  61;
    v[37][6] =  89;
    v[38][6] =  45;
    v[39][6] = 107;

    v[37][7] =  7;
    v[38][7] = 23;
    v[39][7] = 39;
//
//  Check parameters.
//
	if ( dim_num < 1 || DIM_MAX < dim_num )
    {
		CC_Ostringstream os;
		os  << "I4_SOBOL The spatial dimension DIM_NUM" 
			<< dim_num
			<< "should satisfy 2 <= DIM_NUM <=" 
			<< DIM_MAX;
	    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str());
    }

    dim_num_save = dim_num;
//
//  Find the number of bits in ATMOST.
//
    maxcol = i4_bit_hi1 ( atmost );
//
//  Initialize row 1 of V.
//
    for ( j = 0; j < maxcol; j++ )
    {
      v[0][j] = 1;
    }
//
//  Initialize the remaining rows of V.
//
    for ( i = 1; i < dim_num; i++ )
    {
//
//  The bit pattern of the integer POLY(I) gives the form
//  of polynomial I.
//
//  Find the degree of polynomial I from binary encoding.
//
      j = poly[i];
      m = 0;

      while ( true )
      {
        j = j / 2;
        if ( j <= 0 )
        {
          break;
        }
        m = m + 1;
      }
//
//  We expand this bit pattern to separate components
//  of the logical array INCLUD.
//
      j = poly[i];
      for ( k = m-1; 0 <= k; k-- )
      {
        j2 = j / 2;
        includ[k] = ( j != ( 2 * j2 ) );
        j = j2;
      }
//
//  Calculate the remaining elements of row I as explained
//  in Bratley and Fox, section 2.
//
//  Some tricky indexing here.  Did I change it correctly?
//
      for ( j = m; j < maxcol; j++ )
      {
        newv = v[i][j-m];
        l = 1;

        for ( k = 0; k < m; k++ )
        {
          l = 2 * l;

          if ( includ[k] )
          {
            newv = ( newv ^ ( l * v[i][j-k-1] ) );
          }

        }

        v[i][j] = newv;

      }

    }
//
//  Multiply columns of V by appropriate power of 2.
//
    l = 1;
    for ( j = maxcol-2; 0 <= j; j-- )
    {
      l = 2 * l;
      for ( i = 0; i < dim_num; i++ )
      {
        v[i][j] = v[i][j] * l;
      }
    }
//
//  RECIPD is 1/(common denominator of the elements in V).
//
    recipd = 1.0E+00 / ( ( float ) ( 2 * l ) );
  }

  if ( *seed < 0 )
  {
    *seed = 0;
  }

  if ( *seed == 0 )
  {
    l = 1;
    for ( i = 0; i < dim_num; i++ )
    {
      lastq[i] = 0;
    }
  }
  else if ( *seed == seed_save + 1 )
  {
    l = i4_bit_lo0 ( *seed );
  }
  else if ( *seed <= seed_save )
  {
    seed_save = 0;
    l = 1;
    for ( i = 0; i < dim_num; i++ )
    {
      lastq[i] = 0;
    }

    for ( seed_temp = seed_save; seed_temp <= (*seed)-1; seed_temp++ )
    {

      l = i4_bit_lo0 ( seed_temp );

      for ( i = 0; i < dim_num; i++ )
      {
        lastq[i] = ( lastq[i] ^ v[i][l-1] );
      }

    }

    l = i4_bit_lo0 ( *seed );
  }
  else if ( seed_save+1 < *seed )
  {
    for ( seed_temp = seed_save+1; seed_temp <= (*seed)-1; seed_temp++ )
    {

      l = i4_bit_lo0 ( seed_temp );

      for ( i = 0; i < dim_num; i++ )
      {
        lastq[i] = ( lastq[i] ^ v[i][l-1] );
      }

    }

    l = i4_bit_lo0 ( *seed );

  }
//
//  Check that the user is not calling too many times!
//
  if ( maxcol < l )
  {
	  CC_Ostringstream os;
	  os  << "I4_SOBOL Two many calls";
	  throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str());
    }
//
//  Calculate the new components of QUASI.
//  The caret indicates the bitwise exclusive OR.
//
  for ( i = 0; i < dim_num; i++ )
  {
    quasi[i] = ( ( float ) lastq[i] ) * recipd;

    lastq[i] = ( lastq[i] ^ v[i][l-1] );
  }

  seed_save = *seed;
  *seed = *seed + 1;

  return;
# undef MAX_DIM
}


////////////////////////////////////////////////////
///	Class  : ARM_NewSobol
///	Routine: DrawOne
///	Returns: 
///	Action : draw one set of values
////////////////////////////////////////////////////
void ARM_NewSobol::DrawAll()
{
	/// this routine calculates the next random number 
	int seed = itsIndex;
	float r[DIM_MAX];
	int seed_in;
	int seed_out;
	seed_in = seed;
    i4_sobol ( itsDim, &seed, r );
    seed_out = seed;
	for (int k=0;k< itsDim;k++) 
	{
		itsCurrentValues[itsCurrentBucketPos][itsCurrentPos][k]=r[k];
	}
	itsIndex++;
}


////////////////////////////////////////////////////
///	Class  : ARM_NewSobol
///	Routine: Init
///	Returns:  
///	Action : initilisation routine
////////////////////////////////////////////////////
void ARM_NewSobol::Init()
{
	itsIndex=itsSeed;
}

////////////////////////////////////////////////////
///	Class  : ARM_NewSobol
///	Routine: SetDim
///	Returns: 
///	Action : reset the dimension!
////////////////////////////////////////////////////
void ARM_NewSobol::SetDim( size_t dim )
{
	itsDim = dim;
	Init();
}


////////////////////////////////////////////////////
///	Class  : ARM_NewSobol
///	Routine: reset
///	Returns: 
///	Action : reset the generator
////////////////////////////////////////////////////
void ARM_NewSobol::reset( size_t dim, size_t nbOfPoints )
{
	SetDim(dim);
	Init();
}


////////////////////////////////////////////////////
///	Class  : ARM_Sobol
///	Routine: Copy constructor and operator=
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_NewSobol::ARM_NewSobol(const ARM_NewSobol& rhs)
:	ARM_QuasiRandom( rhs), 
	itsIndex(	rhs.itsIndex	),
	itsValues ( rhs.itsValues   ),
	itsSeed		(rhs.itsSeed	)
{}

ARM_NewSobol& ARM_NewSobol::operator=(const ARM_NewSobol& rhs )
{
	if(this!=&rhs)
	{
		ARM_QuasiRandom::operator=( rhs );
		itsIndex	= rhs.itsIndex;
		itsValues   = rhs.itsValues;
		itsSeed		= rhs.itsSeed;
	}
	return *this;
}



////////////////////////////////////////////////////
///	Class  : ARM_NewSobol
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_NewSobol::~ARM_NewSobol()
{}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_NewSobol
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_NewSobol::Clone() const
{
	return new ARM_NewSobol(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_Sobol
///	Routine: toString
///	Returns: 
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_NewSobol::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << "Sobol Sequence with dim " << itsDim;
	return os.str();
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


