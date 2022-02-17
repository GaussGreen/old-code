/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file mersenetwister.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2003
 */


#include "gpnumlib/mersenetwister.h"

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
////////////////////////////////////////////////////
//////////	ARM_FastMersenneTwister ////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
///   \class ARM_FastMersenneTwister
///    
///	  Adapted from A C-program for MT19937, with initialization improved 2002/2/10.
///   Coded by Takuji Nishimura and Makoto Matsumoto.
///   This is a faster version by taking Shawn Cokus's optimization,
///   Matthe Bellew's simplification, Isaku Wada's real version.
///   
///   Before using, initialize the state by using init_genrand(seed) 
///   or init_by_array(init_key, key_length).
///   
///   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
///   All rights reserved.                          
///   
///   Redistribution and use in source and binary forms, with or without
///   modification, are permitted provided that the following conditions
///   are met:
///   
///   1. Redistributions of source code must retain the above copyright
///          notice, this list of conditions and the following disclaimer.
///   
///   2. Redistributions in binary form must reproduce the above copyright
///   notice, this list of conditions and the following disclaimer in the
///      documentation and/or other materials provided with the distribution.
///   
///   3. The names of its contributors may not be used to endorse or promote 
///      products derived from this software without specific prior written 
///      permission.
///
///   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
///   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
///   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
///   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
///   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
///   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
///   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
///   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
///   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
///   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
///   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///
///   
///   Any feedback is very welcome.
///   http://www.math.keio.ac.jp/matumoto/emt.html
///   email: matumoto@math.keio.ac.jp
///
/////////////////////////////////////////////////////////////


/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UMASK 0x80000000UL /* most significant w-r bits */
#define LMASK 0x7fffffffUL /* least significant r bits */
#define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )
#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))


////////////////////////////////////////////////////
///	Class  : ARM_FastMersenneTwister  
///	Routine: constructor and init_genrand
///	Returns: void
///	Action : initializes mt[N] with a seed 
////////////////////////////////////////////////////
ARM_FastMersenneTwister::ARM_FastMersenneTwister ( unsigned long seed)
:	left(1), initf(0), itsSeed(seed), itsState(N)
{
	init_genrand(itsSeed);
}

////////////////////////////////////////////////////
///	Class  : ARM_FastMersenneTwister  
///	Routine: reset
///	Returns: void
///	Action : reset the random number generator
////////////////////////////////////////////////////
void ARM_FastMersenneTwister::reset( size_t dim, const ARM_GP_T_Vector<size_t>& NbOfPoints, size_t factorNb ) 
{
	left = 1;
	initf= 0;
	init_genrand(itsSeed);
}


void ARM_FastMersenneTwister::init_genrand(unsigned long s)
{
    int j;
    itsState[0]= s & 0xffffffffUL;
    for (j=1; j<N; j++) {
        itsState[j] = (1812433253UL * (itsState[j-1] ^ (itsState[j-1] >> 30)) + j); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array itsState[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        itsState[j] &= 0xffffffffUL;  /* for >32 bit machines */
    }
    left = 1; initf = 1;
}



////////////////////////////////////////////////////
///	Class  : ARM_OriginalMersenneTwister  
///	Routine: next_state and DrawOne
///	Returns: double for DrawOne
///	Action : 
////////////////////////////////////////////////////

void ARM_FastMersenneTwister::next_state()
{
// FIXMEFRED: mig.vc8 (22/05/2007 18:13:07):cast
	unsigned long *p=&(*itsState.begin());
    int j;

    /* if init_genrand() has not been called, */
    /* a default initial seed is used         */
    if (initf==0) init_genrand(5489UL);

    left = N;
// FIXMEFRED: mig.vc8 (22/05/2007 18:14:57):cast
	next = &(*itsState.begin());
    
    for (j=N-M+1; --j; p++) 
        *p = p[M] ^ TWIST(p[0], p[1]);

    for (j=M; --j; p++) 
        *p = p[M-N] ^ TWIST(p[0], p[1]);

    *p = p[M-N] ^ TWIST(p[0], itsState[0]);
}

/// generates a random number on [0,1]-real-interval
double ARM_FastMersenneTwister::DrawOne()
{
    unsigned long y;

    if (--left == 0) next_state();
    y = *next++;

    /// Tempering
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (double)y * (1.0/4294967295.0); 
	/// divided by 2^32-1 //// 
}


#undef N
#undef M
#undef MATRIX_A 
#undef UMASK
#undef LMASK 
#undef MIXBITS
#undef TWIST




////////////////////////////////////////////////////
///	Class  : ARM_FastMersenneTwister
///	Routine: Copy constructor and operator=
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_FastMersenneTwister::ARM_FastMersenneTwister(const ARM_FastMersenneTwister& rhs)
:	ARM_UniformGenerator( rhs	), 
	itsSeed(	rhs.itsSeed		),
	itsState(		rhs.itsState		),
	left(		rhs.left		),
	initf(		rhs.initf		),
	next(		rhs.next		)
{}

ARM_FastMersenneTwister& ARM_FastMersenneTwister::operator=(const ARM_FastMersenneTwister& rhs )
{
	if(this!=&rhs)
	{
		ARM_UniformGenerator::operator=( rhs );
		itsSeed	= rhs.itsSeed;
		itsState	= rhs.itsState;
		left	= rhs.left;
		initf	= rhs.initf;
		next	= rhs.next;
	}
	return *this;
}



////////////////////////////////////////////////////
///	Class  : ARM_FastMersenneTwister
///	Routine: Destructor, Clone
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_FastMersenneTwister::~ARM_FastMersenneTwister()
{}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_FastMersenneTwister
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_FastMersenneTwister::Clone() const
{
	return new ARM_FastMersenneTwister(*this);
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
//////////	ARM_FastMersenneTwister ////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
///   \class ARM_OriginalMersenneTwister 
///   
///   Adapted from a C-program for MT19937, with initialization improved 2002/1/26.
///   Coded by Takuji Nishimura and Makoto Matsumoto.
///
///   Before using, initialize the itsState by using init_genrand(seed)  
///
///   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
///   All rights reserved.                          
///
///   Redistribution and use in source and binary forms, with or without
///   modification, are permitted provided that the following conditions
///   are met:
///
///     1. Redistributions of source code must retain the above copyright
///        notice, this list of conditions and the following disclaimer.
///
///     2. Redistributions in binary form must reproduce the above copyright
///        notice, this list of conditions and the following disclaimer in the
///        documentation and/or other materials provided with the distribution.
///
///     3. The names of its contributors may not be used to endorse or promote 
///        products derived from this software without specific prior written 
///        permission.
///
///   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
///   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
///   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
///   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
///   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
///   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
///   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
///   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
///   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
///   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
///   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///
///
///   Any feedback is very welcome.
///   http://www.math.keio.ac.jp/matumoto/emt.html
///   email: matumoto@math.keio.ac.jp
///
/////////////////////////////////////////////////////////////

/// Period parameters
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /// constant vector a
#define UPPER_MASK 0x80000000UL /// most significant w-r bits
#define LOWER_MASK 0x7fffffffUL /// least significant r bits

////////////////////////////////////////////////////
///	Class  : ARM_OriginalMersenneTwister  
///	Routine: constructor and init_genrand
///	Returns: void
///	Action : initializes mt[N] with a seed 
////////////////////////////////////////////////////
ARM_OriginalMersenneTwister ::ARM_OriginalMersenneTwister ( unsigned long seed)
:	
	ARM_UniformGenerator(),
	itsSeed(seed),
	mt(N),
	mti(N+1)
{
	init_genrand(itsSeed);
}


////////////////////////////////////////////////////
///	Class  : ARM_OriginalMersenneTwister  
///	Routine: reset
///	Returns: void
///	Action : reset the random number generator
////////////////////////////////////////////////////
void ARM_OriginalMersenneTwister::reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList,  size_t factorNb )
{
	mti = N+1;
	init_genrand(itsSeed);
}


void ARM_OriginalMersenneTwister::init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) 
	{
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        
		/// See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. 
        /// In the previous versions, MSBs of the seed affect   
        /// only MSBs of the array mt[].                        
        /// 2002/01/09 modified by Makoto Matsumoto             
 
		mt[mti] &= 0xffffffffUL;
        /// for >32 bit machines 
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_OriginalMersenneTwister  
///	Routine: DrawOne
///	Returns: double
///	Action : do a next seed and get next number!
///			based on Mersenne twister random number generator
///			generates a random number on [0,1]-real-interval 
////////////////////////////////////////////////////

double ARM_OriginalMersenneTwister ::DrawOne()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /// mag01[x] = x * MATRIX_A  for x=0,1

    if (mti >= N)
	{ 
		/// generate N words at one time
        int kk;

        if (mti == N+1)				/// if init_genrand() has not been called,
            init_genrand(5489UL);	/// a default initial seed is used

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /// Tempering
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (double)y / 4294967295.0;
}

#undef N
#undef M
#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK




////////////////////////////////////////////////////
///	Class  : ARM_OriginalMersenneTwister
///	Routine: Copy constructor and operator=
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_OriginalMersenneTwister::ARM_OriginalMersenneTwister(const ARM_OriginalMersenneTwister& rhs)
:	ARM_UniformGenerator( rhs	), 
	itsSeed(	rhs.itsSeed		),
	mt(			rhs.mt			),
	mti(		rhs.mti			)
{}


ARM_OriginalMersenneTwister& ARM_OriginalMersenneTwister::operator=(const ARM_OriginalMersenneTwister& rhs )
{
	if(this!=&rhs)
	{
		ARM_UniformGenerator::operator=( rhs );
		itsSeed = rhs.itsSeed;
		mt		= rhs.mt;
		mti		= rhs.mti;	
	}

	return *this;
}



////////////////////////////////////////////////////
///	Class  : ARM_OriginalMersenneTwister
///	Routine: Destructor, Clone
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_OriginalMersenneTwister::~ARM_OriginalMersenneTwister()
{}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_OriginalMersenneTwister
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_OriginalMersenneTwister::Clone() const
{
	return new ARM_OriginalMersenneTwister(*this);
}

CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----