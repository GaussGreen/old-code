/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file sobol.h
 *
 *  \brief General file for the sobol sequence
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _INGPNUMLIB_SOBOL_H
#define _INGPNUMLIB_SOBOL_H

#include "gpbase/port.h"
#include "quasirandom.h"

CC_BEGIN_NAMESPACE( ARM )

static const int MAXBITsobol	= 30;
static const int MAXDIMsobol	= 366;

// degr�s des 366 premiers polynomes primitifs
static const int mdegsobol[]={
	0,
	1 ,
	2,
	3,3,
	4,4,
	5,5,5,5,5,5,
	6,6,6,6,6,6,
	7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
	8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
	9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
	9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
	9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
	10,10,10,10,10,10,10,10,10,10,
	10,10,10,10,10,10,10,10,10,10,
	10,10,10,10,10,10,10,10,10,10,
	10,10,10,10,10,10,10,10,10,10,
	10,10,10,10,10,10,10,10,10,10,
	10,10,10,10,10,10,10,10,10,10,
	11,11,11,11,11,11,11,11,11,11,
	11,11,11,11,11,11,11,11,11,11,
	11,11,11,11,11,11,11,11,11,11,
	11,11,11,11,11,11,11,11,11,11,
	11,11,11,11,11,11,11,11,11,11,
	11,11,11,11,11,11,11,11,11,11,
	11,11,11,11,11,11,11,11,11,11,
	11,11,11,11,11,11,11,11,11,11,
	11,11,11,11,11,11,11,11,11,11,
	11,11,11,11,11,11,11,11,11,11,
	11,11,11,11,11,11,11,11,11,11,
	11,11,11,11,11,11,11,11,11,11,
	11,11,11,11,11,11,11,11,11,11,
	11,11,11,11,11,11,11,11,11,11,
	11,11,11,11,11,11,11,11,11,11,
	11,11,11,11,11,11,11,11,11,11,
	11,11,11,11,11,11,11,11,11,11,
	11,11,11,11,11,11,11,11,11,11,
	11,11,11,11,11,11,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12,12,12,12,12,12,
	12,12,12,12,12
};

// les 366 premiers polyn�mes primitifs sous leur '�criture d�cimale'
static const int ipsobol[]={
	0,
	0,
	1,
	1,2,
	1,4,
	// 5 deg   totale 6
	2,   4,   7,   11,  13,  14,
	// 6 deg   totale 6
	1,   13,  16,  19,  22,  25,
	// 7 deg   totale 18
	1,   4,   7,   8,   14,  19,  21,  28,  31,  32,
	37,  41,  42,  50,  55,  56,  59,  62,
	// 8 deg   totale 16
	14,  21,  22,  38,  47,  49,  50,  52,  56,  67,
	70,  84,  97,  103, 115, 122,            
	// 9 deg   totale 48
	8,   13,  16,  22,  25,  44,  47,  52,  55,  59,
	62,  67,  74,  81,  82,  87,  91,  94,  103, 104,   
	109, 122, 124, 137, 138, 143, 145, 152, 157, 167, 
	173, 176, 181, 182, 185, 191, 194, 199, 218, 220,
	227, 229, 230, 234, 236, 241, 244, 253,
	//degre 10 totale 60
	4,   13,  19,  22,  50,  55,  64,  69,  98,  107,
	115, 121, 127, 134, 140, 145, 152, 158, 161, 171,
	181, 194, 199, 203, 208, 227, 242, 251, 253, 265, 
	266, 274, 283, 289, 295, 301, 316, 319, 324, 346,
	352, 361, 367, 382, 395, 398, 400, 412, 419, 422,
	426, 428, 433, 446, 454, 457, 472, 493, 505, 508,
	//degre 11	totale 176
	2,   11,   21,    22,    35,   49,   50,   56,   61,   70, 
	74,   79,   84,    88,    103,  104,  112,  115,  117,  122,  
	134,  137,  146,   148,   157,  158,  162,  164,  168,  173,  
	185,  186,  191,   193,   199,  213,  214,  220,  227,  236, 
	242,  251,  256,   259,   265,  266,  276,  292,  304,  310,  
	316,  319,  322,   328,   334,  339,  341,  345,  346,  362,   
	367,  372,  375,   376,   381,  385,  388,  392,  409,  415,
	416,  421,  428,   431,   434,  439,  446,  451,  453,  457, 
	458,  471,  475,   478,   484,  493,  494,  499,  502,  517, 
	518,  524,  527,   555,   560,  565,  569,  578,  580,  587, 
	589,  590,  601,   607,   611,  614,  617,  618,  625,  628, 
	635,  641,  647,   654,   659,  662,  672,  675,  682,  684, 
	689,  695,  696,   713,   719,  724,  733,  734,  740,  747,  
	749,  752,  755,   762,   770,  782,  784,  787,  789,  793, 
	796,  803,  805,   810,   815,  824,  829,  830,  832,  841,  
	847,  849,  861,   871,   878,  889,  892,  901,  908,  920,  
	923,  942,  949,   950,   954,  961,  968,  971,  973,  979,  
	982,  986,  998,  1001,  1010,  1012, 
	//degre 12	totale 144 
	41,   52,   61,   62,   76,   104,  117,  131,  143,  145,  
	157,  167,  171,  176,  181,  194,  217,  236,  239,  262,  
	283,  286,  307,  313,  319,  348,  352,  357,  391,  398, 
	400,  412,  415,  422,  440,  460,  465,  468,  515,  536,  
	539,  551,  558,  563,  570,  595,  598,  617,  647,  654,  
	678,  713,  738,  747,  750,  757,  772,  803,  810,  812,  
	850,  862,  906,  908,  929,  930,  954,  964,  982,  985,  
	991,  992,  1067, 1070, 1096, 1099, 1116, 1143, 1165, 1178, 
	1184, 1202, 1213, 1221, 1240, 1246, 1252, 1255, 1267, 1293, 
	1301, 1305, 1332, 1349, 1384, 1392, 1402, 1413, 1417, 1423, 
	1451, 1480, 1491, 1503, 1504, 1513, 1538, 1544, 1547, 1555, 
	1574, 1603, 1615, 1618, 1629, 1634, 1636, 1639, 1657, 1667, 
	1681, 1697, 1704, 1709, 1722, 1730, 1732, 1802, 1804, 1815, 
	1826, 1832, 1843, 1849, 1863, 1905, 1928, 1933, 1939, 1976, 
	1996, 2013, 2014, 2020 
};

// premi�res valeurs de l'algorithme
static int ivsobol[MAXDIMsobol*MAXBITsobol + 1]={
	0,
	1, 1, 1, 1, 1, 1,
	3, 1, 3, 3, 1, 1,
	5, 7, 7, 3, 3, 5,
	15, 11, 5, 15, 13, 9 
};

static const int vsobol[]={
	0,1,3,5,9,17,35,61,131,281,317
};


static const int iwsobol[MAXBITsobol*MAXDIMsobol + 1]={
	0,
	//Impaires de Sobol
	1,		 
	1,3,
	1,1,7,	  //3
	1,3,7,
	1,1,5,3, //5
	1,3,1,1,

	1,1,3,7, 31, //7
	1,3,3,9, 9,
	1,3,7,13,3, //9
	1,1,5,11,27,
	1,3,5,1, 15,//11
	1,1,7,3, 29,
	1,3,7,7, 21,61, //13
	1,1,1,9, 23,37,
	1,3,3,5, 19,33,
	1,1,3,13,11,7, //16

	// points ajoute' ruotao

	1,3,5,15,13,23,
	1,1,1,11,5, 15,   //18

	//degre 7
	1,3,3,7, 17,9,  35, 
	1,1,5,3, 25,41, 119,
	1,3,1,11,29,5,	9,	   //21
	1,1,7,9, 19,1,	55,

	1,3,3,7, 31,19,	73,   //23
	1,1,3,15,7,	27,	59,
	1,3,5,13,9,	47,	93,   //25
	1,1,1,5, 27,59,	15,   
	1,3,7,1, 21,13,	31,  //27
	1,1,5,9, 13,49,	79,
	1,3,5,11,5,	39,	13,
	1,1,1,1, 23,23,	67,   //30
	1,3,7,15,1,	17,	39, 
	1,1,3,7, 15,11,	25,			  //32
	1,3,3,3, 3,	57,	107,
	1,1,1,13,11,31, 23,				  //34
	1,3,1,5, 25,3,	61,

	1,1,5,3, 17,43,	85,		  //36

	//degre 8
	1,3,3,15,29,25,	113, 107,
	1,1,7,7, 19,55,	67,	 25,		  //38

	1,3,7,11,21,21,	125, 43,
	1,1,5,1, 7,	63,	33,	 205,	  //40
	1,3,3,9, 31,53,	31,	 163,
	1,1,1,5, 17,49,	95,	 233,		  //42
	1,3,1,13,15,29,	49,	 29,

	1,1,7,5, 3,	31,	89,	 97,	  //44
	1,3,5,11,23,35,	53,	 187,
	1,1,3,15,19,45,	21,	 53,		  //46
	1,3,1,7, 9,	51,	111, 101,

	1,1,5,3, 11,13,	57,	 147,		  //48
	1,3,7,11,25,55,	121, 79,
	1,1,1,9, 1,	27,	37,	 151,	  //50
	1,3,1,1, 29,15,	77,	 43,
	// tour 5
	1,1,3,11,13,43,	55,	 211,		 //52
};

class ARM_Sobol : public ARM_QuasiRandom
{
public:
	ARM_Sobol(int firstSimulations);
	ARM_Sobol( const ARM_Sobol& rhs );
	ARM_Sobol& operator=(const ARM_Sobol& rhs );
	virtual ~ARM_Sobol();

	virtual void SetDim( size_t dim );
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual void reset( size_t dim, size_t nbOfPoints );

	/// standard ARM support
	virtual ARM_Object* Clone() const;

	/// validation function
	virtual void ValidateResult( double result ) {};

private:
	virtual void DrawAll();
	void newDrawAll();

	void Init();
	void newInit();
	double RandomNb();

	// key variables for the computation of Sobol sequences
	/// index_ represents the number of draws already done
	unsigned long itsIndex; 
	vector<unsigned long> ix;
	vector<unsigned long> iv;
	int * ixsobol;
	int ** iusobol;
	int insobol;
	double facsobol;

	double itsFactor;
	bool itsInitRand;

};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
