// drparrates.cpp: implementation of the drparrates class.
//
//////////////////////////////////////////////////////////////////////

#include "drparrates.h"
#include "drsymboltable.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

const int NUM_PARRATE_PTS = 25;

DRParRates::DRParRates(DRDate& startDate, DRString cmtspdsFile, 
					   DRString yldcrvlcFile)
					   : m_swapParRates(NUM_PARRATE_PTS),
					   m_humpRates(NUM_PARRATE_PTS),
					   m_liborSpreads(NUM_PARRATE_PTS)
{
	LoadCmtspds (cmtspdsFile);
	LoadYldcrvlc (yldcrvlcFile, startDate);
}


void DRParRates::TweakParRate(Tweak tweakPoint, double bpTweak)
{
	bpTweak /= 10000;

	DArray changes(MyTweaks[tweakPoint], m_swapParRates.size());

	changes *= bpTweak / 5.;

	m_swapParRates += changes;
}

void DRParRates::LoadCmtspds (DRString filename)
{
	char temp[20];

	ifstream inFile(CheckPathDelimiters(theSymbolTable.get(filename)).c_str());
	if (!inFile) 
		throw DRException ("Can't open hump file ") << filename;

	inFile >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp;
	int i;
	for (i = 0; i<NUM_PARRATE_PTS && !inFile.eof(); i++) {
		inFile >> temp >> temp >> temp >> temp >> temp >> temp;
		inFile >> m_humpRates[i] >> temp;
	}

	m_humpRates /= 10000;

	if (i!=NUM_PARRATE_PTS) 
		throw DRException("Found end of file prematurely in LoadCmtspds");
}

void DRParRates::LoadYldcrvlc (DRString filename, TDate startDate)
{
	LArray tempDates(NUM_PARRATE_PTS);
	char temp[20];
	double junk;
	DRString input;

	ifstream inFile(CheckPathDelimiters(theSymbolTable.get(filename)).c_str());
	if (!inFile) 
		throw DRException ("Failed to open Yldcrvlc file: ") << filename;

	int i;
	for (i = 0; i < NUM_PARRATE_PTS && !inFile.eof(); i++) {

		do {
			inFile >> temp;
			input = temp;
		} while (!inFile.eof() && input.size() == 0);
		
		input = left_substring(input, -1);
		input = right_substring(input, -1);

		if (input.at(0) == 'O') {
			tempDates[i] = startDate + 1;
		}
		else {
			TDateInterval dateInterval;
			char* inp = const_cast<char*> (input.c_str());
			GtoStringToDateInterval (inp, "", &dateInterval);
			GtoDtFwdAny (startDate, &dateInterval, &tempDates[i]);
		}
	
		inFile >> temp;
		input = temp;
		m_swapParRates[i] = atof (left_substring(input, -1).c_str());

		inFile >> junk;
		inFile >> junk;
		inFile >> m_liborSpreads[i];
	}

	m_swapParRates /= 100;
	m_liborSpreads /= 10000;

	if (i!=NUM_PARRATE_PTS) 
		throw DRException("Found end of file prematurely in ReadYldcrvlc");

	m_rateDates = DRDateList (tempDates, NUM_PARRATE_PTS);;
}

ostream& operator<<(ostream& s, const DRParRates& a) 
{
	s << "Dates:\t" << a.m_rateDates;
	s << "ParRates:\t" << a.m_swapParRates;
	s << "HumpRates:\t" << a.m_humpRates;
	s << "LiborSpreads:\t" << a.m_liborSpreads;
	return s;
}
