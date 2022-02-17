// DRCmdline.h: interface for the DRCmdline class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DRCMDLINE_H__03C3B0B7_6506_11D2_97E0_00C04FD8EB9A__INCLUDED_)
#define AFX_DRCMDLINE_H__03C3B0B7_6506_11D2_97E0_00C04FD8EB9A__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include "drexception.h"
#include "drstring.h"
#include <vector>
#include <map>

class DRCmdLine  
{
public:
	typedef map<string, string, less<string>, MYALLOC(string)> SwitchMap;
	typedef vector<string, MYALLOC(string)> InputVector;

	DRCmdLine(int argc, char** argv);
	virtual ~DRCmdLine() {}

	const SwitchMap& switches () const;
	const InputVector& inputs () const;
private:
	SwitchMap m_switchMap;
	InputVector m_inputVector;
};

inline
const DRCmdLine::SwitchMap& DRCmdLine::switches () const
{return m_switchMap;}

inline
const DRCmdLine::InputVector& DRCmdLine::inputs () const
{return m_inputVector;}

inline
DRCmdLine::DRCmdLine(int argc, char** argv)
{
	for (int i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			string switchName = argv[i] + 1;
			if (i + 1 >= argc)
				throw DRException ("Command line error: unmatched switch: ") << switchName;
			string switchValue = argv[++i];
			m_switchMap[switchName] = switchValue;
		}
		else {
			string input = argv[i];
			m_inputVector.push_back(input);
		}
	}
}

#endif // !defined(AFX_DRCMDLINE_H__03C3B0B7_6506_11D2_97E0_00C04FD8EB9A__INCLUDED_)






