// drinput.h: interface for the drinput class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_KINPUT_H__2F9722F7_ADFE_11D1_80A4_00C04FB91C08__INCLUDED_)
/**@#-*/
#define AFX_KINPUT_H__2F9722F7_ADFE_11D1_80A4_00C04FB91C08__INCLUDED_
/**@#+*/

#include "kstring.h"
#include "kptr.h"

#include FSTREAM_H
#include SSTREAM_H
#include <vector>
#include <map>

class KInput;
class KInputMap;
class KInputVector;


//k This class represents my generic input/output structure.  
//k See web page for details.

/** Generic input strucuture.

  Provides a convenient way to add organization to large inputs and/or
  provide support for polymorphic input (ie the input format depends upon values of the input).

  A detailed description of the syntax can be found 
  <a href = http://dr.ny.jpmorgan/~mbs/docs/mbsio.html>here</a>.

  @see KInputVector
  @see KInputMap
*/

class KInput  
{
public:
 	typedef upper_string String;
	typedef KPtr<String> StringPtr;
	typedef KPtr<KInputVector> KInputVectorPtr;
	typedef KPtr<KInputMap> KInputMapPtr;
	
	/** Creates an empty KInput */
	KInput(); //

	/** Construct from an istream */
	KInput (std::istream& s) {Parse(s);}
	/** Construct from a KString */
	KInput(const String&); // 
	/** Construct from a KInputVector */
	KInput(const KInputVector&); // 
	/** Construct from a KInputMap */
	KInput(const KInputMap&); // 

	/** Opens a file and reads in the KInput

		Current contents of KInput are overwritten.
		*/
	void open (KString filename); // Create from file
	
	/** Gets KInputVector from KInput
		@returns KInputVector
		@exception Fails if KInput is not a KInputVector
		*/
	KInputVector& get_vector(); // Converts to KInputVector
	/** Gets KInputMap from KInput
		@returns KInputMap
		@exception Fails if KInput is not a KInputMap
		*/
	KInputMap& get_map(); // Converts to KInputMap 
	/** Gets String from KInput
		@returns String
		@exception Fails if KInput is not a String
		*/
	String& get_string(); // Converts to String
	
	/** Checks if KInput is a vector */
	bool is_vector(); // 
	/** Checks if KInput is a map */
	bool is_map(); //
	/** Checks if KInput is a KString */
	bool is_string(); //
	/** Checks if KInput is empty */
	bool is_empty(); //
	
	/** Checks equality */
	bool operator==(const KInput&) const;
	
	friend std::ostream& operator<<(std::ostream&, const KInput&); //
	friend std::istream& operator>>(std::istream&, KInput&);

	/**@#-*/
	void print(std::ostream&, int level = 0) const ; // Level is the number of tabs
	/**@#+*/
	
protected:
	enum KInputType {EMPTY, STRING, MAP, VECTOR};
	KInputType m_type;
	
	StringPtr m_stringPtr;
	KInputVectorPtr m_vectorPtr;
	KInputMapPtr m_mapPtr;
	
	KInput Parse(std::istream&);
	void ParseFile(KString fileName);
	void ParseVector (std::istream& inFile);
	void ParseMap (std::istream& inFile);
	void ParseString (std::istream& inFile);
	void ParseQuote (std::istream& inFile);
	
	char PeekNextChar (std::istream& inFile);
	upper_string GetNextString (std::istream& inFile);
};

/** Constructs KInput from a generic type.

	Directs variable to an output stream and attempts to construct 
	input from the stream. 
	@params variable Variable
	@params input KInput to be created
*/

template <class T> inline
void operator>>(const T& variable, KInput& input)
{	
	OSTRINGSTREAM out;
	out << variable;
	KString s = out.str();

#ifdef __SUNPRO_CC
// Solaris uses fucked up conventions in its strstream class.
// It doesn't have a KString class and nothing makes any fucking sense.
	KString temp;
	for (int i = 0; i < out.pcount(); i++) 
		temp += s.at(i);
	s = temp;
#endif

	ISTRINGSTREAM in (s.c_str());
	input = KInput(in);
}

/** Constructs a generic type from a KInput

	Directs input to an output stream and attempts to construct 
	variable from the stream. 
	@params variable Variable to initialize
	@params input KInput
*/
template <class T> inline
void operator>>(const KInput& input, T& variable)
{	
	OSTRINGSTREAM out;
	out << input;
	KString s = out.str();
	ISTRINGSTREAM in (s.c_str());
	in >> variable;
}

inline bool KInput::is_vector() {return (m_type == VECTOR);}
inline bool KInput::is_map() {return (m_type == MAP);}
inline bool KInput::is_string() {return (m_type == STRING);}
inline bool KInput::is_empty() {return (m_type == EMPTY);}

/** Vector of KInputs

  Represents a vector of KInputs (a type of KInput).
  @see KInput
  */
class KInputVector: public KVector(KInput) {
public:
	/** Creates an empty KInputVector */
	KInputVector () {}
	
#ifdef K_MEMBER_TEMPLATES
	template <class T> 
		KInputVector (const KVector(T)& v) {
		for (KVector(T)::const_iterator iter = v.begin(); iter != v.end(); iter++) {
			KInput input;
			(*iter) >> input;
			push_back(input);
		}
	}
	
	template <class T>
		void  toVector(KVector(T)& v) {
		for (const_iterator iter = begin(); iter != end(); iter++) {
			v.clear();
			T t;
			(*iter) >> t;
			v.push_back(t);
		}
	}
#endif
	
	friend std::ostream& operator<<(std::ostream&, const KInputVector&); //
/*@#-*/
	void print(std::ostream&, int level = 0) const ; // Level is the number of tabs
/*@#+*/
};

/** Map of (KInput::String, KInput) pairs

  Represents a map from KInput::String to KInput (a type of KInput).
  @see KInput
  */

class KInputMap : public KMap(KInput::String, KInput) {
public:
	typedef KInput::String String;

	typedef KMap(KInput::String, KInput) super;

	/** Creates an empty KInputMap */
	KInputMap () {}
	/** Gets a KString with default value */
	upper_string get (upper_string name, upper_string defaultValue); // Gets a KString; if doesn't exist, returns default
	/** Gets a KString without default value */
	upper_string get (upper_string name); // Gets a KString; if it doesn't exist, error
	bool is_key (upper_string name){return isKey(name);}
	/** Checks whether key exists */
	bool isKey (upper_string name); // Checks if the key exists
	/** Gets a KInput */
	KInput& get_KInput(upper_string name); // Gets an KInput; if doens't exist, return an error
	
	/** Stores a (String, String) pair */
	KInputMap& store (upper_string name, upper_string value); // Opposite of get

	/** Same as map's operator[], except argument is by value. */
	KInput& operator[](upper_string);
	
	friend std::ostream& operator<<(std::ostream&, const KInputMap&); //
/*@#-*/
	void print(std::ostream&, int level = 0) const ; // Level is the number of tabs
/*@#+*/
	
#ifdef K_MEMBER_TEMPLATES
	template <class T> 
		KInputMap (const KMap(upper_string, T)& m) {
		for (KMap(upper_string, T)::const_iterator iter = m.begin(); iter != m.end(); iter++) {
			KInput input;
			(*iter).second >> input;
			operator[](*iter.first) = push_back(input);
		}
	}
	
	template <class T>
		void  toMap(KMap(upper_string, T)& m) {
		for (const_iterator iter = begin(); iter != end(); iter++) {
			m.clear();
			T t;
			(*iter).second >> t;
			m[(*iter).first] = t;
		}
	}
#endif 
	
};


#endif // !defined(AFX_KINPUT_H__2F9722F7_ADFE_11D1_80A4_00C04FB91C08__INCLUDED_)
