/////////////////////////////////////////////////////////////////////////////
// Name:        tests/factory/main.cpp
// Purpose:     Factory unit test
// Author:      Vadim Zeitlin
// Created:     2005-07-25
// RCS-ID:      $Id: main.cpp,v 1.2 2005/07/27 16:35:03 zeitlin Exp $
// Copyright:   (c) 2003-2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/factory.h"
#include "ito33/cppunit.h"
#include "ito33/string.h"

using std::string;
using ito33::Factory;

static const string KEY_REGISTERED = "foo";
static const string KEY_REGISTERED_SPECIAL = "bar";
static const string KEY_NOT_REGISTERED = "moo";

static const string NAME_SPECIAL = "Special";

struct Product
{
  Product(const string& s = string()) : m_str(s) { }

  string m_str;
};

typedef Factory<string, Product> StringVoidFactory;
typedef Factory<string, Product, string> StringStringFactory;

// use a macro to predefine a product
ITO33_DEFINE_FACTORY_PRODUCT(string, KEY_REGISTERED, Product, Product);

// do it manually
namespace
{
  Product *CreateSpecial()
  {
    return new Product(NAME_SPECIAL);
  }

  StringVoidFactory factorySV(KEY_REGISTERED_SPECIAL, CreateSpecial);


  // for factories with arguments we must always implement creatores manually
  Product *CreateNamed(const string *arg)
  {
    return new Product(*arg);
  }

  StringStringFactory factorySS(KEY_REGISTERED, CreateNamed);
}


// all factories we use must be implemented:

// in the simplest case, this macro taking just key and value types is enough
ITO33_IMPLEMENT_FACTORY(string, Product);

// otherwise a more generic one taking the factory name must be used
ITO33_IMPLEMENT_THE_FACTORY(StringStringFactory);


class FactoryTestCase : public CppUnit::TestCase
{
public:
  FactoryTestCase() { }

private:
  CPPUNIT_TEST_SUITE( FactoryTestCase );
    CPPUNIT_TEST( StringVoid );
    CPPUNIT_TEST( StringString );
    CPPUNIT_TEST( Enumerate );
  CPPUNIT_TEST_SUITE_END();

  void StringVoid();
  void StringString();
  void Enumerate();

  NO_COPY_CLASS(FactoryTestCase);
};

void FactoryTestCase::StringVoid()
{
  Product *product = StringVoidFactory::Create(KEY_REGISTERED);
  CPPUNIT_ASSERT( product );
  CPPUNIT_ASSERT( product->m_str.empty() );

  delete product;

  product = StringVoidFactory::Create(KEY_REGISTERED_SPECIAL);
  CPPUNIT_ASSERT( product );
  CPPUNIT_ASSERT( product->m_str == NAME_SPECIAL );
  delete product;

  product = StringVoidFactory::Create(KEY_NOT_REGISTERED);
  CPPUNIT_ASSERT( !product );
}

void FactoryTestCase::StringString()
{
  static const string NAME_FIXED = "Usual";

  Product *product = StringStringFactory::Create(KEY_REGISTERED, &NAME_FIXED);
  CPPUNIT_ASSERT( product );
  CPPUNIT_ASSERT( product->m_str == NAME_FIXED );
  delete product;

  product = StringStringFactory::Create(KEY_NOT_REGISTERED, &NAME_FIXED);
  CPPUNIT_ASSERT( !product );
}

void FactoryTestCase::Enumerate()
{
  // we register only 2 creators so we should have 2 factories
  StringVoidFactory *factory = StringVoidFactory::GetFirst();
  CPPUNIT_ASSERT( factory );

  factory = factory->GetNext();
  CPPUNIT_ASSERT( factory );

  factory = factory->GetNext();
  CPPUNIT_ASSERT( !factory );
}

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(FactoryTestCase::suite());

  return runner.run("") ? 0 : 1;
}

