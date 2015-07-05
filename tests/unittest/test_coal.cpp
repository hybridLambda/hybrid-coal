/*
 * hybrid-coal is used to compute gene tree probabilities given species network under coalescent process.
 *
 * Copyright (C) 2010 -- 2015 Sha (Joe) Zhu
 *
 * This file is part of hybrid-coal
 *
 * hybrid-coal is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "coal.hpp"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestCoal : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestCoal );
  CPPUNIT_TEST ( testConstructor );
  CPPUNIT_TEST ( testJames1 );
  CPPUNIT_TEST ( testJames2 );
  CPPUNIT_TEST ( testJames3 );
  CPPUNIT_TEST ( testJames4 );
  CPPUNIT_TEST ( testJames5 );
  CPPUNIT_TEST_SUITE_END();
 private:
  string gtstr;
 public:
  void setUp() {
      gtstr = "((A:1,B:1):1,(C:1,((D:1,F:1):1,(E:1,G:1):1):1):1);";
      return;
    }

    void tearDown() {
      return;
    }

    void testConstructor(){
      string spstr = "(((B:0.2,C:0.2):0.2,A:0.2):0.2,((D:0.2,E:0.2):0.2,(F:0.2,G:0.2):0.2):0.2);";
      CPPUNIT_ASSERT_NO_THROW( CoalST sp( spstr ) );
      CPPUNIT_ASSERT_NO_THROW( CoalGT gt( gtstr ) );
    }

    void testJames1 (){
      string spstr1 = "(((B:1,C:1):1,A:1):1,((D:1,E:1):1,(F:1,G:1):1):1);";
      CoalST sp1( spstr1 );
      CoalGT gt1( gtstr );
      CPPUNIT_ASSERT_NO_THROW( gt1.prob_given_sp_tree( sp1 ) ) ;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0001544740, gt1.probability, 0.00000001 );
    }

    void testJames2 (){
      string spstr2 = "(((B:0.5,C:0.5):0.5,A:0.5):0.5,((D:0.5,E:0.5):0.5,(F:0.5,G:0.5):0.5):0.5);";
      CoalST sp2( spstr2 );
      CoalGT gt2( gtstr );
      CPPUNIT_ASSERT_NO_THROW( gt2.prob_given_sp_tree( sp2 ) ) ;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0004767273, gt2.probability, 0.00000001 );
    }

    void testJames3 (){
      string spstr3 = "(((B:0.2,C:0.2):0.2,A:0.2):0.2,((D:0.2,E:0.2):0.2,(F:0.2,G:0.2):0.2):0.2);";
      CoalST sp3( spstr3 );
      CoalGT gt3( gtstr );
      CPPUNIT_ASSERT_NO_THROW( gt3.prob_given_sp_tree( sp3 ) ) ;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000461655, gt3.probability, 0.00000001 );
    }

    void testJames4 (){
      string spstr4 = "(((B:1,C:1):.01,A:1):1,((D:1,E:1):1,(F:1,G:1):1):1);";
      CoalST sp4( spstr4 );
      CoalGT gt4( gtstr );
      CPPUNIT_ASSERT_NO_THROW( gt4.prob_given_sp_tree( sp4 ) ) ;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0004157257, gt4.probability, 0.00000001 );
    }

    void testJames5 (){
      string spstr5 = "(((B:1,C:1):1,A:1):1,((D:1,E:1):1,(F:1,G:1):1):0.01);";
      CoalST sp5( spstr5 );
      CoalGT gt5( gtstr );
      CPPUNIT_ASSERT_NO_THROW( gt5.prob_given_sp_tree( sp5 ) ) ;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0000068945, gt5.probability, 0.00000001 );
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestCoal );
