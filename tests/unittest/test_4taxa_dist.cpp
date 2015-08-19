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

class TestRmTrees4taxa : public CppUnit::TestCase {
  CPPUNIT_TEST_SUITE( TestRmTrees4taxa );
  CPPUNIT_TEST ( testTree1 );
  //CPPUNIT_TEST ( testRemoveSnode1case2 );
  //CPPUNIT_TEST ( testRemoveSnode1case3 );
  //CPPUNIT_TEST ( testRemoveSnode2case1 );
  //CPPUNIT_TEST ( testRemoveSnode2case2 );
  //CPPUNIT_TEST ( testRemoveSnode2case3 );
  //CPPUNIT_TEST ( testRemoveSnode2case4 );
  CPPUNIT_TEST_SUITE_END();
 private:
  string gtstr;
 public:
  void setUp() {
      return;
    }

    void tearDown() {
      return;
    }

    void testTree1(){
        string netStr = "((((B:0,C:0)s1:0)h1#.3:0,A:0)s2:0,(h1#.3:0,D:0)s3:0)r;";
        //./hybrid-coal -sp "(((A:0,B:0):0,C:0):0,D:0);"
        CoalSN spNet (netStr);

        // List of trees, refer to thesis page 96
        gtstr = "(((A,D),C),B);";
        CPPUNIT_ASSERT_NO_THROW( spNet.simplifyNetworks(gtstr) );
        //CPPUNIT_ASSERT_NO_THROW( spNet.DEBUGprintSubNetwork() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)4, spNet.NetStrWizPriorList.size() );

        gtstr = "((A,(C,D)),B);";
        CPPUNIT_ASSERT_NO_THROW( spNet.simplifyNetworks(gtstr) );
        CPPUNIT_ASSERT_EQUAL ( (size_t)4, spNet.NetStrWizPriorList.size() );

        gtstr = "(((A,C),D),B);";
        CPPUNIT_ASSERT_NO_THROW( spNet.simplifyNetworks(gtstr) );
        CPPUNIT_ASSERT_EQUAL ( (size_t)4, spNet.NetStrWizPriorList.size() );

        gtstr = "((A,C),(B,D));";
        CPPUNIT_ASSERT_NO_THROW( spNet.simplifyNetworks(gtstr) );
        CPPUNIT_ASSERT_EQUAL ( (size_t)4, spNet.NetStrWizPriorList.size() );

        gtstr = "(((A,C),B),D);";
        CPPUNIT_ASSERT_NO_THROW( spNet.simplifyNetworks(gtstr) );
        CPPUNIT_ASSERT_EQUAL ( (size_t)4, spNet.NetStrWizPriorList.size() );

        gtstr = "((A,D),(B,C));";
        CPPUNIT_ASSERT_NO_THROW( spNet.simplifyNetworks(gtstr) );
        CPPUNIT_ASSERT_NO_THROW( spNet.DEBUGprintSubNetwork() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)6, spNet.NetStrWizPriorList.size() );

        gtstr = "(A,((B,D),C));";
        CPPUNIT_ASSERT_NO_THROW( spNet.simplifyNetworks(gtstr) );
        CPPUNIT_ASSERT_EQUAL ( (size_t)4, spNet.NetStrWizPriorList.size() );

        gtstr = "(A,(B,(C,D)));";
        CPPUNIT_ASSERT_NO_THROW( spNet.simplifyNetworks(gtstr) );
        CPPUNIT_ASSERT_EQUAL ( (size_t)4, spNet.NetStrWizPriorList.size() );

        gtstr = "(A,((B,C),D));";
        CPPUNIT_ASSERT_NO_THROW( spNet.simplifyNetworks(gtstr) );
        CPPUNIT_ASSERT_EQUAL ( (size_t)6, spNet.NetStrWizPriorList.size() );

        gtstr = "((A,(B,C)),D);";
        CPPUNIT_ASSERT_NO_THROW( spNet.simplifyNetworks(gtstr) );
        CPPUNIT_ASSERT_EQUAL ( (size_t)6, spNet.NetStrWizPriorList.size() );

        gtstr = "(((A,D),B),C);";
        CPPUNIT_ASSERT_NO_THROW( spNet.simplifyNetworks(gtstr) );
        CPPUNIT_ASSERT_EQUAL ( (size_t)4, spNet.NetStrWizPriorList.size() );

        gtstr = "((A,(B,D)),C);";
        CPPUNIT_ASSERT_NO_THROW( spNet.simplifyNetworks(gtstr) );
        CPPUNIT_ASSERT_EQUAL ( (size_t)4, spNet.NetStrWizPriorList.size() );

        gtstr = "(((A,B),D),C);";
        CPPUNIT_ASSERT_NO_THROW( spNet.simplifyNetworks(gtstr) );
        CPPUNIT_ASSERT_EQUAL ( (size_t)4, spNet.NetStrWizPriorList.size() );

        gtstr = "((A,B),(C,D));";
        CPPUNIT_ASSERT_NO_THROW( spNet.simplifyNetworks(gtstr) );
        CPPUNIT_ASSERT_EQUAL ( (size_t)4, spNet.NetStrWizPriorList.size() );

        gtstr = "(((A,B),C),D);";
        CPPUNIT_ASSERT_NO_THROW( spNet.simplifyNetworks(gtstr) );
        CPPUNIT_ASSERT_EQUAL ( (size_t)4, spNet.NetStrWizPriorList.size() );
    }


};

CPPUNIT_TEST_SUITE_REGISTRATION( TestRmTrees4taxa );
