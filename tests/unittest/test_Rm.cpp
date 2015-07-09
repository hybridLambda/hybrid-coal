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

class TestRm : public CppUnit::TestCase {
  CPPUNIT_TEST_SUITE( TestRm );
  CPPUNIT_TEST ( testChooseRmNode1 );
  CPPUNIT_TEST ( testChooseRmNode2 );
  CPPUNIT_TEST ( testExtract_hybrid_para );
  CPPUNIT_TEST ( testRemoveHnodeOneChild );
  CPPUNIT_TEST ( testRemoveHnodeManyChild );
  CPPUNIT_TEST ( testRemoveSnode1 );
  CPPUNIT_TEST ( testRemoveSnode2 );
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

    void testChooseRmNode1(){
      string spstr = "((((B:1,C:1)s1:1)h1#.5:1,A:3)s2:1,(h1#.5:1,D:3)s3:1)r;";
      CPPUNIT_ASSERT_NO_THROW( TmpSN sp( spstr ) );
      TmpSN sp(spstr);
      CPPUNIT_ASSERT_EQUAL (2, sp.toBeRemovedNodeIndex());
    }

    void testChooseRmNode2(){
      string spstr = "(((((B:2.1,C:2.1)s_1:.4)h_1#H1:0.2,D:2.7)s_2:.6,(h_1#H1:.7)h_2#H2:.1)s_3:0.3,(A:3.4,h_2#H2:0.2)s_4:0.2)r;";
      CPPUNIT_ASSERT_NO_THROW( TmpSN sp( spstr ) );
      TmpSN sp(spstr);
      CPPUNIT_ASSERT_EQUAL (2, sp.toBeRemovedNodeIndex());
    }

    void testExtract_hybrid_para (){
        string nodeName = "H13#0.3";
        CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.3, extract_hybrid_para(nodeName), 0.00000001);
    }

    void testRemoveHnodeOneChild(){
        string netStr = "(((A:1)H1#0.3:1,C:1)Int1:1,(H1#0.3:1,D:1)Int2:1)r";
        CoalSN spNet( netStr );
        TmpSN tmpSN ( netStr );
        int rm_node_index = tmpSN.toBeRemovedNodeIndex();
        CPPUNIT_ASSERT_EQUAL (1, rm_node_index);
        CPPUNIT_ASSERT_NO_THROW( spNet.removeHnode ( tmpSN, rm_node_index, spNet.NetStrWizPriorList[0], true, true ) );
        CPPUNIT_ASSERT_NO_THROW( spNet.NetStrWizPriorList.erase(spNet.NetStrWizPriorList.begin() + spNet.currentSubNetworkIndex));
        CPPUNIT_ASSERT_EQUAL ( (size_t)2, spNet.NetStrWizPriorList.size() );
        CPPUNIT_ASSERT_EQUAL ( string("((C:1.000000,A:2.000000)Int1:1.000000,D:2.000000)r;"),
                               spNet.NetStrWizPriorList[0].netStr);
        CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.3, spNet.NetStrWizPriorList[0].prior.omega(), 0.00000001 );
        CPPUNIT_ASSERT_EQUAL ( string("((D:1.000000,A:2.000000)Int2:1.000000,C:2.000000)r;"),
                               spNet.NetStrWizPriorList[1].netStr);
        CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.7, spNet.NetStrWizPriorList[1].prior.omega(), 0.00000001 );
    }

    void testRemoveHnodeManyChild(){
        string netStr = "(((A:1,B:1)H1#0.3:1,C:1)Int1:1,(H1#0.3:1,D:1)Int2:1)r";
        CoalSN spNet( netStr );
        TmpSN tmpSN ( netStr );
        int rm_node_index = tmpSN.toBeRemovedNodeIndex();
        CPPUNIT_ASSERT_EQUAL (2, rm_node_index);
        CPPUNIT_ASSERT_NO_THROW( spNet.removeHnode ( tmpSN, rm_node_index, spNet.NetStrWizPriorList[0], true, true) );
        CPPUNIT_ASSERT_NO_THROW( spNet.NetStrWizPriorList.erase(spNet.NetStrWizPriorList.begin() + spNet.currentSubNetworkIndex));
        CPPUNIT_ASSERT_EQUAL ( (size_t)4, spNet.NetStrWizPriorList.size() );
        CPPUNIT_ASSERT_EQUAL ( string("((C:1.000000,(A:1.000000,B:1.000000)H1L:1.000000)Int1:1.000000,D:2.000000)r;"),
                               spNet.NetStrWizPriorList[0].netStr);
        CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.09, spNet.NetStrWizPriorList[0].prior.omega(), 0.00000001 );
        CPPUNIT_ASSERT_EQUAL ( string("((D:1.000000,(A:1.000000,B:1.000000)H1R:1.000000)Int2:1.000000,C:2.000000)r;"),
                               spNet.NetStrWizPriorList[1].netStr);
        CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.49, spNet.NetStrWizPriorList[1].prior.omega(), 0.00000001 );
        CPPUNIT_ASSERT_EQUAL ( string("((C:1.000000,A:2.000000)Int1:1.000000,(D:1.000000,B:2.000000)Int2:1.000000)r;"),
                               spNet.NetStrWizPriorList[2].netStr);
        CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.21, spNet.NetStrWizPriorList[2].prior.omega(), 0.00000001 );
        CPPUNIT_ASSERT_EQUAL ( string("((C:1.000000,B:2.000000)Int1:1.000000,(D:1.000000,A:2.000000)Int2:1.000000)r;"),
                               spNet.NetStrWizPriorList[3].netStr);
        CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.21, spNet.NetStrWizPriorList[3].prior.omega(), 0.00000001 );

    }

    void testRemoveSnode1(){
        string gtStr = "";
        string netStr = "((((B:1,C:1)s1:1)h1#.3:1,A:3)s2:1,(h1#.3:13,D:3)s3:1)r;";
        CoalSN spNet (netStr);
        TmpSN tmpSN ( netStr );
        int rm_node_index = tmpSN.toBeRemovedNodeIndex();
        CPPUNIT_ASSERT_EQUAL (2, rm_node_index);
        CPPUNIT_ASSERT_NO_THROW( spNet.removeSnode ( gtStr, tmpSN, rm_node_index, spNet.NetStrWizPriorList[0], true, true) );
        CPPUNIT_ASSERT_NO_THROW( spNet.NetStrWizPriorList.erase(spNet.NetStrWizPriorList.begin() + spNet.currentSubNetworkIndex));
        CPPUNIT_ASSERT_EQUAL ( (size_t)2, spNet.NetStrWizPriorList.size() );
    }


    void testRemoveSnode2(){
        string gtStr = "";
        string netStr = "((((B:1,C:1,E:1)s1:1)h1#.3:1,A:3)s2:1,(h1#.3:13,D:3)s3:1)r;";
        CoalSN spNet (netStr);
        TmpSN tmpSN ( netStr );
        int rm_node_index = tmpSN.toBeRemovedNodeIndex();
        CPPUNIT_ASSERT_EQUAL (3, rm_node_index);
        CPPUNIT_ASSERT_NO_THROW( spNet.removeSnode ( gtStr, tmpSN, rm_node_index, spNet.NetStrWizPriorList[0], true, true) );
        CPPUNIT_ASSERT_NO_THROW( spNet.NetStrWizPriorList.erase(spNet.NetStrWizPriorList.begin() + spNet.currentSubNetworkIndex));
        CPPUNIT_ASSERT_EQUAL ( (size_t)5, spNet.NetStrWizPriorList.size() );
    }

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestRm );
