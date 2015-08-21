#!/bin/bash
function test_hybrid-coal {
  echo -n " hybrid-coal $@ "
  echo -n "."

  if [ -f hybrid-coal_dbg ]; then
    # Test using hybrid-coal self-checks
    ./hybrid-coal_dbg $@ > /dev/null
    if [ $? -ne 0 ]; then
      echo ""
      echo "Executing \"./hybrid-coal_dbg $@ \" failed."
      echo "Debug Call: make -mj2 hybrid-coal_dbg && ./hybrid-coal_dbg $@ 2>&1 | less"
      exit 1
    fi
  fi

  ## Test for memory leaks
  valgrind --error-exitcode=1 --leak-check=full -q ./hybrid-coal $@ > /dev/null
  if [ $? -ne 0 ]; then
    echo ""
    echo "Valgrind check of \"./hybrid-coal $@ \" failed."
    exit 1
  fi

  ## Test for distribution
  ./hybrid-coal $@ 2>&1 | grep "Total probability = 1"
  if [ $? -ne 0 ]; then
    echo ""
    echo "Total probability check of \"./hybrid-coal $@ \" failed."
    exit 1
  fi
  echo " done."
}

echo "Testing Examples"
  test_hybrid-coal -sp trees/6_tax_multi_tax_test_hybrid1_topo1 || exit 1
  test_hybrid-coal -sp trees/6_tax_multi_tax_test_hybrid1_topo2 || exit 1
  test_hybrid-coal -sp trees/6_tax_multi_tax_test_hybrid2_topo1 || exit 1
  test_hybrid-coal -sp trees/6_tax_multi_tax_test_hybrid2_topo2 || exit 1
  test_hybrid-coal -sp trees/6_tax_multi_tax_test_hybrid3_topo1 || exit 1
  test_hybrid-coal -sp trees/6_tax_multi_tax_test_hybrid3_topo2 || exit 1
  test_hybrid-coal -sp trees/6_tax_multi_tax_test_hybrid4_topo1 || exit 1
  test_hybrid-coal -sp trees/6_tax_multi_tax_test_hybrid4_topo2 || exit 1
  test_hybrid-coal -sp trees/6_tax_non_binary_sp.tre || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_level2_para || exit 1

  test_hybrid-coal -sp trees/6_tax_sp_nt1_para || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt1_para00_bl06 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt1_para02_bl06 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt1_para05_bl06 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt1_para08_bl06 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt1_para10_bl06 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt1_para_bl0 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt1_para_bl100 || exit 1

  test_hybrid-coal -sp trees/6_tax_sp_nt2_para || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt2_para00_bl06 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt2_para02_bl06 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt2_para05_bl06 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt2_para08_bl06 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt2_para10_bl06 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt2_para2 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt2_para_bl0 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt2_para_bl100 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt3_para || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt3_para00_bl06 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt3_para02_bl06 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt3_para05_bl06 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt3_para08_bl06 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt3_para10_bl06 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt3_para_bl0 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt3_para_bl100 || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt4_para || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt5_para || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt6_para || exit 1
  test_hybrid-coal -sp trees/6_tax_sp_nt7_para || exit 1

########################## Checked
  test_hybrid-coal -sp '((((B:1,C:1)s1:1)h1#.5:1,A:3)s2:1,(h1#.5:1,D:3)s3:1)r;' || exit 1
  test_hybrid-coal -sp trees/3_tax_sp_nt1 || exit 1
  test_hybrid-coal -sp trees/3_tax_sp_nt2 || exit 1
  test_hybrid-coal -sp trees/3_tax_sp_nt3 || exit 1
  test_hybrid-coal -sp trees/3_tax_sp_nt3_bl100 || exit 1
  test_hybrid-coal -sp trees/3_tax_sp_nt4 || exit 1
  test_hybrid-coal -sp trees/3_tax_sp_nt5 || exit 1

  test_hybrid-coal -sp trees/4_tax_sp_nt1 || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt1_para || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt1_para2 || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt1_para3 || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt1_para_diff_bl || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt1_para_flat || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt1_para_skew || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt2 || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt2_para || exit 1

  test_hybrid-coal -sp trees/4_tax_sp_nt3 || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt3_para || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt4 || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt4_para || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt5 || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt5_para || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt6 || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt6_para || exit 1

  test_hybrid-coal -sp trees/5_tax_sp_nt2_para_bl0 || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt2 || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt2_para_bl100 || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt2_para || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt2_para00_bl06 || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt2_para02_bl06 || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt2_para05_bl06 || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt2_para08_bl06 || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt2_para10_bl06 || exit 1

  test_hybrid-coal -sp trees/5_tax_sp_nt3_para_bl0 || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt3_para_bl100 || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt3_para00_bl06 || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt3_para02_bl06 || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt3_para05_bl06 || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt3_para08_bl06 || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt3_para || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt3_para10_bl06 || exit 1


  test_hybrid-coal -sp trees/5_tax_sp_nt1_para || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt1_para00_bl06 || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt1_para_bl0 || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt1_para_bl100 || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt1_para02_bl06 || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt1_para05_bl06 || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt1_para08_bl06 || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt1_para10_bl06 || exit 1

  test_hybrid-coal -sp trees/5_tax_sp_nt4_para || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt5_para || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt6_para || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt7_para || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt8_para || exit 1
  test_hybrid-coal -sp trees/5_tax_sp_nt9_para || exit 1


# special case ...
  #test_hybrid-coal -sp trees/5_tax_sp_nt1 || exit 1
  #test_hybrid-coal -sp trees/5_tax_sp_nt3 || exit 1
  #test_hybrid-coal -sp trees/5_tax_sp_nt6 || exit 1
  #test_hybrid-coal -sp trees/5_tax_sp_nt5 || exit 1
  #test_hybrid-coal -sp trees/5_tax_sp_nt4 || exit 1

echo ""
