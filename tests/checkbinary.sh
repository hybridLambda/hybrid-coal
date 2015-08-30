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

  echo " done."
}

echo "Testing Examples"
  # OK
  test_hybrid-coal -sp trees/7_tax_sp_james02.tre -gt trees/7_tax_gt_james || exit 1
  test_hybrid-coal -sp trees/7_tax_sp_james02.tre || exit 1
  test_hybrid-coal -sp trees/7_tax_sp_james02.tre -gtopo || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt1_para -gtopo || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt1_para -plot || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt1_para -plot -branch || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt1_para -plot -label || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt1_para -dot || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt1_para -print || exit 1
  test_hybrid-coal

  # Not OK
  test_hybrid-coal -sp '((((B:1,C:1)s1:1)h1#.5:1,A:3)s2:1,(h1#.5:1,D:3)s3:1)r;' || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt1_para -gt '(((A,D),C),B);' || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt1_para -gt trees/4_tax_gt4.tre || exit 1
  test_hybrid-coal -sp trees/4_tax_sp_nt1_para -gt trees/4_tax_gt4.tre || exit 1
  #test_hybrid-coal -sp trees/4_tax_sp_nt1_para -gt trees/4_tax_gt4.tre -maple || exit 1
  #test_hybrid-coal -sp trees/4_tax_sp_nt1_para -gt trees/4_tax_gt4.tre -latex || exit 1
  #test_hybrid-coal -sp trees/4_tax_sp_nt1_para -maple -latex || exit 1

echo ""
