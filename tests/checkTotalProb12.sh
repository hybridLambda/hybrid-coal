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
  test_hybrid-coal -sp trees/12_tax_sp_nt1_para02_bl06 || exit 1
echo ""
