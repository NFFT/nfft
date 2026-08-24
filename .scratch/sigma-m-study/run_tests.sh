#!/bin/sh
# Build and run checkall_ng in one precision tree, reporting the CUnit result.
# Usage: run_tests.sh <worktree-root> <build-dir>
ROOT="$1"
TREE="$2"
cmake --build "$ROOT/$TREE" -j8 --target checkall_ng > "$ROOT/$TREE/build.log" 2>&1
if [ $? -ne 0 ]; then
  echo "BUILD FAILED ($TREE)"
  tail -20 "$ROOT/$TREE/build.log"
  exit 1
fi
cd "$ROOT/$TREE/tests" || exit 1
./checkall_ng > run.log 2>&1
rc=$?
fails=$(grep -c "CUNIT_RUN_TEST_FAILURE" CUnitAutomated_ng-Results.xml)
echo "$TREE: checkall_ng exit=$rc cunit_failures=$fails"
if [ "$fails" != "0" ]; then
  grep -B1 -A4 "CUNIT_RUN_TEST_FAILURE" CUnitAutomated_ng-Results.xml | head -40
  exit 1
fi
exit $rc
