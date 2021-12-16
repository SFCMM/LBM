#!/bin/bash

BINPATH=""
SUCCESS=0
FAILED=0

if [[ -z "${SFCMM_LBM_BIN}" ]]; then
  echo "Environment variable SFCMM_LBM_BIN is not set!"
  exit
else
  BINPATH="${SFCMM_LBM_BIN}"
fi

DEFAULT_OMP_THREADS=4
LOG_PATH="run_log_"$(date +%s)

run_test(){
  NAME="$(cut -d'.' -f1 <<<"$1")"

  echo "Running testcase $NAME..."

  if OMP_NUM_THREADS=$DEFAULT_OMP_THREADS $BINPATH "$1" &> "$NAME".stdout; then
    echo "passed"
    SUCCESS=$((SUCCESS+1))
  else
    echo "FAILED"
    FAILED=$((FAILED+1))
  fi

  cp lbm_log ../"$LOG_PATH"/"$NAME"_lbm_log.txt
  cp gridgen_log ../"$LOG_PATH"/"$NAME"_gridgen_log.txt
  cp "$NAME".stdout ../"$LOG_PATH"/"$NAME".stdout
}

mkdir "$LOG_PATH"

./clean.sh
cd couette || exit
run_test couette.json
run_test couette_bnd.json
cd ..

cd poiseuille || exit
run_test poiseuille.json
cd ..

cd poisson || exit
run_test poisson1D.json
run_test poisson2D.json
cd ..

./clean.sh

if [ $FAILED == 0 ]; then
  echo "All $SUCCESS testcases have been successful"
  exit 0
else
  echo "$FAILED testcases have failed"
  exit 1
fi
