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

write_report() {
  touch "$LOG_PATH"/report.csv.txt
  echo "testcase, passed?, time" >>"$LOG_PATH"/report.csv.txt
}

time_test() {
  exec 3>&1 4>&2
  #  time_result=$({ time run_test "$1" 1>&3 2>&4; } 2>&1)
  time_result=$({ time run_test "$1" 1>&3 2>&4; } 2>&1)
  exec 3>&- 4>&-

  #only get user time
  usertime=$(echo -n "$time_result" | grep user)

  #only get numerical value by delimiting the space
  usertime=$(echo $usertime | cut -d' ' -f2)

  #replace , with .
  usertime=$(echo $usertime | tr , .)

  if [ -f ../"$LOG_PATH"/SUCCESS ]; then
    echo "$1" ", passed," "$usertime" >>../"$LOG_PATH"/report.csv.txt
    SUCCESS=$((SUCCESS + 1))
    rm ../"$LOG_PATH"/SUCCESS
  fi
  if [ -f ../"$LOG_PATH"/FAILED ]; then
    echo "$1" ", failed," "$usertime" >>../"$LOG_PATH"/report.csv.txt
    FAILED=$((FAILED + 1))
    rm ../"$LOG_PATH"/FAILED
  fi
}

run_test() {
  NAME="$(cut -d'.' -f1 <<<"$1")"

  echo "Running testcase $NAME..."

  if OMP_NUM_THREADS=$DEFAULT_OMP_THREADS $BINPATH "$1" &>"$NAME".stdout; then
    echo "passed"
    touch ../"$LOG_PATH"/SUCCESS
  else
    echo "FAILED"
    touch ../"$LOG_PATH"/FAILED
  fi

  if [ -f "lbm_log" ]; then
    cp lbm_log ../"$LOG_PATH"/"$NAME"_lbm_log.txt
  fi
  if [ -f "lpt_log" ]; then
    cp lpt_log ../"$LOG_PATH"/"$NAME"_lpt_log.txt
  fi
  if [ -f "gridgen_log" ]; then
    cp gridgen_log ../"$LOG_PATH"/"$NAME"_gridgen_log.txt
  fi
  cp "$NAME".stdout ../"$LOG_PATH"/"$NAME".stdout.txt
}

./clean.sh
mkdir "$LOG_PATH"
write_report

cd couette || exit
time_test couette.json
time_test couette_bnd.json
time_test couette_bnd_bbDirichlet.json
time_test couette_bnd_eq.json
time_test couette_bnd_eq2.json
time_test couette_bnd_NEEM.json
time_test couette_bnd_NEBB.json
time_test couette_bnd_eq_aligned.json
cd ..

cd poiseuille || exit
time_test poiseuille.json
time_test poiseuille_bnd.json
time_test poiseuille_bnd_eq.json
time_test poiseuille_bnd_NEBB.json
time_test poiseuille_bnd_pressure.json
cd ..

cd poisson || exit
time_test poisson1D.json
time_test poisson2D.json
time_test poissonD2Q9.json
cd ..

cd lpt || exit
time_test falling.json
time_test falling1k.json
time_test falling1m.json
time_test falling_implicit.json
time_test falling3d.json
time_test falling3d_terminal.json
time_test falling3d_initial.json
time_test falling_noDrag.json
time_test falling_Buoyancy.json
cd ..

if [ $FAILED == 0 ]; then
  echo "All $SUCCESS testcases have been successful"
  ./clean.sh
  exit 0
else
  echo "$FAILED testcases have failed"
  exit 1
fi
