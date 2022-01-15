#!/bin/bash
# reset the test case directory

# todo: excluding doesn't really work...
find . -name "*.vtp" -type f -delete -not -path "./run_log*/*"
find . -name "*_log" -type f -delete -not -path "./run_log*/*"
find . -name "*.csv" -type f -delete -not -path "./run_log*/*"
find . -name "*.stdout" -type f -delete -not -path "./run_log*/*"