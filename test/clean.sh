#!/bin/bash

find . -name "*.vtp" -type f -delete -not -path "./run_log*/*"
find . -name "*_log" -type f -delete -not -path "./run_log*/*"
find . -name "*.csv" -type f -delete -not -path "./run_log*/*"