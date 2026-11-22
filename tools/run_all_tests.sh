#!/bin/bash
# Phase-2 test driver: Python tests, C++ build + tests, cross-language comparison.
# (cmake is optional; falls back to direct g++.)  MATLAB: run `main('all')` inside MATLAB afterwards,
# then re-run tools/cross_validate.py to include the MATLAB column.
set -e
cd "$(dirname "$0")/.."
python3 tools/generate_cases.py
python3 Python/tests/run_tests.py
mkdir -p build
SRC="Cpp/src/material.cpp Cpp/src/tmm.cpp Cpp/src/case_io.cpp"
g++ -std=c++17 -O2 -Wall -Wextra -ICpp/include $SRC Cpp/tests/test_core.cpp -o build/test_core
g++ -std=c++17 -O2 -Wall -Wextra -ICpp/include $SRC Cpp/src/main.cpp -o build/nfs_run
./build/test_core
python3 tools/cross_validate.py
