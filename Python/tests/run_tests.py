"""Dependency-free runner (pytest was not installable offline): python tests/run_tests.py"""
import sys, os, time, traceback
sys.path.insert(0, os.path.dirname(__file__))
import test_core
tests = [(n, f) for n, f in vars(test_core).items() if n.startswith("test_") and callable(f)]
fail = 0
for n, f in tests:
    try:
        f(); print(f"PASS  {n}")
    except Exception:
        fail += 1; print(f"FAIL  {n}"); traceback.print_exc()
print(f"\n{len(tests)-fail}/{len(tests)} passed")
sys.exit(1 if fail else 0)
