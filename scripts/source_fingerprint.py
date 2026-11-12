#!/usr/bin/env python3
"""Compatibility entry point for the native implementation fingerprint."""
import runpy
from pathlib import Path
runpy.run_path(str(Path(__file__).resolve().parents[1]/'src/pymhd/native/fingerprint.py'),run_name='__main__')
