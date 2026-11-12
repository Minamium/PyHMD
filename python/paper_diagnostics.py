"""Compatibility script; use the installed pymhd package."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from pymhd.diagnostics.paper import *

if __name__=="__main__":
    raise SystemExit(main())
