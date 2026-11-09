"""Compatibility script; use the installed pyhmd package."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from pyhmd.plotting._viz import *
from pyhmd.plotting._viz import _read_table

if __name__=="__main__":
    raise SystemExit(main())
