"""Compatibility entry point for installed diagnostic tools."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from pyhmd.diagnostics.transport import *

if __name__=="__main__":
    raise SystemExit(main())
