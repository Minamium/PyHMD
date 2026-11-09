"""Offline plotting; optional dependencies are loaded only when requested."""
from pathlib import Path

def render_result(result_dir,*,output_dir=None,xlim=None,ylim=None):
    from ._viz import load_snapshot,render_snapshot,render_history,write_vtk
    root=Path(result_dir);out=Path(output_dir) if output_dir is not None else root/'visualization'
    out.mkdir(parents=True,exist_ok=True)
    snapshot=load_snapshot(root)
    window={k:v for k,v in [('xlim',xlim),('ylim',ylim)] if v is not None}
    return [p for p in (render_snapshot(snapshot,out/'overview.png',window=window),render_history(root,out/'history.png'),write_vtk(snapshot,out/'final.vtk')) if p is not None]
