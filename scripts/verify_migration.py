#!/usr/bin/env python3
"""Compare a separately built baseline with the reorganized native engine.

Requires two real executables. A local compiler/platform bitwise gate, not a
cross-compiler precision claim or a paper-reproduction test.
"""
import argparse
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import tempfile

def payload(path):
    data=path.read_bytes();assert data[:16]==b'PyMHD-CHECKPOINT'
    offset=64;n=struct.unpack_from('<Q',data,offset)[0];offset+=8+n
    clock=data[offset:offset+40];offset+=40+16
    arrays=[]
    for width in (10,1):
        n=struct.unpack_from('<Q',data,offset)[0];offset+=8
        size=n*width*8;arrays.append(data[offset:offset+size]);offset+=size
    return clock,arrays

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference',type=Path,required=True);p.add_argument('--candidate',type=Path,default=Path('build/pyhmd_cpu'))
    p.add_argument('--target-grid',action='store_true');p.add_argument('--report',type=Path,required=True)
    args=p.parse_args();root=Path(__file__).resolve().parents[1]
    binaries={'reference':args.reference.resolve(),'candidate':args.candidate.resolve()}
    cases=[('uniform',32,16,.05),('alfven',64,4,.1),('diffusion',64,4,.05),('brio_wu',64,4,.03),('thin_disk',64,32,.1)]
    report={'passed':True,'scope':'same compiler/platform, full CSV and binary U/Psi/clock bits; no GPU or scientific equivalence claim',
            'binary_sha256':{key:hashlib.sha256(binary.read_bytes()).hexdigest() for key,binary in binaries.items()},'cases':[]}
    if args.target_grid: cases.append(('thin_disk',600,128,100))
    with tempfile.TemporaryDirectory(prefix='pyhmd_migration_') as tmp:
        for case,nx,ny,end in cases:
            outputs={}
            for label,binary in binaries.items():
                out=Path(tmp)/f'{case}_{nx}_{label}';outputs[label]=out
                command=[str(binary),'--config',str(root/'configs'/f'{case}.ini'),'--nx1',str(nx),'--nx2',str(ny),
                         '--t-end',str(end),'--output-interval',str(end/2),'--checkpoint-interval',str(end/2),
                         '--max-steps','100' if nx==600 else '100000000','--output',str(out)]
                process=subprocess.run(command,capture_output=True,text=True)
                expected=3 if nx==600 else 0
                assert process.returncode==expected,process.stdout+process.stderr
            left,right=(outputs[k] for k in ('reference','candidate'))
            sa,sb=(json.loads((folder/'status.json').read_text()) for folder in (left,right))
            csv_equal=all(a.read_bytes()==b.read_bytes() for a,b in zip(sorted(left.glob('snapshot_*.csv')),sorted(right.glob('snapshot_*.csv'))))
            assert len(list(left.glob('snapshot_*.csv')))==len(list(right.glob('snapshot_*.csv')))
            ca,cb=(folder/f"checkpoint_step{status['step']}"/'rank000000.bin' for folder,status in ((left,sa),(right,sb)))
            clock_a,arrays_a=payload(ca);clock_b,arrays_b=payload(cb)
            item=dict(case=case,grid=[nx,ny],step=sb['step'],time=sb['time'],csv_bits_equal=csv_equal,
                      U_Psi_bits_equal=arrays_a==arrays_b,clock_bits_equal=clock_a==clock_b)
            report['cases'].append(item)
            report['passed'] &= csv_equal and arrays_a==arrays_b and clock_a==clock_b
            print(json.dumps(item),flush=True)
    args.report.parent.mkdir(parents=True,exist_ok=True);args.report.write_text(json.dumps(report,indent=2)+'\n')
    return 0 if report['passed'] else 1

if __name__=='__main__':raise SystemExit(main())
