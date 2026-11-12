"""Check real checkpoint decoding and detection of isolated accounting corruption."""
import json
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile
import unittest
from pymhd.results.checkpoint import checksum,read_checkpoint,compare_checkpoints,CheckpointError
from pymhd import Config,Engine,RuntimeConfig,ResultOptions
ROOT=Path(__file__).resolve().parents[2];BINARY=ROOT/'build/pymhd_cpu'

@unittest.skipUnless(BINARY.exists(),'build native CPU first')
class CheckpointContracts(unittest.TestCase):
    def test_midplane_observer_is_passive(self):
        with tempfile.TemporaryDirectory(prefix='pymhd_observer_') as temp:
            results=[]
            cfg=Config(problem='thin_disk',grid={'nx1':64,'nx2':32},time={'t_end':.1})
            for enabled in (False,True):
                result=Engine(cfg,runtime=RuntimeConfig(binary=BINARY),result_options=ResultOptions(out_dir=Path(temp)/str(enabled),
                    output_interval=.05,checkpoint_interval=.1,midplane_diagnostics=enabled)).run()
                self.assertTrue(result.complete);results.append(result)
            a,b=[r.out_dir/r.checkpoints[-1]['path'] for r in results]
            self.assertEqual(compare_checkpoints(a,b),[])
            self.assertTrue((results[1].out_dir/'midplane_metrics.csv').is_file())
            self.assertFalse((results[0].out_dir/'midplane_metrics.csv').exists())

    def test_all_state_and_accounting_fields_and_fault_detection(self):
        with tempfile.TemporaryDirectory(prefix='pymhd_cp_contract_') as temp:
            root=Path(temp);out=root/'run'
            subprocess.run([str(BINARY),'--config',str(ROOT/'configs/thin_disk.ini'),'--nx1','16','--nx2','8','--t-end','.01','--checkpoint-interval','.01','--output',str(out)],check=True,stdout=subprocess.DEVNULL)
            state=json.loads((out/'status.json').read_text());checkpoint=out/f'checkpoint_step{state["step"]}'
            shard=read_checkpoint(checkpoint)[0]
            self.assertEqual(set(shard.buffers),{'U','Psi','floors','floor_mass','floor_energy','boundary_energy','thermal','flux_accounting'})
            self.assertEqual(compare_checkpoints(checkpoint,checkpoint),[])
            for field in ('U','thermal','flux_accounting','counters'):
                target=root/field;shutil.copytree(checkpoint,target)
                file=target/'rank000000.bin';raw=bytearray(file.read_bytes())
                if field=='counters':
                    position=raw.find(shard.signature.encode())+len(shard.signature.encode())+len(shard.clock)
                    raw[position:position+8]=struct.pack('<q',shard.rejected_steps+1)
                else:
                    payload=shard.buffers[field];position=raw.find(payload)
                    self.assertGreaterEqual(position,0)
                    value=struct.unpack_from('<d',raw,position)[0]
                    struct.pack_into('<d',raw,position,value+1e-5)
                struct.pack_into('<Q',raw,len(raw)-8,checksum(raw[:-8]));file.write_bytes(raw)
                differences=compare_checkpoints(checkpoint,target)
                self.assertTrue(any(field in entry for entry in differences),differences)
            raw=bytearray((checkpoint/'rank000000.bin').read_bytes());raw[-9]^=1
            bad=root/'bad';shutil.copytree(checkpoint,bad);(bad/'rank000000.bin').write_bytes(raw)
            with self.assertRaises(CheckpointError):read_checkpoint(bad)

if __name__=='__main__':unittest.main()
