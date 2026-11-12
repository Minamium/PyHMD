"""Read-only checkpoint inspection, independent of a native binary's pointers.

Version dispatch is explicit. Unknown future layouts fail instead of guessing
offsets. This module never rewrites signatures or creates restart files.
"""
from dataclasses import dataclass
from pathlib import Path
import struct
import math


class CheckpointError(ValueError):
    pass


def checksum(data):
    value = 14695981039346656037
    for byte in data:
        value = ((value ^ byte) * 1099511628211) & 0xffffffffffffffff
    return value


class Reader:
    def __init__(self, data):
        self.data, self.offset = data, 0

    def take(self, size):
        if size < 0 or self.offset + size > len(self.data):
            raise CheckpointError('truncated checkpoint')
        value = self.data[self.offset:self.offset + size]
        self.offset += size
        return value

    def value(self, format):
        return struct.unpack('<' + format, self.take(struct.calcsize('<' + format)))[0]

    def string(self):
        size = self.value('Q')
        if size > 16384:
            raise CheckpointError('invalid signature/string size')
        return self.take(size).decode('utf-8')

    def block(self, width, expected=None):
        count = self.value('Q')
        if expected is not None and count != expected:
            raise CheckpointError(f'buffer count {count} differs from schema {expected}')
        return self.take(count * width * 8)

    def finish(self):
        position = self.offset
        stored = self.value('Q')
        if self.offset != len(self.data):
            raise CheckpointError('unexpected trailing bytes')
        if stored != checksum(self.data[:position]):
            raise CheckpointError('checkpoint checksum mismatch')
        return stored


@dataclass(frozen=True)
class Checkpoint:
    version: int
    nv: int
    ng: int
    layout: tuple
    signature: str
    clock: bytes
    counters: bytes
    buffers: dict
    file_checksum: int
    transport: dict | None = None

    @property
    def time(self): return struct.unpack('<4dq', self.clock)[0]
    @property
    def step(self): return struct.unpack('<4dq', self.clock)[4]
    @property
    def rejected_steps(self): return struct.unpack('<2q', self.counters)[0]


def _version1(reader):
    if reader.value('I') != 0x01020304:
        raise CheckpointError('unsupported byte order')
    nv, ng = reader.value('I'), reader.value('I')
    if nv != 10 or ng != 3:
        raise CheckpointError('version 1 requires NV=10, NG=3')
    layout = tuple(reader.value('i') for _ in range(8))
    rank, ranks, nx, ny, tx, ty, start, global_nx = layout
    if not (0 <= rank < ranks and nx > 0 and ny > 0 and tx == nx + 2*ng
            and ty == ny + 2*ng and 0 <= start <= global_nx - nx):
        raise CheckpointError('invalid rank/grid layout')
    signature = reader.string()
    clock, counters = reader.take(struct.calcsize('<4dq')), reader.take(struct.calcsize('<2q'))
    time, dt, next_output, interval, step = struct.unpack('<4dq', clock)
    if not (all(math.isfinite(x) for x in (time, dt, next_output, interval))
            and time >= 0 and dt >= 0 and next_output > time and interval > 0 and step >= 0
            and min(struct.unpack('<2q', counters)) >= 0):
        raise CheckpointError('invalid accepted clock/counters')
    cells, vertices = tx*ty, (tx+1)*(ty+1)
    owned = nx*ny
    fields = [('U', nv, cells), ('Psi', 1, vertices), ('floors', 1, owned),
              ('floor_mass', 1, owned), ('floor_energy', 1, owned),
              ('boundary_energy', 1, owned), ('thermal', 8, owned), ('flux_accounting', 4, owned)]
    buffers = {name: reader.block(width, count) for name, width, count in fields}
    return Checkpoint(1, nv, ng, layout, signature, clock, counters, buffers, reader.finish())


READERS = {1: _version1}


def read_shard(path):
    reader = Reader(Path(path).read_bytes())
    if reader.take(16) != b'PyMHD-CHECKPOINT':
        raise CheckpointError('invalid checkpoint magic')
    version = reader.value('I')
    if version not in READERS:
        raise CheckpointError(f'unsupported checkpoint version {version}')
    shard = READERS[version](reader)
    extra = Path(path).with_name('transport_' + Path(path).name)
    if extra.exists():
        transport = Reader(extra.read_bytes())
        if transport.string() != 'PyMHD-TRANSPORT-1' or transport.string() != shard.signature:
            raise CheckpointError('transport schema/signature mismatch')
        if transport.value('Q') != shard.file_checksum or transport.value('d') != shard.time or transport.value('q') != shard.step:
            raise CheckpointError('transport belongs to a different checkpoint generation')
        rank, ranks, nx, ny, tx, ty, start, global_nx = shard.layout
        if tuple(transport.value('i') for _ in range(6)) != (rank, ranks, nx, ny, start, global_nx):
            raise CheckpointError('transport rank/grid mismatch')
        data = {'origin_clock': transport.take(struct.calcsize('<dq'))}
        for name, width, count in [('radial', 8, (nx+1)*ny), ('polar', 8, nx*(ny+1)),
                                    ('vertices', 4, (tx+1)*(ty+1)), ('origin_U', shard.nv, tx*ty),
                                    ('origin_Psi', 1, (tx+1)*(ty+1))]:
            data[name] = transport.block(width, count)
        transport.finish()
        object.__setattr__(shard, 'transport', data)
    # v1 revision7 and phase0 keep the first 26 integer fields unchanged.
    lines = shard.signature.splitlines()
    integers = lines[3 if lines[2]=='fp_contract_off' else 2].split()
    if len(integers) >= 26 and int(integers[25]) == 1 and shard.transport is None:
        raise CheckpointError('required transport sidecar missing')
    return shard


def read_checkpoint(directory):
    directory = Path(directory)
    lines = (directory/'checkpoint.meta').read_text().splitlines()
    name, version, ranks = lines[0].split()
    if name != 'PyMHD_CHECKPOINT' or int(version) not in READERS or int(ranks) <= 0:
        raise CheckpointError('unsupported checkpoint manifest')
    if len(list(directory.glob('rank*.bin'))) != int(ranks):
        raise CheckpointError('checkpoint shard count mismatch')
    shards = [read_shard(directory/f'rank{rank:06}.bin') for rank in range(int(ranks))]
    for rank, shard in enumerate(shards):
        if shard.layout[:2] != (rank, len(shards)) or shard.clock != shards[0].clock:
            raise CheckpointError('inconsistent checkpoint shards')
        tokens = lines[1].split()
        if tokens[:1] != ['time'] or tokens[2:3] != ['step'] or float(tokens[1]) != shard.time or int(tokens[3]) != shard.step:
            raise CheckpointError('manifest clock differs from shards')
    return shards


def compare_checkpoints(left, right):
    a, b = read_checkpoint(left), read_checkpoint(right)
    differences = []
    if len(a) != len(b): return ['rank_count']
    for rank, (x, y) in enumerate(zip(a, b)):
        for name in ('version', 'nv', 'ng', 'layout', 'clock', 'counters', 'buffers', 'transport'):
            if getattr(x, name) != getattr(y, name):
                if name in ('buffers', 'transport') and getattr(x, name) is not None and getattr(y, name) is not None:
                    for key in set(getattr(x, name)) | set(getattr(y, name)):
                        if getattr(x, name).get(key) != getattr(y, name).get(key):
                            differences.append(f'rank{rank}/{name}/{key}')
                else:
                    differences.append(f'rank{rank}/{name}')
    return differences
