"""Hash exactly the bundled native sources, independent of their install path."""
import hashlib
from pathlib import Path

def source_fingerprint(root=None):
    root=Path(root) if root is not None else Path(__file__).resolve().parent
    paths=sorted([root/'Makefile',root/'fingerprint.py',*root.glob('include/**/*.hpp'),*root.glob('src/**/*.cpp')])
    digest=hashlib.sha256()
    for path in paths:
        digest.update(path.relative_to(root).as_posix().encode()+b'\0')
        digest.update(path.read_bytes()+b'\0')
    return digest.hexdigest()

if __name__=='__main__':print(source_fingerprint())
