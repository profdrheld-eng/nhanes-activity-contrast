"""Validate locally downloaded public-use sources before processing."""
import csv
import hashlib
from pathlib import Path
import shutil
import zipfile


def source_records(manifest: Path):
    with manifest.open(newline='') as handle:
        rows = list(csv.DictReader(handle))
    paths = set()
    for row in rows:
        path = Path(row['relative_path'])
        if path.is_absolute() or '..' in path.parts or str(path) in paths:
            raise ValueError('Source paths must be unique, relative and contained')
        if path.name != row['filename']:
            raise ValueError('Source filename and relative path disagree')
        paths.add(str(path))
    return rows


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda:handle.read(8*1024*1024),b''):
            digest.update(chunk)
    return digest.hexdigest()


def verify_file(path: Path, record):
    if path.stat().st_size != int(record['bytes']):
        raise ValueError(f'Source byte size differs: {path.name}')
    if sha256(path) != record['sha256']:
        raise ValueError(f'Source SHA-256 differs: {path.name}')


def extract_xport(archive: Path, destination: Path):
    """Extract one XPORT by streaming to the caller's fixed path, never member paths."""
    with zipfile.ZipFile(archive) as handle:
        members = [m for m in handle.infolist() if not m.is_dir() and m.filename.lower().endswith('.xpt')]
        if len(members) != 1:
            raise ValueError('Expected exactly one XPORT member')
        member = members[0]
        with handle.open(member) as source, destination.open('xb') as output:
            shutil.copyfileobj(source,output,length=8*1024*1024)
        if destination.stat().st_size != member.file_size:
            raise ValueError('Extracted XPORT byte size differs')
