import hashlib
from pathlib import Path
import tempfile
import unittest
import zipfile
from nhanes_activity.sources import verify_file, extract_xport

class SourceTests(unittest.TestCase):
    def test_hash_and_size(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'sample';p.write_bytes(b'abc')
            record={'bytes':'3','sha256':hashlib.sha256(b'abc').hexdigest()}
            verify_file(p,record)
            p.write_bytes(b'abd')
            with self.assertRaises(ValueError):verify_file(p,record)
            p.write_bytes(b'abcd')
            with self.assertRaises(ValueError):verify_file(p,record)

    def test_extraction_and_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as tmp:
            z=Path(tmp)/'source.zip';p=Path(tmp)/'out.xpt'
            with zipfile.ZipFile(z,'w') as out:out.writestr('PAXRAW_C.XPT',b'fixture')
            extract_xport(z,p);self.assertEqual(p.read_bytes(),b'fixture')
            with self.assertRaises(FileExistsError):extract_xport(z,p)

    def test_multiple_members_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            z=Path(tmp)/'source.zip'
            with zipfile.ZipFile(z,'w') as out:
                out.writestr('a.xpt',b'a');out.writestr('b.xpt',b'b')
            with self.assertRaises(ValueError):extract_xport(z,Path(tmp)/'out.xpt')

if __name__=='__main__':unittest.main()
