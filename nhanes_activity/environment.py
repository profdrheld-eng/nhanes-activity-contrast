"""Record imported code versions and detect conflicting distribution metadata."""
import importlib
import importlib.metadata
import re


def package_record(name):
    module = importlib.import_module(name)
    normalize = lambda value: re.sub(r'[-_.]+', '-', value).lower()
    versions = sorted({d.version for d in importlib.metadata.distributions()
                       if normalize(d.metadata.get('Name', '')) == normalize(name)})
    loaded = str(module.__version__)
    return dict(loaded_version=loaded, distribution_versions=versions,
                module_path=module.__file__, metadata_consistent=versions == [loaded])


def analysis_packages():
    records = {name: package_record(name) for name in ('numpy', 'pandas', 'pyreadstat')}
    inconsistent = [name for name, record in records.items() if not record['metadata_consistent']]
    if inconsistent:
        raise RuntimeError('Conflicting package metadata: ' + ', '.join(inconsistent)
                           + '. Use a clean environment; see README.md.')
    return records
