import unittest
from types import SimpleNamespace
from unittest.mock import patch
from nhanes_activity.environment import package_record, analysis_packages


class EnvironmentTests(unittest.TestCase):
    def test_preflight_rejects_conflicting_metadata(self):
        with patch('nhanes_activity.environment.package_record', return_value={'metadata_consistent':False}):
            with self.assertRaisesRegex(RuntimeError, 'clean environment'):
                analysis_packages()

    def test_preflight_accepts_consistent_metadata(self):
        with patch('nhanes_activity.environment.package_record', return_value={'metadata_consistent':True}):
            self.assertEqual(set(analysis_packages()), {'numpy','pandas','pyreadstat'})

    def test_loaded_version_is_not_replaced_by_stale_metadata(self):
        module = SimpleNamespace(__version__='1.2.8', __file__='/synthetic/module.py')
        distributions = [SimpleNamespace(metadata={'Name':'example'},version=v)
                         for v in ('1.3.5','1.2.8')]
        with patch('nhanes_activity.environment.importlib.import_module',return_value=module), \
             patch('nhanes_activity.environment.importlib.metadata.distributions',return_value=distributions):
            r = package_record('example')
        self.assertEqual(r['loaded_version'], '1.2.8')
        self.assertFalse(r['metadata_consistent'])
        self.assertEqual(r['distribution_versions'], ['1.2.8','1.3.5'])

    def test_single_matching_distribution(self):
        module = SimpleNamespace(__version__='2.0', __file__='/synthetic/module.py')
        distribution = SimpleNamespace(metadata={'Name':'example'},version='2.0')
        with patch('nhanes_activity.environment.importlib.import_module',return_value=module), \
             patch('nhanes_activity.environment.importlib.metadata.distributions',return_value=[distribution]):
            self.assertTrue(package_record('example')['metadata_consistent'])
