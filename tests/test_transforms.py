import unittest
from statistics import NormalDist
import numpy as np
from nhanes_activity.transforms import weighted_midrank_normal, weighted_standardize

class TransformTests(unittest.TestCase):
    def test_unequal_weight_ties(self):
        # x=1 has weight 2, x=2 has pooled weight 4, x=4 has weight 4.
        x = [2, 1, 2, 4]; w = [1, 2, 3, 4]
        expected = [NormalDist().inv_cdf(p) for p in [.4, .1, .4, .8]]
        np.testing.assert_allclose(weighted_midrank_normal(x, w), expected, atol=1e-14)

    def test_tail_bounds(self):
        expected = [NormalDist().inv_cdf(.0001), 0, NormalDist().inv_cdf(.9999)]
        np.testing.assert_allclose(weighted_midrank_normal([0,1,2], [1, 1000000, 1]), expected, atol=1e-14)

    def test_population_sd(self):
        z, mean, sd = weighted_standardize([1,3], [1,3])
        self.assertEqual(mean, 2.5)
        self.assertAlmostEqual(sd, np.sqrt(.75))
        self.assertAlmostEqual(np.average(z, weights=[1,3]), 0)
        self.assertAlmostEqual(np.average(z*z, weights=[1,3]), 1)

    def test_invalid_inputs(self):
        for x,w in [([],[]),([1,2],[1]),([[1],[2]],[1,1]),([1,np.nan],[1,1]),
                    ([1,2],[1,0]),([1,2],[1,-1]),([1,2],[1,np.inf])]:
            for fn in [weighted_midrank_normal, weighted_standardize]:
                with self.subTest(x=x,w=w,fn=fn.__name__), self.assertRaises(ValueError): fn(x,w)

    def test_constant_values(self):
        np.testing.assert_array_equal(weighted_midrank_normal([2,2],[1,3]), [0,0])
        with self.assertRaises(ValueError): weighted_standardize([2,2],[1,3])

    def test_order_and_weight_scaling(self):
        x=np.array([4.,1.,2.,2.]); w=np.array([3.,2.,1.,4.]); order=[2,0,3,1]
        z=weighted_midrank_normal(x,w)
        np.testing.assert_allclose(weighted_midrank_normal(x[order],w[order]),z[order],atol=1e-14)
        np.testing.assert_allclose(weighted_midrank_normal(x,w*2),z,atol=1e-14)

if __name__ == '__main__': unittest.main()
