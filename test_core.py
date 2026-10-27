import unittest, tempfile, ast
from pathlib import Path
from collections import Counter
import numpy as np
from tfcds_core import *

def rec(id,p,s,typ='complete'):
    return Record(id,p,s,id,id,typ,1)

def brute(records,k,lo=1,hi=None,drop=False):
    # Independently written string/dictionary oracle, no sparse-statistics code.
    pool=Counter(s.sequence[j:j+k] for s in records for j in range(len(s.sequence)-k+1))
    pops=list(dict.fromkeys(r.population for r in records))
    width=max(len(r.sequence)-k+1 for r in records)
    result=np.zeros((len(pops),width))
    for pi,p in enumerate(pops):
        rr=[r for r in records if r.population==p]
        stop=max(len(r.sequence)-k+1 for r in rr)-(1 if drop else 0)
        for j in range(stop):
            unique={r.sequence[j:j+k] for r in rr if j+k<=len(r.sequence)}
            result[pi,j]=sum(pool[t] for t in unique if pool[t]>=lo and (hi is None or pool[t]<=hi))
    return result

class CoreTests(unittest.TestCase):
    def setUp(self):
        self.rr=[rec('a1','A','AACCGG'),rec('a2','A','AACCAA'),rec('b1','B','AAGGTT'),rec('c1','C','TTCCGG')]
    def test_oracle_many_thresholds(self):
        for k in [1,2,3,6]:
            data=ProfileData(self.rr,k)
            for lo,hi in [(1,None),(2,None),(2,3),(3,3)]:
                for drop in [False,True]:
                    actual,_=data.positional(lower=lo,upper=hi,legacy_end_drop=drop)
                    np.testing.assert_array_equal(actual,brute(self.rr,k,lo,hi,drop))
    def test_subsampling_recomputes_pool_and_presence(self):
        data=ProfileData(self.rr,2);sel=np.array([1,0,1,1])
        actual,_=data.positional(sel)
        np.testing.assert_array_equal(actual,brute([self.rr[i] for i in [0,2,3]],2))
    def test_no_cds_boundary_kmers(self):
        data=ProfileData([rec('a','P','AAA'),rec('b','P','CCC')],3)
        self.assertEqual(set(data.words),{'AAA','CCC'})
        f,_=data.frequency();self.assertAlmostEqual(float(f.sum()),1)
    def test_exact_k_boundary(self):
        data=ProfileData([rec('a','A','AAA'),rec('b','B','CCC')],3)
        np.testing.assert_array_equal(data.positional()[0],[[1],[1]])
        np.testing.assert_array_equal(data.positional(legacy_end_drop=True)[0],[[0],[0]])
    def test_frequency_distance_and_order(self):
        data=ProfileData(self.rr,2);f,_=data.frequency();d=frequency_distance(f)
        np.testing.assert_allclose(d,squareform(pdist(np.sqrt(f.toarray()))),atol=1e-12)
        rev=ProfileData(list(reversed(self.rr)),2,populations=data.populations)
        np.testing.assert_allclose(d,frequency_distance(rev.frequency()[0]),atol=1e-12)
    def test_duplicates_are_local_not_global(self):
        rr=[rec('1','A','AAA'),rec('2','A','AAA'),rec('3','B','AAA')]
        self.assertEqual(len(choose_records(rr,'exact_unique')),2)
    def test_parser_rejects_bad_alphabet_and_keeps_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'a.fa';p.write_text('>A_g_i1.p1 type:complete\nAAA\n>B_g_i1.p1 type:complete\n',encoding='utf-8')
            rr=read_fasta(p);self.assertEqual(len(rr),2);self.assertEqual(rr[1].sequence,'')
            p.write_text('>a\nAAN\n',encoding='utf-8')
            with self.assertRaises(ValueError):read_fasta(p)
    def test_external_pool_weight(self):
        rr=[rec('1','A','AAA'),rec('2','A','AAA'),rec('3','A','CCC'),rec('4','B','AAA'),rec('5','B','CCC')]
        self.assertEqual(ProfileData(rr,3).positional()[0][0,0],5)
        rr.append(rec('6','B','AAA'))
        self.assertEqual(ProfileData(rr,3).positional()[0][0,0],6)
    def test_canonical_explicit(self):
        rr=[rec('1','A','AAA'),rec('2','B','TTT')]
        self.assertEqual(ProfileData(rr,3,True).v,1)
        self.assertEqual(ProfileData(rr,3,False).v,2)
    def test_identical_frequency_vectors_not_false_clusters(self):
        rng=np.random.default_rng(12)
        seq=''.join(rng.choice(list('ACGT'),3000))
        rr=[rec(f'{p}_{i}',p,seq) for p,n in [('A',1),('B',2),('C',4)] for i in range(n)]
        d=frequency_distance(ProfileData(rr,20).frequency()[0])
        self.assertLess(float(d.max()),1e-12)

if __name__=='__main__':unittest.main(verbosity=2)
