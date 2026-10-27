"""Regression tests for clustering parameters and input validation."""
from pathlib import Path
import os,subprocess,sys,tempfile,unittest
import numpy as np
from scipy.spatial.distance import pdist,squareform
from sklearn.metrics import silhouette_score
from tfcds_core import analyze_distance

class ReleaseRegressions(unittest.TestCase):
    def test_requested_K_beyond_diagnostic_scan(self):
        distance=squareform(pdist(np.arange(8,dtype=float)[:,None]**1.5))
        info,labels,_=analyze_distance(distance,6)
        self.assertEqual(len(set(labels)),6)
        self.assertEqual(set(info['silhouette_scan']),{2,3,4,5})
        self.assertAlmostEqual(info['silhouette'],silhouette_score(distance,labels,metric='precomputed'))

    def test_invalid_fixed_K_fails_explicitly(self):
        distance=squareform(pdist(np.arange(8,dtype=float)[:,None]))
        for K in (0,1,8,9,2.5):
            with self.subTest(K=K),self.assertRaisesRegex(ValueError,'K must'):
                analyze_distance(distance,K)

    def test_cli_rejects_population_erased_by_frequency_filter(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); fasta=root/'input.fa'; output=root/'out'
            fasta.write_text('>A_g1\nAAA\n>B_g1\nCCC\n>B_g2\nCCC\n>C_g1\nGGG\n>C_g2\nGGG\n>C_g3\nGGG\n',encoding='utf-8')
            cmd=[sys.executable,str(Path(__file__).with_name('main.py')),'--fasta',str(fasta),'--out',str(output),'--k','3','--lower','2','--clusters','2','--representation','kmer_frequency']
            run=subprocess.run(cmd,capture_output=True,text=True,encoding='utf-8',env={**os.environ,'PYTHONIOENCODING':'utf-8','PYTHONDONTWRITEBYTECODE':'1'})
            self.assertNotEqual(run.returncode,0)
            self.assertIn('Zero-profile population after filtering',run.stderr)
            self.assertFalse(output.exists())

if __name__=='__main__':unittest.main(verbosity=2)
