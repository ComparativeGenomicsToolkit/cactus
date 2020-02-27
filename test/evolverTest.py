import os
import pytest
import unittest
import subprocess
from sonLib.bioio import getTempDirectory

class TestCase(unittest.TestCase):
    def setUp(self):
        self.tempDir = getTempDirectory(os.getcwd())
        unittest.TestCase.setUp(self)

    def tearDown(self):
        unittest.TestCase.tearDown(self)
        os.system("rm -rf %s" % self.tempDir)

    def _job_store(self, binariesMode):
        return os.path.join(self.tempDir, 'js-{}'.format(binariesMode))

    def _out_hal(self, binariesMode):
        return os.path.join(self.tempDir, 'evovler-{}.hal'.format(binariesMode))

    def _run_evolver(self, binariesMode):
        """ Run the full evolver test, putting the jobstore and output in tempDir 
        """
        cmd = ['cactus', self._job_store(binariesMode), './examples/evolverMammals.txt', self._out_hal(binariesMode),
                               '--binariesMode', binariesMode, '--logInfo', '--realTimeLogging', '--workDir', self.tempDir]
        # todo: it'd be nice to have an interface for setting tag to something not latest or commit
        if binariesMode == 'docker':
            cmd += ['--latest']
            
        subprocess.check_call(' '.join(cmd), shell=True)

    def _check_stats(self, halPath, delta_pct=0.25):
        """ Compare halStats otuput of given file to baseline 
        """
        # this is just pasted from a successful run.  it will be used to catch serious regressions
        ground_truth = '''Anc0, 2, 545125, 9, 0, 14471
        Anc1, 2, 569645, 6, 16862, 54019
        simHuman_chr6, 0, 601863, 1, 53463, 0
        mr, 2, 611571, 3, 52338, 55582
        simMouse_chr6, 0, 636262, 1, 56172, 0
        simRat_chr6, 0, 647215, 1, 55687, 0
        Anc2, 2, 577109, 15, 17350, 57130
        simCow_chr6, 0, 602619, 1, 56309, 0
        simDog_chr6, 0, 593897, 1, 55990, 0'''

        # run halStats on the evolver output
        proc = subprocess.Popen(['bin/halStats',  halPath], stdout=subprocess.PIPE)
        output, errors = proc.communicate()
        sts = proc.wait()
        self.assertEqual(sts, 0)
        output_stats = ''
        skip=True
        
        # filter out everything but the csv stats from the halStats output
        for line in output.decode("utf-8").split('\n'):
            toks = [tok.strip() for tok in str(line).split(',')]
            if skip:
                if 'GenomeName' in toks and 'NumChildren' in toks:
                    skip = False                    
            else:
                output_stats += line + '\n'

        # convert csv to dict
        def csv_to_table(s):
            t = {}
            for line in s.split('\n'):
                toks = [tok.strip() for tok in line.strip().split(',')]
                if len(toks) == 6:
                    t[toks[0]] = toks[1:]
            return t

        # make sure the stats are roughly the same
        truth_table = csv_to_table(ground_truth)
        output_table = csv_to_table(output_stats)
        self.assertEqual(len(truth_table), len(output_table))
        for key, val in truth_table.items():
            self.assertTrue(key in output_table)
            oval = output_table[key]
            self.assertEqual(len(val), len(oval))
            for i in range(len(val)):
                is_leaf = key.startswith('sim')
                # exact comparison for NumChildren and, for leaves, Lengh and NumSequences
                exact_comp = i == 0 or (is_leaf and i in [1,2])
                # no comparison for NumSequences in Ancestors as it seems to be all over the place
                no_comp = i == 2 and not is_leaf
                if exact_comp:
                    self.assertEqual(int(oval[i]), int(val[i]))
                elif not no_comp:
                    delta = delta_pct * int(val[i])
                    self.assertGreaterEqual(int(oval[i]), int(val[i]) - delta)
                    self.assertLessEqual(int(oval[i]), int(val[i]) + delta)
            
    def testEvolverLocal(self):
        """ Check that the output of halStats on a hal file produced by running cactus with --binariesMode local is 
        is reasonable
        """
        # run cactus
        self._run_evolver("local")

        # check the output
        self._check_stats(self._out_hal("local"))

    def testEvolverDocker(self):
        """ Check that the output of halStats on a hal file produced by running cactus with --binariesMode docker is 
        is reasonable.  Note: the local image being tested should be set up via CACTUS_DOCKER_ORG (with tag==latest)
        """
        # run cactus
        self._run_evolver("docker")

        # check the output
        self._check_stats(self._out_hal("docker"))



if __name__ == '__main__':
    unittest.main()
