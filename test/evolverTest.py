import os, sys
import pytest
import unittest
import subprocess
from sonLib.bioio import getTempDirectory, popenCatch

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

    def _run_evolver_decomposed(self, binariesMode):
        """ Run the full evolver test, putting the jobstore and output in tempDir
        but instead of doing in in one shot like above, use cactus-prepare, cactus-blast
        and cactus-align to break it into different steps """

        out_dir = os.path.join(self.tempDir, 'output')
        out_seqfile = os.path.join(out_dir, 'evolverMammalsOut.txt')
        in_seqfile = './examples/evolverMammals.txt'
        cmd = ['cactus-prepare', in_seqfile, out_dir, out_seqfile, self._out_hal(binariesMode),
               '--jobStore', self._job_store(binariesMode)]

        job_plan = popenCatch(' '.join(cmd))

        for line in job_plan.split('\n'):
            line = line.strip()
            if len(line) > 0 and not line.startswith('#'):
                subprocess.check_call(line, shell=True)

    def _csvstr_to_table(self, csvstr, header_fields):
        """ Hacky csv parse """
        output_stats = {}
        skip = True
        # filter out everything but the csv stats from the halStats output
        for line in csvstr.split('\n'):
            toks = [tok.strip() for tok in str(line).split(',')]
            if skip:
                if all([header in toks for header in header_fields]):
                    skip = False                    
            elif len(toks) > len(header_fields):
                output_stats[toks[0]] = toks[1:]
        return output_stats

    def _check_stats(self, halPath, delta_pct=0.25):
        """ Compare halStats otuput of given file to baseline 
        """
        # this is just pasted from a successful run.  it will be used to catch serious regressions
        ground_truth = '''GenomeName, NumChildren, Length, NumSequences, NumTopSegments, NumBottomSegments
        Anc0, 2, 545125, 9, 0, 14471
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
        sys.stderr.write("\nComparing stats\n{}\nvs ground truth\n{}\n".format(output.decode("utf-8"), ground_truth))
        output_table = self._csvstr_to_table(output.decode("utf-8"), ['GenomeName', 'NumChildren'])
        truth_table = self._csvstr_to_table(ground_truth, ['GenomeName', 'NumChildren'])        

        # make sure the stats are roughly the same
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

    def _check_coverage(self, halPath, delta_pct=0.05):
        """ Compare halStats otuput of given file to baseline 
        """
        # this is just pasted from a successful run.  it will be used to catch serious regressions
        ground_truth = '''Genome, sitesCovered1Times, sitesCovered2Times
        simMouse_chr6, 1000000, 0
        simHuman_chr6, 721276, 4846
        simCow_chr6, 671888, 7535
        simDog_chr6, 680661, 4527
        simRat_chr6, 902624, 6863'''

        # run halCoverage on the evolver output
        proc = subprocess.Popen(['bin/halCoverage',  halPath, 'simMouse_chr6', '--hdf5InMemory'], stdout=subprocess.PIPE)
        output, errors = proc.communicate()
        sts = proc.wait()
        self.assertEqual(sts, 0)
        sys.stderr.write("\n\nComparing coverage\n{}\nvs ground truth\n{}\n".format(output.decode("utf-8"), ground_truth))
        output_table = self._csvstr_to_table(output.decode("utf-8"), ['Genome', 'sitesCovered1Times'])
        truth_table = self._csvstr_to_table(ground_truth, ['Genome', 'sitesCovered1Times'])

        # make sure the stats are roughly the same
        self.assertEqual(len(truth_table), len(output_table))
        for key, val in truth_table.items():
            self.assertTrue(key in output_table)
            oval = output_table[key]
            self.assertEqual(len(val), len(oval))
            for i in range(len(val)):
                delta = delta_pct * int(val[i])
                self.assertGreaterEqual(int(oval[i]), int(val[i]) - delta)
                self.assertLessEqual(int(oval[i]), int(val[i]) + delta)                    
            
    def testEvolverLocalPrepare(self):
        """ Check that the output of halStats on a hal file produced by running cactus with --binariesMode local is 
        is reasonable
        """
        # run cactus step by step via the plan made by cactus-prepare
        self._run_evolver_decomposed("local")

        # check the output
        self._check_stats(self._out_hal("local"))
        self._check_coverage(self._out_hal("local"))

    def testEvolverDocker(self):
        """ Check that the output of halStats on a hal file produced by running cactus with --binariesMode docker is 
        is reasonable.  Note: the local image being tested should be set up via CACTUS_DOCKER_ORG (with tag==latest)
        """
        # run cactus
        self._run_evolver("docker")

        # check the output
        self._check_stats(self._out_hal("docker"))
        self._check_coverage(self._out_hal("docker"))

if __name__ == '__main__':
    unittest.main()
