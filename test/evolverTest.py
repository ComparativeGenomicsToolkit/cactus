import os, sys
import pytest
import unittest
import subprocess
import shutil
from sonLib.bioio import getTempDirectory, popenCatch, popen

class TestCase(unittest.TestCase):
    def setUp(self):
        self.tempDir = getTempDirectory(os.getcwd())
        unittest.TestCase.setUp(self)
        self.cromwell = False

    def tearDown(self):
        unittest.TestCase.tearDown(self)
        if self.cromwell:
            # try to catch some cromwell root mode directories that above missies
            # see https://github.com/ComparativeGenomicsToolkit/cactus/issues/255
            os.system("docker run --rm -it -v {}:/data evolvertestdocker/cactus:latest rm -rf {}".format(
                os.path.dirname(self.tempDir.rstrip('/')), os.path.join('/data', os.path.basename(self.tempDir))))
        os.system("rm -rf %s || true" % self.tempDir)

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

        sys.stderr.write('Running {}'.format(' '.format(cmd)))
        subprocess.check_call(' '.join(cmd), shell=True)

    def _run_evolver_decomposed(self, name):
        """ Run the full evolver test, putting the jobstore and output in tempDir
        but instead of doing in in one shot like above, use cactus-prepare, cactus-blast
        and cactus-align to break it into different steps """

        out_dir = os.path.join(self.tempDir, 'output')
        out_seqfile = os.path.join(out_dir, 'evolverMammalsOut.txt')
        in_seqfile = './examples/evolverMammals.txt'
        cmd = ['cactus-prepare', in_seqfile, '--outDir', out_dir, '--outSeqFile', out_seqfile, '--outHal', self._out_hal(name),
               '--jobStore', self._job_store(name)]

        job_plan = popenCatch(' '.join(cmd))

        for line in job_plan.split('\n'):
            line = line.strip()
            if len(line) > 0 and not line.startswith('#'):
                # do Anc2 in binariesMode docker to broaden test coverage
                if 'Anc2' in line and line.startswith('cactus-'):
                    line += ' --binariesMode docker --latest'
                sys.stderr.write('Running {}'.format(line))
                subprocess.check_call(line, shell=True)

    def _run_evolver_decomposed_wdl(self, name):
        """ Run the full evolver test, putting the jobstore and output in tempDir
        but instead of doing in in one shot, use cactus-prepare to make a wdl
        script and run that locally through cromwell """

        # hack to use docker to remove filetree in teardown, as cromwell
        # can leave root files hanging around there
        self.cromwell = True
        
        out_seqfile = os.path.join(self.tempDir, 'evolverMammalsOut.txt')
        in_seqfile = './examples/evolverMammals.txt'
        out_wdl = os.path.join(self.tempDir, 'prepared.wdl')
        cmd = ['cactus-prepare', in_seqfile, '--outHal', self._out_hal(name),
               '--jobStore', self._job_store(name), '--wdl',
               '--dockerImage', 'evolvertestdocker/cactus:latest',
               '--preprocessCores', '2',
               '--blastCores', '4',
               '--alignCores', '4']

        # specify an output directory
        cw_optionsfile = os.path.join(self.tempDir, 'options.json')
        cw_output = os.path.join(self.tempDir, 'cromwell-output')
        with open(cw_optionsfile, 'w') as cw_opts:
            cw_opts.write('{\n')
            cw_opts.write('    \"final_workflow_outputs_dir\": \"{}\",\n'.format(cw_output))
            cw_opts.write('    \"use_relative_output_paths\": true\n')
            cw_opts.write('}\n')

        # download cromwell
        subprocess.check_call(['wget', '-q', 'https://github.com/broadinstitute/cromwell/releases/download/49/cromwell-49.jar'],
                              cwd=self.tempDir)

        # run cactus-prepare and write the output to a wdl file
        popen(' '.join(cmd), out_wdl)

        # run the wdl script in cromwell
        subprocess.check_call(['java', '-jar', './cromwell-49.jar', 'run',
                               out_wdl, '-o', cw_optionsfile], cwd=self.tempDir)

        # fish the file out of cromwell
        shutil.copyfile(os.path.join(cw_output, os.path.basename(self._out_hal(name))), self._out_hal(name))

    def _run_evolver_decomposed_toil(self, name, binariesMode):
        """ Run the full evolver test, putting the jobstore and output in tempDir
        but instead of doing with "cactus", use "cactus-prepare-toil" to run
        the toil-in-toil workflow """

        out_dir = os.path.join(self.tempDir, 'output')
        in_seqfile = './examples/evolverMammals.txt'
        cmd = ['cactus-prepare-toil', self._job_store(name), in_seqfile, '--outDir', out_dir, '--outHal', self._out_hal(name),
               "--defaultMemory", "16G", "--defaultDisk", "20G", "--alignCores", "4"]
        subprocess.check_call(cmd)

    def _run_evolver_decomposed_no_outgroup(self, binariesMode):
        """ Run just the mouse-rat alignment.  Inspired by issues arising here
        https://github.com/ComparativeGenomicsToolkit/cactus/pull/216
        https://github.com/ComparativeGenomicsToolkit/cactus/pull/217 """

        out_dir = os.path.join(self.tempDir, 'output')
        out_seqfile = os.path.join(out_dir, 'evolverMammalsOut.txt')
        in_seqfile = os.path.join(self.tempDir, 'evolverMammalsIn.txt')
        with open(in_seqfile, 'w') as inseq:
            inseq.write('(simMouse_chr6:0.084509,simRat_chr6:0.091589);\n')
            inseq.write('simMouse_chr6 http://s3-us-west-2.amazonaws.com/jcarmstr-misc/testRegions/evolverMammals/simMouse.chr6\n')
            inseq.write('simRat_chr6 http://s3-us-west-2.amazonaws.com/jcarmstr-misc/testRegions/evolverMammals/simRat.chr6\n')

        cmd = ['cactus-prepare', in_seqfile, '--outDir', out_dir, '--outSeqFile', out_seqfile, '--outHal', self._out_hal(binariesMode),
               '--jobStore', self._job_store(binariesMode)]
        job_plan = popenCatch(' '.join(cmd))

        for line in job_plan.split('\n'):
            line = line.strip()
            if len(line) > 0 and not line.startswith('#'):
                # todo interface in prepare
                if line.startswith('cactus-'):
                    line += ' --binariesMode {}'.format(binariesMode)
                    if binariesMode == 'docker':
                        line += ' --latest'
                if line.startswith('cactus-align'):
                    #Remove all the id prefixes to pretend the cigars came not cactus-blast
                    subprocess.check_call('sed -i -e \'s/id=[0,1]|//g\' {}/Anc0.cigar*'.format(out_dir), shell=True)
                    line += ' --nonBlastInput'
                sys.stderr.write('Running {}'.format(line))
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

    def _subset_table(self, csvstr, subset={}):
        """ Hack a truth table to a subset """
        if not subset:
            return csvstr
        lines=csvstr.split('\n')
        ret=''
        for i, line in enumerate(lines):
            if i == 0 or line.split(',')[0].strip() in subset:
                ret += line + '\n'
        return ret

    def _check_stats(self, halPath, delta_pct, subset={}):
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

        ground_truth = self._subset_table(ground_truth, subset)

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

    def _check_coverage(self, halPath, delta_pct, subset={}, columns=3):
        """ Compare halStats otuput of given file to baseline
        """
        # this is just pasted from a successful run.  it will be used to catch serious regressions
        ground_truth = '''Genome, sitesCovered1Times, sitesCovered2Times, sitesCovered3Times
        simMouse_chr6, 636262, 24981, 567
        simRat_chr6, 574916, 21476, 1328
        simHuman_chr6, 459323, 3948, 0
        simCow_chr6, 427278, 4610, 0
        simDog_chr6, 433022, 2905, 0'''

        ground_truth = self._subset_table(ground_truth, subset)

        # run halCoverage on the evolver output
        proc = subprocess.Popen(['bin/halStats',  halPath, '--coverage', 'simMouse_chr6', '--hdf5InMemory'], stdout=subprocess.PIPE)
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
            self.assertEqual(len(val[:columns]), len(oval[:columns]))
            for i in range(columns):
                delta = delta_pct * int(val[i])
                self.assertGreaterEqual(int(oval[i]), int(val[i]) - delta)
                self.assertLessEqual(int(oval[i]), int(val[i]) + delta)

    def testEvolverLocal(self):
        """ Check that the output of halStats on a hal file produced by running cactus with --binariesMode local is
        is reasonable
        """
        # run cactus directly, the old school way
        name = "local"
        self._run_evolver(name)

        # check the output
        self._check_stats(self._out_hal(name), delta_pct=0.25)
        self._check_coverage(self._out_hal(name), delta_pct=0.20)

    def testEvolverPrepareWDL(self):

        # run cactus step by step via a WDL workflow made by cactus-prepare
        self._run_evolver_decomposed_wdl("wdl")

        # check the output
        self._check_stats(self._out_hal("wdl"), delta_pct=0.25)
        self._check_coverage(self._out_hal("wdl"), delta_pct=0.20)

    def testEvolverPrepareToil(self):

        # run cactus step by step via toil in toil
        name = "toil-in-toil"
        self._run_evolver_decomposed_toil(name, "docker")

        # check the output
        self._check_stats(self._out_hal(name), delta_pct=0.25)
        self._check_coverage(self._out_hal(name), delta_pct=0.20)                

    def testEvolverDecomposedLocal(self):
        """ Check that the output of halStats on a hal file produced by running cactus with --binariesMode local is
        is reasonable
        """
        # run cactus step by step via the plan made by cactus-prepare
        name = "decomposed"
        self._run_evolver_decomposed(name)

        # check the output
        self._check_stats(self._out_hal(name), delta_pct=0.25)
        self._check_coverage(self._out_hal(name), delta_pct=0.20)

    def testEvolverDocker(self):
        """ Check that the output of halStats on a hal file produced by running cactus with --binariesMode docker is
        is reasonable.  Note: the local image being tested should be set up via CACTUS_DOCKER_ORG (with tag==latest)
        """
        # run cactus
        self._run_evolver("docker")

        # check the output
        self._check_stats(self._out_hal("docker"), delta_pct=0.25)
        self._check_coverage(self._out_hal("docker"), delta_pct=0.20)

    def testEvolverPrepareNoOutgroupDocker(self):

        # run cactus step by step via the plan made by cactus-prepare, hacking to apply --nonBlastInput option to cactus-align
        self._run_evolver_decomposed_no_outgroup("docker")

        # check the output
        self._check_stats(self._out_hal("docker"), delta_pct=2.5, subset=['simMouse_chr6', 'simRat_chr6', 'Anc0'])
        self._check_coverage(self._out_hal("docker"), delta_pct=0.20, subset=['simMouse_chr6', 'simRat_chr6', 'Anc0'], columns=1)

    def testEvolverPrepareNoOutgroupLocal(self):

        # run cactus step by step via the plan made by cactus-prepare, hacking to apply --nonBlastInput option to cactus-align
        self._run_evolver_decomposed_no_outgroup("local")

        # check the output
        self._check_stats(self._out_hal("local"), delta_pct=2.5, subset=['simMouse_chr6', 'simRat_chr6', 'Anc0'])
        self._check_coverage(self._out_hal("local"), delta_pct=0.20, subset=['simMouse_chr6', 'simRat_chr6', 'Anc0'], columns=1)


if __name__ == '__main__':
    unittest.main()
