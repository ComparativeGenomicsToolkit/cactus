import os, sys
import pytest
import unittest
import subprocess
import shutil
from sonLib.bioio import getTempDirectory, popenCatch, popen
import xml.etree.ElementTree as ET
from xml.dom import minidom

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

    def _run_evolver(self, binariesMode, configFile = None, seqFile = './examples/evolverMammals.txt'):
        """ Run the full evolver test, putting the jobstore and output in tempDir
        """
        cmd = ['cactus', self._job_store(binariesMode), seqFile, self._out_hal(binariesMode),
                        '--binariesMode', binariesMode, '--logInfo', '--realTimeLogging', '--workDir', self.tempDir]
        if configFile:
            cmd += ['--configFile', configFile]

        # todo: it'd be nice to have an interface for setting tag to something not latest or commit
        if binariesMode == 'docker':
            cmd += ['--latest']

        sys.stderr.write('Running {}'.format(' '.format(cmd)))
        subprocess.check_call(' '.join(cmd), shell=True)

    def _run_evolver_in_docker(self, seqFile = './examples/evolverMammals.txt'):

        out_hal = self._out_hal('in_docker')
        cmd = ['docker', 'run', '--rm', '-v', '{}:{}'.format(os.path.dirname(out_hal), '/data'),
               '-v', '{}:{}'.format(os.getcwd(), '/workdir'), 'evolvertestdocker/cactus:latest',
               'cactus /data/js /workdir/{} /data/{}'.format(seqFile, os.path.basename(out_hal))]
        sys.stderr.write('Running {}'.format(' '.format(cmd)))
        subprocess.check_call(' '.join(cmd), shell=True)
               
    def _write_primates_seqfile(self, seq_file_path):
        """ create the primates seqfile at given path"""
        with open(seq_file_path, 'w') as seq_file:
            seq_file.write('simHuman\thttps://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/primates/loci1/simHuman.chr6\n')            
            seq_file.write('simChimp\thttps://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/primates/loci1/simChimp.chr6\n')
            seq_file.write('simGorilla\thttps://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/primates/loci1/simGorilla.chr6\n')
            seq_file.write('simOrang\thttps://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/primates/loci1/simOrang.chr6\n')

    def _run_evolver_primates_star(self, binariesMode, configFile = None):
        """ Run cactus on the evolver primates with a star topology
        """
        seq_file_path = os.path.join(self.tempDir, 'primates.txt')
        self._write_primates_seqfile(seq_file_path)
        self._run_evolver(binariesMode, configFile=configFile, seqFile=seq_file_path)

    def _run_evolver_decomposed(self, binariesMode, name):
        """ Run the full evolver test, putting the jobstore and output in tempDir
        but instead of doing in in one shot like above, use cactus-prepare, cactus-blast
        and cactus-align to break it into different steps """

        out_dir = os.path.join(self.tempDir, 'output')
        out_seqfile = os.path.join(out_dir, 'evolverMammalsOut.txt')
        in_seqfile = './examples/evolverMammals.txt'
        cmd = ['cactus-prepare', in_seqfile, '--outDir', out_dir, '--outSeqFile', out_seqfile, '--outHal', self._out_hal(name),
               '--jobStore', self._job_store(name)]
        bm_flag = '--binariesMode {}'.format(binariesMode)
        if binariesMode == "docker":
            bm_flag += ' --latest'
        cmd += ['--cactusOptions', '\"{}\"'.format(bm_flag)]

        job_plan = popenCatch(' '.join(cmd))

        for line in job_plan.split('\n'):
            line = line.strip()
            if len(line) > 0 and not line.startswith('#'):
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
               '--halAppendBatchSize', '2',
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
            inseq.write('simMouse_chr6  https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simMouse.chr6\n')
            inseq.write('simRat_chr6 https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simRat.chr6\n')

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
                sys.stderr.write('Running {}'.format(line))
                subprocess.check_call(line, shell=True)

    def _run_evolver_primates_refmap(self, binariesMode):
        """ primates star test but using refmap pangenome pipeline
        """
        # borrow seqfile from other primates test
        # todo: make a seqfile and add it to the repo
        seq_file_path = os.path.join(self.tempDir, 'primates.txt')
        self._write_primates_seqfile(seq_file_path)

        # do the mapping
        cigar_path = os.path.join(self.tempDir, 'aln.cigar')
        cactus_opts = ['--binariesMode', binariesMode, '--logInfo', '--realTimeLogging', '--workDir', self.tempDir]
        # todo: it'd be nice to have an interface for setting tag to something not latest or commit
        if binariesMode == 'docker':
            cactus_opts += ['--latest']
        
        subprocess.check_call(['cactus-refmap', self._job_store(binariesMode), seq_file_path, "simHuman", cigar_path] + cactus_opts)

        # do the alignment
        subprocess.check_call(['cactus-align', self._job_store(binariesMode), seq_file_path, cigar_path, self._out_hal(binariesMode),
                               '--pangenome'] + cactus_opts)            
                

    def _run_evolver_primates_graphmap(self, binariesMode):
        """ primates star test but using graphmap pangenome pipeline
        """
        # borrow seqfile from other primates test
        # todo: make a seqfile and add it to the repo
        seq_file_path = os.path.join(self.tempDir, 'primates.txt')
        self._write_primates_seqfile(seq_file_path)

        # use the same logic cactus does to get default config
        config_path = 'src/cactus/cactus_progressive_config.xml'
        
        xml_root = ET.parse(config_path).getroot()
        bar_elem = xml_root.find("bar")
        bar_elem.attrib["partialOrderAlignment"] = "1"
        bar_elem.attrib["partialOrderAlignmentMaskFilter"] = "1"

        poa_config_path = os.path.join(self.tempDir, "config.poa.xml")
        with open(poa_config_path, 'w') as poa_config_file:
            xmlString = ET.tostring(xml_root, encoding='unicode')
            xmlString = minidom.parseString(xmlString).toprettyxml()
            poa_config_file.write(xmlString)

        # make an output seqfile for preprocessed sequences
        out_seq_file_path = os.path.join(self.tempDir, 'prepared', 'primates.txt')
        subprocess.check_call(['cactus-prepare', seq_file_path, '--outDir', os.path.dirname(out_seq_file_path)])

        cactus_opts = ['--binariesMode', binariesMode, '--logInfo', '--realTimeLogging', '--workDir', self.tempDir, '--configFile', poa_config_path]
        
        # preprocess with dna-brnn
        subprocess.check_call(['cactus-preprocess', self._job_store(binariesMode), seq_file_path, out_seq_file_path, '--maskAlpha'] + cactus_opts)

        # make the reference graph
        mg_cmd = ['minigraph', '-xggs', '-t', '4']
        with open(seq_file_path, 'r') as seq_file:
            for line in seq_file:
                toks = line.strip().split()
                if len(toks) == 2:
                    name = os.path.basename(toks[1])
                    subprocess.check_call(['wget', toks[1]], cwd=self.tempDir)
                    mg_cmd += [os.path.join(self.tempDir, name)]
        mg_path = os.path.join(self.tempDir, 'refgraph.gfa')
        with open(mg_path, 'w') as mg_file:
            subprocess.check_call(mg_cmd, stdout=mg_file)

        # do the mapping
        paf_path = os.path.join(self.tempDir, 'aln.paf')
        fa_path = os.path.join(self.tempDir, 'refgraph.gfa.fa')
        # todo: it'd be nice to have an interface for setting tag to something not latest or commit
        if binariesMode == 'docker':
            cactus_opts += ['--latest']
        
        subprocess.check_call(['cactus-graphmap', self._job_store(binariesMode), out_seq_file_path, mg_path, paf_path,
                               '--outputFasta', fa_path] + cactus_opts)

        # do the alignment
        subprocess.check_call(['cactus-align', self._job_store(binariesMode), out_seq_file_path, paf_path, self._out_hal(binariesMode),
                               '--pangenome', '--outVG', '--outGFA', '--barMaskFilter', '10000'] + cactus_opts) 

    def _run_yeast_pangenome(self, binariesMode):
        """ yeast pangenome chromosome by chromosome pipeline
        """

        # make a config that replaces stuff like SC4.chrI with id=SC4|chrI
        # also run dna brnn with softmasking
        orig_seq_file_path = './examples/yeastPangenome.txt'
        seq_file_path = os.path.join(self.tempDir, 'pp', os.path.basename(orig_seq_file_path))        
        subprocess.check_call(['cactus-prepare',  orig_seq_file_path, '--outDir', os.path.join(self.tempDir, 'pp')])
        seq_file_fix_path = seq_file_path + '.fix'
        with open(seq_file_fix_path, 'w') as seq_file_fix_file:
            subprocess.check_call(['grep', '-v', 'Anc', seq_file_path], stdout=seq_file_fix_file)
        seq_file_path = seq_file_fix_path
        orig_config_path = 'src/cactus/cactus_progressive_config.xml'
        config_path = os.path.join(self.tempDir, 'config.xml')
        subprocess.check_call('cat {} | sed -e \'s/cutBefore=""/cutBefore="."/g\' | grep -v lastzRepeatMask > {}'.format(
            orig_config_path, config_path), shell=True)
        cactus_opts = ['--binariesMode', binariesMode, '--logInfo', '--realTimeLogging', '--workDir', self.tempDir, '--configFile', config_path, '--maxCores', '8']
        subprocess.check_call(['cactus-preprocess', self._job_store(binariesMode),
                               orig_seq_file_path, seq_file_path, '--maskAlpha', '--minLen', '20000'] + cactus_opts, shell=False)
                               
        # build the minigraph
        mg_cmd = ['minigraph', '-xggs', '-t', '4']
        events = []
        with open(seq_file_path, 'r') as seq_file:
            for line in seq_file:
                toks = line.strip().split()
                if len(toks) == 2:
                    mg_cmd += [toks[1]]
                    events.append(toks[0])
        mg_path = os.path.join(self.tempDir, 'yeast.gfa')
        with open(mg_path, 'w') as mg_file:
            subprocess.check_call(mg_cmd, stdout=mg_file)
        
        # run graphmap in base mode
        paf_path = os.path.join(self.tempDir, 'yeast.paf')
        fa_path = os.path.join(self.tempDir, 'yeast.gfa.fa')        
        subprocess.check_call(['cactus-graphmap', self._job_store(binariesMode), seq_file_path, mg_path, paf_path,
                               '--outputFasta', fa_path, '--base', '--delFilter', '2000000'] + cactus_opts)
            
        # split into chromosomes
        chroms = ['chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII', 'chrIX', 'chrX', 'chrXI', 'chrXIV', 'chrXV']
        gm_out_dir = os.path.join(self.tempDir, 'chroms')
        subprocess.check_call(['cactus-graphmap-split', self._job_store(binariesMode), seq_file_path, mg_path, paf_path,
                               '--refContigs'] + chroms +
                              ['--reference', 'S288C', '--maskFilter', '20000', '--otherContig', 'chrOther', '--outDir', gm_out_dir] + cactus_opts)

        # create the vg graph for each chromosome
        chromfile_path = os.path.join(gm_out_dir, 'chromfile.txt')
        ba_path = os.path.join(self.tempDir, 'batch')
        os.makedirs(ba_path, exist_ok=True)
        subprocess.check_call(['cactus-align-batch', self._job_store(binariesMode), chromfile_path, ba_path, '--alignCores', '1',
                               '--alignOptions', '--pangenome --outVG --barMaskFilter 20000 --realTimeLogging --reference S288C --binariesMode {}'.format(binariesMode)])

        vg_files = [os.path.join(ba_path, c) + '.vg' for c in chroms]
        hal_files = [os.path.join(ba_path, c) + '.hal' for c in chroms]

        # join up the graphs and index for giraffe
        join_path = os.path.join(self.tempDir, 'join')
        # must specify haploid for everything
        rename_opts = ['--rename']
        for event in events:
            rename_opts += ['{}>{}.0'.format(event, event)]
        subprocess.check_call(['cactus-graphmap-join', self._job_store(binariesMode), '--outDir', join_path, '--outName', 'yeast',
                               '--reference', 'S288C', '--vg'] +  vg_files + ['--hal'] + hal_files + ['--gfaffix', '--wlineSep', '.',
                               '--vcf', '--giraffe'] + rename_opts + cactus_opts)

    def _check_yeast_pangenome(self, binariesMode):
        """ yeast pangenome chromosome by chromosome pipeline
        """

        # load up the events
        orig_seq_file_path = './examples/yeastPangenome.txt'
        events = []
        with open(orig_seq_file_path, 'r') as orig_seq_file:
            for line in orig_seq_file:
                toks = line.strip().split()
                if len(toks) == 2:
                    events += [toks[0]]
        assert len(events) == 6

        join_path = os.path.join(self.tempDir, 'join')
        vcf_path = os.path.join(join_path, 'yeast.vcf.gz')

        # check that we have some alts for each sample in the VCF
        vcf_allele_threshold = 40000
        for event in events:
            allele = 0 if event == "S288C" else 1
            proc = subprocess.Popen('bcftools view {} -s {} -a -H | awk \'{{print $10}}\' | grep {} | wc -l'.format(vcf_path, event, allele),
                                    shell=True, stdout=subprocess.PIPE)
            output, errors = proc.communicate()
            sts = proc.wait()
            num_alleles = int(output.strip())
            self.assertGreaterEqual(num_alleles, vcf_allele_threshold)

        # make sure we have about the right sequence counts in the hal
        hal_path = os.path.join(join_path, 'yeast.hal')
        for event in events:
            if event != "_MINIGRAPH_":
                proc = subprocess.Popen('halStats {} --sequenceStats {} | wc -l'.format(hal_path, event), shell=True, stdout=subprocess.PIPE)
                output, errors = proc.communicate()
                sts = proc.wait()
                num_sequences = int(output.strip())
                self.assertGreaterEqual(num_sequences, 11)
                self.assertLessEqual(num_sequences, 15)

        # make sure the vg is sane
        xg_path = os.path.join(join_path, 'yeast.xg')
        proc = subprocess.Popen('vg stats -l {} | awk \'{{print $2}}\''.format(xg_path), shell=True, stdout=subprocess.PIPE)
        output, errors = proc.communicate()
        sts = proc.wait()
        num_bases = int(output.strip())
        self.assertGreaterEqual(num_bases, 10000000)
        self.assertLessEqual(num_bases, 11200000)

        # make sure we have the giraffe indexes
        for giraffe_idx in ['yeast.snarls', 'yeast.gg', 'yeast.dist', 'yeast.min', 'yeast.trans.gz', 'yeast.gbwt']:
            idx_bytes = os.path.getsize(os.path.join(join_path, giraffe_idx))
            self.assertGreaterEqual(idx_bytes, 500000)

        # make sure the chrom splitting stats are somewhat sane
        contig_sizes = {}
        with open(os.path.join(self.tempDir, 'chroms', 'contig_sizes.tsv'), 'r') as sizes_file:
            for line in sizes_file:
                if line.startswith('chr'):
                    toks=line.strip().split()
                    contig_sizes[toks[0]] = toks[1:]
        self.assertEqual(len(contig_sizes), 14)
        for chr,sizes in contig_sizes.items():
            self.assertEqual(len(sizes), 10)
            self.assertGreaterEqual(int(sizes[9]), 200000)
                
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

    def _check_maf_accuracy(self, halPath, delta, dataset="mammals"):
        """ Compare mafComparator output of evolver mammals to baseline
        """
        assert dataset in ('mammals', 'primates')
        # this is just pasted from a successful run.  it will be used to catch serious regressions
        baseline_file = 'test/evolver{}-default.comp.xml'.format(dataset.capitalize())
        # this is downloaded in the Makefile
        ground_truth_file = 'test/{}-truth.maf'.format(dataset)

        # run mafComparator on the evolver output
        subprocess.check_call(['bin/hal2maf', halPath,  halPath + '.maf', '--onlySequenceNames'], shell=False)
        subprocess.check_call(['bin/mafComparator', '--maf1', halPath + '.maf', '--maf2', ground_truth_file, '--samples', '100000000', '--out', halPath + 'comp.xml'])

        # grab the two accuracy values out of the XML
        def parse_mafcomp_output(xml_path):
            xml_root = ET.parse(xml_path).getroot()
            comp_roots = xml_root.findall("homologyTests")
            assert len(comp_roots) == 2            
            first, second = None, None
            for comp_root in comp_roots:
                avg_val = float(comp_root.find("aggregateResults").find("all").attrib["average"])
                if os.path.basename(comp_root.attrib["fileB"]) == os.path.basename(ground_truth_file):
                    first = avg_val
                else:
                    second = avg_val
            assert first is not None and second is not None
            return first, second

        baseline_acc = parse_mafcomp_output(baseline_file)
        acc = parse_mafcomp_output(halPath + 'comp.xml')

        sys.stderr.write("Comparing mafcomp accuracy {},{} to baseline accuracy {},{} with threshold {}\n".format(acc[0], acc[1], baseline_acc[0], baseline_acc[1], delta))

        self.assertGreaterEqual(acc[0] + delta, baseline_acc[0])
        self.assertGreaterEqual(acc[1] + delta, baseline_acc[1])

    def testEvolverLocal(self):
        """ Check that the output of halStats on a hal file produced by running cactus with --binariesMode local is
        is reasonable
        """
        # run cactus directly, the old school way
        name = "local"
        self._run_evolver(name)

        # check the output
        #self._check_stats(self._out_hal(name), delta_pct=0.25)
        #self._check_coverage(self._out_hal(name), delta_pct=0.20)
        self._check_maf_accuracy(self._out_hal(name), delta=0.04)

    def testEvolverPrepareWDL(self):

        # run cactus step by step via a WDL workflow made by cactus-prepare
        self._run_evolver_decomposed_wdl("wdl")

        # check the output
        #self._check_stats(self._out_hal("wdl"), delta_pct=0.25)
        #self._check_coverage(self._out_hal("wdl"), delta_pct=0.20)
        self._check_maf_accuracy(self._out_hal("wdl"), delta=0.04)

    def testEvolverPrepareToil(self):

        # run cactus step by step via toil in toil
        name = "toil-in-toil"
        self._run_evolver_decomposed_toil(name, "docker")

        # check the output
        #self._check_stats(self._out_hal(name), delta_pct=0.25)
        #self._check_coverage(self._out_hal(name), delta_pct=0.20)
        self._check_maf_accuracy(self._out_hal(name), delta=0.04)

    def testEvolverDecomposedLocal(self):
        """ Check that the output of halStats on a hal file produced by running cactus with --binariesMode local is
        is reasonable
        """
        # run cactus step by step via the plan made by cactus-prepare
        name = "decomposed"
        self._run_evolver_decomposed("local", name)

        # check the output
        #self._check_stats(self._out_hal(name), delta_pct=0.25)
        #self._check_coverage(self._out_hal(name), delta_pct=0.20)
        self._check_maf_accuracy(self._out_hal(name), delta=0.04)

    def testEvolverDecomposedDocker(self):
        """ Check that the output of halStats on a hal file produced by running cactus with --binariesMode docker is
        is reasonable
        """
        # run cactus step by step via the plan made by cactus-prepare
        name = "decomposed"
        self._run_evolver_decomposed("docker", name)

        # check the output
        #self._check_stats(self._out_hal(name), delta_pct=0.25)
        #self._check_coverage(self._out_hal(name), delta_pct=0.20)
        self._check_maf_accuracy(self._out_hal(name), delta=0.04)
        
    def testEvolverDocker(self):
        """ Check that the output of halStats on a hal file produced by running cactus with --binariesMode docker is
        is reasonable.  Note: the local image being tested should be set up via CACTUS_DOCKER_ORG (with tag==latest)
        """
        # run cactus
        self._run_evolver("docker")

        # check the output
        #self._check_stats(self._out_hal("docker"), delta_pct=0.25)
        #self._check_coverage(self._out_hal("docker"), delta_pct=0.20)
        self._check_maf_accuracy(self._out_hal("docker"), delta=0.04)

    def testEvolverInDocker(self):
        """ Check that the output of halStats on a hal file produced by running cactus in docker
        is reasonable.  Note: the local image being tested should be set up via CACTUS_DOCKER_ORG (with tag==latest)
        """
        # run cactus
        self._run_evolver_in_docker()

        # check the output
        #self._check_stats(self._out_hal("docker"), delta_pct=0.25)
        #self._check_coverage(self._out_hal("docker"), delta_pct=0.20)
        self._check_maf_accuracy(self._out_hal("in_docker"), delta=0.04)
        
    def testEvolverPrepareNoOutgroupDocker(self):

        # run cactus step by step via the plan made by cactus-prepare
        self._run_evolver_decomposed_no_outgroup("docker")

        # check the output
        self._check_stats(self._out_hal("docker"), delta_pct=2.5, subset=['simMouse_chr6', 'simRat_chr6', 'Anc0'])
        self._check_coverage(self._out_hal("docker"), delta_pct=0.20, subset=['simMouse_chr6', 'simRat_chr6', 'Anc0'], columns=1)

    def testEvolverPrepareNoOutgroupLocal(self):

        # run cactus step by step via the plan made by cactus-prepare
        self._run_evolver_decomposed_no_outgroup("local")

        # check the output
        self._check_stats(self._out_hal("local"), delta_pct=2.5, subset=['simMouse_chr6', 'simRat_chr6', 'Anc0'])
        self._check_coverage(self._out_hal("local"), delta_pct=0.20, subset=['simMouse_chr6', 'simRat_chr6', 'Anc0'], columns=1)

    def testEvolverPOALocal(self):
        """ Check that the output of halStats on a hal file produced by running cactus with --binariesMode local is
        is reasonable... when using POA mode in BAR.
        """
        # use the same logic cactus does to get default config
        config_path = 'src/cactus/cactus_progressive_config.xml'
        
        xml_root = ET.parse(config_path).getroot()
        bar_elem = xml_root.find("bar")
        bar_elem.attrib["partialOrderAlignment"] = "1"

        poa_config_path = os.path.join(self.tempDir, "config.poa.xml")
        with open(poa_config_path, 'w') as poa_config_file:
            xmlString = ET.tostring(xml_root, encoding='unicode')
            xmlString = minidom.parseString(xmlString).toprettyxml()
            poa_config_file.write(xmlString)

        # run cactus directly, the old school way
        name = "local"
        self._run_evolver_primates_star(name, configFile = poa_config_path)

        # check the output
        self._check_maf_accuracy(self._out_hal("local"), delta=0.0025, dataset='primates')

    def testEvolverRefmapLocal(self):
        """ Use the new minimap pangenome pipeline to create an alignment of the primates, then compare with the baseline
        """
        name = "local"
        self._run_evolver_primates_refmap(name)

        # check the output
        self._check_maf_accuracy(self._out_hal("local"), delta=0.01, dataset='primates')

    def testEvolverMinigraphLocal(self):
        """ Use the new minigraph pangenome pipeline to create an alignment of the primates, then compare with the baseline
        """
        name = "local"
        self._run_evolver_primates_graphmap(name)

        # check the output
        # todo: tune config so that delta can be reduced
        self._check_maf_accuracy(self._out_hal("local"), delta=0.025, dataset='primates')

    def testYeastPangenomeLocal(self):
        """ Run pangenome pipeline (including contig splitting!) on yeast dataset """
        name = "local"
        self._run_yeast_pangenome(name)
        
        # check the output
        self._check_yeast_pangenome(name)

if __name__ == '__main__':
    unittest.main()
