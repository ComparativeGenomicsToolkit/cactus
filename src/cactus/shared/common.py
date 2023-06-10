#!/usr/bin/env python3
#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
"""Wrapper functions for assisting in running the various programs of the cactus package.
"""

import os
import sys
import shutil
import subprocess
import logging
import pathlib
import pipes
import uuid
import json
import time
import signal
import hashlib
import tempfile
import math
import threading
import traceback
import errno
import shlex

try:
    import boto3
    import botocore
    has_s3 = True
except:
    has_s3 = False

from urllib.parse import urlparse
from datetime import datetime

from toil.statsAndLogging import logger
from toil.lib.bioio import system
from toil.lib.bioio import getLogLevelString
from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger
from toil.lib.humanize import bytes2human
from toil.lib.threading import cpu_count
from sonLib.bioio import popenCatch
from sonLib.bioio import getTempDirectory

from cactus.shared.version import cactus_commit
    

_log = logging.getLogger(__name__)

subprocess._has_poll = False

def cactus_cpu_count():
    """ try the more cluster-friendly cpu counter before reverting to toil's
    https://github.com/ComparativeGenomicsToolkit/cactus/issues/820
    """
    num_cpus = cpu_count()
    try:
        sched_cpus = len(os.sched_getaffinity(0))
        if sched_cpus and sched_cpus < num_cpus:
            num_cpus = sched_cpus
    except:
        pass
    return num_cpus

def cactus_override_toil_options(options):
    """  Mess with some toil options to create useful defaults. """
    if options.retryCount is None and options.batchSystem.lower() not in ['single_machine', 'singleMachine']:
        # If the user didn't specify a retryCount value, make it 5
        # instead of Toil's default (1).
        options.retryCount = 5

    if options.batchSystem.lower() in ['slurm', 'lsf', 'torque']:
        # disable caching for cluster style batch systems as it seems to
        # lead to weird toil errors?
        # https://github.com/DataBiosphere/toil/issues/4218
        options.disableCaching = True
    
    if not options.realTimeLogging:
        # Too much valuable debugging information to pass up
        logger.info('Enabling realtime logging in Toil')
        options.realTimeLogging = True

def makeURL(path_or_url):
    if urlparse(path_or_url).scheme == '':
        return "file://" + os.path.abspath(path_or_url)
    else:
        return path_or_url

def catFiles(filesToCat, catFile):
    """Cats a bunch of files into one file. Ensures a no more than maxCat files
    are concatenated at each step.
    """
    if len(filesToCat) == 0: #We must handle this case or the cat call will hang waiting for input
        open(catFile, 'w').close()
        return
    maxCat = 25
    system("cat %s > %s" % (" ".join(filesToCat[:maxCat]), catFile))
    filesToCat = filesToCat[maxCat:]
    while len(filesToCat) > 0:
        system("cat %s >> %s" % (" ".join(filesToCat[:maxCat]), catFile))
        filesToCat = filesToCat[maxCat:]

def cactusRootPath():
    """
    function for finding external location
    """
    import cactus
    i = os.path.abspath(cactus.__file__)
    return os.path.split(i)[0]

def getLogLevelString2(logLevelString):
    """Gets the log level string for the binary
    """
    if logLevelString == None:
        return getLogLevelString()
    return logLevelString

def getOptionalAttrib(node, attribName, typeFn=None, default=None, errorIfNotPresent=False):
    """Get an optional attrib, or default if not set or node is None
    """
    if node != None and attribName in node.attrib:
        if typeFn != None:
            if typeFn == bool:
                aname = node.attrib[attribName].lower()
                if aname == 'false':
                    return False
                elif aname == 'true':
                    return True
                else:
                    return bool(int(node.attrib[attribName]))
            return typeFn(node.attrib[attribName])
        return node.attrib[attribName]
    if errorIfNotPresent:
        raise RuntimeError("Could not find attribute %s in %s node" % (attribName, node))
    return default

def findRequiredNode(configNode, nodeName):
    """Retrieve an xml node, complain if it's not there."""
    nodes = configNode.findall(nodeName)
    if nodes == None:
        raise RuntimeError("Could not find any nodes with name %s in %s node" % (nodeName, configNode))
    assert len(nodes) == 1, "More than 1 node for %s in config XML" % nodeName
    return nodes[0]

#############################################
#############################################
#All the following provide command line wrappers
#for core programs in the cactus pipeline.
#############################################
#############################################


def _fn(toilDir,
      logLevel=None, retryCount=0,
      batchSystem="single_machine",
      rescueJobFrequency=None,
      buildAvgs=False,
      buildHal=False,
      buildFasta=False,
      toilStats=False,
      maxThreads=None,
      maxCpus=None,
      defaultMemory=None,
      logFile=None):
    logLevel = getLogLevelString2(logLevel)
    args = [toilDir, "--logLevel", logLevel]
    if buildAvgs:
        args += ["--buildAvgs"]
    if buildHal:
        args += ["--buildHal"]
    if buildFasta:
        args += ["--buildFasta"]
    #Jobtree args
    if batchSystem is not None:
        args += ["--batchSystem", batchSystem]
    if retryCount is not None:
        args += ["--retryCount", str(retryCount)]
    if rescueJobFrequency is not None:
        args += ["--rescueJobFrequency", str(rescueJobFrequency)]
    if toilStats:
        args += ["--stats"]
    if maxThreads is not None:
        args += ["--maxThreads", str(maxThreads)]
    if maxCpus is not None:
        args += ["--maxCpus", str(maxCpus)]
    if defaultMemory is not None:
        args += ["--defaultMemory", str(defaultMemory)]
    if logFile is not None:
        args += ["--logFile", logFile]
    return args

def runCactusWorkflow(experimentFile,
                      toilDir,
                      logLevel=None, retryCount=0,
                      batchSystem="single_machine",
                      rescueJobFrequency=None,
                      skipAlignments=False,
                      buildAvgs=False,
                      buildHal=False,
                      buildFasta=False,
                      toilStats=False,
                      maxThreads=None,
                      maxCpus=None,
                      defaultMemory=None,
                      logFile=None,
                      intermediateResultsUrl=None,
                      extraToilArgumentsString=""):
    args = ["--experiment", experimentFile] + _fn(toilDir,
                      logLevel, retryCount, batchSystem, rescueJobFrequency,
                      buildAvgs, buildHal, buildFasta, toilStats, maxThreads, maxCpus, defaultMemory, logFile)
    if intermediateResultsUrl is not None:
        args += ["--intermediateResultsUrl", intermediateResultsUrl]

    import cactus.pipeline.cactus_workflow as cactus_workflow
    cactus_workflow.runCactusWorkflow(args)
    logger.info("Ran the cactus workflow okay")

def runCactusProgressive(seqFile,
                         configFile,
                         toilDir,
                         logLevel=None, retryCount=0,
                         batchSystem="single_machine",
                         rescueJobFrequency=None,
                         skipAlignments=False,
                         buildHal=True,
                         buildAvgs=False,
                         toilStats=False,
                         maxCpus=None):
    opts = Job.Runner.getDefaultOptions(toilDir)
    opts.batchSystem = batchSystem if batchSystem is not None else opts.batchSystem
    opts.logLevel = logLevel if logLevel is not None else opts.logLevel
    opts.maxCores = maxCpus if maxCpus is not None else opts.maxCores
    # Used for tests
    opts.scale = 0.1
    opts.retryCount = retryCount if retryCount is not None else opts.retryCount
    # This *shouldn't* be necessary, but it looks like the toil
    # deadlock-detection still has issues.
    opts.deadlockWait = 3600

    opts.buildHal = buildHal
    opts.buildAvgs = buildAvgs
    opts.buildFasta = True
    if toilStats:
        opts.stats = True
    opts.seqFile = seqFile
    opts.configFile = configFile
    opts.database = 'kyoto_tycoon'
    opts.root = None
    opts.outputHal = '/dev/null'
    opts.intermediateResultsUrl = None
    from cactus.progressive.cactus_progressive import runCactusProgressive as runRealCactusProgressive
    runRealCactusProgressive(opts)

def runToilStats(toil, outputFile):
    system("toil stats %s --outputFile %s" % (toil, outputFile))
    logger.info("Ran the job-tree stats command apparently okay")

def runGetChunks(sequenceFiles, chunksDir, chunkSize, overlapSize, work_dir=None):
    chunks = cactus_call(work_dir=work_dir,
                         check_output=True,
                         parameters=["faffy", "chunk",
                                     "--logLevel", getLogLevelString(),
                                     "--chunkSize", str(chunkSize),
                                     "--overlap", str(overlapSize),
                                     "--dir", chunksDir] + sequenceFiles)
    return [chunk for chunk in chunks.split("\n") if chunk != ""]

def pullCactusImage():
    """Ensure that the cactus Docker image is pulled."""
    if os.environ.get('CACTUS_DOCKER_MODE') == "0":
        return
    if os.environ.get('CACTUS_USE_LOCAL_IMAGE', 0) == "1":
        return
    image = getDockerImage()
    call = ["docker", "pull", image]
    process = subprocess.Popen(call, stdout=subprocess.PIPE,
                                 stderr=sys.stderr, bufsize=-1)
    output, _ = process.communicate()
    if process.returncode != 0:
        raise RuntimeError("Command %s failed with output: %s" % (call, output))

def getDockerOrg():
    """Get where we should find the cactus containers."""
    if "CACTUS_DOCKER_ORG" in os.environ:
        return os.environ["CACTUS_DOCKER_ORG"]
    else:
        return "quay.io/comparative-genomics-toolkit"

def getDockerTag():
    """Get what docker tag we should use for the cactus image
    (either forced to be latest or the current cactus commit)."""
    if 'CACTUS_DOCKER_TAG' in os.environ:
        return os.environ['CACTUS_DOCKER_TAG']
    elif 'CACTUS_USE_LATEST' in os.environ:
        return "latest"    
    else:
        return cactus_commit

def getDockerImage():
    """Get fully specified Docker image name."""
    return "%s/cactus:%s" % (getDockerOrg(), getDockerTag())

def getDockerRelease(gpu=False):
    """Get the most recent docker release."""
    r = "quay.io/comparative-genomics-toolkit/cactus:v2.5.2"
    if gpu:
        r += "-gpu"
    return r

def maxMemUsageOfContainer(containerInfo):
    """Return the max RSS usage (in bytes) of a container, or None if something failed."""
    if containerInfo['id'] is None:
        # Try to get the internal container ID from the docker name
        try:
            id = popenCatch("docker inspect -f '{{.Id}}' %s" % containerInfo['name']).strip()
            containerInfo['id'] = id
        except:
            # Not yet running
            return None
    # Try to check for the maximum memory usage ever used by that
    # container, in a few different possible locations depending on
    # the distribution
    possibleLocations = ["/sys/fs/cgroup/memory/docker/%s/memory.max_usage_in_bytes",
                         "/sys/fs/cgroup/memory/system.slice.docker-%s.scope/memory.max_usage_in_bytes"]
    possibleLocations = [s % containerInfo['id'] for s in possibleLocations]
    for location in possibleLocations:
        try:
            with open(location) as f:
                return int(f.read())
        except IOError:
            # Not at this location, or sysfs isn't mounted
            continue
    return None

# send a time/date stamped message to the realtime logger, truncating it
# if it's too long (so it's less likely to be dropped)
def cactus_realtime_log(msg, max_len = 1500, log_debug=False):
    if len(msg) > max_len:
        msg = msg[:max_len-207] + " <...> " + msg[-200:]
    if not log_debug:
        RealtimeLogger.info("{}: {}".format(datetime.now(), msg))
    else:
        RealtimeLogger.debug("{}: {}".format(datetime.now(), msg))
        
def setupBinaries(options):
    """Ensure that Cactus's C/C++ components are ready to run, and set up the environment."""
    if options.latest:
        os.environ["CACTUS_USE_LATEST"] = "1"
    if options.binariesMode is not None:
        # Mode is specified on command line
        mode = options.binariesMode
    else:
        # Might be specified through the environment, or not, in which
        mode = os.environ.get("CACTUS_BINARIES_MODE")

    def verify_docker():
        # If running without Docker, verify that we can find the Cactus executables
        from distutils.spawn import find_executable
        if find_executable('docker') is None:
            raise RuntimeError("The `docker` executable wasn't found on the "
                               "system. Please install Docker if possible, or "
                               "use --binariesMode local and add cactus's bin "
                               "directory to your PATH.")
    def verify_local():
        from distutils.spawn import find_executable
        if find_executable('cactus_consolidated') is None:
            raise RuntimeError("Cactus isn't using Docker, but it can't find "
                               "the Cactus binaries. Please add Cactus's bin "
                               "directory to your PATH (and run `make` in the "
                               "Cactus directory if you haven't already).")

    if mode is None:
        # there is no mode set, we use local if it's available, otherwise default to docker
        try:
            verify_local()
            mode = "local"
        except:
            verify_docker()
            mode = "docker"
    elif mode == "docker":
        verify_docker()
    elif mode == "local":
        verify_local()
    else:
        assert mode == "singularity"
        jobStoreType, locator = Toil.parseLocator(options.jobStore)
        if jobStoreType == "file":
            # if not using a local jobStore, then don't set the `SINGULARITY_CACHEDIR`
            # in this case, the image will be downloaded on each call
            if options.containerImage:
                imgPath = os.path.abspath(options.containerImage)
                os.environ["CACTUS_USE_LOCAL_SINGULARITY_IMG"] = "1"
            else:
                # When SINGULARITY_CACHEDIR is set, singularity will refuse to store images in the current directory
                if 'SINGULARITY_CACHEDIR' in os.environ:
                    imgPath = os.path.join(os.environ['SINGULARITY_CACHEDIR'], "cactus.img")
                else:
                    imgPath = os.path.join(os.path.abspath(locator), "cactus.img")
            os.environ["CACTUS_SINGULARITY_IMG"] = imgPath
            if options.workDir:
                os.environ['SINGULARITY_TMPDIR'] = os.path.abspath(options.workDir)
            elif 'TMPDIR' in os.environ:
                os.environ['SINGULARITY_TMPDIR'] = os.environ['TMPDIR']
            else:
                os.environ['SINGULARITY_TMPDIR'] = os.path.abspath(locator)

    os.environ["CACTUS_BINARIES_MODE"] = mode

    if options.workDir:
        os.environ['TMPDIR'] = os.path.abspath(options.workDir)

def importSingularityImage(options):
    """Import the Singularity image from Docker if using Singularity."""
    mode = os.environ.get("CACTUS_BINARIES_MODE", "docker")
    localImage = os.environ.get("CACTUS_USE_LOCAL_SINGULARITY_IMG", "0")
    if mode == "singularity" and Toil.parseLocator(options.jobStore)[0] == "file":
        imgPath = os.environ["CACTUS_SINGULARITY_IMG"]
        # If not using local image, pull the docker image
        if localImage == "0":
            # Singularity will complain if the image file already exists. Remove it.
            try:
                os.remove(imgPath)
            except OSError:
                # File doesn't exist
                pass
            # Singularity 2.4 broke the functionality that let --name
            # point to a path instead of a name in the CWD. So we change
            # to the proper directory manually, then change back after the
            # image is pulled.
            # NOTE: singularity writes images in the current directory only
            #       when SINGULARITY_CACHEDIR is not set
            oldCWD = os.getcwd()
            os.chdir(os.path.dirname(imgPath))
            # --size is deprecated starting in 2.4, but is needed for 2.3 support. Keeping it in for now.
            try:
                subprocess.check_call(["singularity", "pull", "--size", "2000", "--name", os.path.basename(imgPath),
                                       "docker://" + getDockerImage()], stderr=subprocess.PIPE)
            except subprocess.CalledProcessError:
                # Call failed, try without --size, required for singularity 3+
                subprocess.check_call(["singularity", "pull", "--name", os.path.basename(imgPath),
                                       "docker://" + getDockerImage()])
            os.chdir(oldCWD)
        else:
            logger.info("Using pre-built singularity image: '{}'".format(imgPath))

def singularityCommand(tool=None,
                       work_dir=None,
                       parameters=None,
                       port=None,
                       file_store=None,
                       gpus=None,
                       cpus=None):

    if parameters is None:
        parameters = []
    if work_dir is None:
        work_dir = os.getcwd()

    base_singularity_call = ['singularity', '--silent', 'exec']

    # Mount workdir as /mnt and work in there.
    # Hope the image actually has a /mnt available.
    # Otherwise this silently doesn't mount.
    # But with -u (user namespaces) we have no luck pointing in-container
    # home at anything other than our real home (like something under /var
    # where Toil puts things).
    # Note that we target Singularity 3+.
    base_singularity_call += ['-u', '-B', '{}:{}'.format(os.path.abspath(work_dir), '/mnt'), '--pwd', '/mnt']
    if gpus:
        base_singularity_call += ['--nv']
    #todo: it seems like this would be useful (eg to hopefully limit squashfs resources) but it doesn't work
    #      on our cluster (crytpic cgroups errors)
    #if cpus:
    #    base_singularity_call += ['--cpus', str(cpus)]
    
    if "CACTUS_SINGULARITY_IMG" in os.environ:
        # old logic: just run a local image
        # (this was toggled by only setting CACTUS_SINGULARITY_IMG when using a local jobstore in cactus_progressive.py)
        base_singularity_call += [os.environ['CACTUS_SINGULARITY_IMG']] 
        base_singularity_call.extend(parameters)
        return base_singularity_call
    else:
        # workaround for kubernetes toil: explicitly make a local image
        # (see https://github.com/vgteam/toil-vg/blob/master/src/toil_vg/singularity.py)

        # Problem: Multiple Singularity downloads sharing the same cache directory will
        # not work correctly. See https://github.com/sylabs/singularity/issues/3634
        # and https://github.com/sylabs/singularity/issues/4555.

        # As a workaround, we have out own cache which we manage ourselves.
        home_dir = str(pathlib.Path.home())
        default_singularity_dir = os.path.join(home_dir, '.singularity')
        cache_dir = os.path.join(os.environ.get('SINGULARITY_CACHEDIR',  default_singularity_dir), 'toil')
        os.makedirs(cache_dir, exist_ok=True)

        # hack to transform back to docker image
        if tool == 'cactus':
            tool = getDockerImage()
        # not a url or local file? try it as a Docker specifier
        if not tool.startswith('/') and '://' not in tool:
            tool = 'docker://' + tool

        # What name in the cache dir do we want?
        # We cache everything as sandbox directories and not .sif files because, as
        # laid out in https://github.com/sylabs/singularity/issues/4617, there
        # isn't a way to run from a .sif file and have write permissions on system
        # directories in the container, because the .sif build process makes
        # everything owned by root inside the image. Since some toil-vg containers
        # (like the R one) want to touch system files (to install R packages at
        # runtime), we do it this way to act more like Docker.
        #
        # Also, only sandbox directories work with user namespaces, and only user
        # namespaces work inside unprivileged Docker containers like the Toil
        # appliance.
        sandbox_dirname = os.path.join(cache_dir, '{}.sandbox'.format(hashlib.sha256(tool.encode('utf-8')).hexdigest()))

        if not os.path.exists(sandbox_dirname):
            # We atomically drop the sandbox at that name when we get it

            # Make a temp directory to be the sandbox
            temp_sandbox_dirname = tempfile.mkdtemp(dir=cache_dir)

            # Download with a fresh cache to a sandbox
            download_env = os.environ.copy()
            download_env['SINGULARITY_CACHEDIR'] = file_store.getLocalTempDir() if file_store else tempfile.mkdtemp(dir=work_dir)
            build_cmd = ['singularity', 'build', '-s', '-F', temp_sandbox_dirname, tool]

            cactus_realtime_log("Running the command: \"{}\"".format(' '.join(build_cmd)))
            start_time = time.time()
            subprocess.check_call(build_cmd, env=download_env)
            run_time = time.time() - start_time
            cactus_realtime_log("Successfully ran the command: \"{}\" in {} seconds".format(' '.join(build_cmd), run_time))

            # Clean up the Singularity cache since it is single use
            shutil.rmtree(download_env['SINGULARITY_CACHEDIR'])

            try:
                # This may happen repeatedly but it is atomic
                os.rename(temp_sandbox_dirname, sandbox_dirname)
            except OSError as e:
                if e.errno in (errno.EEXIST, errno.ENOTEMPTY):
                    # Can't rename a directory over another
                    # Make sure someone else has made the directory
                    assert os.path.exists(sandbox_dirname)
                    # Remove our redundant copy
                    shutil.rmtree(temp_sandbox_dirname)
                else:
                    raise

            # TODO: we could save some downloading by having one process download
            # and the others wait, but then we would need a real fnctl locking
            # system here.
        return base_singularity_call + [sandbox_dirname] + parameters


def dockerCommand(tool=None,
                  work_dir=None,
                  parameters=None,
                  rm=True,
                  port=None,
                  dockstore=None,
                  entrypoint=None,
                  gpus=None,
                  cpus=None):
    # This is really dumb, but we have to work around an intersection
    # between two bugs: one in CoreOS where /etc/resolv.conf is
    # sometimes missing temporarily, and one in Docker where it
    # refuses to start without /etc/resolv.conf.
    while not os.path.exists('/etc/resolv.conf'):
        pass

    base_docker_call = ['docker', 'run',
                        '--interactive',
                        '--net=host',
                        '--log-driver=none',
                        '-u', '%s:%s' % (os.getuid(), os.getgid()),
                        '-v', '{}:/data'.format(os.path.abspath(work_dir))]
    if gpus:
        if 'SLURM_JOB_GPUS' in os.environ:
            # this allows slurm to identify which gpus are free
            base_docker_call += ['--gpus', '"device={}"'.format(os.environ['SLURM_JOB_GPUS'])]            
        else:                    
            base_docker_call += ['--gpus', str(gpus)]
    if cpus:
        base_docker_call += ['--cpus', str(cpus)]

    if entrypoint is not None:
        base_docker_call += ['--entrypoint', entrypoint]
    else:
        base_docker_call += ['--entrypoint', '/opt/cactus/wrapper.sh']

    if port is not None:
        base_docker_call += ["-p", "%d:%d" % (port, port)]

    containerInfo = { 'name': str(uuid.uuid4()), 'id': None }
    base_docker_call.extend(['--name', containerInfo['name']])
    if rm:
        base_docker_call.append('--rm')

    docker_tag = getDockerTag()
    tool = "%s/%s:%s" % (dockstore, tool, docker_tag)
    call = base_docker_call + [tool] + parameters
    return call, containerInfo

def prepareWorkDir(work_dir, parameters):
    if not work_dir:
        # Make sure all the paths we're accessing are in the same directory
        files = [par for par in parameters if os.path.isfile(par)]
        folders = [par for par in parameters if os.path.isdir(par)]
        work_dirs = set([os.path.dirname(fileName) for fileName in files] + [os.path.dirname(folder) for folder in folders])
        _log.info("Work dirs: %s" % work_dirs)
        if len(work_dirs) > 1:
            work_dir = os.path.commonprefix(list(work_dirs))
        elif len(work_dirs) == 1:
            work_dir = work_dirs.pop()

    #If there are no input files, or if their MRCA is '' (when working
    #with relative paths), just set the current directory as the work
    #dir
    if work_dir is None or work_dir == '':
        work_dir = os.getcwd()
    _log.info("Docker work dir: %s" % work_dir)

    #We'll mount the work_dir containing the paths as /data in the container,
    #so set all the paths to their basenames. The container will access them at
    #/data/<path>
    def adjustPath(path, wd):
        # Hack to relativize paths that are not provided as a
        # single argument (i.e. multiple paths that are
        # space-separated and quoted)
        if wd != '.':
            if not wd.endswith('/'):
                wd = wd + '/'
            return path.replace(wd, '')
        else:
            return path

    if work_dir and os.environ.get('CACTUS_DOCKER_MODE', 1) != "0":
        parameters = [adjustPath(par, work_dir) for par in parameters]
    return work_dir, parameters

def cactus_call(tool=None,
                work_dir=None,
                parameters=None,
                rm=True,
                check_output=False,
                infile=None,
                outfile=None,
                outappend=False,
                stdin_string=None,
                server=False,
                shell=False,
                port=None,
                check_result=False,
                dockstore=None,
                soft_timeout=None,
                job_name=None,
                features=None,
                fileStore=None,
                returnStdErr=False,
                realtimeStderrPrefix=None,
                gpus=None,
                cpus=None):
    mode = os.environ.get("CACTUS_BINARIES_MODE", "docker")
    if dockstore is None:
        dockstore = getDockerOrg()
    if parameters is None:
        parameters = []
    if tool is None:
        tool = "cactus"
    
    entrypoint = None
    if (len(parameters) > 0) and isinstance(parameters[0], list):
        # We have a list of lists, which is the convention for commands piped into one another.
        flattened = [i for sublist in parameters for i in sublist]
        chain_params = [' '.join(p) for p in [list(map(pipes.quote, q)) for q in parameters]]
        parameters = ['bash', '-c', 'set -eo pipefail && ' + ' | '.join(chain_params)]
        if mode == "docker":
            # We want to shell into bash directly rather than going
            # through the default cactus entrypoint.
            entrypoint = '/bin/bash'
            parameters = parameters[1:]
            work_dir, _ = prepareWorkDir(work_dir, flattened)

    if mode in ("docker", "singularity"):
        work_dir, parameters = prepareWorkDir(work_dir, parameters)

    if mode == "docker":
        call, containerInfo = dockerCommand(tool=tool,
                                            work_dir=work_dir,
                                            parameters=parameters,
                                            rm=rm,
                                            port=port,
                                            dockstore=dockstore,
                                            entrypoint=entrypoint,
                                            gpus=gpus, cpus=cpus)
    elif mode == "singularity":
        call = singularityCommand(tool=tool, work_dir=work_dir,
                                  parameters=parameters, port=port, file_store=fileStore,
                                  gpus=gpus, cpus=cpus)
    else:
        assert mode == "local"
        call = parameters

    if stdin_string:
        stdinFileHandle = subprocess.PIPE
    elif infile:
        stdinFileHandle = open(infile, 'r')
    else:
        stdinFileHandle = subprocess.DEVNULL
    stdoutFileHandle = None
    if outfile:
        stdoutFileHandle = open(outfile, 'a' if outappend else 'w')
    if check_output:
        stdoutFileHandle = subprocess.PIPE

    _log.info("Running the command %s" % call)
    rt_message = 'Running the command: \"{}\"'.format(' '.join(call))
    if features:
        rt_message += ' (features={})'.format(features)
    cactus_realtime_log(rt_message, log_debug = 'ktremotemgr' in call)

    # hack to keep track of memory usage for single machine
    time_v = os.environ.get("CACTUS_LOG_MEMORY") is not None and 'ktserver' not in call and 'redis-server' not in call

    # use /usr/bin/time -v to get peak memory usage
    if time_v:
        if not shell:
            shell = True
            call = ' '.join(shlex.quote(t) for t in call)
        call = '/usr/bin/time -f "CACTUS-LOGGED-MEMORY-IN-KB: %M" {}'.format(call)

    # optionally pipe stderr (but only if realtime logging enabled)
    # note the check below if realtime logging is enabled is rather hacky
    pid = None
    if realtimeStderrPrefix and RealtimeLogger.getLogger().level < logging.CRITICAL:
        # Make our pipe
        rfd, wfd = os.pipe()
        rfile = os.fdopen(rfd, 'rb', 0)
        wfile = os.fdopen(wfd, 'wb', 0)
        # And a different pipe for the memory log
        mlrfd, mlwfd = os.pipe()
        mlrfile = os.fdopen(mlrfd, 'rb', 0)
        mlwfile = os.fdopen(mlwfd, 'wb', 0)
                
        # Fork our child process (pid == 0) to catch stderr and log it
        pid = os.fork()
        if pid == 0:
            wfile.close()
            mlrfile.close()
            mem_log_line = ''
            while 1:
                data = rfile.readline()
                if not data:
                    break
                line = data.strip().decode()
                if 'CACTUS-LOGGED-MEMORY-IN-KB:' in line:
                    mem_log_line = line
                else:
                    RealtimeLogger.info('{}: {}'.format(realtimeStderrPrefix, line))                    
            rfile.close()
            mlwfile.write(mem_log_line.encode())
            mlwfile.close()            
            os._exit(0)
        else:
            assert pid > 0
            # main process carries on, but sending stderr to the pipe
            rfile.close()
            mlwfile.close()
            # note that only call_directly below actually does anything with errfile at the moment
            errfile = wfile
    else:
        errfile = subprocess.PIPE

    # hack to force consolidated to keep tmp files local (required to run at all in some singularity setups where /tmp is not writable)
    sub_env = None
    if (shell and 'cactus_consolidated' in call) or (not shell and any(['cactus_consolidated' in c for c in call])):
        sub_env = os.environ.copy()
        sub_env['TMPDIR']='.'

    process = subprocess.Popen(call, shell=shell, encoding=None,
                               stdin=stdinFileHandle, stdout=stdoutFileHandle,
                               stderr=errfile,
                               bufsize=-1, cwd=work_dir, env=sub_env)

    if server:
        return process

    memUsage = 0
    first_run = True
    start_time = time.time()
    output = stderr = None  # used later to report errors
    while True:
        try:
            # Wait a bit to see if the process is done
            output, stderr = process.communicate(stdin_string.encode() if first_run and stdin_string else None, timeout=10)
        except subprocess.TimeoutExpired:
            if mode == "docker":
                # Every so often, check the memory usage of the container
                updatedMemUsage = maxMemUsageOfContainer(containerInfo)
                if updatedMemUsage is not None:
                    assert memUsage <= updatedMemUsage, "memory.max_usage_in_bytes should never decrease"
                    memUsage = updatedMemUsage
            first_run = False
            if soft_timeout is not None and time.time() - start_time > soft_timeout:
                # Soft timeout has been triggered. Just return early.
                process.send_signal(signal.SIGINT)
                return None
        else:
            break
    if mode == "docker" and job_name is not None and features is not None and fileStore is not None:
        # Log a datapoint for the memory usage for these features.
        fileStore.logToMaster("Max memory used for job %s (tool %s) "
                              "on JSON features %s: %s" % (job_name, parameters[0],
                                                           json.dumps(features), memUsage))

    mem_log_line = None
    if pid and pid > 0:
        # It's not enough that the forked process is exited, it must be waited for or it's
        # deemed a zombine process:
        # https://medium.com/@BeNitinAgarwal/an-init-system-inside-the-docker-container-3821ee233f4b
        # and Toil will complain forever about it:
        # https://github.com/ComparativeGenomicsToolkit/cactus/issues/610#issuecomment-1015759593
        wfile.close()
        mem_log_line = mlrfile.readline().strip().decode()
        mlrfile.close()        
        os.wait()

    if stderr is not None:
        stderr = stderr.decode()
    if output is not None:
        output = output.decode()
        
    if process.returncode == 0:
        run_time = time.time() - start_time
        if time_v:
            call = call[len('/usr/bin/time -f "CACTUS-LOGGED-MEMORY-IN-KB: %M" '):]
        rt_message = "Successfully ran: \"{}\"".format(' '.join(call) if not shell else call)
        if features:
            rt_message += ' (features={})'.format(features)
        rt_message += " in {} seconds".format(round(run_time, 4))
        if time_v:
            if stderr:
                for line in stderr.split('\n'):
                    if 'CACTUS-LOGGED-MEMORY-IN-KB:' in line:
                        mem_log_line = line
                        break
            if mem_log_line:
                rt_message += ' and {} memory'.format(bytes2human(int(mem_log_line.split()[-1]) * 1024))
        cactus_realtime_log(rt_message, log_debug = 'ktremotemgr' in call)

    if check_result:
        return process.returncode

    if process.returncode != 0:
        out = "stderr={}".format(stderr) if stderr else ''

        sigill_msg = ''
        if abs(process.returncode) == 4:
            # this (rightfully) trips up many users, try to make logging a little more clear
            # https://github.com/ComparativeGenomicsToolkit/cactus/issues/574
            sigill_msg  = '\n\n\n'
            sigill_msg += '***********************************************************************************\n'
            sigill_msg += '***********************************************************************************\n'
            sigill_msg += 'ERROR: Your CPU is incompatible with the Cactus binaries.\n'
            sigill_msg += '       This is often due to it being too old to support AVX2 instructions.\n'
            sigill_msg += '       Please see the release notes for more details:\n'
            sigill_msg += '       https://github.com/ComparativeGenomicsToolkit/cactus/releases\n'
            sigill_msg += '***********************************************************************************\n'
            sigill_msg += '***********************************************************************************\n\n\n'
        
        if process.returncode > 0:
            raise RuntimeError("{}Command {} exited {}: {}".format(sigill_msg, call, process.returncode, out))
        else:
            raise RuntimeError("{}Command {} signaled {}: {}".format(sigill_msg, call, signal.Signals(-process.returncode).name, out))

    if check_output:
        return (output, stderr) if returnStdErr else output

    if returnStdErr:
        return stderr

class RunAsFollowOn(Job):
    def __init__(self, job, *args, **kwargs):
        Job.__init__(self, cores=0.1, memory=100000000, preemptable=True)
        self._args = args
        self._kwargs = kwargs
        self.job = job

    def run(self, fileStore):
        return self.addFollowOn(self.job(*self._args, **self._kwargs)).rv()

class RoundedJob(Job):
    """Thin wrapper around Toil.Job to round up resource requirements.

    Rounding is useful to make Toil's Mesos scheduler more
    efficient--it runs a process that is O(n log n) in the number of
    different resource requirements for every offer received, so
    thousands of slightly different requirements will slow down the
    leader and the workflow.
    """
    # Default rounding amount: 100 MiB
    roundingAmount = 100*1024*1024
    def __init__(self, memory=None, cores=None, disk=None, preemptable=None,
                 unitName=None, checkpoint=False, accelerators=None):
        if memory is not None:
            memory = self.roundUp(memory)
        if disk is not None:
            # hack: we may need extra space to cook up a singularity image on the fly
            #       so we add it (1.5G) here.
            # todo: only do this when needed
            disk = 1500*1024*1024 + self.roundUp(disk)
        super(RoundedJob, self).__init__(memory=memory, cores=cores, disk=disk,
                                         preemptable=preemptable, unitName=unitName,
                                         checkpoint=checkpoint, accelerators=accelerators)

    def roundUp(self, bytesRequirement):
        """
        Round the amount up to the next self.roundingAmount.

        >>> j = RoundedJob()
        >>> j.roundingAmount = 100000000
        >>> j.roundUp(1000)
        10000000
        >>> j.roundUp(200000000)
        200000000
        >>> j.roundUp(200000001)
        300000000
        """
        if bytesRequirement % self.roundingAmount == 0:
            return bytesRequirement
        return (bytesRequirement // self.roundingAmount + 1) * self.roundingAmount

    def _runner(self, *args, jobStore=None, fileStore=None, **kwargs):
        # We aren't supposed to override this. Toil can change the signature at
        # any time. But Toil has passed all the arguments by name=value for a
        # while, so we are pretty safe fetching out just the jobStore and
        # fileStore like this.
        if jobStore.config.workDir is not None:
            os.environ['TMPDIR'] = fileStore.getLocalTempDir()

        super(RoundedJob, self)._runner(*args, jobStore=jobStore,
                                        fileStore=fileStore, **kwargs)

def readGlobalFileWithoutCache(fileStore, jobStoreID):
    """Reads a jobStoreID into a file and returns it, without touching
    the cache.

    Works around toil issue #1532.
    """
    f = fileStore.getLocalTempFile()
    fileStore.jobStore.readFile(jobStoreID, f)
    return f

class ChildTreeJob(RoundedJob):
    """Spreads the child-job initialization work among multiple jobs.

    Jobs with many children can often be a bottleneck (because they
    are written serially into the jobStore in a consistent-write
    fashion). Subclasses of this job will automatically spread out
    that work amongst a tree of jobs, increasing the total work done
    slightly, but reducing the wall-clock time taken dramatically.
    """
    def __init__(self, memory=None, cores=None, disk=None, preemptable=None,
                 unitName=None, checkpoint=False, maxChildrenPerJob=20):
        self.queuedChildJobs = []
        self.maxChildrenPerJob = maxChildrenPerJob
        super(ChildTreeJob, self).__init__(memory=memory, cores=cores, disk=disk,
                                           preemptable=preemptable, unitName=unitName,
                                           checkpoint=checkpoint)

    def addChild(self, job):
        self.queuedChildJobs.append(job)
        return job

    def _run(self, *args, **kwargs):
        # We really shouldn't be overriding _run, but we need to to hook in
        # some code after the derived class's run() but before it saves all its
        # child relationships. So we handle our arguments with tweezers,
        # because they could be anything.

        # Pass them along to the real _run, which will call back to our derived
        # class's run()
        ret = super(ChildTreeJob, self)._run(*args, **kwargs)

        # Now we can do our actual work.
        if len(self.queuedChildJobs) <= self.maxChildrenPerJob:
            # The number of children is small enough that we can just
            # add them directly.
            for childJob in self.queuedChildJobs:
                super(ChildTreeJob, self).addChild(childJob)
        else:
            # Too many children, so we have to build a tree to avoid
            # bottlenecking on consistently serializing all the jobs.

            # compute the number of levels (after root) of our job tree
            num_levels = math.floor(math.log(len(self.queuedChildJobs), self.maxChildrenPerJob))

            # fill out all the internal nodes of the tree, where the root is self
            # they will be empty RoundedJobs
            prev_level = [self]
            level = []
            for i in range(num_levels):
                for parent_job in prev_level:
                    # with this check, we allow a partial split of the last level
                    # to account for rounding
                    if len(level) * self.maxChildrenPerJob < len(self.queuedChildJobs):
                        for j in range(self.maxChildrenPerJob):
                            child_job = RoundedJob()
                            if parent_job is self:
                                super(ChildTreeJob, self).addChild(child_job)
                            else:
                                parent_job.addChild(child_job)
                            level.append(child_job)
                    else:
                        level.append(parent_job)

                prev_level = level
                level = []

            # add the leaves.  these will be the jobs in self.queuedChildJobs
            leaves_added = 0
            for parent_job in prev_level:
                num_children = min(len(self.queuedChildJobs) - leaves_added, self.maxChildrenPerJob)
                for j in range(num_children):
                    if parent_job is self:
                        super(ChildTreeJob, self).addChild(self.queuedChildJobs[leaves_added])
                    else:
                        parent_job.addChild(self.queuedChildJobs[leaves_added])
                    leaves_added += 1
            assert leaves_added == len(self.queuedChildJobs)

        return ret

def dumpStacksHandler(signal, frame):
    """Signal handler to print the stacks of all threads to stderr"""
    fh = sys.stderr
    print("###### stack traces {} ######".format(datetime.now().isoformat()), file=fh)
    id2name = dict([(th.ident, th.name) for th in threading.enumerate()])
    for threadId, stack in sys._current_frames().items():
        print("# Thread: {}({})".format(id2name.get(threadId,""), threadId), file=fh)
        traceback.print_stack(f=stack, file=fh)
    print("\n", file=fh)
    fh.flush()

def enableDumpStack(sig=signal.SIGUSR1):
    """enable dumping stacks when the specified signal is received"""
    signal.signal(sig, dumpStacksHandler)

def unzip_gzs(job, input_paths, input_ids):
    """ go through a list of files and unzip any that end with .gz and return a list 
    of updated ids.  files that don't end in .gz are just passed through.  relying on the extension
    is pretty fragile but better than nothing """
    unzipped_ids = []
    for input_path, input_id in zip(input_paths, input_ids):
        if input_path.endswith('.gz'):
            unzip_job = job.addChildJobFn(unzip_gz, input_path, input_id, disk=10*input_id.size)
            unzipped_ids.append(unzip_job.rv())
        else:
            unzipped_ids.append(input_id)
    return unzipped_ids

def unzip_gz(job, input_path, input_id):
    """ unzip a single file """
    work_dir = job.fileStore.getLocalTempDir()
    assert input_path.endswith('.gz')
    fa_path = os.path.join(work_dir, os.path.basename(input_path))
    job.fileStore.readGlobalFile(input_id, fa_path, mutable=True)
    cactus_call(parameters=['gzip', '-fd', os.path.basename(fa_path)], work_dir=work_dir)
    return job.fileStore.writeGlobalFile(fa_path[:-3])

def zip_gzs(job, input_paths, input_ids, list_elems = None):
    """ zip up some files.  the input_ids can be a list of lists.  if it is, then list_elems
    can be used to only zip a subset (leaving everything else) on each list."""
    zipped_ids = []
    for input_path, input_list in zip(input_paths, input_ids):
        if input_path.endswith('.gz'):
            if isinstance(input_list, list) or isinstance(input_list, tuple):
                output_list = []
                for i, elem in enumerate(input_list):
                    if not list_elems or i in list_elems:
                        output_list.append(job.addChildJobFn(zip_gz, input_path, elem, disk=2*elem.size).rv())
                    else:
                        output_list.append(elem)
                zipped_ids.append(output_list)
            else:
                zipped_ids.append(job.addChildJobFn(zip_gz, input_path, input_list, disk=2*input_list.size).rv())
        else:
            zipped_ids.append(input_list)
    return zipped_ids
    
def zip_gz(job, input_path, input_id):
    """ zip a single file """
    work_dir = job.fileStore.getLocalTempDir()
    fa_path = os.path.join(work_dir, os.path.basename(input_path))
    if fa_path.endswith('.gz'):
        fa_path = fa_path[:-3]
    job.fileStore.readGlobalFile(input_id, fa_path, mutable=True)
    cactus_call(parameters=['gzip', '-f', os.path.basename(fa_path)], work_dir=work_dir)
    return job.fileStore.writeGlobalFile(fa_path + '.gz')

def get_aws_region(full_path):
    """ parse aws:region:url  to just get region (toil surely has better way to do this but in rush)"""
    if full_path.startswith('aws:'):
        return full_path.split(':')[1]
    else:
        return None

def write_s3(local_path, s3_path, region=None):
    """ cribbed from toil-vg.  more convenient just to throw hal output on s3
    than pass it as a promise all the way back to the start job to export it locally """
    assert s3_path.startswith('s3://')
    bucket_name, name_prefix = s3_path[5:].split("/", 1)
    botocore_session = botocore.session.get_session()
    botocore_session.get_component('credential_provider').get_provider('assume-role').cache = botocore.credentials.JSONFileCache()
    boto3_session = boto3.Session(botocore_session=botocore_session)

    # Connect to the s3 bucket service where we keep everything
    s3 = boto3_session.client('s3')
    try:
        s3.head_bucket(Bucket=bucket_name)
    except:
        if region:
            s3.create_bucket(Bucket=bucket_name, CreateBucketConfiguration={'LocationConstraint':region})
        else:
            s3.create_bucket(Bucket=bucket_name)

    s3.upload_file(local_path, bucket_name, name_prefix)

def get_faidx_subpath_rename_cmd():
    """
    transform chr1:10-15 (1-based inclusive) into chr1_sub_9_15 (0-based end open)
    this is a format that contains no special characters in order to make assembly hubs
    happy.  But it does require conversion going into vg which wants chr[9-15] and
    hal2vg is updated to do this autmatically
    """
    return ['sed', '-e', 's/\([^:]*\):\([0-9]*\)-\([0-9]*\)/echo "\\1_sub_$((\\2-1))_\\3"/e']

