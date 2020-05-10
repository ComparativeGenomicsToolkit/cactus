
"""Tests the database servers using DbServerService
"""
import unittest
import os
from toil.common import Toil
from toil.job import Job

from cactus.shared.experimentWrapper import DbElemWrapper
from cactus.shared.common import cactus_call
from cactus.pipeline.dbServerToil import DbServerService
from cactus.pipeline.ktserverControl import KtServer
from cactus.pipeline.redisServerControl import RedisServer

import xml.etree.ElementTree as ET

KT_CONF_STRING = '<st_kv_database_conf type="kyoto_tycoon"><kyoto_tycoon in_memory="1" port="{port}" snapshot="0"/></st_kv_database_conf>'
REDIS_CONF_STRING = '<st_kv_database_conf type="redis"><redis in_memory="1" port="{port}" snapshot="0"/></st_kv_database_conf>'


class DbTestJob(Job):
    def __init__(self, dbElem, testFunc, snapshotPath=None, memory="100M", cores=1, disk="100M"):
        self.dbElem = dbElem
        self.testFunc = testFunc
        self.snapshotPath = snapshotPath
        Job.__init__(self, memory=memory, cores=cores, disk=disk, preemptable=False)

    def run(self, fileStore):
        snapshotID = None
        if self.snapshotPath is not None:
            snapshotID = fileStore.importFile("file://" + self.snapshotPath)
        dbService = DbServerService(dbElem=self.dbElem, isSecondary=False,
                                    existingSnapshotID=snapshotID,
                                    memory="100M", cores=1, disk="100M")

        s = self.addService(dbService)
        xmlString, snapshotExportID = s.rv(0), s.rv(1)
        self.addChildFn(self.testFunc, xmlString,
                        memory="100M", cores=1, disk="100M")
        return snapshotExportID

    def addService(self, service):
        """Works around toil issue #1695, returning a Job rather than a Promise."""
        super(DbTestJob, self).addService(service)
        return self._services[-1]


def runAndTestDbServerService(dbElem, testFunc, inputSnapshotPath=None, outputSnapshotPath=None):
    dbTestJob = DbTestJob(dbElem, testFunc, inputSnapshotPath)
    options = Job.Runner.getDefaultOptions("./testDbServerWorkflow")
    options.logLevel = "INFO"
    options.clean = "always"

    with Toil(options) as toil:
        snapshotID = toil.start(dbTestJob)
        if outputSnapshotPath is not None:
            toil.exportFile(snapshotID, "file://" + outputSnapshotPath)


def _testKtSetFlag(xmlString):
    """ set FLAG to 1 in a Kyoto_Tycoon database (the first server)"""
    dbElem = DbElemWrapper(ET.fromstring(xmlString))
    ktServer = KtServer(dbElem)
    cactus_call(parameters=['ktremotemgr', 'set'] + ktServer.getRemoteParams() + ['FLAG', '1'])


def _testKtGetFlag(xmlString):
    """ get and check FLAG in a Kyoto_Tycoon database (the second server)"""
    dbElem = DbElemWrapper(ET.fromstring(xmlString))
    ktServer = KtServer(dbElem)
    flag = cactus_call(parameters=['ktremotemgr', 'get'] + ktServer.getRemoteParams() + ['FLAG'], check_output=True)
    assert flag.strip() == '1'


def _testRedisSetFlag(xmlString):
    """ set FLAG to 1 in a Redis database (the first server)"""
    dbElem = DbElemWrapper(ET.fromstring(xmlString))
    redisServer = RedisServer(dbElem)
    base_call = ['redis-cli'] + redisServer.getRemoteParams()
    out = cactus_call(parameters=base_call + ['set', 'FLAG', '1'], check_output=True)
    assert out.strip() == "OK"


def _testRedisGetFlag(xmlString):
    """ get and check FLAG in a Redis database (the second server)"""
    dbElem = DbElemWrapper(ET.fromstring(xmlString))
    redisServer = RedisServer(dbElem)
    base_call = ['redis-cli'] + redisServer.getRemoteParams()
    flag = cactus_call(parameters=base_call + ['get', 'FLAG'], check_output=True)
    assert flag.strip() == '1'


class TestCase(unittest.TestCase):
    def setUp(self) -> None:
        """ make a file for saving snapshots """
        os.mkdir("./testSnapshot")
        self.snapshotTempPath = os.path.abspath("./testSnapshot/snapshot.temp")
        open(self.snapshotTempPath, "w").close()

    #@unittest.skip("skip KT!")
    def testKtServerService(self):
        initialDbElem = DbElemWrapper(ET.fromstring(KT_CONF_STRING))
        # The first database saves the snapshot after termination
        runAndTestDbServerService(dbElem=initialDbElem, testFunc=_testKtSetFlag,
                                  outputSnapshotPath=self.snapshotTempPath)
        # The second database loads the snapshot
        runAndTestDbServerService(dbElem=initialDbElem, testFunc=_testKtGetFlag,
                                  inputSnapshotPath=self.snapshotTempPath)

    #@unittest.skip("skip Redis!")
    def testRedisServerService(self):
        initialDbElem = DbElemWrapper(ET.fromstring(REDIS_CONF_STRING))
        # The first database saves the snapshot after termination
        runAndTestDbServerService(dbElem=initialDbElem, testFunc=_testRedisSetFlag,
                                  outputSnapshotPath=self.snapshotTempPath)
        # The second database loads the snapshot
        runAndTestDbServerService(dbElem=initialDbElem, testFunc=_testRedisGetFlag,
                                  inputSnapshotPath=self.snapshotTempPath)

    def tearDown(self) -> None:
        os.remove(self.snapshotTempPath)
        os.rmdir("./testSnapshot")


if __name__ == '__main__':
    unittest.main()
