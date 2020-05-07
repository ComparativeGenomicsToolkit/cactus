import unittest
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
    def __init__(self, dbElem, testFunc, memory="100M", cores=1, disk="100M"):
        self.dbElem = dbElem
        self.testFunc = testFunc
        Job.__init__(self, memory=memory, cores=cores, disk=disk, preemptable=False)

    def run(self, fileStore):
        dbService = DbServerService(dbElem=self.dbElem, isSecondary=False,
                                    memory="100M", cores=1, disk="100M")
        xmlString = self.addService(dbService).rv(0)
        self.addChildFn(self.testFunc, xmlString,
                        memory="100M", cores=1, disk="100M")

    def addService(self, service):
        """Works around toil issue #1695, returning a Job rather than a Promise."""
        super(DbTestJob, self).addService(service)
        return self._services[-1]


def runAndTestDbServerService(dbElem, testFunc):
    dbTestJob = DbTestJob(dbElem, testFunc)
    options = Job.Runner.getDefaultOptions("./testDbServerWorkflow")
    options.logLevel = "INFO"
    options.clean = "always"

    with Toil(options) as toil:
        toil.start(dbTestJob)


def _testKt(xmlString):
    dbElem = DbElemWrapper(ET.fromstring(xmlString))
    ktServer = KtServer(dbElem)
    cactus_call(parameters=['ktremotetest', 'order'] + ktServer.getRemoteParams() + ['10'])


def _testRedis(xmlString):
    dbElem = DbElemWrapper(ET.fromstring(xmlString))
    redisServer = RedisServer(dbElem)
    pingOut = cactus_call(parameters=['redis-cli'] + redisServer.getRemoteParams() + ['PING'],
                          check_output=True)
    assert pingOut.strip() == "PONG"


class TestCase(unittest.TestCase):

    #@unittest.skip("skip KT!")
    def testKtServerService(self):
        initialDbElem = DbElemWrapper(ET.fromstring(KT_CONF_STRING))
        runAndTestDbServerService(initialDbElem, _testKt)

    #@unittest.skip("skip Redis!")
    def testRedisServerService(self):
        initialDbElem = DbElemWrapper(ET.fromstring(REDIS_CONF_STRING))
        runAndTestDbServerService(initialDbElem, _testRedis)


if __name__ == '__main__':
    unittest.main()