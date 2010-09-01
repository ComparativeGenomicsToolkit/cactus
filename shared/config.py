import xml.etree.ElementTree as ET
import os

from sonLib.bioio import system

MYSQL_DEFINED = 0

def initialiseMysqlConf(host, port, user, password, database):
    pass

class CactusWorkflowExperiment:
    def __init__(self, sequences, newickTreeString, tempDir='./', configFile=None):
        self.experiment = ET.Element("cactus_workflow_experiment")
        self.kvDatabaseName = "cactusDisk" #Needs to be unique
        assert os.path.exists(tempDir) #Check it exists
        #Do the database first
        database = ET.SubElement(self.experiment, "cactus_disk")
        if MYSQL_DEFINED:
            pass
        else:
            databaseConf = ET.SubElement(database, "st_kv_database_conf")
            databaseConf.attrib["type"] = "tokyo_cabinet"
            tokyoCabinet = ET.SubElement(databaseConf, "tokyo_cabinet")
            self.databaseFile = os.path.join(tempDir, self.kvDatabaseName)
            assert not os.path.exists(self.databaseFile)
            tokyoCabinet.attrib["database_dir"] = self.databaseFile
        #Now add in the user stuff..
        self.experiment.attrib["sequences"] = " ".join(sequences)
        self.experiment.attrib["species_tree"] = newickTreeString
        self.experiment.attrib["config"] = "default"
        if configFile != None:
            self.experiment.attrib["config"] = configFile
    
    def writeExperimentFile(self, experimentFile):
        """Writes a description of the config into this directory.
        """
        fileHandle = open(experimentFile, 'w')
        ET.ElementTree(self.experiment).write(fileHandle)
        fileHandle.close()
        
    def cleanupDatabase(self):
        """Removes the database that was created.
        """
        if MYSQL_DEFINED:
            pass
        else:
            system("rm -rf %s" % self.databaseFile)
        
    def getDatabaseString(self):
        """Gets the database conf string.
        """
        return ET.tostring(self.experiment.find("cactus_disk").find("st_kv_database_conf"))
