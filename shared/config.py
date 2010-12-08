import xml.etree.ElementTree as ET
import os
import sys

from sonLib.bioio import system
from sonLib.bioio import getRandomAlphaNumericString

def checkDatabaseConf(databaseConf):
    """Function checks the database conf is as expected and creates useful exceptions
    if not"""
    dataString = ET.tostring(databaseConf)
    if databaseConf.tag != "st_kv_database_conf":
        raise RuntimeError("The database conf string is improperly formatted: %s" % dataString)
    if not databaseConf.attrib.has_key("type"):
        raise RuntimeError("The database conf string does not have a type attrib: %s" % dataString)
    typeString = databaseConf.attrib["type"]
    if typeString == "mysql":
        mysql = databaseConf.find("mysql")
        if mysql == None:
            raise RuntimeError("Database conf is of type mysql but there is no nested mysql tag: %s" % dataString)
        if set(mysql.attrib.keys()) != set(("host", "port", "password", "user", "database_name")):
            raise RuntimeError("Mysql tag is improperly formatted: %s" % dataString)
    if typeString == "postgresql":
        postgresql = databaseConf.find("postgresql")
        if postgresql == None:
            raise RuntimeError("Database conf is of type postgresql but there is no nested postgresql tag: %s" % dataString)
        if set(mysql.attrib.keys()) != set(("host", "password", "user", "database_name")):
            raise RuntimeError("Postgresql tag is improperly formatted: %s" % dataString)
    elif typeString == "tokyo_cabinet":
        tokyoCabinet = databaseConf.find("tokyo_cabinet")
        if tokyoCabinet == None:
            raise RuntimeError("Database conf is of type tokyo cabinet but there is no nested tokyo cabinet tag: %s" % dataString)
        if not tokyoCabinet.attrib.has_key("database_dir"):
            raise RuntimeError("The tokyo cabinet tag has no database_dir tag: %s" % dataString)
    else:
        raise RuntimeError("Unrecognised database type in conf string: %s" % typeString)

class CactusWorkflowExperiment:
    """Object used for generating cactus workflow experiment config files,
    using the inputs to generate valid config strings and files.
    """
    def __init__(self, sequences, newickTreeString, outputDir=None, databaseName=None, databaseConf=None, configFile=None):
        self.experiment = ET.Element("cactus_workflow_experiment")
        if databaseName == None:
            self.databaseName = "cactusDisk_%s_%i" % (getRandomAlphaNumericString(), os.getpid()) #Needs to be unique
        else:
            self.databaseName = databaseName
        self.databaseFile = None
        #Do the database first
        database = ET.SubElement(self.experiment, "cactus_disk")
        self.mysql = 0
        self.postgresql = 0
        if databaseConf != None:
            checkDatabaseConf(databaseConf)
            databaseConf = ET.fromstring(ET.tostring(databaseConf))
            checkDatabaseConf(databaseConf) #This is just a redundant check
            database.append(databaseConf)
            if databaseConf.attrib["type"] == "mysql":
                self.mysql = 1
                #Add a table name:
                databaseConf.find("mysql").attrib["table_name"] = self.databaseName
            elif databaseConf.attrib["type"] == "postgresql":
                self.postgresql = 1
                #Add a table name:
                databaseConf.find("postgresql").attrib["table_name"] = self.databaseName
            else:
                assert databaseConf.attrib["type"] == "tokyo_cabinet"
                tokyoCabinet = databaseConf.find("tokyo_cabinet")
                tokyoCabinet.attrib["database_dir"] = os.path.join(tokyoCabinet.attrib["database_dir"], self.databaseName)
                self.databaseFile = tokyoCabinet.attrib["database_dir"]
                assert not os.path.exists(self.databaseFile)
        else:
            databaseConf = ET.SubElement(database, "st_kv_database_conf")
            databaseConf.attrib["type"] = "tokyo_cabinet"
            tokyoCabinet = ET.SubElement(databaseConf, "tokyo_cabinet")
            assert outputDir != None
            assert os.path.exists(outputDir) #Check it exists
            self.databaseFile = os.path.join(outputDir, self.databaseName)
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
        if self.mysql:
            assert self.databaseFile == None
            i = self.experiment.find("cactus_disk").find("st_kv_database_conf").find("mysql").attrib
            #Connect to MYSQL and remove the table..
            system("mysql --host=%s --port=%s --user=%s --password='%s' --execute='drop table %s.%s'" \
                   % (i["host"], i["port"], i["user"], i["password"], i["database_name"], i["table_name"]))
        elif self.postgresql:
            assert self.databaseFile == None
            i = self.experiment.find("cactus_disk").find("st_kv_database_conf").find("mysql").attrib
            #Connect to MYSQL and remove the table..
            system("psql --host=%s --user=%s --dbname=%s --command='drop table %s'" \
                   % (i["host"], i["user"], i["database_name"], i["table_name"]))
        else:
            assert self.databaseFile != None
            system("rm -rf %s" % self.databaseFile)
        
    def getDatabaseString(self):
        """Gets the database conf string.
        """
        return ET.tostring(self.experiment.find("cactus_disk").find("st_kv_database_conf"))
