#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Fixes paths, which tend to be absolute in the project and experiment xml files 
to facilitate copying the database. (ex run this program after copying the db in order to sync
the xml files with the new path).  old xml files are backed up (with .old added to name)
""" 
import os
import sys
from optparse import OptionParser
import xml.etree.ElementTree as ET
import copy
import socket

from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.progressive.experimentWrapper import ExperimentWrapper
from cactus.progressive.configWrapper import ConfigWrapper
from cactus.progressive.outgroup import GreedyOutgroup
from sonLib.bioio import system

def updateProject(path):
    mcProj = MultiCactusProject()
    mcProj.readXML(path)
    basePath, name = os.path.split(path)
    
    for name,oldPath in mcProj.expMap.items():
        fileName = os.path.basename(oldPath)
        dirName = os.path.dirname(oldPath).rpartition('/')[2] 
        newPath = os.path.join(basePath, dirName, fileName)
        
        if not os.path.isfile(newPath):
            raise RuntimeError("Experiment file %s not found\n" % newPath)
        
        mcProj.expMap[name] = newPath   
        
        exp = ExperimentWrapper(ET.parse(newPath).getroot())
        
        oldDbDir = exp.getDbDir()
        if oldDbDir is not None:
            dbDirName = oldDbDir[oldDbDir.find(name):]
            newDbDir = os.path.join(basePath, dbDirName)
            exp.setDbDir(newDbDir)
        
        oldRefPath = exp.getReferencePath()
        if oldRefPath is not None:
            refName = oldRefPath[oldRefPath.find(name):]
            newRefPath = os.path.join(basePath, refName)
            exp.setReferencePath(newRefPath)
    
        oldHalPath = exp.getHALPath()
        if oldHalPath is not None:
            halName = oldHalPath[oldHalPath.find(name):]
            newHalPath = os.path.join(basePath, halName)
            exp.setHALPath(newHalPath)
        
        oldHostName = exp.getDbHost()
        if oldHostName is not None:
            newHostName = socket.gethostname()
            exp.setDbHost(newHostName)
        
        system("cp %s %s.old" %(newPath, newPath))
        exp.writeXML(newPath)
    
    mcProj.writeXML(path)

def main():
    usage = "usage: %prog <project xml path>"
    description = "Update the _project.xml and _experiment.xml files to use current directory structure"
    parser = OptionParser(usage=usage, description=description)
    
    options, args = parser.parse_args()
    
    if len(args) != 1:
        parser.print_help()
        raise RuntimeError("Wrong number of arguments")
    
    path = os.path.abspath(args[0])
    
    if not os.path.isfile(path):
        raise RuntimeError("Project file %s not found\n" % path)
    
    system("cp %s %s.old" %(path, path))
    updateProject(path)
    
    return 0
    

if __name__ == '__main__':    
    main()
