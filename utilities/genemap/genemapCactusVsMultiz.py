#!/usr/bin/env python

import os, sys, re
from optparse import OptionParser
import xml.etree.ElementTree as ET

class Gene(list):
    def __init__(self, name, exonCount, chr, start, end):
        self.name = name
        self.exonCount = exonCount
        self.chr = chr
        self.start = start
        self.end = end
    def setStatus(self, status):
        self.status = status
    def setMissedSpecies(self, missedSpecies):
        self.missedSpecies = missedSpecies


class Species(list):
    def __init__(self, name, exonCount):
        self.name = name
        self.exonCount = exonCount
    def setStatus(self, status):
        self.status = status
    def setAlnBases(self, count):
        self.alnBases = count


class Exon:
    def __init__(self, line):
        #self.description = line      
        items = line.strip().split('\t')
        if len(items) < 20:
            sys.stderr.write("Wrong format: multizLine has less than 16 fields.\n%s\n" %(line))
        self.id = items[7]
        self.gene = items[0]
        self.chr = items[1]
        self.genestart = items[2]
        self.geneend = items[3]
        self.start = items[8]
        self.end = items[9]
        self.exonCount = int(items[6])
        self.codonBases = int(items[11])
        self.insBases = int(items[13])
        self.delBases = int(items[14])
        self.indelBases = self.insBases - self.delBases
        #self.indelBases = int(items[13]) - int(items[14])
        self.species = items[10]
        self.alnBases = int(items[12]) #number of aligned bases
        self.missingBases = int(items[15]) + int(items[16])
        self.numContigs = int(items[17])
        self.orderCons = items[19]
        maxIndel = self.codonBases*0.1
        #self.gapBases = int(items[16])
        #self.numRanges = int(items[16])
        if self.missingBases > 0:#if has any gap at all, marked it as having 'missing data'
            self.status = "MD" #Missing Data
        elif self.alnBases == 0:
            self.status = "NA" #Not Applicaple 
        #elif items[17] == "yes": 
        elif self.numContigs > 1 or self.insBases > maxIndel or self.delBases > maxIndel or self.indelBases%3 != 0 :
            self.status = "N" #Incomplete / non-conserved
        else:
            self.status = "C" #Complete/ conserved


def mapId2name(file):
    #mapped refseq id to gene name
    f = open(file, "r")
    id2name = {}
    name2id = {}
    for line in f.readlines():
        items = line.strip().split('\t')
        if len(line) < 2:
            continue
        else:
            id2name[items[1]] = items[0]
            name2id[items[0]] = items[1]
    f.close()
    return (id2name, name2id)

def getGene(genes, name):
    for gene in genes:
        if gene.name == name:
            return gene
    return None

def getSpecies(gene, name):
    for species in gene:
        if species.name == name:
            return species
    return None

def getExon(species, id):
    for exon in species:
        if exon.id == id:
            return exon
    return None

def multizgene_getExon(gene, id):
    for species in gene:
        exon = getExon(species, id)
        if exon:
            return exon
    return None

#==== cactus
def cactus_getGene(genes, name):
    for gene in genes:
        if gene.get("name") == name:
            return gene
    return None

def cactus_getSpecies(gene, name):
    species = gene.findall("species")
    for spc in species:
        if spc.get("name") == name:
            return spc
    return None

def cactus_getSpcExonStatus(spc, id):
    exons = spc.getiterator("exon")
    #spcStatus = spc.get("status")
    mdexons = spc.get("missingdataexons").split(',')

    hasExon = False
    for exon in exons:
        if exon.get("id") == id:
            hasExon = True
            if exon.get("status") == "1":
                return "C"
    if not hasExon: #so NA could include cases where the exons are skipped because of MissingData
        status = "NA"
    elif id in mdexons:
        status = 'MD'
    else:
        status = "N"
    return status 

def readMultiz(file, id2name):
    f = open(file, "r")
    genes = []
    for line in f.readlines():
        if re.search('#', line): #comment line
            continue
        exon = Exon(line)
        if exon.gene not in id2name:
            continue

        genename = id2name[exon.gene]
        gene = getGene(genes, genename)
        if not gene:#add gene to list if it is not in the list yet
            gene = Gene(genename, exon.exonCount, exon.chr, exon.genestart, exon.geneend)
            genes.append(gene)

        species = getSpecies(gene, exon.species)
        if not species:#add species to gene if it is not in gene yet
            species = Species(exon.species, exon.exonCount)
            gene.append(species)

        if not getExon(species, exon.id):
            species.append(exon)#add exon to species
    return genes

def setMultizSpeciesStatus(species):
    #A species is marked as having 'complete' homolog when all the exons of the reference present and 
    #are conserved (status = yes) and [the order and orientation] of the ref. exons are conserved
    #
    #A species is marked as having 'partial' homolog when exons that have missing data are allowed, and
    #(what about exons that get skipped because of msising data?)
    #
    #A species is marked as having incomplete homolog when it has some Non-conserved exon (status = no)
    #
    #If there is no homolog at all mapped to the gene, the species is marked as 'NA'
    status = "C"
    alnBases = 0
    hasMissingData = False
    for exon in species:
        alnBases += exon.alnBases
        if exon.status == "NA" or exon.status == "N" or exon.orderCons == "no":#species status is no if any of its exon is no
            status = "N"
            break
        #if exon.alnBases == 0 and exon.gapBases == 0: #totally miss the exon
        #    status = "N"
        #    break
        if exon.status == "MD": #has missing data exon
            hasMissingData = True

    if alnBases == 0:
        status = "NA"
    elif status != "N" and hasMissingData: #
        status = "MD"
    species.setStatus(status)
    #sys.stderr.write("Getting spc status: %s\t%s\t%d\n" %(species.name, species.status, species.alnBases)) 
    return

def setMultizGeneStatus(genes):
    #a gene is OK (yes) if it presents in all species with status Yes
    #return 'no' if it presents in all species but broken in at least one
    #'missed' if it does not aligned in at least one species, but not all
    #'NA' if it does not present in any of the species

    for gene in genes:
        status = "C"
        numNAspecies = 0
        hasMDspc = False
        for species in gene:
            setMultizSpeciesStatus(species)
        for species in gene:
            #sys.stderr.write("%s\t%s\t%s\n" %(gene.name, species.name, species.status))
            if species.status == "NA":#gene is absent in current species
                numNAspecies += 1
                hasMDspc = True #TESTING...
                continue
            elif species.status == "N":#incomplete (non-conserved) homolog
            #if species.status == "N":#incomplete (non-conserved) homolog
                status = "N"
                break
            elif species.status == "MD": #species with partially complete homolog
                hasMDspc = True

        if numNAspecies == len(gene):
            gene.setStatus('NA')
        elif species.status != "N" and hasMDspc:
        #if species.status != "N" and hasMDspc:
            gene.setStatus("MD")
        else:
            gene.setStatus(status)
        #gene.setMissedSpecies(missedSpecies)
        #sys.stderr.write("%s\t%s\n\n" %(gene.name, gene.status))
    return

def cactus_getSpcExonString(species):
    exons = species.getiterator("exon")
    exonStr = ""
    for exon in exons:
        exonStr += "%s_%s_%s," %(exon.get("id"), exon.get("status"), exon.get("strand"))
    return exonStr

def cactus_getSpcStatus2(species, exonCount):
    #assert(exonCount > 0)
    exonStr = cactus_getSpcExonString(species)
    #sys.stderr.write("species: %s\n" % species.get("name"))
    #sys.stderr.write("exonStr: [%s]\n" %exonStr)
    strands = ["\+", "-"]
    for strand in strands:
        forwardStr = "\D*0_1_%s.*" %(strand) #first exon
        for id in range(1, exonCount):
            forwardStr += "\D+%d_1_%s.*" %(id, strand)
        #sys.stderr.write("forwardStr: [%s]\n" %forwardStr)
        
        backwardStr = "\D*%d_1_%s.*" %(exonCount -1, strand) #last exon
        for id in range(exonCount -2, -1, -1):
            backwardStr += "\D+%d_1_%s.*" %(id, strand)
        #sys.stderr.write("backwardStr: [%s]\n" %backwardStr)

        forwardPattern = re.compile(forwardStr)
        fm = forwardPattern.search(exonStr)
        backwardPattern = re.compile(backwardStr)
        bm = backwardPattern.search(exonStr)
        #if fm:
        #    sys.stderr.write("Match forward String\n")
        #if bm:
        #    sys.stderr.write("Match backward String\n")
        if fm or bm:
            return "C"

    #Missingdata cases:
    mdExonsStr = species.get("missingdataexons")
    mdExons = mdExonsStr.split(',')
    for strand in strands:
        first = "0"
        if first in mdExons:
            forwardStr = "\D*0_[01]_%s.*" %(strand) #first exon
        else:
            forwardStr = "\D*0_1_%s.*" %(strand) #first exon
        
        for id in range(1, exonCount):
            if str(id) in mdExons:
                forwardStr += "\D+%d_[01]_%s.*" %(id, strand)
            else:
                forwardStr += "\D+%d_1_%s.*" %(id, strand)
        #sys.stderr.write("forwardStr: [%s]\n" %forwardStr)
        
        if str(exonCount -1) in mdExons:
            backwardStr = "\D*%d_[01]_%s.*" %(exonCount -1, strand) #last exon
        else:
            backwardStr = "\D*%d_1_%s.*" %(exonCount -1, strand) #last exon
        for id in range(exonCount -2, -1, -1):
            if str(id) in mdExons:
                backwardStr += "\D+%d_[01]_%s.*" %(id, strand)
            else:
                backwardStr += "\D+%d_1_%s.*" %(id, strand)
        #sys.stderr.write("backwardStr: [%s]\n" %backwardStr)

        forwardPattern = re.compile(forwardStr)
        fm = forwardPattern.search(exonStr)
        backwardPattern = re.compile(backwardStr)
        bm = backwardPattern.search(exonStr)
        if fm or bm:
            return "MD"

    return "N"

def cactus_getGeneStatus2(gene):
    species = gene.findall("species")
    exonCount = int(gene.get("exonCount"))
    hasMDspecies = False
    for spc in species:
        spcStatus = spc.get("status")
        if spcStatus == "N":
            spcStatus = cactus_getSpcStatus2(spc, exonCount)
        if spcStatus == "MD":
            hasMDspecies = True
        
        if spcStatus == "N":
            return "N"
    if hasMDspecies:
        return "MD"
    else:
        return 'C'


def getCactusExonStatus(exonid, gene):
    #return 'C' if input exon present in all species and perfectly conserved 
    #"N" when exon is broken (non-conserved) or absent in at least one species
    #'MD': exon marked as having missing data
    #return 'NA' if exon does not present in any species
    species = gene.findall("species")
    numAlnSpecies = 0
    hasMDspc = 0
    for spc in species:
        spcStatus = spc.get("status")
        mdExonsStr = spc.get("missingdataexons")
        mdExons = mdExonsStr.split(',')
        #sys.stderr.write("*%s*\n" %(mdExonsStr))

        #if spcStatus != "incomplete" and spcStatus != "NA":#species is ok ==> all exon is ok (exon that ovelapped with missing data)
        #    continue
        hasexon = "no"
        exons = spc.getiterator("exon")
        exonstatus = "NA"
        for exon in exons:
            if exon.get("id") == exonid: #found exon in current species!
                hasexon = "yes"
                #if exon.get("status") == "1":
                if exon.get("status") == "1":
                    exonstatus = 'C'
                    break

        if hasexon == "yes": #found exon in current species
            numAlnSpecies += 1 
            if exonstatus != 'C':#didn't find any complete exon
                if exonid in mdExons: #but found exon that was marked as having missing data
                    hasMDspc += 1
                else:#exon is non-conserved in current species --> exon is non-conserved in all
                    return "N"
        #elif spcStatus == "MD":#doesn't have the exon, but species status = P
        else:#doesn't have the exon, but species status = P
            hasMDspc += 1

    if numAlnSpecies == 0:#does not present in any species
        return "NA"
    elif numAlnSpecies != len(species): #absent in at least one species
        return "N"
    elif hasMDspc > 0:#only have exon marked missing data in some species or exon absent in some species
        return "MD"
    else:#present in all species
        return "C"

def getMultizExonStatus(exonid, gene):
    #return no if any of the involved species has an exonStatus "no"
    if not gene:
        return "NA"
    numAlnSpecies = 0
    hasMDspc = 0
    for species in gene:
        exon = getExon(species, exonid)
        if exon.alnBases > 0:
        #if exon:
            numAlnSpecies += 1
            #if exon.status == "NA":
            #    continue
            if exon.status == "N":
                return "N"
            elif exon.status == "MD":
                hasMDspc += 1
    if numAlnSpecies == 0:
        return 'NA'
    elif numAlnSpecies < len(gene):
        return 'N'
    elif hasMDspc > 0:
    #elif hasMDspc > 0 or numAlnSpecies < len(gene):
        return 'MD'
    else:
        return "C"

def cactusVSmultiz_perSpecies(multizGenes, cactusTree, name2id, allowDup, fh):
    cactusRoot = cactusTree.getroot()
    cGenes = cactusRoot.findall("gene") #cactus genes

    fh.write("Chr\tStart\tEnd\tSpecies\tGeneName\tGeneAccession\tCactusStatus\tMultizStatus\n")
    
    for cgene in cGenes:
        species = cgene.findall("species")
        genename = cgene.get("name")
        exonCount = int(cgene.get("exonCount"))
        mgene = getGene(multizGenes, genename) #corresponding multiz gene
        
        for cspc in species:
            cstatus = cspc.get("status")
            #relax the definition of cactus complete homolog by allowing for duplicated pieces of exons in between exons:
            if allowDup and cstatus == "N":
                cstatus = cactus_getSpcStatus2(cspc, exonCount)
            
            if not mgene:
                mstatus = "NA"
            else:
                mspc = getSpecies(mgene, cspc.get("name"))
                if not mspc:
                    mstatus = "NA"
                else:
                    mstatus = mspc.status

            fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(mgene.chr, mgene.start, mgene.end, cspc.get("name"), genename, name2id[genename], cstatus, mstatus))
        
    return


def cactusVSmultiz(multizGenes, cactusTree, name2id, allowDup):
    cactusRoot = cactusTree.getroot()
    cGenes = cactusRoot.findall("gene") #cactus genes
    exonTab = {}#first column is cactus, second is multiz
    geneTab = {}#first column is cactus, second is multiz

    sys.stdout.write("Chr\tStart\tEnd\tExonId\tGeneName\tGeneAccession\tCactusStatus\tMultizStatus\n")
    for cgene in cGenes:
        genename = cgene.get("name")

        mgene = getGene(multizGenes, genename) #corresponding multiz gene
        cgstatus = cgene.find("status").text

        #relax the definition of cactus complete homolog by allowing for duplicated pieces of exons in between exons:
        if allowDup and cgstatus == "N":
            cgstatus = cactus_getGeneStatus2(cgene)

        if not mgene:
            mgstatus = "NA"
        else:
            mgstatus = mgene.status
       
        #Missing data is consider non-conserved if one of the methods (multiz or cactus) has a conserved homolog
        if cgstatus == "MD" and mgstatus == "C":
            cgstatus = "N"
        if mgstatus == "MD" and cgstatus == "C":
            mgstatus = "N"
        #NA (not-aligned) is considered as non-conserved
        if cgstatus == "NA":
            cgstatus = "N"
        if mgstatus == "NA":
            mgstatus = "N"

        #comparing the gene status
        gStatus = cgstatus + mgstatus
        if gStatus not in geneTab:
            geneTab[gStatus] = 1
        else:
            geneTab[gStatus] += 1
       
        sys.stdout.write("%s\t%s\t%s\tgene\t%s\t%s\t%s\t%s\n" %(mgene.chr, mgene.start, mgene.end, genename, name2id[genename], cgstatus, mgstatus))
        
        #comparing exon status:
        exonCount = int(cgene.get("exonCount"))
        for id in range(exonCount):
            cexonStatus = getCactusExonStatus(str(id), cgene)
            mexonStatus = getMultizExonStatus(str(id), mgene)
            
            #Missing data is consider non-conserved if one of the methods (multiz or cactus) has a conserved homolog
            if cexonStatus == "MD" and mexonStatus == "C":
                cexonStatus = "N"
            if mexonStatus == "MD" and cexonStatus == "C":
                mexonStatus = "N"
            
            if cexonStatus == "NA":
                cexonStatus = "N"
            if mexonStatus == "NA":
                mexonStatus = "N"

            exonStatus = cexonStatus + mexonStatus
            
            mexon =  multizgene_getExon(mgene, str(id))
            if mexon:
                sys.stdout.write("%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\n" %(mexon.chr, mexon.start, mexon.end, id, genename, name2id[genename], cexonStatus, mexonStatus))
            else:
                sys.stdout.write("\t\t\t%d\t%s\t%s\t%s\t%s\n" %(id, genename, name2id[genename], cexonStatus, mexonStatus))
            
            if exonStatus not in exonTab:
                exonTab[exonStatus] = 1
            else:
                exonTab[exonStatus] += 1

    return (exonTab, geneTab)

#================ get cactus gene/exon status of species of interest ========
def getMultizGeneSpecies(genes, genename, species):
    for gene in genes:
        if gene.name == genename: #skip genes with missing data if missingData is not allowed
            for spc in gene:
                if spc.name == species:
                    return spc
    return None

def cactusSpcExonStatus(spc):
    mdExonsStr = spc.get("missingdataexons")
    mdExons = mdExonsStr.split(',')
    
    exon2status = {}
    exons = spc.getiterator("exon")
    for exon in exons:
        exonid = exon.get("id")
        status = exon.get("status")

        if exonid not in exon2status:
            exon2status[exonid] = status
        if (status == "1"):
            exon2status[exonid] = "C"
    
    for id in exon2status:
        if exon2status[id] == "0":
            if id in mdExons:
                exon2status[id] = "MD"
            else:
                exon2status[id] = "N"
    return exon2status

#================ multiz species to genes info ===========
#def multizSpc2genes(genes, cactusTree, mdgenes, allSpcGenes, strict, allowDup):
#def multizSpc2genes(genes, cactusTree, allSpcGenes, strict, allowDup):
def multizSpc2genes(genes, cactusTree, allowDup):
    spc2total = {} #key = species, val = {total number of genes this species covered}
    spc2cons = {} #key = species, val = {number of conserved genes}
    root = cactusTree.getroot()
    cgenes = root.findall("gene")
    
    for gene in genes:
        cgene = cactus_getGene(cgenes, gene.name)
        for spc in gene:
            cspc = cactus_getSpecies(cgene, spc.name)
            cstatus = cspc.get("status")
            if allowDup and cstatus == "N":
                #cstatus = cactus_getGeneStatus2(cgene)
                cstatus = cactus_getSpcStatus2(cspc, spc.exonCount)

            status = spc.status
            if (status == "MD" and cstatus == "C") or status == "NA":
                status = "N"

            if status != "C" and cstatus == "MD":
                continue
            
            if status == "C" or status == "N":
                if spc.name not in spc2total:
                    spc2total[spc.name] = 1
                    spc2cons[spc.name] = 0 #initialize spc2cons
                else:
                    spc2total[spc.name] += 1

                if status == "C":#if conserved, report
                    spc2cons[spc.name] +=1

    return (spc2total, spc2cons)


#def multizSpcExonCount(spc, gene, cactusSpc, allSpcExons, strict):
def multizSpcExonCount(spc, cactusSpc):
    total = 0
    cons = 0
   
    cactusExon2status = cactusSpcExonStatus(cactusSpc)
    for exon in spc:
        status = exon.status
        if status == "MD" and (exon.id in cactusExon2status) and cactusExon2status[exon.id] == "C":
            status = "N"
        if status == "NA":
            status = "N"
        if status != "C" and (exon.id in cactusExon2status) and cactusExon2status[exon.id] == "MD":
            continue
        if status == "C" or status == "N":
            total +=1
            if status == "C":
                cons +=1
    return (total, cons)

#def multizSpc2exons(genes, cactusTree, mdgenes, allSpcExons, strict):
def multizSpc2exons(genes, cactusTree):
    spc2total = {}
    spc2cons = {}
    
    root = cactusTree.getroot()
    cgenes = root.findall("gene")
    
    for gene in genes:
        #spc2exon2status = {}
        cgene = cactus_getGene(cgenes, gene.name)
        for spc in gene:
            cspc = cactus_getSpecies(cgene, spc.name)
            (total, cons) = multizSpcExonCount(spc, cspc)
            #spc2exon2status[spc.name] = exon2status
            if spc.name not in spc2total:
                spc2total[spc.name] = total
                spc2cons[spc.name] = cons
            else:
                spc2total[spc.name] += total
                spc2cons[spc.name] += cons

        #for i in range(gene.exonCount):
        #    if spc2exon2status["panTro2"][str(i)] != "C" and spc2exon2status["rheMac2"][str(i)] == "C":
        #        spc = getSpecies(gene, "rheMac2")  
        #        mexon = getExon(spc, str(i))
        #        sys.stdout.write("%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\n" %(mexon.chr, mexon.start, mexon.end, i, gene.name, mexon.gene, spc2exon2status["panTro2"][str(i)], spc2exon2status["rheMac2"][str(i)]))

    return (spc2total, spc2cons)

#================== cactus species to genes info ==============
#def cactusSpc2genes(tree, multizGenes, allowMissingData, allSpcGenes, strict, allowDup):
def cactusSpc2genes(tree, multizgenes, allowDup):
    spc2total = {}
    spc2cons = {}
    root = tree.getroot()
    genes = root.findall("gene")
    
    if len(genes) <= 0:
        return None

    #initialized the hashes:
    species = genes[0].findall("species")
    for spc in species:
        spc2total[spc.get("name")] = 0
        spc2cons[spc.get("name")] = 0
    
    for gene in genes:
        mgene = getGene(multizgenes, gene.get("name"))
        
        species = gene.findall("species")
        for spc in species:
            spcname = spc.get("name")
            status = spc.get("status")
            if allowDup and status == "N":
                status = cactus_getSpcStatus2(spc, mgene.exonCount)

            mspc = getSpecies(mgene, spcname) #multiz species
            if status != "C" and mspc.status == "MD":
                continue

            if (status == "MD" and mspc.status == "C") or status == "NA":
                status = "N"

            if status == "C" or status == "N":
                #if status == "N":
                #    if allowDup:
                #        status = cactus_getSpcStatus2(spc, int(gene.get("exonCount")))
                
                #if status == "MD" and mspc.status == "C":
                #    status = "N"

                if status == "C" or status == "N":
                    spc2total[spcname] += 1
                if status == "C":
                    spc2cons[spcname] += 1

    return (spc2total, spc2cons)

def cactusSpcExonCount(spc, multizspc):
    total = 0
    cons = 0
    exon2status = {}
    exons = spc.getiterator("exon")
    
    mdExonsStr = spc.get("missingdataexons")
    mdExons = mdExonsStr.split(',')
    
    total = 0
    for exon in exons:
        exonid = exon.get("id")
        status = exon.get("status")

        if exonid not in exon2status:
            exon2status[exonid] = status
        if (status == "1"):
            exon2status[exonid] = 'C'
    
    for id in range(multizspc.exonCount):
        if str(id) not in exon2status:
            exon2status[str(id)] = "N"

    for id in exon2status:
        mexon = getExon(multizspc, id)
        cstatus = exon2status[id]
        if cstatus == "0":
            if (id in mdExons) and mexon.status != "C":
                cstatus = "MD"
            else:
                cstatus = "N"
        if mexon.status == "MD" and cstatus != "C":
            continue
        if cstatus == "C" or cstatus == "N":#exon is complete
            total +=1
            if cstatus == "C":
                cons += 1
    #sys.stderr.write("spc %s, total %d, broken %d\n" %(spc.get("name"), total, broken))
    return (total, cons)

    
#def cactusSpc2exons(tree, multizGenes, allowMissingData, allSpcExons, strict):
def cactusSpc2exons(tree, multizgenes):
    spc2total = {}
    spc2cons = {}
    root = tree.getroot()
    genes = root.findall("gene")

    if len(genes) <= 0:
        return None

    #initialized the hashes:
    species = genes[0].findall("species")
    for spc in species:
        spc2total[spc.get("name")] = 0
        spc2cons[spc.get("name")] = 0
    
    for gene in genes:
        species = gene.findall("species")
        mgene = getGene(multizgenes, gene.get("name"))
        for spc in species:
            spcname = spc.get("name")
            mspc = getSpecies(mgene, spcname)
            (total, cons) = cactusSpcExonCount(spc, mspc) 
            spc2total[spcname] += total
            spc2cons[spcname] += cons

    return (spc2total, spc2cons)

def printSpc2Genes(f, title, totalM, consM, totalC, consC):
    #f.write("%s\n" %(title))
    for spc in totalM:
        #sys.stderr.write("species: %s\n" %(spc))
        tmultiz = totalM[spc]
        cmultiz = consM[spc]
        if tmultiz == 0:
            sys.stderr.write("species %s, total (%s) of multiz is 0\n" %(spc, title))
        #percentM = float(bmultiz)/tmultiz*100.0
        percentM = float(cmultiz)/tmultiz*100.0
        tcactus = totalC[spc]
        ccactus = consC[spc]
        if tcactus == 0:
            sys.stderr.write("species %s, total (%s) of cactus is 0\n" %(spc, title))
        #percentC = float(bcactus)/tcactus*100
        percentC = float(ccactus)/tcactus*100
        #f.write("%s\t%s\t%d\t%d\t%.2f\t%d\t%d\t%.2f\n" %(title, spc, bcactus, tcactus, percentC, bmultiz, tmultiz, percentM))
        f.write("%s\t%s\t%d\t%d\t%.2f\t%d\t%d\t%.2f\n" %(title, spc, ccactus, tcactus, percentC, cmultiz, tmultiz, percentM))
    f.write("#\n")
    return


#def species2genes(f, multizGenes, cactusTree, allowMissingData, allSpcGenes, allSpcExons, strict, allowDup):
def species2genes(f, multizGenes, cactusTree, allowDup):

    #speciesToGenes:
    (gTotalMultiz, gConsMultiz) = multizSpc2genes(multizGenes, cactusTree, allowDup)
    (gTotalCactus, gConsCactus) = cactusSpc2genes(cactusTree, multizGenes, allowDup)
    printSpc2Genes(f, "spcGeneTab", gTotalMultiz, gConsMultiz, gTotalCactus, gConsCactus)
    
    (eTotalCactus, eConsCactus) = cactusSpc2exons(cactusTree, multizGenes)
    (eTotalMultiz, eConsMultiz) = multizSpc2exons(multizGenes, cactusTree)
    printSpc2Genes(f, "spcExonTab", eTotalMultiz, eConsMultiz, eTotalCactus, eConsCactus)
    return

def printTab(f, title, tab, strict):
    #f.write("%s\n" %(title))
    totalCount = 0
    selectedKeys = ["CC", "NN", "CN", "NC"]
    for key in tab:
        if not strict or key in selectedKeys:
            totalCount += int(tab[key])

    for key in tab:
        if not strict or key in selectedKeys:
            percent = float(tab[key])/totalCount*100
            f.write("%s\t%s\t%d\t%.2f\n" %(title, key, tab[key], percent))
            
    f.write("#\n")

def dupStats(tree):
    root = tree.getroot()
    genes = root.findall("gene")
    total = {"gene":0, "geneexon": 0, "partialdupgene":0, "exon":0, "intron":0}
    dup = {"gene":0, "geneexon": 0, "partialdupgene":0, "exon":0, "intron": 0}

    for gene in genes:
        numDupExon = 0
        genedup = False
        geneexondup = False
        total["gene"] += 1
        total["geneexon"] += 1
        total["partialdupgene"] += 1
        exonCount = int(gene.get("exonCount"))
        total["exon"] += exonCount
        total["intron"] += exonCount -1
        intronStr = gene.find("intronDups").text
        if intronStr:
            intronDups = intronStr.split(',')
        else:
            intronDups = []

        exonStr = gene.find("exonDups").text
        if exonStr:
            exonDups = exonStr.strip(',').split(',')
        else:
            exonDups = []
        
        if len(exonDups) != exonCount:
            sys.stderr.write("ERR\n")
        for i in range(exonCount):
            if exonDups[i] == "1":
                genedup = True
                geneexondup = True
                dup["exon"] += 1
                numDupExon +=1
        for i in range(exonCount -1):
            if intronDups[i] == "1":
                genedup = True
                dup["intron"] +=1
        if genedup:
            dup["gene"] +=1
        if geneexondup:
            dup["geneexon"] +=1
            if numDupExon < exonCount:
                dup["partialdupgene"] +=1

    for key in total:
        t = total[key]
        d = dup[key]
        if t == 0:
            sys.stderr.write("total %s 0\n" %(key))
        percent = float(d)/t*100.0
        sys.stdout.write("%s\t%d/%d\t%.2f%%\n" %(key, d, t, percent))
    return



def usage():
    sys.stderr.write("Usage: cactusVSmultiz.py geneNameToRefseqIdFile multizMarkFile cactusXMLFile\n")
    sys.exit(2)

def main():
    if len(sys.argv) < 5:
        sys.stderr.write("Too few input arguments.\n")
        usage()
    usg = "usage: %prog [options] txId2nameFile consEvidFile cactusHomolog.xml outputFile"
    parser = OptionParser(usage=usg)
    parser.add_option("-a", "--perSpecies", dest="perSpcOut", \
                      help="Print the status per species", default="perSpeciesCmp")
    #parser.add_option("-m", "--missingData", dest="missingData", action="store_true", \
    #                  help="If specified, will allow for missingData. Default is false.", default=False)
    #parser.add_option("-s", "--strict", dest="strict", action="store_true", \
    #                  help="If specified, only look at genes/exons that present in all species in both multiz and cactus. Default is false.", default=False)
    parser.add_option("-d", "--dup", dest="dup", action="store_true", \
                      help="If specified, relax the definition of cactus gene conservation and allow duplicated pieces of exons in between exons. Default is false.", default=False)
    (options, args) = parser.parse_args()

    (id2name, name2id) = mapId2name(args[0])

    multizGenes = readMultiz(args[1], id2name)
    setMultizGeneStatus(multizGenes)

    cactusTree = ET.parse(args[2])
    (exonTab, geneTab) = cactusVSmultiz(multizGenes, cactusTree, name2id, options.dup)
    f = open(args[3], "w")
    
    printTab(f, "exonTab", exonTab, True)
    printTab(f, "geneTab", geneTab, True)
    
    printTab(f, "exonTabSupplement", exonTab, False)
    printTab(f, "geneTabSupplement", geneTab, False)

    species2genes(f, multizGenes, cactusTree, options.dup)
    f.close()

    #====== print out the difference in gene conservation status per species =========
    fh = open(options.perSpcOut, "w")
    cactusVSmultiz_perSpecies(multizGenes, cactusTree, name2id, options.dup, fh)
    fh.close()
    
    #
    dupStats(cactusTree)

    #TEMP:
    #multizSpc2exons(multizGenes, [], options.strict)

if __name__ == "__main__":
    main()




