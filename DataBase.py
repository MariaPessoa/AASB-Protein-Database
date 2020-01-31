# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 15:40:08 2019

@author: liamo
"""
from re import finditer
from Protein import Protein       
from SubstMatrix import SubstMatrix
from AlignSeq import AlignSeq
from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.ExPASy import ScanProsite
from upgma import UPGMA
from MultipleAlign import MultipleAlign
from StackQueue import Queue
import os
from Bio.Align.Applications import ClustalwCommandline


class DataBase:
    def __init__(self):
        self.dseqs = {}
        self.alphabet = "ACDEFGHIKLMNPQRSTVWY_"
        self.sm = None
        self.matdist = []
        self.email = None
        self.tree = None
        self.domains = {}
        self.history = Queue()
        
    def __str__(self):
        res = ""
        for pt in self.dseqs.values():
            res += "\n" + pt.seq
        return res
    
    def __len__(self):
        return len(self.dseqs)
    
    def set_email(self):
        self.email = str(input("Entrez email: "))
        return None
    
    def AddSeq(self,pseq):
        pt = Protein(pseq)
        if pt.isvalid(): #check validity
            for Ptn in self.dseqs.values(): #check repeat sequences
                if pt.seq == Ptn.seq: return False 
            self.dseqs[len(self.dseqs)+1] = pt #1-based integers
            return True
        return False
    
    def LoadFile(self, filename, filetype="text", sep=None):
        f = open(filename, "r")
        lines = f.readlines()
        if sep is not None: lines = lines.split(sep)
        if filetype == "fasta":
            title = None
            temp = []
            for line in lines:
                if line.startswith(">"):
                    if title is None:
                        title = line
                    else: 
                        title = line
                        pt = Protein("".join(temp))
                        if pt.isvalid: #only adds valid seqs
                            self.dseqs[len(self.dseqs)+1] = pt
                        temp = []
                else: 
                    temp.append(line.strip("\n"))
            if len(temp) != 0:
                pt = Protein("".join(temp))
                if pt.isvalid: #only adds valid seqs
                    self.dseqs[len(self.dseqs)+1] = pt
            self.history.enqueue(("LoadFile",filename,"fasta",sep))
        elif filetype == "text":
            for line in lines:
                pt = Protein(line) 
                if pt.isvalid(pt.seq): self.dseqs[len(self.dseqs)+1] = line
        else: return False
        return True
    
    def SaveFile(self,filename="Proteins.txt"): #fasta format
        with open(filename,"w") as f:
            for k,v in self.dseqs.items():
                f.write(f">SEQ{k}\n")
                f.write(f"{v.seq}\n")
        return None
    
    def SMFromFile(self, filename, sep="\t"):
        self.sm = SubstMatrix()
        self.alphabet = self.sm.loadFromFile(filename,sep)
        if sep == "\t":
            self.history.enqueue(("SMFromFile",filename))
        else: self.history.enqueue(("SMFromFile",filename,sep))
        return None
    
    def FindInDB(self,dbID): 
        if dbID in self.dseqs.keys(): 
            print("Sequence size:",len(self.dseqs[dbID]))
            print("Sequence (First 1k bp)):",
                  self.dseqs[dbID].seq[min(1000,len(self.dseqs[dbID].seq))])
            if input("Update Sequence? [y/n]: ") == "y":
                new_seq = input("New sequence (key number maintained): ")
                pt = Protein(new_seq)
                if pt.isvalid: self.dseqs[dbID] = pt
                else: print("Invalid sequence")
            if len(self.dseqs[dbID].otherDBs.values()) > 0:
                print(self.dseqs[dbID].otherDBs)
            else: print("No external DB IDs")
            if input("Update external DB IDs? [y/n]: ") == "y":
                new_db = input("Database: ")
                new_id = input("Database ID: ")
                if len(new_db) > 2 and len(new_id) > 2:
                    self.dseqs[dbID].addotherDBs(new_db,new_id)
                else: print("Invalid input")
        else: 
            print("Invalid ID")
            return False
    
    def SimilarInDB(self,query, g=-2): 
        if self.sm is None:
            align = AlignSeq(self.sm, g)
            maxseq = -999
            for k in self.dseqs.keys():
                if align.smithWaterman(query,self.dseqs[k].seq) > maxseq:
                    maxseq = align.smithWaterman(query,self.dseqs[k].seq)
                    key = k
            return (key,maxseq)
        else:
            print("No SubsMatrix loaded")
            return False
            
    def RegexInDB(self,query,seqid=None): 
        res = []
        if seqid is None:
            for i in self.dseqs.keys():
                mos = finditer(query, self.dseqs[i].seq) 
                for x in mos:
                    res.append((i,x.span()[0]))
        else:
            mos = finditer(query, self.dseqs[seqid].seq) 
            for x in mos:
                res.append(x.span()[0])
        self.history.enqueue(("RegexInDB",query,seqid))
        return res

    def FreqInDB(self,query,seqid=None): 
        res = {}
        if seqid is None: 
            for i in self.dseqs.keys():
                mos = finditer(query, self.dseqs[i].seq) 
                for x in mos:
                    if i not in res.keys(): res[i] = 1
                    else: res[i] += 1
        else:
            if seqid in self.dseqs.keys():
                mos = finditer(query, self.dseqs[seqid].seq) 
                res[seqid] = 0
                for x in mos: res[seqid] += 1
            else: 
                print("Invalid ID")
                return False
        print("Number of occurrences: ",len(res))
        return res

    def ProteinBlast(self,seqid,program="blastp",blastdb="nr",evalue=0.05,
              max_number=10, filename="blast.xml"):
        if self.email is None: self.set_email()
        Entrez.email = self.email
        if seqid in self.dseqs.keys(): 
            print("Running blast...")
            result_handle = NCBIWWW.qblast(program,blastdb, self.dseqs[seqid].seq,
                                   hitlist_size=max_number, expect=evalue)
            print("Blast done. Parsing...")
            with open(filename, "w") as out_handle:
                out_handle.write(result_handle.read())
        else: 
            print("Invalid ID")
            return False
                
    def BlastParse(self,filename):
        with open(filename) as file:
            blast_record = NCBIXML.read(file)
            #acessions = []
            allrec = {}
            for align in blast_record.alignments:
                #acessions.append(align.accession)
                allrec[align.accession] = align.hit_def
        return allrec 
    
    def ProteinFetch(self,query): #FALTA ADICIONAR AS otherDBs/source
        if type(query) == "list" and len(query) > 1:
            query = " ".join(query) 
        elif type(query) == "list" and len(query) == 1:
            query = query[0]
        fetch_handle = Entrez.efetch(db="nucleotide", rettype="fasta",
                                     retmode="text", id=query)
        protein_record = fetch_handle.read()
        with open("BlastFetch.txt", "w") as out_handle:
            out_handle.write(protein_record)
        self.LoadFile("BlastFetch.txt","fasta")
        return None
    
    def SearchDomains(self,seqid): 
        if seqid in self.dseqs.keys():
            fastaseq = ">SEQ1\n" + str(self.dseqs[seqid].seq)
            result_handle = ScanProsite.scan(seq=fastaseq)
            result = ScanProsite.read(result_handle)
            if seqid in self.domains.keys():
                self.domains[seqid].append(result)
            else: self.domains[seqid] = [result]
            return result
        else: 
            print("Invalid ID")
            return False
    
    def DBupgma(self):
        lseq = []
        for k in self.dseqs.keys(): lseq.append(self.dseqs[k].seq)
        if self.sm is not None:
            alseq = AlignSeq(self.sm, -1)
            up  = UPGMA(lseq, alseq)
            up.criaMatDists()
            self.history.enqueue(("DBupgma"))
            self.tree = up.run()
            return self.tree
        else:
            print("No Substitution Matrix")
            return False
    
    def DBmultiAlign(self,seqids,g=-2):
        if self.sm is None: 
            print("No Substitution Matrix")
            return False
        else:
            if self.tree is None: self.DBupgma()
            final = []
            order = self.tree.getCluster() #0-based
            for s in range(len(order)): 
                if (order[s]+1) in seqids: #1-based
                    final.append(self.dseqs[order[s]+1].seq)
            return MultipleAlign(final,AlignSeq(self.sm,g)).alignConsensus()
                    
    def ExportHistory(self,filename="history.txt"):
        self.SaveFile()
        temp = Queue()
        with open(filename,"w") as file:
            for k,v in self.domains.items():
                if v != []:
                    file.write(f"SearchDomains,{k},{v}")
            while self.history.isEmpty() is False: 
                temp.enqueue(self.history.front())
                file.write(str(self.history.dequeue()) + "\n")
        while temp.isEmpty() is False:
             self.history.enqueue(temp.dequeue())
        return None
    
    def ImportHistory(self,filename="history.txt"):
        self.LoadFile("Proteins.txt",filetype="fasta")
        with open(filename) as file: filelist = file.readlines()
        for line in filelist:
            line = line.strip("\n").strip("(").strip(")")
            temp = line.split(",")
            for i in range(len(temp)): temp[i] = temp[i].replace("'","").strip(" ")
            print(temp)
            if "SMFromFile" in line:   
                if len(temp) == 3: self.SMFromFile(temp[1],temp[2][1:])
                else: self.SMFromFile(temp[1])
                print("Substitution Matrix loaded")
            elif "SearchDomains" in line: 
                self.domains[temp[1]] = temp[2]
            elif "DBupgma" in line: 
                self.DBupgma()
                print("UPGMA evolutionary tree loaded")
            else: pass
        print("Done")
        return None    
    
    def Clustal(self,clustal_path,file):
        #Please fill in the software location correctly Ex:
        #clustalw_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"   
        clustalw_exe = clustal_path
        clustalw_cline = ClustalwCommandline(clustalw_exe, infile=file)
        assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
        stdout, stderr = clustalw_cline()

    