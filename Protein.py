# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 11:05:57 2019

@author: liamo
"""


class Protein:
    def __init__(self, seq):
        self.alfabet = "ACDEFGHIKLMNPQRSTVWY_"
        self.seq = seq.upper()      
        self.otherDBs = {}
        self.source = None
    
    def __eq__(self,obj): 
        if isinstance(obj, Protein):
            return self.seq == obj.seq
        else:
            return False
    
    def __str__(self):
        return self.seq
    
    def __getitem__(self,index):
        return self.seq[index]
    
    def __len__(self): 
        return len(self.seq)
    
    def __getslice__(self, i, j):
        return self.seq[i:j]
        
    def isvalid(self):
        res = True
        i = 0
        while i < len(self.seq) and res:
            if self.seq[i] not in self.alfabet: res = False
            else: i += 1
        return res
    
    def addotherDBs(self,database,dbID):
        if database not in self.otherDBs.keys():
            self.otherDBs[database] = [dbID]
        else: self.otherDBs[database].append(dbID)
        return None
    
    def setsource(self,organism):
        self.source = organism
        return None
