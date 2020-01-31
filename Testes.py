
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 15:56:13 2019

@author: liamo
"""   
from DataBase import DataBase
pts = ["ACDEFGHI","DEFGHIKL","FGHIKLMN"]


def test_AddSeq():
    db = DataBase()
    db.AddSeq("MEIDKFVKEEDIPFEYGVVRERDNAVSWSRYL") #Valid
    db.AddSeq("MEIDKFVKEE111VRERDNAVSWSRYL") #Invalid 
    db.AddSeq("MEIDKFVKEEDIPFEYGVVRERDNAVSWSRYL") #Repeated
    for i in db.dseqs.values(): print(i.seq)
    
def test_BlastParse():
    #ProteinBlast > BlastParse > BlastFetch pipeline
    db = DataBase() 
    #db.ProteinBlast(1)
    ANs = db.BlastParse("Blast_mitf_on.xml") #dict = {ANs : description}
    for i in ANs.keys(): print(i,"\n",ANs[i],end="\n\n")
    db.ProteinFetch(list(ANs.keys())[:3])     #choose results
    print("Sequences added: ",len(db.dseqs))

def test_Clustal():
    db = DataBase()
    for i in pts: db.AddSeq(i)
    db.SMFromFile("blosum62.mat")
    db.Clustal(path,"testfasta.txt")

def test_DBupgma():
    db = DataBase() 
    db.LoadFile("testfasta.txt", filetype="fasta")
    db.SMFromFile("blosum62.mat")
    db.DBupgma()
    db.tree.printtree()
    
def test_DBmultiAlign():
    db = DataBase() 
    for i in pts: db.AddSeq(i)
    db.SMFromFile("blosum62.mat")
    db.DBupgma()
    print(db.DBmultiAlign([1,2]))
    
def test_FindInDB(): 
    db = DataBase() 
    db.LoadFile("testfasta.txt", filetype="fasta")
    db.FindInDB(2)
  
def test_FreqInDB():
    global pts
    db = DataBase()
    for i in pts: db.AddSeq(i)
    print(db.FreqInDB("FG"))
    print(db.FreqInDB("FG",1))
    #print(db.FreqInDB("*")) #gives error
    print(db.FreqInDB("O"))

def test_history():
    global pts
    db = DataBase()
    psmotif = "AAGSGGAAGQAASAAAGAGKGLAA"
    db.AddSeq(psmotif)
    for i in pts: db.AddSeq(i)
    db.SMFromFile("blosum62.mat")
    db.DBupgma()
    #db.SearchDomains(1) #slow but works
    db.ExportHistory()
    db2 = DataBase()
    db2.ImportHistory()
    print(db2.domains)
    
def test_LoadFile(): 
    db = DataBase() 
    print(db.LoadFile("testfasta.txt", filetype="fasta"))
    for i in db.dseqs.values(): print(i.seq,end="\n\n")
    
def test_RegexInDB():
    global pts
    db = DataBase()
    for i in pts: db.AddSeq(i)
    print(db.RegexInDB("FGHI")) #list of tuples (key in db, index in seq)
    print(db.RegexInDB("FGHI",1)) #list of indexes in seq
    print(db.RegexInDB("E*FGHI"))

def test_SaveFile():
    global pts
    db = DataBase()
    for i in pts: db.AddSeq(i)
    db.SaveFile()
    
def test_SearchDomains():
    psmotif = "AAGSGGAAGQAASAAAGAGKGLAA"
    db = DataBase() 
    db.AddSeq(psmotif)
    db.SearchDomains(1)
    print(db.domains)
    
def test_SimilarInDB():
    global pts
    db = DataBase()
    for i in pts: db.AddSeq(i)
    db.SMFromFile("blosum62.mat")
    print(db.SimilarInDB("ACDEFGHI"))
    
def test_SMFromFile():
    db = DataBase()
    db.SMFromFile("blosum62.mat")
    print(db.sm)
    print(db.alphabet)
    
def menu():
    try: 
        db = DataBase()
        db.ImportHistory()
    except: print("Unable to import into database")
    
    while True:
        print("""
        ======== Funções ============
        0. Sair
        1.DataBase.AddSeq
        2.DataBase Blast
        3.DataBase.DBupgma
        4.DataBase.DBmultiAlign
        5.DataBase.FindInDB
        6.DataBase.FreqInDB
        7.DataBase history
        8.DataBase.LoadFile
        9.DataBase.RegexInDB
        10.DataBase.SaveFile
        11.DataBase.SearchDomains
        12.DataBase.SimilarInDB
        13.DataBase.SMFromFile
        """)
        op = input("Escolha uma opção: ")
        if op == "1":
            test_AddSeq()
        elif op == "2":
            test_BlastParse()
        elif op == "3":
            test_DBupgma()
        elif op == "4":
            test_DBmultiAlign()
        elif op == "5":
            test_FindInDB()
        elif op == "6":
            test_FreqInDB()
        elif op == "7":
            test_history()
        elif op == "8":
            test_LoadFile()
        elif op == "9":
            test_RegexInDB()
        elif op == "10":
            test_SaveFile()
        elif op == "11":
            test_SearchDomains()
        elif op == "12":
            test_SimilarInDB()
        elif op == "13":
            test_SMFromFile()
        elif op == "0":
            db.ExportHistory()
            break
        else:
            print("\n Invalid option.")

if __name__ == "__main__":
    menu()
    