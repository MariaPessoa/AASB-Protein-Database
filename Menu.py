#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 18:35:05 2020

@author: mariajoao
"""

from DataBase import DataBase

def menu():
    try: 
        db = DataBase()
        db.ImportHistory()
    except: print("Unable to import into database")
    
    while True:
        print("""
        ======== Funções ============
        0. Sair
        1. Add sequence
        2. Add sequence(s) from file
        3. Search in the Data Base
        4. Search in the Data Base for a similar sequence
        5. Calculate frequency of a symbol ou sub-sequence in the Data Base
        6. Protein Blast
        7. Print Protein Blast results
        8. Search for Domains
        9. UPGMA
        10.Consensus
        11.Import History
        """)
        op = input("Escolha uma opção: ")
        if op == "1":
            db.AddSeq(input("Sequence: "))
        elif op == "2":
            db.LoadFile(input("File name: "),input("File type (Default = txt): "))
        elif op == "3":
            db.FindInDB(input("ID: "))
        elif op == "4":
            db.SimilarInDB(input("Query: "),input("g (Default = -2): "))
        elif op == "5":
            db.FreqInDB(input("Query: "),input("ID (Default = None): "))
        elif op == "6":
            db.ProteinBlast(input("ID: "), input("Program (Default = blastp): "),
                            input("Blast DB (Default = nr): "),
                            input("evalue (Default = 0.05): "),
                            input("Max number (Default = 10): "),
                            input("Filename (Default = 'blast.xml'): "))
        elif op == "7":
            db.BlastParse(input("Filename: "))
        elif op == "8":
            db.SearchDomains(input("ID: "))
        elif op == "9":
            db.DBupgma()
        elif op == "10":
            db.DBmultiAlign(input("ID's: "),input("g (Default = -2): "))
        elif op == "11":
            db.ImportHistory(input("Filename (Default ='history.txt'): "))
        elif op == "0":
            db.ExportHistory()
            break
        else:
            print("\n Invalid option.")

if __name__ == "__main__":
    menu()
