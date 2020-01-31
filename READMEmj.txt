Method descriptions
Requires the Bioython package installed and all the python files provided.

####Protein Class####
The Protein class is meant to store information in an organized way, so that the main DataBase class can access it to perform the tasks required. However, besides the built-in sequence validity check (isvalid method) it also allows the user to manually add external database IDs (addotherDBs method) and set the source organism of a sequence (setsource method). This last method is not used in the DataBase class, even in the BLAST method pipeline.

####DataBase Class####
Methods:
AddSeq - The AddSeq method allows manual insertion of one sequence, which is passed as the argument. The function checks the sequences' validity by comparison with the "alphabet" class attribute using the Protein class' isvalid method, adding it do the class attribute dictionary "dseqs" with an auto-incremented 1-based ID as the key and a Protein class object with the sequence as the value. The method also checks for repeat sequences and returns a Boolean value.

LoadFile - The LoadFile method loads sequences from files of 2 formats, text and FASTA. The method receives 3 arguments, the file name to open, the type of file ("text" or "fasta", default is "text") and an optional separator argument (no default). As this method does not call the AddSeq method, it allows repeat sequences to be added to the database, but does check for sequence validity. The method returns a Boolean value.

SaveFile - The SaveFile method saves the sequences in the database (dictionary class attribute "dseqs") in fasta format. It receives the name of the file as an optional argument. The default file name is "Proteins.txt".

SMFromFile - The SMFromFile method creates a substitution matrix using the file name passed as an argument and saves it in the class attribute "sm". The method has an optional separator argument, with the default set to "tab" ("\t"). 

FindInDB - If the database ID passed as the argument is valid (that is, exists in the database) the FindInDB method prints the sequence size and the first 1000bp (or the entire sequence if it's smaller). Then the function offers the user the option to update the sequence or to add external database IDs.

SimilarInDB - The function searches for the sequence most similiar to the query passed as the argument, using SmithWaterman's local alignment algorithm. The function returns the database ID and the score of the best alignment (in a tuple).

RegexInDB - The function takes two arguments, a regular expression to search for and an optional database ID to select a particular sequence. If no ID is passed, the function searches the entire database. The function returns a list of tuples with the database ID and start position if searching the entire database, or a list with the start position(s) in the sequence corresponding to the ID passed.

FreqInDB - The FreqInDb function takes the query sequence or expression and calculates its frequency in the database (number of total occurrences). If no database ID is passed, the function searches in the entire database. The funtion prints the frequency and returns a dictionary with database IDs as keys and the number of occurrences of the query per sequence as values.

The ProteinBlast, BlastParse and BlastFetch functions are designed to be used in pipeline, but can also work independently. The purpose of this split into 3 functions was to permit separate functionality and allow the user to use previous BLAST results (saved in an XML file), for example.

ProteinBlast - The function receives the database ID of the sequence to use in the BLAST, the program (default is "blastp"), the external database to search in (default is "nr" or non redundant proteins), the e-value (default is 0.05), the max number of results and a file name (an XML file) for the saved BLAST results. All arguments except the database ID are optional.
The function does not return anything, instead saving all results in the specified XML file or under the default name ("blast.xml").

BlastParse - The BlastParse method takes only one argument, the name of the XML file with the BLAST results. It makes use f Biopython's XML file parser and returns a dictionary with the accession numbers obtained as keys and a small discription as values. This is meant to allow the user to choose which sequences to save with a little more information than simply their accession numbers, minimizing the number of requests to the Entrez database.

BlastFetch - The BlastFetch method takes a list, with one or more elements, as an argument. This list should contain valid accession numbers, which BlastFetch will attempt to retrieve. Currently, the function does not support specification of start and stop positions for a single sequence. The function fetches the sequence in FASTA format using Biopython's Entrez module, saves them in an output file ("BlastFetch.txt") and then loads them into the database using the LoadFile method.

SearchDomains - The SearchDomains method takes only one argument, a valid database ID. The method uses Biopython's ExPASy module to scan Prosite's database and returns the motifs found (or an empty list) after adding said results to the class attribute dictionary "domains".

DBupgma - The DBupgma method takes no arguments. The method requires, however, a substitution matrix to be loaded beforehand. Using the UPGMA algorithm, the method creates a guide tree and returns it after adding it to the class attribute "tree" (BinaryTree object).

DBmultiAlign - The DBmultiAlign method takes 2 arguments, a list with the database IDs (the sequences to align) and the gap penalty value (default is -2). The method requires a substitution matrix to be loaded beforehand and  returns the alingment consensus.

ExportHistory - The ExportHistory method takes only one optional argument, the name of the file to export to (default is "history.txt"). The method calls the SaveFile method, saving all the sequences in the database and then manually saves the information stored via the class attribute Queue "history". The method does not return anything and the information stored in the history Queue is not lost.

!!! The ImportHistory has only been used with the files created by ExportHistory and does not check if said files exist when called. !!!

ImportHistory - The ImportHistory method takes only on optional argument, the name of a file resulting from the ExportHistory method. The method first tries to open the "Proteins.txt" file resulting from the SaveFile method, to load the sequences into the database. Then the method loads the information stored in the export file (substitution matrix, search domains, guide tree) passed as argument (default is "history.txt"). This method does not return anything.

Clustal - The Clustal method takes two arguments, the cluster path and the name of the file with the sequences. The method does the alignment of the sequences using an external program, CLUSTAL.


