if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
#Installieren von Biostrings
library(Biostrings)
#einschalten von Biostrings
s1 <- "GAATC"
s2 <- "CATAC"
#definieren von zu alignierenden Paketen
sigma <- nucleotideSubstitutionMatrix(match = 10, mismatch = -5, baseOnly = TRUE)
#sigma ist eine Matrix für die Scores des alignments
sigma["T","C"] <- 0
sigma["C","T"] <- 0
sigma["A","G"] <- 0
sigma["G","A"] <- 0
#in den Vorgaben der Vorlesung wurden manche mismatch-Werte anders bewertet, dementsprechend werden die Werte hier geändert
sigma
#sigma anzeigen
alignment <- pairwiseAlignment(s1, s2, substitutionMatrix = sigma, gapOpening = -4, gapExtension = -4, scoreOnly = FALSE)
#jetzt wird ein alignment zwischen den beiden Sequenzen erzeugt
alignment
#Anzeigen des alignment mit Score
#alignment wird angezeigt

BiocManager::install("DECIPHER")
#Installieren von Decipher
library(DECIPHER) 
# festlegen des Verzeichnisses
# fas <- "C:/Users/Andrea Nino/Desktop/ (1)/Studium/Biologie/Bioinformatik/Übung 8/FOXA1_seq.fasta" 
# Dateipfadvariable für die verwendete Datei(andere Möglichkeit für eigene Daten)
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
#Dateipfadvariable für die verwendete Datei, mit system.file wird eine Beispieldatei aufgerufen
dna <- readDNAStringSet(fas)
#in dna werden die Sequenzen aus der Datei eingelsesen
dna
#Sequenzen werden angezeigt

AA <- AlignTranslation(dna, type="AAStringSet") 
#Sequenzen werden als Proteine übersetzt bewertet aneinandergelegt (aligned) Ausgabe sind Proteinsequenzen
BrowseSeqs(AA, highlight=1) 
#Proteinsequenzen werden als html angezeigt
DNA <- AlignSeqs(dna) 
BrowseSeqs(DNA, highlight=1) 
#DNSSequenzen der Ausgangsdatei werden wieder bewertet aneinandergelegt (aligned) Ausgabe sind DNS Sequenzen
DNA <- AlignTranslation(dna) 
BrowseSeqs(DNA, highlight=1) 
#Sequenzen werden auf Basis der Proteinsequenzen aligniert
writeXStringSet(DNA, file="C:\Users\Andrea Nino\Desktop\ (1)\Studium\Biologie\Bioinformatik\Übung 8")
#Ergebnis wird in eine fasta-Datei geschrieben

#folgende Schritte sind nur für Sequenzbündel, die doppelte Sequenzen beinhalten, also für diese Aufgabe unwichtig

u_dna <- unique(dna) 
#DNS doppelte DNS Sequenzen werden entnommen
index <- match(dna, u_dna) 
#Index für doppelte Sequenzen
U_DNA <- AlignSeqs(u_dna) 
#Sequenzen werden aligniert
DNA <- U_DNA[index]
#DNA wird auf einzigartige Sequenzen beschränkt
names(DNA) <- names(dna) 
#Namen von ursprünglicher Datei werden übernommen
DNA
#Ausgabe der Alignments
