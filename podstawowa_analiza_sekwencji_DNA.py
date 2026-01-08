# -*- coding: utf-8 -*-
"""
Created on Sun Mar 16 15:48:09 2025

@author: magor
"""

# Otwórz plik z sekwencjami DNA
from Bio import SeqIO

for seq_records in SeqIO.parse("sequences.txt", "fasta"):
    print("\nPlik gotowy do pracy")

try:                                            # sprawdzenie czy plik istnieje
    with open("sequences.txt", "r") as plik:
        zawartosc = plik.read()
except FileNotFoundError:
    print("Nie znaleziono pliku z sekwencjami DNA")


# Zapisanie sekwencji z pliku jako lista
lista_sekwencji = []

for seq_records in SeqIO.parse("sequences.txt", "fasta"):
    lista_sekwencji.append(seq_records.seq)


# Zliczanie poszczególnych nukleotydów
def liczenie_nukleotydow(sekwencja):
    A = sekwencja.count("A")
    C = sekwencja.count("C")
    G = sekwencja.count("G")
    T = sekwencja.count("T")
    GC = sekwencja.count("GC")
    return (A, C, G, T, GC)

# Wywołanie funkcji dla każdej sekwencji z pliku
for i, seq in enumerate(lista_sekwencji):
    print(f"\nPrzetwarzam sekwencję: {i+1}: {seq},\n\nID sekwencji: {seq_records.id}") 
    print("\nDługosć sekwencji:", len(seq))
    
    a, c, g, t, gc = liczenie_nukleotydow(seq) 
    
    print("\nZliczone nukleotydy:")
    print(f"A: {a}")
    print(f"C: {c}")
    print(f"G: {g}")
    print(f"T: {t}")
    print(f"\nIlosć wystąpień par GC: {gc}")
    
# Procentowe wystąpienie danego nukleotydu
    procent_a = (a / len(seq)) * 100
    procent_c = (c / len(seq)) * 100
    procent_g = (g / len(seq)) * 100
    procent_t = (t / len(seq)) * 100
    procent_par_gc = gc / len(seq) * 100

    print(f"\nZawartosć procentowa nukleotydów\nA: {procent_a:.2f}%\nC: {procent_c:.2f}%\nG: {procent_g:.2f}%\nT: {procent_t:.2f}%")
    print(f"\nZawartosć procentowa par GC:\n{procent_par_gc:.2f}%")

# Pierwsze wystąpienie pary GC
    lokalizacja_sekwencji_gc = seq.find("GC")   
    print("\nIndeks pierwszego wystąpienia pary GC: ", lokalizacja_sekwencji_gc)


# Zapisanie wyników do Data Frame
import pandas as pd

# Utworzenie pustej listy z uzyskanymi danymi (wynikami)
dane = []

# Przetwarzanie sekwencji i zapis do listy
for i, seq_record in enumerate(SeqIO.parse("sequences.txt", "fasta")):
    seq = str(seq_record.seq)
    a, c, g, t, gc = liczenie_nukleotydow(seq)

    procent_a = (a / len(seq_records)) * 100
    procent_c = (c / len(seq_records)) * 100
    procent_g = (g / len(seq_records)) * 100
    procent_t = (t / len(seq_records)) * 100
    procent_par_gc = (gc / len(seq_records)) * 100
                      
    lokalizacja_sekwencji_gc = seq.find("GC")   

# Dodanie danych do listy
    dane.append([seq_record.id, len(seq_records), a, c, g, t, gc, procent_a, procent_c, procent_g, procent_t, procent_par_gc])

# Tworzenie DataFrame
df = pd.DataFrame(dane, columns=["ID", "Lenght", "A", "C", "G", "T", "GC", "A%", "C%", "G%", "T%", "GC%"])

# Zapis do pliku CSV
df.to_csv("wyniki.csv", index=False)
print("Wyniki zapisane do 'wyniki.csv'")


# Wizualizacja wyników
import matplotlib.pyplot as plt

# Wczytujemy dane z pliku CSV
df = pd.read_csv("wyniki.csv")

# Wykres słupkowy dla procentowej zawartości A, C, G, T
plt.figure(figsize=(10, 6))
plt.bar(df["ID"], df["A%"], color="blue", label="A")
plt.bar(df["ID"], df["C%"], color="green", label="C")
plt.bar(df["ID"], df["G%"], color="red", label="G")
plt.bar(df["ID"], df["T%"], color="orange", label="T")

# Dostosowanie wykresu
plt.title("Procentowa zawartość nukleotydów A, C, G, T")
plt.xlabel("ID sekwencji")
plt.ylabel("%")
plt.xticks(rotation=90)
plt.legend(title="Nukleotydy")
plt.tight_layout()

plt.show()


# Wykres słupkowy dla procentowej zawartości par GC
plt.figure(figsize=(10, 6))
plt.bar(df["ID"], df["GC%"], color="purple", label="GC")

# Dostosowanie wykresu
plt.title("Procentowa zawartość par GC")
plt.xlabel("Zawartosć GC")
plt.ylabel("%")
plt.xticks(rotation=90)
plt.legend(title="Nukleotydy")
plt.tight_layout()

plt.show()
