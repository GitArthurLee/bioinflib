#! /usr/bin/env python

Transcribe = {
    "a": "a", "A": "A",
    "t": "u", "T": "U",
    "g": "g", "G": "G",
    "c": "c", "C": "C",
    "u": "t", "U": "T"
}

ComplementRNA = {
    "a": "t", "A": "T",
    "u": "a", "U": "A",
    "g": "c", "G": "C",
    "c": "g", "C": "G",
    "n": "n", "N": "N"
}

ComplementDNA = {
    "a": "t", "A": "T",
    "t": "a", "T": "A",
    "g": "c", "G": "C",
    "c": "g", "C": "G",
    "n": "n", "N": "N"
}


def reverse(seq):
    return seq[::-1]


def transcribe(seq):
    result = ''.join([Transcribe[i] for i in seq])
    return result


def complement_RNA(seq):
    result = ''.join([ComplementRNA[i] for i in seq])
    return result


def complement_DNA(seq):
    result = ''.join([ComplementDNA[i] for i in seq])
    return result


def reverse_complement_dna(seq):
    return reverse(complement_DNA(seq))


def reverse_complement_rna(seq):
    return reverse(complement_RNA(seq))


def is_dna(seq):
    unique_symbol = set(seq)
    nucleotides = set('ATGCatgc')
    return unique_symbol <= nucleotides


def is_rna(seq):
    unique_symbol = set(seq)
    nucleotides = set('AUGCaugc')
    return unique_symbol <= nucleotides


def ssDNA_MW(seq):
    return len(seq)*330


def ssRNA_MW(seq):
    return len(seq)*340


def GC_content(seq):
    content = 0
    for nucl in seq:
        if nucl == 'G' or nucl == 'g' or nucl == 'C' or nucl == 'c':
            content += 1
    return round(content/len(seq)*100)


def Tm_primer(seq):
    content_AT = 0
    content_GC = 0
    for nucl in seq:
        if nucl == 'A' or nucl == 'a' or nucl == 'T' or nucl == 't':
            content_AT += 1
    for nucl in seq:
        if nucl == 'G' or nucl == 'g' or nucl == 'C' or nucl == 'c':
            content_GC += 1
    return content_AT*2 + content_GC*4
