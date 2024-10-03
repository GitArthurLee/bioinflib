#! /usr/bin/env python

def gc_content(seq):
    content = 0
    for nucl in seq:
        if nucl == 'G' or nucl == 'g' or nucl == 'C' or nucl == 'c':
            content += 1
    return (content / len(seq)) * 100


def average_quality(quality_string):
    total_quality = 0
    for symbol in quality_string:
        quality = ord(symbol) - 33
        total_quality += quality
    return total_quality / len(quality_string)
