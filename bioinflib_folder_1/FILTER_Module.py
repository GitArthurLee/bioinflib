def find_gc_content(seq: str) -> float:
    content = 0
    for nucl in seq:
        if nucl == 'G' or nucl == 'g' or nucl == 'C' or nucl == 'c':
            content += 1
    return (content / len(seq)) * 100


def find_average_quality(quality_string: str) -> float:
    total_quality = sum(ord(char) - 33 for char in quality_string)
    return total_quality / len(quality_string)
