import re

def tryptic_pep(position, sequence):
    """Trypsin digestion, allowing 2 missed cleavages"""
    m = re.finditer(r"K|R", sequence)
    temp = [match.start() for match in m]
    temp.append(position - 1)
    temp = sorted(set(temp))
    idx = temp.index(position - 1)

    if idx - 2 < 0 and idx + 2 >= len(temp):
        pep_cut = sequence
        new_pos = position
    elif idx - 2 < 0:
        pep_cut = sequence[0 : temp[idx + 2] + 1]
        new_pos = position
    elif idx + 2 >= len(temp):
        pep_cut = sequence[temp[idx - 2] :]
        new_pos = position - temp[idx - 2]
    else:
        pep_cut = sequence[temp[idx - 2] : temp[idx + 2] + 1]
        new_pos = position - temp[idx - 2]

    return [new_pos, pep_cut]

def lysC_pep(position, sequence):
    """Lys-C digestion, allowing 2 missed cleavages"""
    m = re.finditer(r"K", sequence)
    temp = [match.start() for match in m]
    temp.append(position - 1)
    temp = sorted(set(temp))
    idx = temp.index(position - 1)

    if idx - 2 < 0 and idx + 2 >= len(temp):
        pep_cut = sequence
        new_pos = position
    elif idx - 2 < 0:
        pep_cut = sequence[0 : temp[idx + 2] + 1]
        new_pos = position
    elif idx + 2 >= len(temp):
        pep_cut = sequence[temp[idx - 2] :]
        new_pos = position - temp[idx - 2]
    else:
        pep_cut = sequence[temp[idx - 2] : temp[idx + 2] + 1]
        new_pos = position - temp[idx - 2]

    return [new_pos, pep_cut]
