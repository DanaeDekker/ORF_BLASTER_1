import re

def loadfile(bestand):
    """leest het bestand1 in en maakt hier een lijst van
    :param bestand-str
    :return: filelijst-list with header and sequence
    """
    orflist = []
    headerlist = []
    header = ""
    line_x = ""
    with open(bestand) as file:
        for line in file:
            line = line.rstrip()
            if ">" in line:
                if line_x:
                    orflist.append(line_x)
                    headerlist.append(header)
                    line_x = ""
                header = line
            else:
                line_x += line.upper()
        orflist.append(line_x)
        headerlist.append(header)
    return orflist, headerlist


def getorfnumber(header):
    regex_match= re.search("ORF[1-9]{0,3}",header)
    print(regex_match.group(0))
    



if __name__ == '__main__':
    bestand = "orfs_nuc.fasta"
    orflist, headerlist = loadfile(bestand)
    getorfnumber(headerlist[2])
    