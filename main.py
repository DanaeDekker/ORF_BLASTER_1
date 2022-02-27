from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
import re
import time


def read_file(file_name):
    """ Read the fasta-file.

    :param file_name: string - name of the fasta-file
    :return headers: list - contains the headers of the sequences
    :return seqs: list - contains the sequences
    """
    # Define the variables.
    headers = []
    seqs = []
    header = ""
    line_x = ""

    # Get the data from the file.
    with open(file_name) as file_open:
        for line in file_open:
            line = line.rstrip()
            if ">" in line:
                if line_x:
                    seqs.append(line_x)
                    headers.append(header)
                    line_x = ""
                header = line
            else:
                line_x += line.upper()
        seqs.append(line_x)
        headers.append(header)

    # Return the headers and sequences.
    return headers, seqs


def blast(headers, seqs):
    """ Function to BLAST the sequences.

    :param headers: list - contains the headers
    :param seqs: list - contains the sequences
    """
    # Start the BLAST for each sequence in the list.
    for i in range(len(seqs)):
        # regex_match = re.search("ORF[0-9]{1,3}", headers[i])
        # orf = regex_match.group()
        print(headers[i].replace(">", ""))
        result_handle = NCBIWWW.qblast("blastp",
                                       "swissprot",
                                       seqs[i],
                                       expect=0.01,
                                       word_size=6,
                                       hitlist_size=10,
                                       entrez_query="txid4891[ORGN]")

        # Write the results in an xml-file.
        with open("{}.xml".format(headers[i].replace(">", "")), "w") \
                as out_handle:
            out_handle.write(result_handle.read())

        # Close the result handle.
        result_handle.close()
        print("Sleep...")
        time.sleep(300.0)


def parse_xml(header, seq):
    """ Parse the xml-file.

    :param header: string - header of the sequence
    :param seq: string - sequence
    """
    # Define the variables.
    values = []
    rows = []

    # Open the xml-file and read it using Biopython's parse-function.
    # regex_match = re.search("ORF[0-9]{1,3}", header)
    # orf = regex_match.group(0)
    result_handle = open("{}.xml".format(header.replace(">", "")))
    blast_record = NCBIXML.read(result_handle)

    # Add the alignment-data to a list. Each list contains ten lists
    # with alignments, or less if less alignments were found.
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            values.append(alignment.title)
            values.append(alignment.accession)
            values.append(hsp.align_length)
            values.append(hsp.expect)
            values.append(hsp.score)
            values.append(hsp.identities)
            values.append(hsp.positives)
            values.append(hsp.gaps)
            values.append(hsp.sbjct)
            values.append(hsp.query)
            values.append(hsp.match)

        if values:
            # Find the scientific name of the species in the title.
            title = str(values[0]).split("[")
            species = title[len(title) - 1]
            scientific = species.replace("]", "")

            # Call the functions to get the neccesary data.
            taxid = get_taxid(scientific)
            protein = get_protein(values[0])

            # Append the results.
            row = str("augustus" + "\t" + str(header.replace(">", ""))
                      + "\t" + str(seq) + "\t" + str(taxid) +
                      "\t" + str(scientific) + "\t" + str(values[1]) +
                      "\t" + str(protein) + "\t" + str(values[2]) +
                      "\t" + str(values[3]) + "\t" + str(values[4]) +
                      "\t" + str(values[5]) + "\t" + str(values[6]) +
                      "\t" + str(values[7]) + "\t" + str(values[8]) +
                      "\t" + str(values[9]) + "\t" + str(values[10]))

            # Append the rows.
            rows.append(row)

            # Empty the list with values.
            values = []

    # Return the lists.
    return rows


def get_taxid(scientific):
    """ Get the taxid and lineage of an organism.

    :param scientific: string - scientific name of the organism
    :return taxid_lineage: string - contains the taxid, domain, genus
    and species of the organism
    """
    # Extract the taxonomy from the Entrez-database.
    handle = Entrez.esearch(db="Taxonomy", term=scientific)
    record = Entrez.read(handle)
    taxid = record["IdList"][0]

    # Return the string.
    return taxid


def get_protein(title):
    """ Get the accession and the name of a protein.

    :param title: string - contains the name of the protein and the
    scientific name of the organism
    :return protein: string - contains the name of the protein
    """
    # Find the name of the protein.
    full = re.search(r"Full=(.*)", title)
    match = full.group(0)
    if ";" in match:
        protein = match.split(";")[0]
    else:
        protein = match.split("[")[0]
    protein = protein.replace("Full=", "")

    # Return the string.
    return protein


if __name__ == '__main__':
    # Define the variables.
    file_name = "augustus_output.txt"

    # Read the txt-file with the headers, sequences and quality-codes.
    headers, seqs = read_file(file_name)

    # Use BLAST on the sequences.
    # blast(headers, seqs)

    # Write a header for the tsv-file.
    with open("results.tsv", "w") as results_open:
        results_open.write("tool\tgene/orf\tseq\ttaxid\torganism"
                           "\taccession\tprotein\tlength\texpect\tscore"
                           "\tidentities\tpositives\tgaps\tsubject"
                           "\tquery\tmatch\n")

    # Parse the xml-files and append the results to a tsv-file.
    with open("results.tsv", "a") as results_open:
        for n in range(len(headers)):
            results = parse_xml(headers[n], seqs[n])
            for result in results:
                if result:
                    results_open.write(result + "\n")
