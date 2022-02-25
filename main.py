from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
import re


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
        regex_match = re.search("ORF[1-9]{0,3}", headers[i])
        orf = regex_match.group(0)
        print(orf)
        result_handle = NCBIWWW.qblast("blastx",
                                       "swissprot",
                                       seqs[i],
                                       expect=0.01,
                                       word_size=6,
                                       hitlist_size=10,
                                       entrez_query=27300,
                                       genetic_code=12)

        # Write the results in an xml-file.
        with open("%s.xml" % orf, "w") as out_handle:
            out_handle.write(result_handle.read())

        # Close the result handle.
        result_handle.close()


def parse_xml(header, seq):
    """ Parse the xml-file.

    :return results: list - contains the results per alignment
    :return organisms: list - contains the organism per alignment
    :return proteins: list - contains the protein per alignment
    """
    # Define the variables.
    values = []
    rows = []

    # Open the xml-file and read it using Biopython's parse-function.
    regex_match = re.search("ORF[1-9]{0,3}", header)
    orf = regex_match.group(0)
    result_handle = open("{}.xml".format(orf))
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
            # Find the name of the species in the alignment-title.
            species = re.search(r"(?<=\[).+?(?=\])", values[0])
            scientific = str(species.group())

            # Call the functions to get the neccesary data.
            taxid = get_organism(scientific)
            protein = get_protein(values[0], scientific)

            # Append the results.
            row = str(str(orf) + "\t" + str(seq) + "\t" + str(taxid) +
                      "\t" + str(scientific) + "\t" + str(values[1]) +
                      "\t" + str(protein) + "\t" + str(values[3]) +
                      "\t" + str(values[2]) + "\t" + str(values[4]) +
                      "\t" + str(values[5]) + "\t" + str(values[6]) +
                      "\t" + str(values[7]) + "\t" + str(values[8]) +
                      "\t" + str(values[9]) + "\t" + str(values[10]))

            # Append the rows.
            rows.append(row)

            # Empty the list with values.
            values = []

    # Return the lists.
    return rows


def get_organism(scientific):
    """ Get the taxid and lineage of an organism.

    :param scientific: string - scientific name of the organism
    :return taxid_lineage: string - contains the taxid, domain, genus
    and species of the organism
    """
    # Extract the taxonomy from the Entrez-database.
    handle = Entrez.esearch(db="Taxonomy", term=scientific)
    record = Entrez.read(handle)
    taxid = record["IdList"][0]
    handle = Entrez.efetch(db="Taxonomy", id=taxid, retmode="xml")
    records = Entrez.read(handle)
    taxonomy = records[0]["Lineage"]
    lineage = taxonomy.split("; ")

    # Split the scientific name into a list of genus and species.
    scientific_list = scientific.split(" ")

    # Sometimes, only the family or genus is in the name. We want a
    # scientific name with
    while len(scientific_list) < 2:
        scientific_list.append("unknown")

    # Create a string that contains the taxid, domain, genus en species,
    # separated by tabs.
    taxid_lineage = str(taxid + "\t" + lineage[1] + "\t" +
                        scientific_list[0] + "\t" + scientific_list[1])

    # Return the string.
    return taxid, taxid_lineage


def get_protein(title, scientific):
    """ Get the accession and the name of a protein.

    :param title: string - contains the name of the protein and the
    scientific name of the organism
    :param scientific: string - scientific name of the organism
    :return accession_protein: string - contains the accession code
    and the name of the protein
    """
    # Remove the scientific name from the title.
    scientific_removed = title.replace("[" + scientific + "]", "")

    # If there are multiple species in the title, a > is used, so split
    # the string there so only one species is taken.
    multispecies = scientific_removed.split(" >")

    # In front of the protein name is the accession code. As we want to
    # get them seperately, split on the |.
    protein = multispecies[0].split("| ")

    # Return the string.
    return protein


if __name__ == '__main__':
    # Define the variables.
    file_name = "orfs_nuc.fa"

    # Read the txt-file with the headers, sequences and quality-codes.
    headers, seqs = read_file(file_name)

    # Use BLAST on the sequences.
    blast(headers, seqs)

    # Parse the xml-files. If results are found, append the lists to
    # another list to combine them.
    with open("results.tsv", "w") as results_open:
        results_open.write("orf\tseq\ttaxid\taccession\tprotein\tlength"
                           "\texpect\tscore\tidentities\tpositives\t"
                           "gaps\tsubject\tquery\tmatch\n")
        for i in range(len(headers)):
            results = parse_xml(headers[i], seqs[i])
            for result in results:
                if result:
                    results_open.write(result + "\n")