def parse_gff(gfffile_name, txtfile_name):
    """ Parses the headers and sequences from a gff-file.

    :param gfffile_name: string - name of the gff-file
    :param txtfile_name: string - name of the txt-file with the parsed
    headers and sequences, almost in a fasta-format
    """
    # Define the variables.
    read = False

    # Parse the headers and sequences from the gff-file.
    with open(txtfile_name, "w") as txtfile_open:
        with open(gfffile_name, "r") as gfffile_open:
            for line in gfffile_open:
                # This line indicates the header.
                if line.startswith("# start"):
                    gene = ">" + line.strip().split(" ")[3]
                    seq = ""
                # This line indicates the start of the sequence.
                elif line.startswith("# protein sequence = "):
                    if "]" in line.strip():
                        seq += line.strip().split("[")[1].replace("]",
                                                                  "")
                    else:
                        read = True
                        seq += line.strip().split("[")[1]
                # This line indicates indicates the continuation of
                # the sequence.
                elif read:
                    if not line.strip().endswith("]"):
                        seq += line.strip().split(" ")[1]
                    else:
                        seq += line.strip().split(" ")[1].replace("]",
                                                                  "")
                        read = False
                # This line indicates the end of the sequence.
                elif line.startswith("# end"):
                    txtfile_open.write(gene + "\n")
                    txtfile_open.write(seq + "\n")


if __name__ == '__main__':
    # Define the variables.
    input_file = "augustus_output.gff"
    output_file = "augustus_output.txt"

    # Call the functions.
    parse_gff(input_file, output_file)
