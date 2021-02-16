"""Sequence Modules used in president."""
import tempfile
import os


def to_valid_upper(dna_seq):
    """
    Transform input sequence to valid output sequence by removing all non-ACGT characters with Ns.

    Parameters
    ----------
    dna_seq, str
              DNA sequence

    Returns
    -------
    str, fixed DNA sequence
    """
    valid_nt = set(list("acgtACGT"))
    return "".join([i.upper() if i in valid_nt else "N" for i in dna_seq])


def preprocess(multifasta, suffix):
    """Remove white spaces from IDs in a FASTA file.

       Removes trailing spaces and newlines.

       Converts sequences to upper case.

       Concatenates n fasta files received as *.fasta

       Replaces - with . for pblat.

    Parameters
    ----------
    fasta : str/[str]
        FASTA file(s) to preprocess.
    suffix : str
        Suffix for the temporary output file.

    Returns
    -------
    processed FASTA file.

    """
    # Create temp file to save output
    _, output = tempfile.mkstemp(suffix=suffix)

    # If input is a string then convert it to a list
    if isinstance(multifasta, str):
        multifasta = [multifasta]

    with open(output, "w") as fout:
        for fasta in multifasta:
            print(f"Preprocessing file: {fasta} -> {output}")
            with open(fasta, "r") as fin:
                for line in fin:
                    # removes trailing spaces and linebreaks
                    line = line.strip()
                    # replaces whitespaces with a placeholder %space% to avoid errors
                    # needs to be undone later in the code
                    if len(line) > 0 and line.startswith('>'):
                        # Append filename to sequence ID
                        header = str(line + ":" + os.path.basename(fasta)).replace(" ", "%space%")
                        fout.write(header + "\n")
                    else:
                        # convert sequence to upper case
                        line = line.upper()
                        # replaces gaps
                        line = line.replace('-', '.')
                        fout.write(line + "\n")
    return output
