"""Sequence Modules used in president."""
import tempfile


def preprocess(fasta):
    """Remove white spaces from IDs in a FASTA file.

       Removes trailing spaces and newlines.

       Converts sequences to upper case.

       Replaces - with . for pblat.

    Parameters
    ----------
    fasta : str
        FASTA file to remove white spaces.

    Returns
    -------
    processed FASTA file.

    """
    _, output = tempfile.mkstemp()
    with open(output, "w") as fout:
        with open(fasta, "r") as fin:
            for line in fin:
                # removes trailing spaces and linebreaks
                line = line.strip()
                # replaces whitespaces with a placeholder %space% to avoid errors
                # needs to be undone later in the code
                if len(line) > 0 and line.startswith('>'):
                    fout.write(line.replace(" ", "%space%") + "\n")
                else:
                    # convert sequence to upper case
                    line = line.upper()
                    # replaces gaps
                    line = line.replace('-', '.')
                    fout.write(line + "\n")
    return output
