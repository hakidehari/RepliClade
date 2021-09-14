import os


class ReplicladeSettings:

    ALIGNMENTS_PATH_MAIN = (
        os.getcwd()
        + os.path.sep
        + "filesystem"
        + os.path.sep
        + "alignment"
        + os.path.sep
    )
    DNA_PATH = (
        os.getcwd() + os.path.sep + "filesystem" + os.path.sep + "DNA" + os.path.sep
    )
    RESULTS_PATH = os.getcwd() + os.path.sep + "filesystem" + os.path.sep + "results" + os.path.sep

    EXECUTABLES_PATH = ALIGNMENTS_PATH_MAIN + "executables" + os.path.sep
    ALIGNMENTS_PATH = ALIGNMENTS_PATH_MAIN + "align" + os.path.sep
    BLAST_PATH = ALIGNMENTS_PATH_MAIN + "blast" + os.path.sep
    FASTAS_PATH = ALIGNMENTS_PATH_MAIN + "fastas" + os.path.sep
