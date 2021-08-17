"""
step00: for general purpose within all step
"""
import hashlib
import hmac
import os
import pickle
import tarfile
import tempfile
import typing
import numpy

key = bytes("asdf", "UTF-8")
matplotlib_parameters = {"font.size": 50, "axes.labelsize": 50, "axes.titlesize": 75, "xtick.labelsize": 50, "ytick.labelsize": 50, "font.family": "serif", "legend.fontsize": 30, "legend.title_fontsize": 30, "figure.dpi": 300}
derivations = ("Accuracy", "Balanced_Accuracy", "Sensitivity", "Specificity", "Precision")
numeric_columns = {"Gestational Week", "Weight", "Mother Age", "Hospitalized Day", "Apgar Score", "Weight gain"}
markers = ["o", "v", "^", "<", ">", "8", "s", "p", "*", "h", "H", "D", "d"]
detailed_PTB = ("Extremely PTB", "Very PTB", "Late PTB", "Normal")


def file_list(path: str) -> typing.List[str]:
    """
    file_list: return a list of files in path
    """
    return list(filter(lambda x: os.path.isfile(x), list(map(lambda x: os.path.join(path, x), os.listdir(path)))))


def directory_list(path: str) -> typing.List[str]:
    """
    directory_list: return a list of directories in path
    """
    return list(filter(lambda x: os.path.isdir(x), list(map(lambda x: os.path.join(path, x), os.listdir(path)))))


def _make_hmac(message: bytes) -> bytes:
    """
    make_hmac: return a HMAC
    """
    return hmac.new(key, message, hashlib.sha512).digest()


def make_pickle(path: str, data: typing.Any) -> None:
    """
    make_pickle: create a pickle
    """
    if not path.endswith(".tar.gz"):
        raise ValueError("Path must end with .tar.gz")

    pkl = pickle.dumps(data, protocol=pickle.HIGHEST_PROTOCOL)
    key = _make_hmac(pkl)

    with tempfile.TemporaryDirectory() as tmp_dir:
        with open(os.path.join(tmp_dir, "data.pkl"), "wb") as f:
            f.write(pkl)
        with open(os.path.join(tmp_dir, "key.txt"), "wb") as f:
            f.write(key)

        with tarfile.open(path, "w:gz") as tar:
            tar.add(os.path.join(tmp_dir, "data.pkl"), arcname="data.pkl")
            tar.add(os.path.join(tmp_dir, "key.txt"), arcname="key.txt")


def read_pickle(path: str) -> typing.Any:
    """
    read_pickle: read a pickle file
    """
    if not path.endswith(".tar.gz"):
        raise ValueError("Path must end with .tar.gz")
    if not tarfile.is_tarfile(path):
        raise ValueError("Path cannot be read as a tar file")

    with tempfile.TemporaryDirectory() as tmp_dir:
        with tarfile.open(path, "r:gz") as tar:
            tar.extractall(tmp_dir)

        with open(os.path.join(tmp_dir, "data.pkl"), "rb") as f:
            pkl = f.read()
        with open(os.path.join(tmp_dir, "key.txt"), "rb") as f:
            key = f.read()

    if not hmac.compare_digest(_make_hmac(pkl), key):
        raise ValueError("Data is not valid")

    return pickle.loads(pkl)


def consistency_taxonomy(taxonomy: str) -> str:
    """
    consistency_taxonomy: make taxonomy information with consistency
    """
    if taxonomy == "Unassigned":
        return taxonomy
    else:
        return ";".join(list(filter(lambda x: len(x) > 3, list(map(lambda x: x.strip().replace("[", "").replace("]", ""), taxonomy.split(";")))))[-3:])


def simplified_taxonomy(taxonomy: str) -> str:
    """
    simplified_taxonomy: simplified taxonomy information for file name
    """
    if taxonomy == "Unassigned":
        return taxonomy
    else:
        return " ".join(list(filter(None, list(map(lambda x: x.strip().replace("[", "").replace("]", "")[3:], taxonomy.split(";")))))[-3:])


def aggregate_confusion_matrix(confusion_matrix: numpy.ndarray, derivation: str = "") -> float:
    """
    aggregate_confusion_matrix: derivations from confusion matrix
    """

    assert (derivation in derivations)
    assert confusion_matrix.shape == (2, 2)

    TP, FP, FN, TN = confusion_matrix[0][0], confusion_matrix[0][1], confusion_matrix[1][0], confusion_matrix[1][1]
    # assert TP and FP and FN and TN

    if derivation == "Sensitivity":
        return TP / (TP + FN)
    elif derivation == "Specificity":
        return TN / (TN + FP)
    elif derivation == "Precision":
        return TP / (TP + FP)
    elif derivation == "Accuracy":
        return (TP + TN) / (TP + TN + FP + FN)
    elif derivation == "Balanced_Accuracy":
        return TP / 2 * (TP + FN) + TN / 2 * (TN + FP)
    else:
        raise Exception("Something went wrong!!")


def star(p: float) -> str:
    if 0 <= p <= 1.00e-04:
        return "****"
    elif 1.00e-04 < p <= 1.00e-03:
        return "***"
    elif 1.00e-03 < p <= 1.00e-02:
        return "**"
    elif 1.00e-02 < p <= 5.00e-02:
        return "*"
    elif 5.00e-02 < p <= 1.00e+00:
        return "ns"
    else:
        raise ValueError("Something went wrong: {0} !!".format(p))
