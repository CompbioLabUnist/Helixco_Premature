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
import matplotlib.patches
import matplotlib.transforms
import numpy

key = bytes("asdf", "UTF-8")
small = 10 ** 3
big = 10 ** 6

matplotlib_parameters = {"font.size": 50, "axes.labelsize": 50, "axes.titlesize": 75, "xtick.labelsize": 50, "ytick.labelsize": 50, "font.family": "sans-serif", "legend.fontsize": 30, "legend.title_fontsize": 30, "figure.dpi": 500, "pdf.fonttype": 42, "ps.fonttype": 42}
markers = ["o", "v", "^", "<", ">", "8", "s", "p", "*", "h", "H", "D", "d"]
derivations = ("Accuracy", "Balanced Accuracy", "Sensitivity", "Specificity", "Precision")
taxonomics_ranks = ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

numeric_columns = {"Gestational Week", "Weight", "Mother Age", "Hospitalized Day", "Apgar Score", "Weight gain", "Cholesterol", "TG", "HDL", "LDL", "Glucose", "WBC", "Hb", "Hct", "ESR", "hsCRP", "AST", "ALT", "SBP", "DBP"}
detailed_PTB = ("Early PTB", "Late PTB", "Normal")
PTB_colors = {"Early PTB": "tab:red", "Late PTB": "tab:pink", "Normal": "w"}

selected_sites = ("M", "C", "V")
selected_long_sites = ("Mouth", "Cervix", "Vagina")
selected_sites_dict = dict(zip(selected_sites, selected_long_sites))

pdist_list = ["braycurtis", "euclidean", "hamming", "jaccard"]


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


def remove_preceding_underscores(x: str) -> str:
    x = x.strip()
    if x.startswith("__"):
        return x[2:]
    elif x[1:3] == "__":
        return x[3:]
    else:
        return x


def consistency_taxonomy(taxonomy: str, number: int = 2) -> str:
    """
    consistency_taxonomy: make taxonomy information with consistency
    """
    if taxonomy == "Unassigned":
        return taxonomy
    else:
        return ";".join(list(filter(lambda x: len(x) > 3, list(map(lambda x: remove_preceding_underscores(x), taxonomy.split(";")))))[-number:]).replace("_", " ")


def aggregate_confusion_matrix(confusion_matrix: numpy.ndarray, derivation: str = "") -> float:
    """
    aggregate_confusion_matrix: derivations from confusion matrix
    """

    assert (derivation in derivations)
    assert confusion_matrix.shape == (2, 2)

    TP, FP, FN, TN = confusion_matrix[0][0], confusion_matrix[0][1], confusion_matrix[1][0], confusion_matrix[1][1]

    if derivation == "Sensitivity":
        return TP / (TP + FN)
    elif derivation == "Specificity":
        return TN / (TN + FP)
    elif derivation == "Precision":
        return TP / (TP + FP)
    elif derivation == "Accuracy":
        return (TP + TN) / (TP + TN + FP + FN)
    elif derivation == "Balanced Accuracy":
        return TP / (2 * (TP + FN)) + TN / (2 * (TN + FP))
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


def confidence_ellipse(x, y, ax, n_std=2.0, facecolor="none", **kwargs):
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = numpy.cov(x, y)
    pearson = cov[0, 1] / numpy.sqrt(cov[0, 0] * cov[1, 1])
    ell_radius_x = numpy.sqrt(1 + pearson)
    ell_radius_y = numpy.sqrt(1 - pearson)
    ellipse = matplotlib.patches.Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2, facecolor=facecolor, **kwargs)

    scale_x = numpy.sqrt(cov[0, 0]) * n_std
    mean_x = numpy.mean(x)

    scale_y = numpy.sqrt(cov[1, 1]) * n_std
    mean_y = numpy.mean(y)

    transf = matplotlib.transforms.Affine2D().rotate_deg(45).scale(scale_x, scale_y).translate(mean_x, mean_y)
    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


def filtering_taxonomy(taxonomy: str) -> bool:
    taxonomy_list = taxonomy.split("; ")
    if len(taxonomy_list) < 6:
        return False
    if taxonomy_list[5] in ("g__", "__"):
        return False
    else:
        return True


def select_taxonomy(taxonomy: str) -> str:
    taxonomy_list = taxonomy.split("; ")
    return remove_preceding_underscores(taxonomy_list[1]).replace("_", " ")


def simplified_taxonomy(taxonomy: str) -> str:
    taxonomy_list = taxonomy.split("; ")
    if (len(taxonomy_list) == 6):
        return remove_preceding_underscores(taxonomy_list[-1]).replace("_", " ") + " spp."
    elif (taxonomy_list[-1] == "__") or (taxonomy_list[-1] == "s__"):
        return remove_preceding_underscores(taxonomy_list[-2]).replace("_", " ") + " spp."
    else:
        return remove_preceding_underscores(taxonomy_list[-2])[0] + ". " + remove_preceding_underscores(taxonomy_list[-1]).replace("_", " ")


if __name__ == "__main__":
    taxo = "k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Prevotellaceae; g__Prevotella; __"
    print(consistency_taxonomy(taxo, 1))
