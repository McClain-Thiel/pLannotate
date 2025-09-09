import os
import subprocess
import sys
import ssl
from datetime import date
from importlib.resources import files
from pathlib import Path
from tempfile import NamedTemporaryFile
import hashlib
import tarfile
import urllib.request
import urllib.error
import shutil

import pandas as pd
import yaml
from platformdirs import user_cache_dir
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from plannotate import __version__ as plannotate_version

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

valid_genbank_exts = [".gbk", ".gb", ".gbf", ".gbff"]
valid_fasta_exts = [".fa", ".fasta", ".fas", ".fna"]
MAX_PLAS_SIZE = 50000

# Environment variable names for configuring database caching
CACHE_ENV = "PLANNOTATE_DB_DIR"
AUTO_ENV = "PLANNOTATE_AUTO_DOWNLOAD"
SKIP_ENV = "PLANNOTATE_SKIP_DB_DOWNLOAD"

DB_ARCHIVE = "BLAST_dbs.tar.gz"
DB_CHECKSUM = (
    "34c7bacb1c73bd75129e16990653f73b3eba7e3cdb3816a55d3989a7601f2137"
)
DB_URL_TEMPLATE = (
    "https://github.com/mmcguffi/pLannotate/releases/download/"
    "v1.2.0/BLAST_dbs.tar.gz"
)


def get_cache_root():
    """Return the root directory for cached databases."""
    return Path(os.environ.get(CACHE_ENV, user_cache_dir("pLannotate")))


def set_db_cache_dir(path):
    """Set a custom cache directory for databases."""
    os.environ[CACHE_ENV] = str(path)


def _confirm_download():
    """Determine whether the user has agreed to a database download."""
    if os.environ.get(SKIP_ENV):
        return False
    if os.environ.get(AUTO_ENV):
        return True
    if not sys.stdin.isatty():
        return False
    reply = input(
        "pLannotate requires additional databases (~900MB). Download now? [y/N]: "
    )
    return reply.strip().lower() in {"y", "yes"}


def download_db(cache_root=None, url=None, checksum=DB_CHECKSUM, force=False):
    """Download and extract the pLannotate databases.

    Parameters
    ----------
    cache_root: str or Path, optional
        Directory in which to place the downloaded data. Defaults to the
        user cache directory.
    url: str, optional
        URL from which to download the databases archive. Defaults to the
        official release matching the installed pLannotate version.
    checksum: str, optional
        Expected SHA256 checksum of the archive. The download will be removed
        and an exception raised if the checksum does not match.
    force: bool, optional
        If True, re-download even if the databases already exist.
    """

    cache_root = Path(cache_root) if cache_root else get_cache_root()
    db_dir = cache_root / "BLAST_dbs"
    if db_dir.exists() and not force:
        return str(db_dir)

    cache_root.mkdir(parents=True, exist_ok=True)
    if not _confirm_download():
        raise RuntimeError("Database download declined or skipped.")

    url = url or DB_URL_TEMPLATE
    archive_path = cache_root / DB_ARCHIVE

    try:
        print(f"Downloading databases from {url}")
        
        # Get file size for progress bar
        with urllib.request.urlopen(url) as response:
            total_size = int(response.headers.get('Content-Length', 0))
        
        # Download with progress bar
        with urllib.request.urlopen(url) as response, open(archive_path, "wb") as fh:
            if total_size > 0:
                with tqdm(total=total_size, unit='B', unit_scale=True, desc="Downloading") as pbar:
                    while True:
                        chunk = response.read(8192)
                        if not chunk:
                            break
                        fh.write(chunk)
                        pbar.update(len(chunk))
            else:
                # Fallback without progress bar if size unknown
                shutil.copyfileobj(response, fh)
                
    except urllib.error.URLError as e:
        if "CERTIFICATE_VERIFY_FAILED" in str(e):
            python_version = f"{sys.version_info.major}.{sys.version_info.minor}"
            error_msg = (
                f"SSL Certificate verification failed. This is a system configuration issue.\n"
                f"To fix this, run the following command in your terminal:\n\n"
                f"    open '/Applications/Python {python_version}/Install Certificates.command'\n\n"
                f"Or try running:\n"
                f"    /Applications/Python\\ {python_version}/Install\\ Certificates.command\n\n"
                f"This will install the necessary SSL certificates for Python to download files.\n"
                f"Original error: {e}"
            )
            raise RuntimeError(error_msg) from e
        else:
            raise RuntimeError(f"Failed to download databases: {e}") from e
    except Exception as e:
        raise RuntimeError(f"Unexpected error during download: {e}") from e

    if checksum:
        print("Verifying download integrity...")
        digest = hashlib.sha256()
        file_size = archive_path.stat().st_size
        with open(archive_path, "rb") as fh:
            with tqdm(total=file_size, unit='B', unit_scale=True, desc="Verifying") as pbar:
                for chunk in iter(lambda: fh.read(8192), b""):
                    digest.update(chunk)
                    pbar.update(len(chunk))
        if digest.hexdigest() != checksum:
            archive_path.unlink(missing_ok=True)
            raise ValueError("Checksum mismatch for downloaded databases.")
        print("✓ Download integrity verified")

    print("Extracting databases...")
    with tarfile.open(archive_path, "r:gz") as tar:
        members = tar.getmembers()
        with tqdm(total=len(members), desc="Extracting") as pbar:
            for member in members:
                tar.extract(member, cache_root)
                pbar.update(1)
    archive_path.unlink(missing_ok=True)
    print(f"✓ Databases successfully installed to: {db_dir}")

    return str(db_dir)


def get_db_dir(download=True):
    """Return the directory containing the databases, downloading if needed."""
    db_dir = get_cache_root() / "BLAST_dbs"
    if download and not os.environ.get(SKIP_ENV) and not db_dir.exists():
        download_db()
    return db_dir

DF_COLS = [
    "sseqid",
    "qstart",
    "qend",
    "sstart",
    "send",
    "sframe",
    "score",
    "evalue",
    "qseq",
    "length",
    "slen",
    "pident",
    "qlen",
    "db",
    "Feature",
    "Description",
    "Type",
    "priority",
    "percmatch",
    "abs percmatch",
    "pi_permatch",
    "wiggle",
    "wstart",
    "wend",
    "kind",
    "qstart_dup",
    "qend_dup",
    "fragment",
]


def get_resource(group, name):
    return str(files(__package__) / f"data/{group}/{name}")


def get_image(name):
    return get_resource("images", name)


def get_template(name):
    return get_resource("templates", name)


def get_example_fastas():
    return get_resource("fastas", "")


def get_yaml_path():
    return get_resource("data", "databases.yml")


# def get_yaml(blast_database_loc):
#     return parse_yaml(,blast_database_loc)


def get_details(name):
    return get_resource("data", name)


def get_name_ext(file_loc):
    base = os.path.basename(file_loc)
    name = os.path.splitext(base)[0]
    ext = os.path.splitext(base)[1]
    return name, ext


def validate_file(file, ext, max_length=MAX_PLAS_SIZE):
    if ext in valid_fasta_exts:
        # This catches errors on file uploads via Biopython
        temp_fileloc = NamedTemporaryFile()
        record = list(SeqIO.parse(file, "fasta"))
        try:
            record[0].annotations["molecule_type"] = "DNA"
        except IndexError:
            error = (
                "Malformed fasta file --> please submit a fasta file in standard format"
            )
            raise ValueError(error)
        SeqIO.write(record, temp_fileloc.name, "fasta")
        record = list(SeqIO.parse(temp_fileloc.name, "fasta"))
        temp_fileloc.close()

        if len(record) != 1:
            error = "FASTA file contains many entries --> please submit a single FASTA file."
            raise ValueError(error)

    elif ext in valid_genbank_exts:
        temp_fileloc = NamedTemporaryFile()
        try:
            record = list(SeqIO.parse(file, "gb"))[0]
        except IndexError:
            error = "Malformed Genbank file --> please submit a Genbank file in standard format"
            raise ValueError(error)
        # submitted_gbk = record # for combining -- not current imlementated
        SeqIO.write(record, temp_fileloc.name, "fasta")
        record = list(SeqIO.parse(temp_fileloc.name, "fasta"))
        temp_fileloc.close()

    else:
        error = "must be a FASTA or GenBank file"
        raise ValueError(error)

    if len(record) != 1:
        error = (
            "FASTA file contains many entries --> please submit a single FASTA file."
        )
        raise ValueError(error)

    inSeq = str(record[0].seq)

    validate_sequence(inSeq, max_length)

    return inSeq


def validate_sequence(inSeq, max_length=MAX_PLAS_SIZE):
    IUPAC = "GATCRYWSMKHBVDNgatcrywsmkhbvdn"
    if not set(inSeq).issubset(IUPAC):
        error = "Sequence contains invalid characters -- must be ATCG and/or valid IUPAC nucleotide ambiguity code"
        raise ValueError(error)

    if len(inSeq) > max_length:
        error = f"Are you sure this is an engineered plasmid? Entry size is too large -- must be {max_length} bases or less."
        raise ValueError(error)


def get_gbk(inDf, inSeq, is_linear=False, record=None):
    record = get_seq_record(inDf, inSeq, is_linear, record)

    # converts gbk into straight text
    outfileloc = NamedTemporaryFile()
    with open(outfileloc.name, "w") as handle:
        record.annotations["molecule_type"] = "DNA"
        SeqIO.write(record, handle, "genbank")
    with open(outfileloc.name) as handle:
        record = handle.read()
    outfileloc.close()

    return record


def get_seq_record(inDf, inSeq, is_linear=False, record=None):
    # this could be passed a more annotated df
    inDf = inDf.reset_index(drop=True)

    if inDf.empty:
        inDf = pd.DataFrame(columns=DF_COLS)

    def FeatureLocation_smart(r):
        # creates compound locations if needed
        if r.qend > r.qstart:
            return FeatureLocation(r.qstart, r.qend, r.sframe)
        elif r.qstart > r.qend:
            first = FeatureLocation(r.qstart, r.qlen, r.sframe)
            second = FeatureLocation(0, r.qend, r.sframe)
            if r.sframe == 1 or r.sframe == 0:
                return first + second
            elif r.sframe == -1:
                return second + first

    # adds a FeatureLocation object so it can be used in gbk construction
    inDf["feat loc"] = inDf.apply(FeatureLocation_smart, axis=1)

    # make a record if one is not provided
    if record is None:
        record = SeqRecord(seq=Seq(inSeq), name="plasmid")

    record.annotations["data_file_division"] = "SYN"

    if "comment" not in record.annotations:
        record.annotations["comment"] = (
            f"Annotated with pLannotate v{plannotate_version}"
        )
    else:
        record.annotations["comment"] = (
            f"Annotated with pLannotate v{plannotate_version}. {record.annotations['comment']}"
        )

    if "date" not in record.annotations:
        record.annotations["date"] = date.today().strftime("%d-%b-%Y").upper()

    if "accession" not in record.annotations:
        record.annotations["accession"] = "."

    if "version" not in record.annotations:
        record.annotations["version"] = "."

    if is_linear:
        record.annotations["topology"] = "linear"
    else:
        record.annotations["topology"] = "circular"

    # this adds "(fragment)" to the end of a feature name
    # if it is a fragment. Maybe a better way show this data in the gbk
    # for downstream analysis, though this may suffice. change type to
    # non-canonical `fragment`?
    def append_frag(row):
        if row["fragment"] is True:
            return f"{row['Feature']} (fragment)"
        else:
            return f"{row['Feature']}"

    inDf["Feature"] = inDf.apply(lambda x: append_frag(x), axis=1)

    inDf["Type"] = inDf["Type"].str.replace("origin of replication", "rep_origin")
    for index in inDf.index:
        record.features.append(
            SeqFeature(
                inDf.loc[index]["feat loc"],
                type=inDf.loc[index]["Type"],  # maybe change 'Type'
                qualifiers={
                    "note": ["pLannotate"],
                    "label": [inDf.loc[index]["Feature"]],
                    "database": [inDf.loc[index]["db"]],
                    "identity": [str(round(inDf.loc[index]["pident"], 1))],
                    "match_length": [str(round(inDf.loc[index]["percmatch"], 1))],
                    "fragment": [str(inDf.loc[index]["fragment"])],
                    "other": [inDf.loc[index]["Type"]],
                },
            )
        )  # maybe change 'Type'

    return record


def get_clean_csv_df(recordDf):
    # change sseqid to something more legible
    columns = [
        "sseqid",
        "qstart",
        "qend",
        "sframe",
        "pident",
        "slen",
        "length",
        "abs percmatch",
        "fragment",
        "db",
        "Feature",
        "Type",
        "Description",
        "qseq",
    ]
    cleaned = recordDf[columns]
    replacements = {
        "qstart": "start location",
        "qend": "end location",
        "sframe": "strand",
        "pident": "percent identity",
        "slen": "full length of feature in db",
        "qseq": "sequence",
        "length": "length of found feature",
        "abs percmatch": "percent match length",
        "db": "database",
    }
    cleaned = cleaned.rename(columns=replacements)
    return cleaned


# parse yaml file
# def parse_yaml(file_name):
#     with open(file_name, 'r') as f:
#         dbs = yaml.load(f, Loader = yaml.SafeLoader)

#     for db in dbs.keys():
#         method = dbs[db]['method']
#         try:
#             parameters = " ".join(dbs[db]['parameters'])
#         except KeyError:
#             parameters = ""
#         details = dbs[db]['details']
#         #print(f'{method} {parameters} {details}')
#         return method, parameters, details

#         print()


def get_yaml(yaml_file_loc, db_dir=None):
    """Load database configuration from YAML and resolve data locations."""

    with open(yaml_file_loc, "r") as f:
        dbs = yaml.load(f, Loader=yaml.SafeLoader)

    cache_dir = Path(db_dir) if db_dir else get_db_dir()

    for db in dbs.keys():
        try:
            dbs[db]["parameters"] = " ".join(dbs[db]["parameters"])
        except KeyError:
            dbs[db]["parameters"] = ""

        location = dbs[db].get("location", "Default")
        if location == "Default":
            dbs[db]["db_loc"] = str(cache_dir / db)
        else:
            dbs[db]["db_loc"] = os.path.join(location, db)

    return dbs



