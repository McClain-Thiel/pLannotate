import os

# Avoid network/database downloads during tests
os.environ.setdefault("PLANNOTATE_SKIP_DB_DOWNLOAD", "1")
