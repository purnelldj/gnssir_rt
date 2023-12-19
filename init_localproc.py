import shutil
from pathlib import Path

# make a directory called 'localproc'
Path("localproc").mkdir(parents=True, exist_ok=True)

# unzip sjdlr.zip to localproc
shutil.unpack_archive("tests/testdata/sjdlr.zip", "localproc/sjdlr", "zip")
