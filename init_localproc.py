import shutil
from pathlib import Path

# make a directory called 'localproc'
Path("localproc").mkdir(parents=True, exist_ok=True)

# now copy data from testdata to localproc
shutil.copytree("tests/testdata/sjdlr", "localproc/sjdlr", dirs_exist_ok=True)
shutil.copytree("tests/testdata/rv3s", "localproc/rv3s", dirs_exist_ok=True)
