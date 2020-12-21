import os
import sys
import pprint
import difflib
import inspect
import re
from subprocess import Popen, PIPE

bindir = "../build"

regressiondir = "data/regression"
tmpdir = "tmp"

def cat(cmd):
    head(cmd, -1)

def head(cmd, lines=4):
    # out = subprocess.check_output(f"{bindir}/{cmd}",shell=True)
    cmd2 = f"{bindir}/{cmd}".split()
    # process = Popen(cmd2, stdout=PIPE, stderr=PIPE, shell=True)
    # out = process.communicate()[0]
    p = Popen(cmd2, stdout=PIPE, stderr=PIPE)
    output = p.communicate()
    out = output[0]
    if len(out) == 0:
        # if stdout is empty fetch stderr
        out = output[1]
    # out=subprocess.check_output(cmd2, universal_newlines=True)
    header = out.decode().expandtabs(tabsize=8).split("\n")[0:lines]
    header = ['>' if i=='' else i for i in header]
    print("\n".join(header))

def run_stdout(cmd, ext = None):
    os.makedirs(tmpdir,exist_ok=True)
    curframe = inspect.currentframe()
    # pp = pprint.PrettyPrinter(indent=4)
    # pp.pprint(inspect.getouterframes(curframe))
    # print("------------n")
    calframe = inspect.getouterframes(curframe, 1)
    p = re.compile('\[([0-9])\]')
    index = p.findall(calframe[1].filename)[0]

    name = calframe[1].filename[0:-4]
    if "doctest" in name:
        name = name[9:-3]
    else:
        name = name[18:-1]

    name += "_"+index
    if ext:
        name += "." + ext

    tmpfn = tmpdir + "/" + name
    os.system(f"{bindir}/{cmd} > {tmpfn}")
    cmpfn = regressiondir+"/"+name
    sys.stdout.writelines(difflib.unified_diff(open(cmpfn).readlines(),open(tmpfn).readlines(),cmpfn,tmpfn,n=1))
    print(f"output in <a href=\"../data/regression/{name}\">{name}</a>")
