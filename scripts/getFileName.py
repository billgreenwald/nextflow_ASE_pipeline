import sys
filePath=sys.argv[1]
fn=filePath.split("/")[-1]
fn=fn.split(".")[:-1]
fn=".".join(fn)
print fn