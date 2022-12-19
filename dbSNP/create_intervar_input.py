import sys

for line in open(sys.argv[1]):
    if line.startswith("1 rs"):
        info = line.split()
        out = info[2] + "\t" + info[3] + "\t" + info[4] + "\t" + info[5].split("/")[0] + "\t" + info[5].split("/")[1]
        print(out)
