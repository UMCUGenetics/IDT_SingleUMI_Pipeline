import gzip
import sys

filestr = str(sys.argv[1])
filestrsplit = filestr.split(".")[0]
f = gzip.open(f'{filestrsplit}.fastq.gz', "rt")
g = open(str(sys.argv[2]),"w")

nt = int(sys.argv[3])

for i, line in enumerate(f):
  if line.startswith("@"):
    split = line.split()[1]
    split = split.split(":")[3]

    umi = split[8:17]
    nextline = next(f)
    read = nextline[0:nt]
    
    g.write(line)
    g.write(f'{umi}{read}\n')
    g.write(next(f))

    nextline = next(f)
    g.write(f'{nextline[0:9+nt]}\n')
