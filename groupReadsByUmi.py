import pysam
import sys

bamfile = sys.argv[1]

output = sys.argv[2]
print(output)
samfile = pysam.AlignmentFile(bamfile, "rb", check_sq = False)

umis = {}

for read in samfile.fetch(until_eof=True):
  if read.is_read1:
    umi = read.get_tag("RX")
    readname = read.query_name
    try:
      umis[umi].append(readname)
    except KeyError:
      umis[umi] = []
      umis[umi].append(readname)
  else:
    continue

groups = {}

for i, key in enumerate(umis.keys()):
  for j in umis[key]:
    groups[j] = i
  
samfile = pysam.AlignmentFile(bamfile, "rb", check_sq = False)
outputfile = pysam.AlignmentFile(output, "wb", template = samfile)
for read in samfile.fetch(until_eof=True):
  readname = read.query_name
  group = groups[readname]
  tags = read.get_tags()
  tags.append(("MI",str(group)))
  read.set_tags(tags)
  outputfile.write(read)
  
