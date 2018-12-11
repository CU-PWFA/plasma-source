def format(infile, outfile):
	with open(infile, 'r') as f, open(outfile, 'w') as f0:
		for line in f:
			f0.write(line.replace(',',''))