#!/usr/bin/env python

# Original work by...
# Rajarshi Guha <rguha@indiana.edu>
# 08/26/07

# RDkit Hack by Chris Arthur 1/11/2015

import sys, getopt

from rdkit import Chem

max_path_len = 9
verbose = False
scale_type = 'raw'
error_file = 'error.txt'

# Schneiders SMARTS's patterns
#cats_smarts = {
#    'D' : ['[#6H]', '[#7H,#7H2]'],
#    'A' : ['[#6]', '[#7H0]'],
#    'P' : ['[*+]', '[#7H2]'],
#    'N' : ['[*-]', '[C&$(C(=O)O),P&$(P(=O)),S&$(S(=O)O)]'],
#    'L' : ['[Cl,Br,I]', '[S;D2;$(S(C)(C))]']
#}

# donor and acceptor SMARTS taken from 
# http://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html#H_BOND
#
# I think the L patterns for the carbon adjacent only to carbons covers
# all possible cases. I have not taken into account primary carbon
# atoms (i.e., methyl groups)
cats_smarts = {
    'D' : ['[!$([#6,H0,-,-2,-3])]'],
    'A' : ['[!$([#6,F,Cl,Br,I,o,s,nX3,#7v5,#15v5,#16v4,#16v6,*+1,*+2,*+3])]'],
    'P' : ['[*+]', '[#7H2]'],
    'N' : ['[*-]', '[C&$(C(=O)O),P&$(P(=O)),S&$(S(=O)O)]'],
    'L' : ['[Cl,Br,I]', '[S;D2;$(S(C)(C))]',
	   '[C;D2;$(C(=C)(=C))]', '[C;D3;$(C(=C)(C)(C))]', '[C;D4;$(C(C)(C)(C)(C))]',
	   '[C;D3;H1;$(C(C)(C)(C))]', '[C;D2;H2;$(C(C)(C))]',]
}
    
cats_desc = ['DD','AD', 'DP', 'DN', 'DL',
	     'AA', 'AP', 'AN', 'AL',
	     'PP', 'NP', 'LP',
	     'NN', 'LN', 'LL']


def getPcoreGroups(mol, smarts):
    """
    Given a molecule it assigns PPP's to individual atoms

    The return value is a list of length number of atoms in
    the input molecule. The i'th element of the list contains
    a list of PPP labels that were identified for the i'th atom
    """

    ret = ['' for x in range(0,mol.GetNumAtoms())]

    labels = smarts.keys()
    for label in labels:
	patterns = smarts[label]
	for pattern in patterns:
         
	    patt = Chem.MolFromSmarts(pattern)
	    matched = False
	    for matchbase in mol.GetSubstructMatches(patt, uniquify=True):
		for idx in matchbase:
		    if ret[idx] == '': ret[idx] = [label]
		    else: 
			tmp = ret[idx]
			tmp.append(label)
			ret[idx] = tmp
		matched = True
	    if matched: break
    return ret

def _getZeroMatrix(r,c):
    return [[0 for x in range(0,c)] for x in range(0,r)]

def getAdjacencyMatrix(mol):
    """
    Generates an adjacency matrix for a molecule. Note that this
    will be a matrix with 0 along the diagonals
    """
    n = mol.GetNumAtoms()
    admat = _getZeroMatrix(n,n)
    for bond in mol.GetBonds():
	bgn_idx = bond.GetBeginAtomIdx()
	end_idx = bond.GetEndAtomIdx()
	admat[bgn_idx][end_idx] = 1
	admat[end_idx][bgn_idx] = 1
    return admat

def getTopologicalDistanceMatrix(admat):
    """
    Generates the topological distance matrix given
    an adjacency matrix
    """
    n = len(admat)
    d = _getZeroMatrix(n,n)
    for i in range(0,n):
	for j in range(0,n):
	    if admat[i][j] == 0: d[i][j] = 99999999;
	    else: d[i][j] = 1
    for i in range(0,n): d[i][i] = 0
    
    for k in range(0,n):
	for i in range(0,n):
	    for j in range(0,n):
		if d[i][k]+d[k][j] < d[i][j]: d[i][j] = d[i][k]+d[k][j]
    return d

def getPPPMatrix(admat, ppp_labels):
    pppm = {}
    n = len(admat)
    for i in range(0,n):
	ppp_i = ppp_labels[i]
	if ppp_i == '': continue
	for j in range(0,n):
	    ppp_j = ppp_labels[j]
	    if ppp_j == '': continue
	    pairs = []
	    for x in ppp_i:
		for y in ppp_j:
		    if (x,y) not in pairs and (y,x) not in pairs:
			## make sure to add the labels in increasing
			## lexicographical order
			if x < y: tmp = (x,y)
			else: tmp = (y,x)
			pairs.append(tmp)
	    pppm[i,j] = pairs
    return pppm
	    
def usage():
    print """
    cats2d.py [options]

    Generates CATS2D descriptors which are topological pharmacophore
    fingerprints.

    -e, --error    Error log file. Default is error.txt
    -h, --help     Print this message
    -i, --in       Name of the input SMILES file
    -l, --len      Maximum path length (default is 9)
    -o, --out      Name of the output file (default is cats2d.txt)
    -s, --scale    What type of scaling to use. The possible values
                   are: raw, num and occ. The default is raw and occ
                   is not implemented at the moment
    -v, --verbose  Verbose output
    """

if __name__ == '__main__':
    
    ifilename = None
    ofilename = 'cats2d.txt'

    if len(sys.argv) == 1:
	usage()
	sys.exit(-1)

    try:
	opt,args = getopt.getopt(sys.argv[1:], 'o:i:l:s:e:hv', \
				 ['help', 'in=', 'out=', 'verbose', 
				  'len=', 'scale=', 'error='])
    except getopt.GetoptError:
	usage()
	sys.exit(-1)

    for o,a in opt:
	if o in ('-e', '--error'):
	    error_file = a
	if o in ('-s', '--scale'):
	    if a not in ['raw', 'num', 'occ']:
		usage()
		sys.exit(-1)
	    else: scale_type = a
	if o in ('-h', '--help'):
	    usage()
	    sys.exit(-1)
	if o in ('-i', '--in'):
	    ifilename = a
	if o in ('-o', '--out'):
	    ofilename = a
	if o in ('-v', '--verbose'):
	    verbose = True
	if o in ('-l', '--len'):
	    max_path_len = int(a)
    
    if max_path_len < 1:
	print 'Using max_path_len = 9'
	max_path_len = 9

    if not ifilename:
	print 'Must specify an input file'
	usage()
	sys.exit(-1)

    ofile = open(ofilename, 'w')
    for label in cats_desc:
	for i in range(0, max_path_len+1):
	    ofile.write('%s.%d ' % (label,i))
    ofile.write('\n')

    
    suppl = Chem.SmilesMolSupplier(ifilename)


    nmol = 0
    for mol in suppl:
	natom = mol.GetNumAtoms()
 
	ppp_labels = getPcoreGroups(mol, cats_smarts)
	admat = getAdjacencyMatrix(mol)
	tdistmat = getTopologicalDistanceMatrix(admat)
	pppmat = getPPPMatrix(tdistmat, ppp_labels)

	# get the occurence of each of the PPP's
	ppp_count = dict(zip(['D','N', 'A', 'P', 'L'], [0]*5))
	for label in ppp_labels:
	    for ppp in label:
		ppp_count[ppp] = ppp_count[ppp] + 1

	# lets calculate the CATS2D raw descriptor
	desc = [ [0 for x in range(0,max_path_len+1)] for x in range(0,15)]
	for x,y in pppmat.keys():
	    labels = pppmat[x,y]
	    dist = tdistmat[x][y]
	    if dist > max_path_len: continue
	    for pair in labels:
		id = '%s%s' % (pair[0],pair[1])
		idx = cats_desc.index(id)
		vals = desc[idx]
		vals[dist] += 1
		desc[idx] = vals
	
	if scale_type == 'num':
	    for row in range(0, len(desc)):
		for col in range(0, len(desc[0])):
		    desc[row][col] = float(desc[row][col]) / natom
	elif scale_type == 'occ':
	    #  get the scaling factors
	    facs = [0]*len(cats_desc)
	    count = 0
	    for ppp in cats_desc:
		facs[count] = ppp_count[ppp[0]] + ppp_count[ppp[1]]
		count += 1

	    # each row in desc corresponds to a PPP pair
	    # so the scale factor is constant over cols of a row
	    count = 0
	    for i in range(0, len(desc)):
		if facs[i] == 0: continue
		for j in range(0, len(desc[0])):
		    desc[i][j] = desc[i][j] / float(facs[i])

#	ofile.write('%s ' % (mol.GetTitle()))
	for row in desc:
	    for col in row: ofile.write('%3.5f ' % (col))
	ofile.write('\n')

	nmol += 1
	if nmol % 100 == 0 and verbose:
	    sys.stdout.write('\rProcessed %d molecules' % (nmol))
	    sys.stdout.flush()


    ofile.close()



	
