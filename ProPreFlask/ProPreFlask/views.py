"""
Routes and views for the flask application.
"""

# Important modules for our ProPre
from ProPreFlask.build_seq import build_seq
from io import StringIO
from Bio import SeqIO
from pymol import cmd
from flask import jsonify

import dnachisel
from dnachisel.biotools import reverse_translate
from Bio.PDB import PDBParser

from itertools import product
from collections import deque


from datetime import datetime
from flask import render_template

import os
from flask import Flask, flash, request, redirect, url_for
from werkzeug.utils import secure_filename

from flask import send_file
from ProPreFlask import app

@app.route('/')
@app.route('/home')
def home():
    """Renders the home page."""
    return render_template(
        'index.html',
        title='Home Page',
        year=datetime.now().year,
    )

@app.route('/contact')
def contact():
    """Renders the contact page."""
    return render_template(
        'contact.html',
        title='Contact',
        year=datetime.now().year,
        message='Your contact page.'
    )

@app.route('/about')
def about():
    """Renders the about page."""
    return render_template(
        'about.html',
        title='About',
        year=datetime.now().year,
        message='Your application description page.'
    )


'''
Old code for protein processing:
'''

def translate_rna(s):
    '''
    Returns the AminoAcids Sequence from RNA Codons.
    '''
    codon2aa = {"AAA":"K", "AAC":"N", "AAG":"K", "AAU":"N", 
                "ACA":"T", "ACC":"T", "ACG":"T", "ACU":"T", 
                "AGA":"R", "AGC":"S", "AGG":"R", "AGU":"S", 
                "AUA":"I", "AUC":"I", "AUG":"M", "AUU":"I", 

                "CAA":"Q", "CAC":"H", "CAG":"Q", "CAU":"H", 
                "CCA":"P", "CCC":"P", "CCG":"P", "CCU":"P", 
                "CGA":"R", "CGC":"R", "CGG":"R", "CGU":"R", 
                "CUA":"L", "CUC":"L", "CUG":"L", "CUU":"L", 

                "GAA":"E", "GAC":"D", "GAG":"E", "GAU":"D", 
                "GCA":"A", "GCC":"A", "GCG":"A", "GCU":"A", 
                "GGA":"G", "GGC":"G", "GGG":"G", "GGU":"G", 
                "GUA":"V", "GUC":"V", "GUG":"V", "GUU":"V", 

                "UAA":"_", "UAC":"Y", "UAG":"_", "UAU":"T", 
                "UCA":"S", "UCC":"S", "UCG":"S", "UCU":"S", 
                "UGA":"_", "UGC":"C", "UGG":"W", "UGU":"C", 
                "UUA":"L", "UUC":"F", "UUG":"L", "UUU":"F"}

    l = [codon2aa.get(s[n:n+3], 'X') for n in range(0, len(s), 3)]
    return "".join(l)


def translate(seq):
    '''
    Returns the AminoAcids Sequence from DNA Codons.
    '''
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein = ""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein+= table[codon]
    return protein

def read_seq(inputfile):
    '''
    Returns the sequence of DNA Codons without line endings.
    '''
    with open(inputfile, "r") as f:
        seq = f.read()
    seq = seq.replace("\n", "")
    seq = seq.replace("\r", "")
    return seq

def dna2rna(seq):
    '''
    Returns the RNA version of a DNA.
    '''
    temp = ''
    temp = str(seq).replace('T', 'U')
    return temp


def readfasta(path=''):
    '''
    Returns all sequences from a fasta file.
    '''
    Reads = str()
    with open(path, 'r') as lines:
        for line in lines:
            if not line.startswith('>'):
                line = line.replace("\n", "")
                line = line.replace("\r", "")
                Reads += line[1:]
    return Reads

'''
reads = readfasta()
 
for read in reads:
    rna = dna2rna(read)
    print(translate_rna(rna))
'''


document = '_document'
filename = 'P69905.fasta'

# Different types of helix visualization.

#_type = '3/10 helix'
_type = 'helix'
#_type = 'parallel'
#_type = 'polypro'
#_type = 'antiparallel'

import flask

@app.route('/GetAminoSeq', methods=['PUT'])
def GetAminoSeq():
    #>My fasta file with headers
    #AGAWDFAWDAWDAD.......
    seq = flask.request.json
    seq = '>Alpha\n' + seq + '\n' # Virtual file handle from sequence string.
    for seq_record in SeqIO.parse(StringIO(seq), "fasta"):
        build_seq(seq_record.seq, _type) # The actual sequence.
        cmd.cmd.select(document,"all") # Determining document type.
        cmd.cmd.save(f"temp{str(seq)[-15:-1]}.pdb", document, -1, 'pdb') # Protein Data Bank file format.
        #cmd.cmd.save(f"temp{str(seq)[-15:-1]}.fasta", document, -1, 'fasta')
        cmd.cmd.delete("all") # Delete unnecessary.
    return send_file(f'../temp{str(seq)[-15:-1]}.pdb') # Flask method to return files as content response result.


@app.route('/PDB2SEQ', methods=['POST'])
def GetSequenceFromPDB():
    # Get file from request
    uploaded_file = request.files['file']
    if uploaded_file.filename != '':
        uploaded_file.save(uploaded_file.filename)
    # You can use a dict to convert three letter code to one letter code
    d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    # Just an example input pdb
    record = f'{uploaded_file.filename}'
    # run parser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', record)    
    # iterate each model, chain, and residue
    # printing out the sequence for each chain

    totalFileTest = str()
    for model in structure:
        for chain in model:
            seq = []
            for residue in chain:
                try:
                    seq.append(d3to1[residue.resname])
                except:
                    pass
            totalFileTest = '>some_header\n' + ''.join(seq)
    myfileName = get_file_name_no_extension(uploaded_file)
    with open(f'{myfileName}.fasta', 'w') as f:
        f.write(totalFileTest)

    recordNew = dnachisel.load_record(f"{myfileName}.fasta")
    hhhh = reverse_translate(str(recordNew.seq))
    return str(hhhh)

def get_file_name_no_extension(file) -> str:
    filename_parts = file.filename.split('.')
    if (len(filename_parts) > 2):
        name_parts = [filename_parts[i] for i in range(len(filename_parts) - 1)]
        return ".".join(name_parts)
    return filename_parts[0]



@app.route('/GetAminoFas', methods=['POST'])
def GetAminoFas():
    #>My fasta file from response files
    uploaded_file = request.files['file']
    if uploaded_file.filename != '':
        uploaded_file.save(uploaded_file.filename)
        for seq_record in SeqIO.parse(f"../{uploaded_file.filename}", "fasta"):
            build_seq(seq_record.seq, _type) # The actual sequence.
            cmd.cmd.select(document,"all") # Determining document type.
            cmd.cmd.save(uploaded_file.filename, document, -1, 'pdb') # Protein Data Bank file format.
            cmd.cmd.delete(uploaded_file.filename) # Delete unnecessary.
    return send_file(f'../{uploaded_file.filename}') # Flask method to return files as content response result.


@app.route('/GetAlignScore/<string:seq>/<string:diseaseName>')
def GetAlignSeq(seq='', diseaseName=''):
    testSeq = readfasta(app.root_path + '/static/fastaDiseases/' + diseaseName)
    #testSeq = 'ATGGTTGCTGCCAACCCAGAAAAG'
    tempScore = list()
    # Old slow brute force method.
    #alignments = list(all_alignments(testSeq, seq))
    # New needleman_wunsch method.
    myAlignment = print_alignment(testSeq, seq, needleman_wunsch(testSeq, seq))
    return jsonify(
            Score = align_score_nw(myAlignment),
            Alignment = myAlignment
        )

def align_score_nw(alignment = [], matchScr = 1, gabScr = 0, misScr = -1):
    '''
    Returns the score of an alignment in X%, based on score schema used.
    '''
    totalScore = int()
    lengthSeq = len(alignment[0])
    for i in range(0, lengthSeq):
        if alignment[0][i] == alignment[1][i]:
            totalScore+=matchScr
        if alignment[1][i] == '-':
            totalScore+=gabScr
        elif alignment[0][i] != alignment[1][i]:
            totalScore+=misScr
    perfectScore = lengthSeq * matchScr
    realScore = float((totalScore / perfectScore) * 100)
    return realScore

def print_alignment(x, y, alignment):
    line1 = "".join(
        "-" if i is None else x[i] for i, _ in alignment
    )
    line2 = "".join(
        "-" if j is None else y[j] for _, j in alignment
    )
    listRet = list()
    listRet.append(line1)
    listRet.append(line2)
    return listRet


def needleman_wunsch(x, y):
    """Run the Needleman-Wunsch algorithm on two sequences.

    x, y -- sequences.

    Code based on pseudocode in Section 3 of:

    Naveed, Tahir; Siddiqui, Imitaz Saeed; Ahmed, Shaftab.
    "Parallel Needleman-Wunsch Algorithm for Grid." n.d.
    https://upload.wikimedia.org/wikipedia/en/c/c4/ParallelNeedlemanAlgorithm.pdf
    """
    N, M = len(x), len(y)
    s = lambda a, b: int(a == b)

    DIAG = -1, -1
    LEFT = -1, 0
    UP = 0, -1

    # Create tables F and Ptr
    F = {}
    Ptr = {}

    F[-1, -1] = 0
    for i in range(N):
        F[i, -1] = -i
    for j in range(M):
        F[-1, j] = -j

    option_Ptr = DIAG, LEFT, UP
    for i, j in product(range(N), range(M)):
        option_F = (
            F[i - 1, j - 1] + s(x[i], y[j]),
            F[i - 1, j] - 1,
            F[i, j - 1] - 1,
        )
        F[i, j], Ptr[i, j] = max(zip(option_F, option_Ptr))

    # Work backwards from (N - 1, M - 1) to (0, 0)
    # to find the best alignment.
    alignment = deque()
    i, j = N - 1, M - 1
    while i >= 0 and j >= 0:
        direction = Ptr[i, j]
        if direction == DIAG:
            element = i, j
        elif direction == LEFT:
            element = i, None
        elif direction == UP:
            element = None, j
        alignment.appendleft(element)
        di, dj = direction
        i, j = i + di, j + dj
    while i >= 0:
        alignment.appendleft((i, None))
        i -= 1
    while j >= 0:
        alignment.appendleft((None, j))
        j -= 1

    return list(alignment)


def alignment_score(x, y, alignment):
    """Score an alignment.

    x, y -- sequences.
    alignment -- an alignment of x and y.

    different is 0, gap is 1 & same are 2.
    """
    score_gap = +1
    score_same = +2
    score_different = 0

    score = 0
    for i, j in alignment:
        if (i is None) or (j is None):
            score += score_gap
        elif x[i] == y[j]:
            score += score_same
        elif x[i] != y[j]:
            score += score_different

    return score

def all_alignments(x, y):
    """Return an iterable of all alignments of two
    sequences.

    x, y -- Sequences.
    """

    def F(x, y):
        """A helper function that recursively builds the
        alignments.

        x, y -- Sequence indices for the original x and y.
        """
        if len(x) == 0 and len(y) == 0:
            yield deque()

        scenarios = []
        if len(x) > 0 and len(y) > 0:
            scenarios.append((x[0], x[1:], y[0], y[1:]))
        if len(x) > 0:
            scenarios.append((x[0], x[1:], None, y))
        if len(y) > 0:
            scenarios.append((None, x, y[0], y[1:]))

        # NOTE: "xh" and "xt" stand for "x-head" and "x-tail",
        # with "head" being the front of the sequence, and
        # "tail" being the rest of the sequence. Similarly for
        # "yh" and "yt".
        for xh, xt, yh, yt in scenarios:
            for alignment in F(xt, yt):
                alignment.appendleft((xh, yh))
                yield alignment

    alignments = F(range(len(x)), range(len(y)))
    return map(list, alignments)