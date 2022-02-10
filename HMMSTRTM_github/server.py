# server.py
# HMMSTR server
__author__ = 'TLB'

from flask import Flask, flash, redirect, render_template, request, session, url_for
import os
from werkzeug.utils import secure_filename
import matplotlib.pyplot as plt
app = Flask(__name__)

# do some tracking of the desired HMM file
HMM_file = './HMMSTR/HMMSTR_100.hmm'
HMM_file_prefix = HMM_file[:HMM_file.rfind('.')]
gamma_dir = HMM_file_prefix + 'GAMMA/'
from hmm import *
HMM = hmm(HMM_file)
HMM.read()
n = HMM.n

html_tail = '<a href="/">HOME</a> &nbsp;&nbsp;&nbsp; <a href="/shutdown">Close Server</a>'


# enusre proper file uplaod
ALLOWED_EXTENSIONS = '.pdb'
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/')
def index():
    return render_template('home.html')

def get_seq_len(seqfile):
    f = open(seqfile)
    i = 0
    for l in f:
        if i == 2:
            nres = len(l) - 15
        i += 1
    return nres

def seq2html(seqfile):
    s = "" 
    f = open(seqfile)
    line_length = 0 
    i = 0
    name = ""
    for l in f:
        if i == 0:
            name = l.strip()
        elif i == 2:
            line_length = len(l) - 15 # 14 is the number of spaces in pre-formating e.g. 'AA SEQUENCE  '
        i += 1
        s += l
    s += '</pre></tt>'
    html =  "<tt><pre>\n"
    html += " " * 14
    for l in range(0,int(line_length/10)+1):
        spaces = 10
        strl = str(l)
        spaces -= len(strl)
        html += strl
        html += " " * spaces
    html += "\n"
    html += " " * 14
    for l in range(0, line_length):
        '<a href="/shutdown">Close Server</a>'
        html += '<a href="/gamma-pie/'+name+'/'+str(l)+'">'+str(l%10)+'</a>'
    html += "\n"
    html += s
    return html

@app.route('/run_pdb', methods=['POST', 'GET'])
def run_pdb():
    if request.method == 'POST':
        pdb_file = request.files['pdbfile']
        if pdb_file.filename == '':
            flash('No selected file.')
            return redirect(request.url)
        if pdb_file and allowed_file(pdb_file.filename):
            filename = secure_filename(pdb_file.filename)
            pdb_file.save(os.path.join('./tmp', filename))
        pdb_name = pdb_file.filename[:pdb_file.filename.rfind('.')]
        chain = request.form['chain']
        os.system('./RUNPDB.csh ' + HMM_file + ' '  + pdb_name + ' ' + chain)
        seqfile = HMM_file_prefix+pdb_name+'.seq'
        html_string = """ <head>
        <link rel="stylesheet" type="text/css" href="static/test.css">
        </head>"""
        html_string += seq2html(seqfile)
        html_string += html_tail
        return html_string
    else:
        html_string = '<a href="/">HOME</a>'
        return html_string
    
@app.route('/gamma-pie/<string:name>/<int:titr>')
def gamma_pie(name, titr):
    # get the gamma pie chart for a given time
    filename = gamma_dir+name
    seqfile = HMM_file_prefix + name + '.seq'
    nres = get_seq_len(seqfile)
    s = './xreadgamma ' + HMM_file + ' ' + filename + ' ' + str(titr) + ' ' + str(nres) + ' > ' + filename + str(titr)
    print(s)
    os.system(s)
    ftitr = filename + str(titr)
    s = ' python3 gamma_pie.py -gf ' + ftitr + " -s True -mf " + HMM_file
    print(s)
    os.system(s) # create gamma pie plot .jpg file
    print(name+str(titr)+'.jpg')
    html = '<h3> GAMMA VOTING PIE CHARTS FOR RESIDUE ' + str(titr) + ' OF ' + name.upper() + ' </h3><br>'
    html += '<img src=../../static/'+name+str(titr)+'.jpg alt="Gamma piechart" width="450" height="333">'
    html += '<img src=../../static/'+name+str(titr)+'rama.jpg alt="Gamma-voted Ramachandran prediction piechart" width="450" height="333">'
    html += '<img src=../../static/'+name+str(titr)+'ss.jpg alt="Gamma-voted Secondary Structure prediction piechart" width="450" height="333">'
    html += html_tail
    return html


@app.route('/inject', methods=['POST'])
def inject():
    sequence_lines = request.form['sequence'].splitlines()
    sequence = ""
    for l in sequence_lines:
        sequence += l.strip()
    code = request.form['code']
    chain = request.form['chain']
    if request.form.get('RUN_BLAST'):
        # create tmp/*.fasta
        fasta_f = './tmp/'+code+chain+'.fasta'
        f = open(fasta_f, 'w')
        f.write('>'+code+':'+chain+'\n')
        f.write(sequence)
        f.close()
        os.system('./BLAST_RUN.csh ' + HMM_file + ' ' + code + ' ' + chain)
    else:
        # inject
        os.system('./xinject ' + HMM_file + ' ' + sequence + ' ' + code + ' ' + chain)
    # seq2html
    seqfile = HMM_file_prefix+code+'.seq'
    html_string = """ <head>
    <link rel="stylesheet" type="text/css" href="static/test.css">
    </head>"""
    html_string += seq2html(seqfile)
    html_string += html_tail
    return html_string

def shutdown_server():
    func = request.environ.get('werkzeug.server.shutdown')
    if func is None:
        raise RuntimeError('Not running with the Werkzeug Server')
    func()
    
@app.get('/shutdown')
def shutdown():
    shutdown_server()
    return 'Server shutting down...'

if __name__ == '__main__':
    app.secret_key = 'super secret key'
    app.config['SESSION_TYPE'] = 'filesystem'
    app.debug = True
    app.run(host = '0.0.0.0', port = 4534)
