# server.py
# HMMSTR server
__author__ = 'TLB'

from flask import Flask, flash, redirect, render_template, request, session, url_for
import os
from werkzeug.utils import secure_filename

app = Flask(__name__)

# enusre proper file uplaod
ALLOWED_EXTENSIONS = '.pdb'
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/')
def index():
    return render_template('test.html')

def seq2html(seqfile):
    s = "<tt><pre>\n"
    f = open(seqfile)
    for l in f:
        s += l 
    s += '</pre></tt>'
    return s

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
        os.system('./RUNPDB.csh ./HMMSTR/HMMSTR_100.hmm ' + pdb_name + ' ' + chain)
        seqfile = './HMMSTR/HMMSTR_100'+pdb_name+'.seq'
        html_string = """ <head>
        <link rel="stylesheet" type="text/css" href="static/test.css">
        </head>"""
        html_string += seq2html(seqfile)
        html_string += '<a href="/">HOME</a> <br>'
        html_string += '<a href="/shutdown">Close Server</a>'
        return html_string
    else:
        html_string = '<a href="/">HOME</a>'
        return html_string
    


@app.route('/inject', methods=['POST'])
def inject():
    sequence_lines = request.form['sequence'].splitlines()
    sequence = ""
    for l in sequence_lines:
        sequence += l.strip()
    code = request.form['code']
    chain = request.form['chain']
    # inject
    os.system('./xinject ./HMMSTR/HMMSTR_100.hmm ' + sequence + ' ' + code + ' ' + chain)
    # seq2html
    seqfile = './HMMSTR/HMMSTR_100'+code+'.seq'
    html_string = """ <head>
    <link rel="stylesheet" type="text/css" href="static/test.css">
    </head>"""
    html_string += seq2html(seqfile)
    html_string += '<a href="/">HOME</a> <br>'
    html_string += '<a href="/shutdown">Close Server</a>'
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
    app.run(host = '0.0.0.0', port = 4481)
