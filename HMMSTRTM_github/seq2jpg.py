# seq2jpg.py
import argparse
import subprocess
import os

def get_seq_len(seqfile):
    f = open(seqfile)
    i = 0
    for l in f:
        if i == 2:
            nres = len(l) - 15
        i += 1
    return nres

def seq2html(seqfile, t):
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
        try:
            s += l[:14+t-1] + '<mark>' + l[14+t-1] + '</mark>' + l[14+t:]
        except:
            s += l 
    s += '</pre></tt></h1>'
    html =  "<h1><tt><pre>\n"
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
        if l == t -1:
            html += '<mark>'
        html += str(l%10)
        if l == t -1:
            html += '</mark>'
    html += "\n"
    html += s
    return html


def parse_args():
    parser = argparse.ArgumentParser(description = 'gamma_heat.py')
    parser.add_argument('-sf', dest = 'seq_file', default = "./HMMSTR/HMMSTR_1001cis.seq")
    parser.add_argument('-od', dest = 'outdir', default = "/home/tiburon/Warehouse/Summer2021/HMMSTRTM_github/HMMSTR/HMMSTR_100GAMMA/1cisdir/")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    seq_file = args.seq_file
    outdir = args.outdir
    
    nres = get_seq_len(seq_file)
    for titr in range(1, nres+1):
        html = seq2html(seq_file, titr)
        pref = ""
        if titr < 10:
            pref = "00"
        elif titr < 100:
            pref = "0"
        html_file = outdir + 'html_'+pref+str(titr)+'.html'
        jpg_file = outdir + 'html_'+pref+str(titr)+'.jpg'
        f = open(html_file, 'w')
        f.write(html)
        #s = 'wkhtmltoimage -f jpg --width 7416 ' + html_file + ' ' + jpg_file +'\n' 
        s = 'wkhtmltoimage -f jpg '+html_file+ ' ' + jpg_file 
        print(s)
        #res = os.system(s.strip())
        #print(res)
        #list_s = s.split(' ')
        #result = subprocess.run(list_s)
        #print(result)

