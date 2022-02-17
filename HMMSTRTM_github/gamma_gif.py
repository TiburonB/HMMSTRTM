# gamma_gif.py
# generate a gif of all time points of a single gamma drct file

from hmm import *
from hmm_viz import *
import os
from os import listdir
from os.path import isdir
import argparse
from gamma_heat import * 

def get_seq_len(seqfile):
    f = open(seqfile)
    i = 0
    for l in f:
        if i == 2:
            nres = len(l) - 15
        i += 1
    return nres

def parse_args():
    parser = argparse.ArgumentParser(description = 'gamma_gif.py')
    parser.add_argument('-mf', dest = 'model_file', default = "./HMMSTR/HMMSTR_100.hmm")
    parser.add_argument('-gf', dest = 'gamma_file', default = './HMMSTR/HMMSTR_100GAMMA/1cis') 
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    model_file = args.model_file
    gamma_file = args.gamma_file

    HMM = hmm(model_file)
    HMM.read()
    name = gamma_file[gamma_file.rfind('/')+1:]  
    gamma_dir = gamma_file + 'dir'
    if not(isdir(gamma_dir)):
        print('mkdir ' + gamma_dir)
        os.system('mkdir ' + gamma_dir)
    
    model_file_prefix = model_file[:model_file.rfind('.')]
    seqfile = model_file_prefix + name + '.seq'
    nres = get_seq_len(seqfile) # get nres ( chain length ) 
    for titr in range(1, nres):
        out_file = gamma_dir + '/' + str(titr)
        s = './xreadgamma ' + model_file + ' ' + gamma_file + ' ' + str(titr) + ' ' + str(nres) + ' > ' + out_file
        print(s)
        os.system(s)
    jpg_files = []
    # run seq2jpg
    # os.system('python3 seq2jpg.py -sf ' + seqfile + ' -od ' + gamma_dir + '/')
    # this doesn't work, need to gather commands and run mannually ( annoying... ) 

    for titr in range(1, nres):
        jpg_file = gamma_dir + '/heat_'
        stitrprefix = ""
        if titr < 10:
            stitrprefix = '00'
        elif titr < 100:
            stitrprefix += '0'
        jpg_file = stitrprefix + str(titr)
        heat_file = gamma_dir + '/heat_'+jpg_file+'.jpg'
        # create gamma heat map
        s = 'python3 gamma_heat.py -gf ' + gamma_dir+'/'+str(titr) + ' -s False -mf ' + model_file + ' -of ' + heat_file[:len(heat_file)-4]
        print(s)
        os.system(s)
        # create gamma pie charts
        s = 'python3 gamma_pie.py -gf ' + gamma_dir +'/'+str(titr) + " -s False -mf " + model_file + ' -od ' + gamma_dir+'/'
        print(s)
        os.system(s)
        
        # stitch pie-chart images together
        gamma_pie_file = gamma_dir +'/gamma_pie_'+stitrprefix+str(titr)+'.jpg'
        gamma_rama_file = gamma_dir +'/gamma_rama_'+stitrprefix+str(titr)+'.jpg'
        gamma_ss_file = gamma_dir +'/gamma_ss_'+stitrprefix+str(titr)+'.jpg'
        # original =  width = 450, height = 999
        # to match heat map, height = 702 
        scaling_factor = float(4135/999)
        width = int(450 * scaling_factor)
        height = int(333 * scaling_factor)
        gamma_pie_outfile = gamma_dir +'/gamma_piecharts'+stitrprefix+str(titr)+'.jpg'

        s = 'montage ' + gamma_pie_file + ' ' + gamma_rama_file + ' ' + gamma_ss_file + ' -tile 1x3 -geometry '+str(width)+'x'+str(height)+'+1+1 ' + gamma_pie_outfile
        print(s)
        os.system(s)
        
        # add time-labeled sequence to gamma heat map image
        htmljpg_file = gamma_dir + '/html_'+jpg_file+'.jpg'
        s = "convert -composite -gravity NorthEast " + heat_file + ' ' + htmljpg_file + ' ' + heat_file
        print(s)
        os.system(s)

        # combine pie charts and gamma heat
        gheatf = gamma_dir + '/gheat_'+stitrprefix+str(titr)+'.jpg' 
        s = 'montage  ' + heat_file + ' ' + gamma_pie_outfile + ' -tile 2x1 -geometry 3708x4135+1+1 ' + gheatf
        print(s)
        os.system(s)


       # jpg_files.append(jpg_file + ".jpg")
    s = 'mencoder "mf://'+gamma_dir+'/gheat*.jpg" -mf fps=6 -o '+gamma_dir+'/'+name+'.avi -ovc lavc -lavcopts vcodec=mjpeg:vbitrate=3000:vhq'
    print(s)
    os.system(s)
    #print(jpg_files)    

