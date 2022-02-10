# gamma_pie.py
import matplotlib.pyplot as plt
import argparse
from hmm import *

def parse_args():
    parser = argparse.ArgumentParser(description = 'gamma_pie.py')
    parser.add_argument('-gf', dest = 'gamma_file', default = "")
    parser.add_argument('-mf', dest = 'model_file', default = "")
    parser.add_argument('-s', dest = 'server', default = "False")
    args = parser.parse_args()
    return args

def str2bool(s):
    if s in ["True", 'T', "t", 'true']:
        return True
    return False

if __name__ == '__main__':
    args = parse_args()
    ftitr = args.gamma_file
    modelfile = args.model_file
    server = str2bool(args.server)
    f = open(ftitr)
    

    gamma = []
    
    for l in f:
        gamma = [ float(x) for x in l.split() ]
    sum_gamma = sum(gamma)
    for i in range(0,len(gamma)):
        gamma[i] = gamma[i] / sum_gamma * 100
    
    fig, ax1 = plt.subplots()
    labels = []
    for i,e in enumerate(gamma):
        if e < 2:
            labels.append('')
        else:
            labels.append(str(float('%.1g'%e)) + ', ' + str(i))
            
    ax1.pie(gamma, labels = labels,startangle=90)
    ax1.set_title('Gamma Distribution by State')
    ax1.axis('equal')
    #plt.show()
    name = ftitr[ftitr.rfind('/')+1:]
    if server: # server expects image to reside in static directory
        save_path = './tmp/'+name+'.jpg'
        print(save_path)
        plt.savefig(save_path)
    else:
        plt.savefig(ftitr+'.jpg')

    HMMSTR = hmm(modelfile)
    HMMSTR.read()
    rama = [ 0 for i in range(0, len(HMMSTR.emissions["RAMA"])) ]
    if HMMSTR.n != len(gamma):
        print("model file doesn't match gamma input ( unequal n ) .")
        raise
    for i in range(0, HMMSTR.n):
        for k in range(0, len(rama)):
            rama[k] += gamma[i] * HMMSTR.node_data[i]["RAMA"][k]
    RAMA_labels = HMMSTR.emissions["RAMA"]
    fig, ax1 = plt.subplots()
    labels = []
    # reweight rama distribution according to reweight values
    # character(len=mrama),parameter :: rama_string="HGBEdbeLlxc"
    for i in range(0, len(rama)):
        if i < 2:
            rama[i] = rama[i] * HMMSTR.rama_reweight[0]
        elif i == 2 or i == 5:
            rama[i] = rama[i] * HMMSTR.rama_reweight[1]
        elif i == 3 or i == 4 or i == 6:
            rama[i] = rama[i] * HMMSTR.rama_reweight[2]
        elif i == 7 or i == 8:
            rama[i] = rama[i] * HMMSTR.rama_reweight[3]
        elif i == 9:
            rama[i] = rama[i] * HMMSTR.rama_reweight[4]
        else:
            rama[i] = rama[i] * HMMSTR.rama_reweight[5]
    sum_rama = sum(rama)
    for i in range(0, len(rama)):
        rama[i] = (rama[i] / sum_rama) * 100
    print(rama)
    for i,e in enumerate(rama):
        labels.append(str(float('%.1g'%e)) + ', ' + RAMA_labels[i])
    ax1.pie(rama, labels = labels, startangle = 90)
    ax1.set_title('Gamma-voted Ramachandran Distribution.')
    ax1.axis('equal')
    if server:
        save_path = './tmp/'+name+'rama.jpg'
        print(save_path)
        plt.savefig(save_path)
    else:
        plt.savefig(ftitr+'rama.jpg')

    # make ss pie chart
    ss = [ 0 for i in range(0, len(HMMSTR.emissions["SS"]))]
    for i in range(0, HMMSTR.n):
        for k in range(0, len(ss)):
            ss[k] += gamma[i] * HMMSTR.node_data[i]["SS"][k]
    SS_labels = HMMSTR.emissions["SS"] # HEGST_
    fig, ax1 = plt.subplots()
    labels = []
    for i in range(0, len(ss)):
        if i == 0 or i == 2:
            ss[i] = ss[i] * HMMSTR.ss_reweight[0]
        elif i >= 3:
            ss[i] = ss[i] * HMMSTR.ss_reweight[1]
        else:
            ss[i] = ss[i] * HMMSTR.ss_reweight[2]
    sum_ss = sum(ss)
    for i in range(0, len(ss)):
        ss[i] = (ss[i] / sum_ss) * 100
    print(ss)
    for i , e in enumerate(ss):
        labels.append(str(float('%.1g'%e)) + ', ' +SS_labels[i])
    ax1.pie(ss, labels = labels, startangle = 90)
    ax1.set_title('Gamma-voted Secondary Structure Distribution.')
    ax1.axis('equal')
    if server:
        save_path = './tmp/'+name+'ss.jpg'
        print(save_path)
        plt.savefig(save_path)
    else:
        plt.savefig(ftitr+'ss.jpg')
    
        


