# datamine.py
# 1/3/22
# interact with drct_io.f95 
import os
import pickle
import argparse
from os import  listdir
from os.path import isfile
from tplot import view_np, view_prob_line_plot
from tqdm import tqdm

# return a 1d array of lengths of consequtive emission of datafact in data_dict
def get_period(data_dict, datafact, emission):
    data = [ k[datafact] for k in data_dict ]
    periods = []
    for seq in tqdm(data):
        start = -1
        if seq.find(emission) == -1:
            continue
        for i in range(0, len(seq)):
            if seq[i] == emission and seq[i-1] != emission:
                start = i
            elif start != -1 and seq[i] != emission:
                periods.append(i-start)
                start = -1
    return periods

    # get probability distribution of lth TMH polar residue being i residues from END of the helix
def polar_res(l, length, tm_type = "H", data_file = './tm_training.drct'): 
    data = datamine(data_file)
    print('finding distribution of ' + str(l) + "th polar residue in helix.")
    polar = [0 for i in range(0, length)] # length = 27
    for datum in tqdm(data):
        start = -1
        polar_this_TM = 0
        lth_polar_ind = -1
        if datum['TMSEQ'].find('H') != -1: # sequence contains TM
            for i in range(0, len(datum["AASEQ"])):
                if datum["TMSEQ"][i] == "H" and datum["TMSEQ"][i-1] != "H":
                    start = i
                    lth_polar_ind = - 1
                    polar_this_TM = 0
                elif start != -1 and datum["TMSEQ"][i] != "H":
                    if lth_polar_ind != -1 and i-lth_polar_ind < len(polar):
                        polar[i-lth_polar_ind] += 1
                    lth_polar_ind = - 1
                    start = -1
                    polar_this_TM = 0 
                elif datum["TMSEQ"][i] != tm_type or start == -1: 
                    continue
                if ( i - start ) < length and ( i-start ) > 7: 
                    if datum["AASEQ"][i] in ["R", "K", "D", "E", "Q"]: 
                        polar_this_TM += 1
                        if polar_this_TM == l:
                            lth_polar_ind = i
    return polar


# return list that represents the probability the first polar residue in a TMH is in position x
# array is 1d of length = max TM periodicity ( for H, this is ~30).
# in position i is the probability of position i in a PDBTM-labeled TMH being a polar residue
def get_tm_polar_res(data, length, polar_number): # polar number specifies which # polar residue in a TMH to observe the distribution of
    # if polar_number = 1 then we report data on the first observed polar residue in the TMH
    prob_polar = [(0,0) for i in range(0, length)]
    for datum in tqdm(data):
        start = -1
        polar_this_TM = 0 
        if datum["TMSEQ"].find("H") != -1: # sequence contains TMH
            for i in range(0, len(datum["AASEQ"])):
                if datum["TMSEQ"][i] == "H" and datum["TMSEQ"][i-1] != "H":
                    start = i
                    polar_this_TM = 0 
                elif start != -1 and datum["TMSEQ"][i] != "H":
                    start = -1
                    polar_this_TM = 0 
                    continue
                elif datum["TMSEQ"][i] != "H" or start == -1:
                    continue
                if ( i - start ) < length and ( i - start ) > 8: 
                    if datum["AASEQ"][i] in ["R", "K", "D", "E", "Q"]: 
                        if polar_this_TM == polar_number -1:
                            prob_polar[i-start] = (prob_polar[i-start][0] + 1, prob_polar[i-start][1])
                        polar_this_TM += 1
                    else:
                        if polar_this_TM == polar_number -1:
                            prob_polar[i-start] = (prob_polar[i-start][0], prob_polar[i-start][1] + 1)
    
    # polar_probability_distribution = [ i[0]/(i[0]+i[1]) for i in prob_polar ]
    ppd = []
    for i in prob_polar:
        d = i[0] + i[1]
        if d == 0:
            ppd.append(0)
        else:
            ppd.append(i[0]/d)
    return ppd


def parse(data):
    gen_data = ( k for k in data)
    data_dict = []
    while True:
        try:
            datum = {}
            datum["AASEQ"] = next(gen_data)
            datum["TMSEQ"] = next(gen_data)
            datum["SSSEQ"] = next(gen_data)
            datum["RAMASEQ"] = next(gen_data)
            cc = next(gen_data)
            datum["CODE"] = cc[:4]
            datum["CHAIN"] = cc[5]
            data_dict.append(datum)
        except:
            return data_dict

def datamine(data_file, code = "null") : # returns datadict
    # compile drct_io.f95 #
    os.system('gfortran drct_io.f95 -o dio.out')
    # ################### #  

    # run the thing
    s = './dio.out ' + data_file + ' ' + code
    print(s)
    os.system(s)
    # parse
    print("parsing")
    data_f = open('dio.txt', 'r')
    data = [ l.strip() for l in data_f ] 
    data_dict = parse(data)
    return data_dict

def parse_args():
    parser = argparse.ArgumentParser(description = 'datamine.py')
    parser.add_argument('-df', dest = 'data_file', default = "./tm_training.drct")
    parser.add_argument('-code', dest = 'code', default = "null")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    data_file = args.data_file
    code = args.code

    if (isfile(data_file) == False):
        print(data_file + " not found.")

    data_dict = datamine(data_file, code = code)
    # polar_probability distributions:
    for polar_number in range(1, 5):
        ppd = get_tm_polar_res(data_dict, 30, polar_number)
        print(ppd)
        view_prob_line_plot(ppd, str(polar_number), [str(i+1) for i in range(0,30)], prediction_name = "TMH_PPD", path = './')
    raise

    # TMH periodicity
    print("getting periods")
    periods = get_period(data_dict, "TMSEQ", "H")
    pickle.dump(periods, open('./TMH_periods.pickle', 'wb'))
    # visualize
    print("visualize")
    view_np(periods, "TMSEQ H periodicity", './')
