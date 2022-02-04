# tplot.py
# TLB 8/5/2021

import matplotlib.pyplot as plt
import scipy
# histogram
def view_np(data, name, path, nf = 100, QUIET = True, x = None):
    plt.clf()
    plt.style.use('seaborn-white')
    plt.style.use('seaborn-white')
    bins = nf # CONSTRAINT_STYLES[name]['bins']
    rnge = ( min(data), max(data) )#CONSTRAINT_STYLES[name]['range']
    density = False #CONSTRAINT_STYLES[name]['density']
    density = True
    plt.hist(data, bins, density = density, alpha = 1, range = rnge)
    if x is None:
        plt.xlabel(str(name),fontsize=16,fontweight='bold')
    else:
        plt.xlabel(str(x),fontsize=16,fontweight='bold')
    if density:
        plt.ylabel('% occurence', fontsize=16,fontweight='bold')
    else:
        plt.ylabel('# of times occured')
    plt.title(str(name),fontsize=24,fontweight='bold')
    #print(len(data))
    plt.grid(True)
    print(path + name + '.png')
    plt.savefig(path+name+'.png')
    if not(QUIET):
        plt.show()

# T 10/12/21
def view_TM_ROC(data, name, path, nf, bin_splits):
    plt.clf()
    plt.style.use('seaborn-white')
    bins = nf # CONSTRAINT_STYLES[name]['bins']
    rnge = ( min(data), max(data) )#CONSTRAINT_STYLES[name]['range']
    density = False #CONSTRAINT_STYLES[name]['density']


    plt.hist(data, bins, density = density, alpha = 0.75, range = rnge)
    plt.xlabel(str(name))
    if density:
        plt.ylabel('% occurence')
    else:
        plt.ylabel('# of times occured')
    #if bin_splits is not None:
    #    for x in bin_splits:
    #        plt.axvline(x, color = 'black', linewidth = 1, linestyle = '--')
    plt.title('Histogram of ' + str(name))
    #print(len(data))
    plt.grid(True)
    print(path + name + '.png')
    plt.savefig(path+name+'.png')
    plt.show()

# TLB 10/11/21 TMPRED PLOT
def view_prob_line_plot(data, name, data_c, prediction_name = "TMPRED", path = './HMMSTRTM/TMPRED_'):
    plt.clf()
    x_axis = [ i for i,f in enumerate(data) ] 
    plt.style.use('seaborn-white')
    plt.plot(x_axis,data, color = 'blue') # prediction ( TM Probability )
    plt.xticks(range(len(data_c)), data_c, fontsize = 16) # ground truth characters
    plt.title(prediction_name, fontsize = 24, fontweight = 'bold')
    plt.savefig(path+name+'.png')
    plt.show()

# TLB 10/11/21 TMPRED_ADVANCED PLOT
def view_TMPRED_plot(PTM, PINSIDE, POUTSIDE, path, name, data_c=None):
    plt.clf()
    fig,ax = plt.subplots()
    x_axis = [ i*2 for i,f in enumerate(PTM) ]
    plt.style.use('seaborn-white')
    plt.plot(x_axis,PTM, color = 'red') # prediction ( TM Probability )
    plt.plot(x_axis,PINSIDE, color = 'blue') # prediction ( INSIDE Probability )
    plt.plot(x_axis,POUTSIDE, color = 'pink') # prediction ( OUTSIDE Probability )
    d = scipy.zeros(len(x_axis))
    #ax.fill_between(x_axis, PTM, where=PTM<=d, interpolate = True, color = 'red')
    if data_c is not None:
        plt.xticks(x_axis, data_c, fontsize=5) # ground truth characters
        GTTM = []
        for c in data_c:
            if c == "U":
                GTTM.append(0.5)
            elif c != "1" and c != "2":
                GTTM.append(1)
            else:
                GTTM.append(0)
        plt.plot(x_axis,GTTM, color = 'green') # prediction ( OUTSIDE Probability )
    plt.title(str(name))
    plt.ylabel('Probability')
    plt.savefig(path + 'TMRPRED_'+name+'.png')
    fig = plt.gcf()
    plt.show()
    #fig.show()

# TLB 10/11/21 TMPRED_ADVANCED PLOT
def view_TMPRED_plot_more(PTM_B, PTM_H, PTM_C, PTM_I, PTM_L, PTM_F, P_U, PINSIDE, POUTSIDE, path, name, data_c):
    plt.clf()
    fig,ax = plt.subplots()
    x_axis = [ i*2 for i,f in enumerate(PTM_B) ]
    plt.style.use('seaborn-white')
    plt.plot(x_axis,PTM_B, color = 'teal') # prediction ( B-TM Probability )
    plt.plot(x_axis,PTM_H, color = 'purple') # prediction ( H-TM Probability )
    plt.plot(x_axis,PTM_C, color = 'cyan') # prediction ( C-TM Probability )
    plt.plot(x_axis,PTM_I, color = 'red') # prediction ( I-TM Probability )
    plt.plot(x_axis,PTM_L, color = 'yellow') # prediction ( L-TM Probability )
    plt.plot(x_axis,PTM_F, color = 'orange') # prediction ( F-TM Probability )
    plt.plot(x_axis,P_U, color = 'gold') # prediction ( F-TM Probability )
    plt.plot(x_axis,PINSIDE, color = 'blue') # prediction ( INSIDE Probability )
    plt.plot(x_axis,POUTSIDE, color = 'pink') # prediction ( OUTSIDE Probability )
    GTTM = []
    for c in data_c:
        if c == "U":
            GTTM.append(0.5)
        elif c != "1" and c != "2":
            GTTM.append(1)
        else:
            GTTM.append(0)
    plt.plot(x_axis,GTTM, color = 'green') # prediction ( OUTSIDE Probability )
    d = scipy.zeros(len(x_axis))
    #ax.fill_between(x_axis, PTM, where=PTM<=d, interpolate = True, color = 'red')

    #plt.xticks(x_axis, data_c, fontsize=5) # ground truth characters
    x_axis = [ i*50 for i in range(0,int(len(PTM_B)/50)+1) ]
    plt.xticks([x*2 for x in x_axis],[str(x) for x in x_axis] , fontsize=10) # ground truth characters
    plt.title(str(name))
    plt.ylabel('Probability')
    plt.savefig(path + 'TMRPRED_'+name+'.png')
    fig = plt.gcf()
    plt.show()
    #fig.show()

# TLB 10/11/21 TMPRED_ADVANCED PLOT
def view_RAMAPRED_plot(PHG,PBb,PEed,Pc,Px,PlL, path, name, data_c):
    plt.clf()
    fig,ax = plt.subplots()
    x_axis = [ i*2 for i,f in enumerate(PHG) ]
    plt.style.use('seaborn-white')
    plt.plot(x_axis,PHG, color = 'red') # prediction ( TM Probability )
    plt.plot(x_axis,PBb, color = 'blue') # prediction ( INSIDE Probability )
    plt.plot(x_axis,PEed, color = 'green') # prediction ( INSIDE Probability )
    plt.plot(x_axis,Pc, color = 'yellow') # prediction ( INSIDE Probability )
    plt.plot(x_axis,Px, color = 'black') # prediction ( OUTSIDE Probability )
    plt.plot(x_axis,PlL, color = 'orange') # prediction ( OUTSIDE Probability )
    rama_bin_chars = ["H","B","E","c","x","l"]
    for i in range(0, len(PHG)):
        l = [PHG[i],PBb[i],PEed[i],Pc[i],Px[i],PlL[i]]
        k = l.index(max(l))
        data_c[i] += '\n' + rama_bin_chars[k]
        
    plt.xticks(x_axis, data_c, fontsize=5) # ground truth characters
    plt.title(str(name))
    plt.ylabel('Probability')
    plt.savefig(path + 'RAMARPRED_'+name+'.png')
    fig = plt.gcf()
    plt.show()
    #fig.show()



