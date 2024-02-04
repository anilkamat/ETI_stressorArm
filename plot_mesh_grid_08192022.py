import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from matplotlib.pyplot import figure

def plot_3d_mesh(z,ModelName,m,config2,best_model_, connections):
    model_para_grid = config2['model_para_grid']
    x = model_para_grid['c']
    y =  model_para_grid['gamma']
    kernel = best_model_[ModelName+'kernel']
    c = best_model_[ModelName+'c']
    gamma = best_model_[ModelName+'gamma']
    #print(f'kernel: {kernel}, c : {x}, gamma: {y}')
    figure(figsize=(10, 8), dpi=80)
    ax =sns.heatmap(z, annot=True,cmap='BrBG')#cmap='BrBG'
    ax.set_xticklabels(x, rotation=45,horizontalalignment='right')
    ax.set_yticklabels(y, rotation=45)
    plt.xlabel('C')
    plt.ylabel('gamma')
    if m == 'all_C':
        ttl= ModelName,',c',c,',gamma',gamma,',kernel',kernel #,',best_test_acc',str('{:.3f}'.format(config2['acc_test'])
        figName = ModelName+'all_C'
    elif m =='single_C_drop':
        ttl= ModelName,'connection#',connections,',c',c,',gamma',gamma,',kernel',kernel    #,',best_test_acc',str('{:.3f}'.format(config2['acc_test'])
        figName = ModelName+'conn#'+str(connections)
    else:
        ttl= ModelName,'set#',m,',c',c,',gamma',gamma,',kernel',kernel  #,',best_test_acc',str('{:.3f}'.format(config2['acc_test'])
        figName = ModelName+'set#'+str(m)
    plt.title(ttl)
    plt.savefig(os.path.join(r'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Results\SVM\normalETI_succ_Vs_unsucss_Day3_T1_5',figName))
    #plt.show()
    plt.close()