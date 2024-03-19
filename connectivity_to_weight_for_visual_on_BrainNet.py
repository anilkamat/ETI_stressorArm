from email import header
import pandas as pd
import numpy as np
import os


""" The code has been used to find the contribution of each connectivity to the accuracy of kSVM model
.It converts the cummulative accuracy into individual contribution and writes them to the csv and .edge file. 
.edge file is used in BrainNet for visualization of connectivity onto brain. Just the file path can be changed for
different experiment/comparision"""

config = {
    'days':['D1','D3'],#'D3'],
    'trials':'T1_3'
}
connection_GC = [[' ','LPFC-->RPFC','LPFC-->LPMC','LPFC-->RPMC','LPFC-->SMA'], 
    ['RPFC-->LPFC',' ','RPFC-->LPMC','RPFC-->RPMC','RPFC-->SMA'],
    ['LPMC-->LPFC','LPMC-->RPFC',' ','LPMC-->RPMC','LPMC-->SMA'],
    ['RPMC-->LPFC','RPMC-->RPFC','RPMC-->LPMC',' ','RPMC-->SMA'],
    ['SMA-->LPFC','SMA-->RPFC','SMA-->LPMC','SMA-->RPMC',' ']]
connection_GC = connection_GC
connection_GC = list(map(list, zip(*connection_GC)))
# import .txt file 
filename = 'top_connectivity_and_cumulative_accuracies.txt'
savefilename = 'connectivity_matrix_for_BrainNetVisual.edge'
csvsavefilename = 'connectivity_matrix_for_BrainNetVisual.csv'
# cumulative accuracy of each of the connectivity calculated from kSVM with feature selection.
num_region = 5
#for day in config['days']:
for i in range(1):
    # connectivity_acc = pd.read_csv(os.path.join(r'C:\Users\kamata3\Work\Brain\Brain_network\Stressor_arm\Results\SVM\normalETIvsSA',day+config['trials']+'_SA',filename), sep=' ', header=None)
    # savepath = os.path.join(r'C:\Users\kamata3\Work\Brain\Brain_network\Stressor_arm\Results\SVM\normalETIvsSA',day+config['trials']+'_SA',savefilename)
    # connectivity_acc = pd.read_csv(os.path.join(r'C:\Users\kamata3\Work\Brain\Brain_network\Stressor_arm\Results\SVM\ETI_succ_vs_unsucc',day+config['trials']+'_MA',filename), sep=' ', header=None)
    # savepath = os.path.join(r'C:\Users\kamata3\Work\Brain\Brain_network\Stressor_arm\Results\SVM\ETI_succ_vs_unsucc',day+config['trials']+'_MA',savefilename)
    # savepath2 = os.path.join(r'C:\Users\kamata3\Work\Brain\Brain_network\Stressor_arm\Results\SVM\ETI_succ_vs_unsucc',day+config['trials']+'_MA',csvsavefilename)
    
    # connectivity_acc = pd.read_csv(os.path.join(r'C:\Users\kamata3\Work\Brain\Brain_network\Stressor_arm\Results\SVM\SVM_from_computer_34\normalETI_vs_reten_succ','D4T1_3_Succ_vsNormal'+day+config['trials'],filename), sep=' ', header=None)
    # savepath = os.path.join(r'C:\Users\kamata3\Work\Brain\Brain_network\Stressor_arm\Results\SVM\SVM_from_computer_34\normalETI_vs_reten_succ','D4T1_3_Succ_vsNormal'+day+config['trials'],savefilename)
    # savepath2 = os.path.join(r'C:\Users\kamata3\Work\Brain\Brain_network\Stressor_arm\Results\SVM\SVM_from_computer_34\normalETI_vs_reten_succ','D4T1_3_Succ_vsNormal'+day+config['trials'],csvsavefilename)
    
    # connectivity_acc = pd.read_csv(os.path.join(r'C:\Users\kamata3\Work\Brain\Brain_network\Stressor_arm\Results\SVM\SVM_from_computer_34\normalReten_vs_stressorReten_succ', filename), sep=' ', header=None)
    # savepath = os.path.join(r'C:\Users\kamata3\Work\Brain\Brain_network\Stressor_arm\Results\SVM\SVM_from_computer_34\normalReten_vs_stressorReten_succ',savefilename)
    # savepath2 = os.path.join(r'C:\Users\kamata3\Work\Brain\Brain_network\Stressor_arm\Results\SVM\SVM_from_computer_34\normalReten_vs_stressorReten_succ',csvsavefilename)

    connectivity_acc = pd.read_csv(os.path.join(r'C:\Users\kamata3\Work\Brain\Brain_network\Stressor_arm\Results\SVM\normalETIvsSA\D1T1_5_Succ', filename), sep=' ', header=None)
    savepath = os.path.join(r'C:\Users\kamata3\Work\Brain\Brain_network\Stressor_arm\Results\SVM\normalETIvsSA\D1T1_5_Succ',savefilename)
    savepath2 = os.path.join(r'C:\Users\kamata3\Work\Brain\Brain_network\Stressor_arm\Results\SVM\normalETIvsSA\D1T1_5_Succ',csvsavefilename)
   
    num_conn = connectivity_acc.shape[0]
    for i in range(num_conn-1,0,-1): 
        connectivity_acc.iloc[i,1] = connectivity_acc.iloc[i,1]-connectivity_acc.iloc[i-1,1]

    temp_mat = np.zeros((num_region,num_region))
    for i in range(num_conn):
        for j in range(num_region):
            for k in range(num_region):
                # print(connectivity_acc.iloc[i,0])
                # print(connection_GC[j][k])
                if connectivity_acc.iloc[i,0] == connection_GC[j][k]:
                    temp_mat[j,k] = connectivity_acc.iloc[i,1]
    print('temp_mat: ',temp_mat)
    temp_mat =  pd.DataFrame(temp_mat)

    temp_mat.to_csv(savepath, header=None,sep = ' ', index=None)

    #temp_mat =  pd.DataFrame(temp_mat)
    connectivity_acc.to_csv(savepath2, header=None,sep = ' ', index=None)
    #write to a txt file with .edge extension.




