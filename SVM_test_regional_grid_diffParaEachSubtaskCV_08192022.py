from email import header
from math import gamma
from operator import index
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import time
import os
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score, plot_confusion_matrix
from sklearn.svm import SVC
# from Rank_variables_SVM_accuracy import rankVariable
from grid_search_parameters_subtask_08192022 import grid_search
from plot_mesh_grid_08192022 import plot_3d_mesh
import warnings
warnings.filterwarnings('ignore')
s_time = time.time()
config={
  'PATH' :r'..\Results\SVM\normalETI_succ_Vs_unsucss_Day3_T1_5',  
  #'PATH2' :r'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\EEG_fNIRS_paper_Brain_informatics\fNIRS_codes_results\Results\Results_Jun2_region_level\SVM_results',
  'trials':['T1_5'],
  'Days' : ['D_2'],#,'Tr_2','Tr_3']
    'model_para_grid':{'c' :np.arange(0.1,10,0.5),  # 'model_para_grid':{'c' :np.arange(0.1,10,0.5),# if change value, change in the grid search fxn as well
    'gamma':np.arange(2,10,0.5),   #'gamma':np.arange(0.1,10,0.5), 
    'kernel_' : ['rbf']}
}

# """ READ THE  .MAT FILE FROM THE MATLAB
# import scipy.io
# # mat = scipy.io.loadmat(r'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Results\SA_normal.mat')
# mat2 = mat['SA_normal']['Day_1'].item() # stores days
# trials = mat2[0][0] # stores trials
# # trials[9] # get the dataset
# TR = trials[9].item() # get the dataset
# print('TR',TR[0][:].shape)"""

filenames = []
Data_all = []           # all subjects, all trial, all HTA
for D in config['Days']:
    for Tr in config['trials']:
        #file_ = 'SA_normal_'+str(D)+str(Tr)+'.csv'
        #filename = os.path.join(config['PATH'], file_)
        #filename = r'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Datasets\Mannequin_Day1_alltrials.csv'
        filename = r'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Results\connectivities\ETI_normal\fNIRS_GC_D3T1_5_succ_vs_unsucc.csv'

        file_ = 'fNIRS_GC_D3T1_5_succ_vs_unsucc.csv'
        df = pd.read_csv(filename, header=0, index_col=False,sep=',')
        df = df[(df!=0).all(1)]
        df = df.sample(frac=1) # shuffle the rows
        connections = df.columns[0:20]
        
        config_trial={'filename': file_,     # config of current subtasks
        'Day':D,
        'trial':Tr,
        'connections': connections}
        print('Analysis of fullconnectivity')
        print('shape of df:',df.shape)
        para_model_allConnectivity,best_model_allConnectivity =grid_search(df,config_trial,config,Type='all_connectivity')  # grid search for best parameters 
        #pickle.dump(para_model_allConnectivity, open('para_model_allC'+'.sav', 'wb'))
        #pickle.dump(best_model_allConnectivity, open('best_model_allC'+'.sav', 'wb'))
        ModelName = 'all_connectivity_'+D+Tr
        acc_allC = best_model_allConnectivity[ModelName+'test_acc_CV']
        m= 'all_C'
        plot_3d_mesh(para_model_allConnectivity[str(ModelName)], str(ModelName),m,config,best_model_allConnectivity, connections)
        print('baseline accuracy: ',acc_allC)
        df_write = pd.DataFrame({'acc':acc_allC},index=['1']) #,index=False,header=0
        set_number = 0
        fname = D+Tr+'.xlsx'
        xls_fileName = os.path.join(config['PATH'],fname)
        df_write.to_excel(xls_fileName,sheet_name='set_'+str(set_number))

        print('drop_connectivity analysis')
        dict = {}        
        ModelName = 'drop_connectivity_'+D+Tr
        m = 'single_C_drop'
        i = 0
        for connection in connections:
            data = df[[connection,'Class']]
            para_model_,best_model_ =grid_search(data,config_trial,config,Type='drop_connectivity')  # grid search for best parameters 
            dict[connection] = best_model_[ModelName+'test_acc_CV']
            plot_3d_mesh(para_model_[str(ModelName)], str(ModelName),m,config,best_model_, i)
            i+=1
        sort_orders = sorted(dict.items(), key=lambda x: x[1], reverse=True)
        selected_conn = sort_orders[0:1]

        def refine_connectivity(set1,set3,df,set1_top_acc,set_number): # refines the discering connectivities.
            dict2={}
            for connection in set3:
                set1 = set1+[connection]  # add a new connection from set3 to set1
                #print('set1 initial:',set1)
                data = df[set1+['Class']]
                para_model_,best_model_ =grid_search(data,config_trial,config,Type='drop_connectivity')  # grid search for best parameters 
                acc_set1 = best_model_[ModelName+'test_acc_CV']
                if acc_set1 > set1_top_acc:  # if acc increased
                    #print('accuracy increased! keep searching')
                    set_number += 1
                    set1_top_acc = acc_set1  # increased accuracy
                    dict2[tuple(set1)] = best_model_[ModelName+'test_acc_CV']
                    #print(dict2)
                    df_write = pd.DataFrame({'set':set1,'acc':set1_top_acc})
                    with pd.ExcelWriter(xls_fileName, mode='a') as writer:  
                        df_write.to_excel(writer,sheet_name='set_'+str(set_number))
                    plot_3d_mesh(para_model_[str(ModelName)], str(ModelName),set_number,config,best_model_,connections)
                    #print('testCV_acc: ',best_model_[ModelName+'test_acc_CV'])
                    #print('connectivity set1: ',set1)
                    set3 = set3.drop(connection)
                    break
                else: 
                    #print(f'accuracy decreased! for {connection} stop searching')
                    set3 = set3.drop(connection)
                    set1.remove(connection) # remove because it was added above to construct the dataset for grid_search. 
                    #print(f'set1.type {type(set1)}, set1: {set1}')
                    #print(f'len(set3): {len(set3)} connectivity set3: {set3}')
                    break
            print(f'set1: {set1}, set2: {set3}')
            return set1, set3,set1_top_acc,dict2,set_number

        set1 = [selected_conn[0][0]] # always selects the top performing connectivity
        set_number =1
        df_write = pd.DataFrame(selected_conn)
        with pd.ExcelWriter(xls_fileName, mode='a') as writer:  
            df_write.to_excel(writer,sheet_name='set_'+str(set_number))
        set3 = connections.drop(set1)  # set after dropping the set1 connetivity
        set1_top_acc = selected_conn[0][1]
        print('set#1: ',selected_conn)
        flag = 1 # continue searching if flage = 1
        dict2 = selected_conn
        while(flag):
            # set2 = connections
            if len(set3) >=1:
                set1,set3,set1_top_acc,dict2,set_number= refine_connectivity(set1,set3,df,set1_top_acc,set_number)
            else: 
                flag = 0
        # df3 = pd.DataFrame(dict2)
        # df3 = df3.T
        # FullFilename2 = os.path.join(config['PATH2'],str(trial),str(ST)+'.csv')
        # df3.to_csv(FullFilename2)

        print(f'final set1: {set1}, acc: {set1_top_acc} and set3: {set3}')
print(f'Execution time :,{(time.time()-s_time)/60} minutes')

def plots(acc, connection, plt_title):
    n = np.arange(len(connection))
    fig = plt.figure(figsize = (8, 4))
    # creating the bar plot
    plt.bar( connection, acc, color ='maroon', width = 0.2)
    plt.xlabel("Leave one variables out")
    plt.ylabel("SVM Accuracy")
    plt.title(plt_title)
    ext2 ='.png'
    plt_title = plt_title +ext2
    plt.xticks(rotation=45)
    for i in range(len(acc)):
        plt.annotate(str(acc[i]), xy=(n[i],acc[i]), ha='center', va='bottom')
    #dd_value_label(acc)
    plt.tight_layout()
    #plt.savefig(plt_title)
    #files.download(plt_title) 
    plt.show()
    plt.close('all')