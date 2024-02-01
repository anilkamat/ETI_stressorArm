from ast import Mod
import numpy as np
import matplotlib.pyplot as plt
from plot_mesh_grid_08192022 import plot_3d_mesh
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score, plot_confusion_matrix
from sklearn.svm import SVC
from sklearn.model_selection import LeaveOneOut, KFold
import pandas as pd
import pickle
import os

def grid_search(data,config_subtask,config,Type = 'dropped_connectivity'):
    def SVM(X_train,y_train):
        svclassifier.fit(X_train, y_train)
        #y_pred_test = svclassifier.predict(X_test)
        y_pred_train = svclassifier.predict(X_train)
        #print(y_pred)
        #print(confusion_matrix(y_test,y_pred_test))
        #print(classification_report(y_test,y_pred_test))
        #   plot_confusion_matrix(svclassifier, X_test, y_test)
        #   plt.show()
        # on train data
        #print(confusion_matrix(y_train,y_pred_train))
        metrics = classification_report(y_train,y_pred_train)
        #print(metrics)
        acc_train = accuracy_score(y_train,y_pred_train)
        #acc_test = accuracy_score(y_test,y_pred_test)
        #print('The accuracy',acc)
        return acc_train

    def SVM_test(X_test,y_test): # for test datasets
        y_pred_test = svclassifier.predict(X_test)
        #print(y_pred)
        #print(confusion_matrix(y_test,y_pred_test))
        #print(classification_report(y_test,y_pred_test))
        #plot_confusion_matrix(svclassifier, X_test, y_test)  
        acc_test = accuracy_score(y_test,y_pred_test)
        #print('The accuracy',acc)
        return acc_test
    
    # model_para_grid={'c' :np.arange(0.1,10,2),  # dictionary within a dictionary
    #     'gamma':np.arange(0.1,10,2), 
    #     'kernel' : ['linear','rbf']}
    # model_para_grid = {
    # 'c' :[0.001, 0.01, 0.1, 1, 10, 100,1000,10000],
    # 'gamma':[0.001, 0.01, 0.1,1,10,100,1000],
    # 'kernel' : ['linear','rbf']
    # }
    model_para_grid = config['model_para_grid']
    para_model_ = {}
    best_model_ = {}
    test_acc_dropped_ = np.zeros((20))
    column_ = ['Zero']
    Day = config_subtask['Day']
    trial= config_subtask['trial']#[0:1]
    filename = config_subtask['filename']
    connections = config_subtask['connections']
    # print('fnames: ',fnames)
    # print('connections: ',connections)
    # print('trials: ',trials)
    # print('subtasks: ',subtasks)
    if Type == 'all_connectivity':
        ModelName = 'all_connectivity_'+Day+trial
        para_model_[ModelName]=np.zeros(((len(model_para_grid['gamma'])), len(model_para_grid['c'])))
        best_model_[ModelName+'c'] = 0.0
        best_model_[ModelName+'gamma'] = 0.0
        best_model_[ModelName+'kernel'] = ''
        best_model_[ModelName+'test_acc'] = 0.0
        best_model_[ModelName+'test_acc_CV'] =0.0
    else:
        for m in range(20):  # for feature/connectivity dropped datasets.
            ModelName = 'drop_connectivity_'+Day+trial
            #print('modelname:',ModelName)
            para_model_[ModelName]=np.zeros(((len(model_para_grid['gamma'])), len(model_para_grid['c'])))
            best_model_[ModelName+'c'] = 0.0
            best_model_[ModelName+'gamma'] = 0.0
            best_model_[ModelName+'kernel'] = ''
            best_model_[ModelName+'test_acc'] = 0.0
            best_model_[ModelName+'test_acc_CV'] =0.0
    test_acc = []
    for kernel_ in model_para_grid['kernel_']:
        n = 0
        for gamma_ in model_para_grid['gamma']:
            o = 0
            for c_ in model_para_grid['c']:
                if kernel_ == 'linear':
                    svclassifier = SVC(kernel=kernel_,random_state=1)#random_state=1,) # ,
                else:
                    svclassifier = SVC(kernel=kernel_, gamma=gamma_, C=c_,random_state=1)# gamma=100, C=0.1,random_state=1,) # ,
                temp_ = np.zeros((20))

                if Type=='all_connectivity':
                    X = data.iloc[:,:-1]
                    Y = data.iloc[:,-1]
                    mean = X.mean(axis = 0)
                    std = X.std( axis = 0)
                    X = (X-mean)/std
                    #print('normalized Xtrain :', X_train)
                    acc_test_CV = []
                    # loo = LeaveOneOut()
                    # loo.get_n_splits(X)
                    # for train_indx, test_indx in loo.split(X):
                    #     #print(f'train_indx: {train_indx} test_indx: {test_indx}')
                    #     x_train,x_test = X.iloc[train_indx],X.iloc[test_indx]
                    #     y_train,y_test = Y.iloc[train_indx],Y.iloc[test_indx]
                    #     #print(f'x_train: {x_train}')
                    #     acc_train = SVM(x_train,y_train)
                    #     acc_test = SVM_test(x_test,y_test)
                    #     acc_test_CV.append(acc_test)
                    # acc_test_CV_mean = np.mean(acc_test_CV)
                    # kfold
                    kf = KFold(n_splits=10)
                    for train_indx, test_indx in kf.split(X):
                        #print(f'train_indx: {train_indx} test_indx: {test_indx}')
                        x_train,x_test = X.iloc[train_indx],X.iloc[test_indx]
                        y_train,y_test = Y.iloc[train_indx],Y.iloc[test_indx]
                        acc_train = SVM(x_train,y_train)
                        acc_test = SVM_test(x_test,y_test)
                        acc_test_CV.append(acc_test)
                    acc_test_CV_mean = np.mean(acc_test_CV)


                    #acc_test_CV_std = np.std(acc_test_CV)
                    #print(f'Train acc:{acc_train}')
                    ModelName = 'all_connectivity_'+Day+trial
                    #print('model_name:',ModelName)
                    #para_model_T1_[ModelName[:-4]] = np.zeros((20,7,8))   # (connection, gamma, c)
                    if (acc_test_CV_mean > best_model_[ModelName+'test_acc_CV']) :
                        #pickle.dump(svclassifier, open(ModelName+'.sav', 'wb'))     #save the best model so far
                        best_model_[ModelName+'c'] = c_
                        best_model_[ModelName+'gamma'] = gamma_
                        best_model_[ModelName+'kernel']= kernel_
                        best_model_[ModelName+'test_acc'] = acc_test
                        best_model_[ModelName+'test_acc_CV'] =acc_test_CV_mean
                        print('best test acc:', acc_test)
                    # store all the parameters and corrosponding accuracies
                    para_model_[ModelName][n,o] = acc_test_CV_mean
                    test_acc=best_model_[ModelName+'test_acc_CV']

                if Type=='drop_connectivity':
                    X = data.iloc[:,:-1]
                    Y = data.iloc[:,-1]
                    mean = X.mean(axis = 0)
                    std = X.std( axis = 0)
                    X = (X-mean)/std
                    #print('normalized Xtrain :', X_train)
                    loo = LeaveOneOut()
                    loo.get_n_splits(X)
                    acc_test_CV = []
                    for train_indx, test_indx in loo.split(X):
                        x_train,x_test = X.iloc[train_indx],X.iloc[test_indx]
                        y_train,y_test = Y.iloc[train_indx],Y.iloc[test_indx]
                        acc_train = SVM(x_train,y_train)
                        acc_test = SVM_test(x_test,y_test)
                        acc_test_CV.append(acc_test)
                    acc_test_CV_mean = np.mean(acc_test_CV)
                    #acc_test_CV_std = np.std(acc_test_CV)
                    #print(f'Train acc:{acc_train}')
                    ModelName = 'drop_connectivity_'+Day+trial
                    #print('model_name:',ModelName)
                    #para_model_T1_[ModelName[:-4]] = np.zeros((20,7,8))   # (connection, gamma, c)
                    if (acc_test_CV_mean > best_model_[ModelName+'test_acc_CV']) :
                        pickle.dump(svclassifier, open(ModelName+'.sav', 'wb'))     #save the best model so far
                        best_model_[ModelName+'c'] = c_
                        best_model_[ModelName+'gamma'] = gamma_
                        best_model_[ModelName+'kernel']= kernel_
                        best_model_[ModelName+'test_acc'] = acc_test
                        best_model_[ModelName+'test_acc_CV'] =acc_test_CV_mean
                    # store all the parameters and corrosponding accuracies
                    para_model_[ModelName][n,o] = acc_test_CV_mean
                    test_acc=best_model_[ModelName+'test_acc_CV']
                    #print('best Mean_test acc:', test_acc)
                    

                o+=1
            n +=1 

    # if Type == 'dropped_connectivity':
    #     ModelName = 'drop_connectivity'+trial+subtask
    #     config2={
    #         'c': best_model_[ModelName+'c'],
    #         'gamma':best_model_[ModelName+'gamma'],
    #         'kernel_': best_model_[ModelName+'kernel'],
    #         'acc_test':best_model_[ModelName+'test_acc'],
    #         'figType': 'C_dropped'
    #     }
    #     plot_3d_mesh(para_model_[ModelName], str(ModelName),m,config2, connections)
    # else:
    #     ModelName = 'all_connectivity'+trial+subtask
    #     m = 'None'
    #     config2={
    #         'c': best_model_[ModelName+'c'],
    #         'gamma':best_model_[ModelName+'gamma'],
    #         'kernel_': best_model_[ModelName+'kernel'],
    #         'acc_test':best_model_[ModelName+'test_acc'],
    #         'figType': 'C_all'
    #     }
    # plot_3d_mesh(para_model_[str(ModelName)], str(ModelName),m,config2, connections)
    #plot_3d_mesh(para_model_T2_, connections,'second trial')
    #plot_3d_mesh(para_model_T3_, connections,'third trial')
    return para_model_,best_model_

    