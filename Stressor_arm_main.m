% ###### setup: import all the directories and packages ###########
clc;clear; close all; warning off;
addpath 'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\GCCA'             % toolbox directory
addpath 'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\GCCA\utilities'   % toolbox utility directory
addpath 'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\bsmart'           % toolbox utility directory
addpath(genpath('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Homer3-master\Homer3-master'));   % Need this directory to load the .mat data files
addpath 'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\toolbox_original\mvgc_v1.0\demo' % toolbox directory
addpath (genpath ('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Datasets\fNIRS\Day1\homerOutput')); % dataset directory
addpath 'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Datasets\Trialwise_txts_noSS';
addpath(genpath('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Codes\CStrAinBP_20090913'));
connections = cell(5,5);
Cnames = {'LPFC','RPFC','LPMC','RPMC','SMA'};
for i =1:5
    for j = 1:5
        connections{j,i} = sprintf('%s --> %s ',Cnames{i},Cnames{j});
    end
end
conn_temp = reshape(connections,25,1);
connection_GC = {'LPFC-->RPFC','LPFC-->LPMC','LPFC-->RPMC','LPFC-->SMA',...
    'RPFC-->LPFC','RPFC-->LPMC','RPFC-->RPMC','RPFC-->SMA',...
    'LPMC-->LPFC','LPMC-->RPFC','LPMC-->RPMC','LPMC-->SMA',...
    'RPMC-->LPFC','RPMC-->RPFC','RPMC-->LPMC','RPMC-->SMA',...
    'SMA-->LPFC','SMA-->RPFC','SMA-->LPMC','SMA-->RPMC'};
connection_func = {'LPFC_RPFC','LPFC_LPMC','LPFC_RPMC','LPFC_SMA',...
    'RPFC_LPMC','RPFC_RPMC','RPFC_SMA','LPMC_RPMC','LPMC_SMA','RPMC_SMA'};                                                                                            % phases of HTA
LPFC = []; RPFC = []; LSMA = []; RSMA = []; SMA = []; LPMC = []; RPMC = [];Index_temp = [];

%%
nG = 1; %3       % number of groups
nD = 3; %3;      % total no.of days (3)
nT = 10; %5;     % number of trials
% Data_Reg( any(isnan(Data_Reg),2),:) = [];
freq = [0.01:0.001:0.07];           % Frequency range for GC calculation
r = 1;
sub_trial_name = sprintf('ETI_StressorArm');

MyFolderInfo = dir('..\Datasets\Trialwise_txts_noSS');
fileList = {MyFolderInfo.name}';
fileList = fileList(3:end);
nF = numel(fileList);
cF = 0;  % counts number of subjct matched in the folder
stressor_performance = readcell('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Datasets\StressorArm_Performance_Updated_v2.xlsx');
tic

for d = 1:nD                    % loop over days
    day = sprintf('Day_%d',d);
    for t = 1:nT                % loop over Trials
        subject = sprintf('_D%02d_%02d',d,t);
        trial = sprintf('Tr_%d',t);
        W = sprintf('Trail_%d',t);
        s = 1;
        us = 1;
        for f= 1:nF
            file = cell2mat(fileList(f));
            if strcmp(file(4:10),subject)==1
                Index = find(contains(stressor_performance(:,4),file(1:end-4)));
                class = cell2mat(stressor_performance(Index,12));
                
                fprintf('F %s:::S %s :::C %d \n',file, subject, class);
                cF = cF+1;
                %fprintf('%s found \n',file)
                subj_data = importdata(file); % load signals
                data_Reg_LS = get_LS_channelData(subj_data);  % regional long seperation signals
                
                % Granger causality
                NL = size(data_Reg_LS,2);
                [bic,aic] = cca_find_model_order(data_Reg_LS,2,9);
                mo = min(bic,aic);                                      % selection of model order
                [GW,COH,pp,waut,cons]= cca_pwcausal(data_Reg_LS,1,NL,mo,5,freq, 1);   % in fNIRS the low model order was selected to make the algorithm work
                idx1 = find(freq == 0.01);                                  %mean in the neurophysiology frequency band.
                idx2 = find(freq == 0.07);
                A_temp = [];
                GC_temp = mean(GW(:,:,idx1:idx2),3);
                
                sz = numel(GC_temp);
                %reshape(GC_fqmean.ExpFLS.(sub_trial_name).(trial).(W),1,sz)
                if class ==1 % NOTE: 1: successful and -1: unsuccessful
                    ETI_StressorArm.(day).(trial).Succ(s,:)  = [reshape(GC_temp,1,sz) class];   % SA_(stressor arm) -> sub x connection x window
                    s = s+1;
                elseif class ==-1 
                    %[reshape(GC_temp,1,sz) class]
                    ETI_StressorArm.(day).(trial).Unsucc(us,:)  = [reshape(GC_temp,1,sz) class];   % SA_(stressor arm) -> sub x connection x window
                    us = us+1;
                end
            end
        end
    end
end
toc
fprintf('%d Number of subjects in data folder \n',cF)

for i = 1:size(fieldnames(ETI_StressorArm),1)
    d = sprintf('Day_%d',i);
    for j =1:size(fieldnames(ETI_StressorArm.(d)))
        trial_temp = sprintf('Tr_%d',j);
        ETI_StressorArm.(d).(trial_temp).Succ(:,[1 7 13 19 25],:) = [];       % Remove all the diagonals
        ETI_StressorArm.(d).(trial_temp).Unsucc(:,[1 7 13 19 25],:) = [];       % Remove all the diagonals
        %WCOH_Exp.(trial_temp).Exp(:,[1 7 13 19 25],:) = [];       % Remove all the diagonals
    end
end
save('..\Results\ETI_StressorArm.mat','ETI_StressorArm')  % saves the connectivity to the local current directory
%%  ############ Box plot #############
fpath = '..\Results\Plots_ETIstressor';
% plot for the successful subjects
for i = 1:size(fieldnames(ETI_StressorArm),1)
    d = sprintf('Day_%d',i);
    for j =1:size(fieldnames(ETI_StressorArm.(d)))
        t = sprintf('Tr_%d',j);
        close all;
        f = figure('visible','off');
        boxplot(ETI_StressorArm.(d).(t).Succ(:,1:end-1),'Labels',connection_GC, 'Colors','b','jitter',0, 'symbol', ''); % 'whisker', inf,  'PlotStyle','Compact'
        set(gca,'FontSize',10,'XTickLabelRotation',45)
        xlabel('Connection')
        ylabel('Granger Causality')
        ylim([-0.01 3.5])
        ttl = sprintf('Stressor arm successful: Day-%d, Trial-%d',i,j);
        title(ttl);
        baseFileName = sprintf('stressorETI_succ_D%d_T%d.png',i,j);
        fullFileName = fullfile(fpath, baseFileName);
        saveas(f,fullFileName)
    end
end
% plot for the successful subjects
for i = 1:size(fieldnames(ETI_StressorArm),1)
    d = sprintf('Day_%d',i);
    for j =1:size(fieldnames(ETI_StressorArm.(d)))
        t = sprintf('Tr_%d',j);
        close all;
        f = figure('visible','off');
        boxplot(ETI_StressorArm.(d).(t).Unsucc(:,1:end-1),'Labels',connection_GC, 'Colors','b','jitter',0, 'symbol', ''); % 'whisker', inf,  'PlotStyle','Compact'
        set(gca,'FontSize',10,'XTickLabelRotation',45)
        xlabel('Connection')
        ylabel('Granger Causality')
        ylim([-0.01 3.5])
        ttl = sprintf('Stressor arm unsuccessful: Day-%d, Trial-%d',i,j);
        title(ttl);
        baseFileName = sprintf('stressorETI_unsucc_D%d_T%d.png',i,j);
        fullFileName = fullfile(fpath, baseFileName);
        saveas(f,fullFileName)
    end
end
%% ############# Plot the means of all the trails ##############
% plot the boxplot successful
q = 1; % grand mean of all the trials
for con= 1:numel(connection_GC)  % loop over connectivity
    close all
    p =1; % counter of total number of trials in entire experiment
    for i=1:size(fieldnames(ETI_StressorArm),1)
        d = sprintf('Day_%d',i);
        for j=1:size(fieldnames(ETI_StressorArm.(d)),1)
            t = sprintf('Tr_%d',j);
            M(p) = mean(ETI_StressorArm.(d).(t).Succ(:,con));
            STD(p) = std(ETI_StressorArm.(d).(t).Succ(:,con));
            GM(q) = M(p);
            p = p+1;
            q = q+1;
        end
    end
    f = figure('visible','off');
    errorbar(M,STD,"-s","MarkerSize",10,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90])
    ylim([-0.15 2.5])
    title('Connectivity over trials: successful')
    legend(connection_GC(con))
    xlabel('Trial')
    ylabel('Connectivity')
    vline2([10, 20])
    text([5 15 25], 2*[1 1 1], {'Day 1', 'Day 2', 'Day 3'})
    baseFileName = sprintf('Conn%d_successful.png',con);
    fullFileName = fullfile(fpath, baseFileName);
    saveas(f,fullFileName)
end
% plot the boxplot unsuccessful
q = 1; % grand mean of all the trials
for con= 1:numel(connection_GC)  % loop over connectivity
    close all
    p =1; % counter of total number of trials in entire experiment
    for i=1:size(fieldnames(ETI_StressorArm),1)
        d = sprintf('Day_%d',i);
        for j=1:size(fieldnames(ETI_StressorArm.(d)),1)
            t = sprintf('Tr_%d',j);
            M(p) = mean(ETI_StressorArm.(d).(t).Unsucc(:,con));
            STD(p) = std(ETI_StressorArm.(d).(t).Unsucc(:,con));
            GM(q) = M(p);
            p = p+1;
            q = q+1;
        end
    end
    f = figure('visible','off');
    errorbar(M,STD,"-s","MarkerSize",10,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90])
    ylim([-0.15 2.5])
    title('Connectivity over trials: unsuccessful')
    legend(connection_GC(con))
    xlabel('Trial')
    ylabel('Connectivity')
    vline2([10, 20])
    text([5 15 25], 2*[1 1 1], {'Day 1', 'Day 2', 'Day 3'})
    baseFileName = sprintf('Conn%d_unsuccessful.png',con);
    fullFileName = fullfile(fpath, baseFileName);
    saveas(f,fullFileName)
end
% plot histogram of all the connectivities of all the trials
%figure(6)
f = figure(); %('visible','off');
histogram(GM,30)
title('Distribution: connectivity')
xlabel('Granger causality')
ylabel('samples')

%%
for d = 1:3  % loop over days
    D = sprintf('Day_%d',d);
    for tr =1:10  % loop ove trials
        Tr = sprintf('Tr_%d',tr);
        fileName_ = sprintf('C:\\Users\\_Kamat_\\Desktop\\RPI\\ResearchWork\\Papers_\\Effective_Connectivity\\Stressor_arm\\Results\\connectivities\\ETI_stressor\\fNIRS_GC_D%dT%d.csv',d,tr);
        %temp = [B_Exp.(Tr).Exp(:,:,i) ones(size(B_Exp.(Tr).Exp(:,:,i),1),1); B_Nov.(Tr).Nov(:,:,i) -ones(size(B_Nov.(Tr).Nov(:,:,i),1),1)];
        %T = array2table(temp);
        %temp_WCOH = [WCOH_Exp.(Tr).Exp(:,:,i) ;
        GC_data = [ETI_StressorArm.(D).(Tr).Succ ; ETI_StressorArm.(D).(Tr).Unsucc];
        T = array2table(GC_data);
        T.Properties.VariableNames(1:21) = [connection_GC 'Class'];
        %writetable(T,fileName_B_exp)
        writetable(T,fileName_)
    end
end
fprintf('writing to connectivites done!')

%%  writing connectivity from T1-T5 
for d = 1:3  % loop over days
    GC_data1 = [];
    GC_data2 = [];
    GC_data = [];
    D = sprintf('Day_%d',d);
    fileName_ = sprintf('C:\\Users\\_Kamat_\\Desktop\\RPI\\ResearchWork\\Papers_\\Effective_Connectivity\\Stressor_arm\\Results\\connectivities\\ETI_stressor\\fNIRS_GC_D%dT1_5_succ_vs_unsucc.csv',d);
    for tr =1:5  % loop ove trials
        Tr = sprintf('Tr_%d',tr);
        GC_data1 = [GC_data1 ; ETI_StressorArm.(D).(Tr).Succ];
        GC_data2 = [GC_data2 ; ETI_StressorArm.(D).(Tr).Unsucc];
    end
    GC_data = [GC_data1; GC_data2];
    T = array2table(GC_data);
    T.Properties.VariableNames(1:21) = [connection_GC 'Class'];
    %writetable(T,fileName_B_exp)
    writetable(T,fileName_)
end
fprintf('writing to connectivites done!')
