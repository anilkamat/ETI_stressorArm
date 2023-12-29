% ###### setup: import all the directories and packages ###########
clc; clear all; close all; warning off;
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

nG = 1; %3       % number of groups
nD = 4; %3;      % total no.of days, 3 triaining days, 4th day is retention and transfer
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
single_sub_name = 'S20'; % Name of a single subject for GC analysis
tic
for d = 1:nD                    % loop over days
    day = sprintf('Day_%d',d);
    for t = 1:nT                % loop over Trials
        subject = sprintf('_D%02d_%02d',d,t);
        trial = sprintf('Tr_%d',t);
        W = sprintf('Trail_%d',t);
        s = 1;  % counter for successfual trials
        us = 1; % counter for Unsuccessfual trials
        as = 1; % counter for all subjects
        for f= 1:nF
            file = cell2mat(fileList(f));
            if strcmp(file(4:10),subject)==1  && strcmp(file(1:3),single_sub_name)
                Index = find(contains(stressor_performance(:,4),file(1:end-4)));
                %fprintf('Index %d, %s, %s',Index, file(1:end-4));
                class = cell2mat(stressor_performance(Index,12));
                %fprintf('F %s:::S %s :::C %d \n',file, subject, class);
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
                ETI_StressorArm.(day).(trial).AllCohort(as,:)  = [reshape(GC_temp,1,sz) class];   % SA_(stressor arm) -> sub x connection x window
                as = as+1; 
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
fn = fieldnames(ETI_StressorArm);
p = 1;
exp_count = [];     % trigger for training
exp_count2 = [];    % trigger for retention and transfer
Tt = 0;             % total number of trial counter
for i = 1:size(fieldnames(ETI_StressorArm),1)  % loop over days
    d = cell2str(fn(i));
    if strcmp(fn(p),d)==1  % check if the field (here day1,...) exists, for some sub it doesn't exists
        fn(p)
        p = p+1;
        fnn = fieldnames(ETI_StressorArm.(d));
        q = 1;  % counter
        for j =1:size(fieldnames(ETI_StressorArm.(d)))  %loop over trials
            trial_temp = sprintf('Tr_%d',j);
            if strcmp(fnn(q),trial_temp)==1
                q = q+1;
                if isfield(ETI_StressorArm.(d).(trial_temp),'Succ')== 1
                    ETI_StressorArm.(d).(trial_temp).Succ(:,[1 7 13 19 25],:) = [];       % Remove all the diagonals
                    ETI_StressorArm.(d).(trial_temp).AllCohort(:,[1 7 13 19 25],:) = [];
                    Tt = Tt+1;
                elseif isfield(ETI_StressorArm.(d).(trial_temp),'Unsucc')== 1
                    ETI_StressorArm.(d).(trial_temp).Unsucc(:,[1 7 13 19 25],:) = [];       % Remove all the diagonals
                    ETI_StressorArm.(d).(trial_temp).AllCohort(:,[1 7 13 19 25],:) = [];
                    Tt = Tt+1;
                end
                if strcmp(d,'Day_4') && strcmp(trial_temp, 'Tr_3')==1
                    exp_count2(1) = Tt;
                end
            end
        end
    end
    exp_count(i) = Tt;
end
exC = [exp_count, exp_count2];  % combine the experiment trigger. 
exC = sort(exC);
save('..\Results\GC_indiviual_subjects\Connectivity_stressor\S20sETI_StressorArm.mat','ETI_StressorArm');  % saves the connectivity to the local current directory

% for k=1:20 % loop over connectivities
%     m = 1;
%     for i = 1:size(fieldnames(ETI_StressorArm),1) 
%         d = sprintf('Day_%d',i);
%         for j =1:size(fieldnames(ETI_StressorArm.(d)))
%             t = sprintf('Tr_%d',j);
%             GC_indiv(k,m)= ETI_StressorArm.(d).(t).AllCohort(:,[k],:);
%             label_indiv(m) = ETI_StressorArm.(d).(t).AllCohort(:,end);
%             m = m+1;
%         end
%     end
% end

for k=1:20 % loop over connectivities
    m = 1;
    p = 1;
    for i = 1:size(fieldnames(ETI_StressorArm),1) % loop over days
        d = cell2str(fn(i));
        if strcmp(fn(p),d)==1  % check if the field (day1,...) exists, for some sub it doesn't exists
            p = p+1;
            fnn = fieldnames(ETI_StressorArm.(d));
            q = 1;  % counter
            for j =1:size(fieldnames(ETI_StressorArm.(d)))
                t = sprintf('Tr_%d',j);
                if strcmp(fnn(q),t)==1 % check if the trial exists
                    q = q+1;
                    GC_indiv(k,m)= ETI_StressorArm.(d).(t).AllCohort(:,[k],:);
                    label_indiv(m) = ETI_StressorArm.(d).(t).AllCohort(:,end);
                    m = m+1;
                end
            end
        end
    end
end

clc; close all;
% LOAD THE GC FROM TRAINING TRIALS
fpath = '..\Results\\GC_indiviual_subjects\Plots_ETI_stressor';
path_dir = sprintf('C:\\Users\\_Kamat_\\Desktop\\RPI\\ResearchWork\\Papers_\\Effective_Connectivity\\Stressor_arm\\Results\\GC_indiviual_subjects\\Connectivity_stressor');

original_files = dir([path_dir, '/*.mat']);
for k=1:length(original_files)
    filename = [path_dir '/' original_files(k).name];
    ETI_GC(k) = load(filename);
end
%%
subjID = {'S01','S02','S03','S04','S05','S06','S07',...
          'S08','S09','S10','S12','S13','S14','S15',...
          'S16','S17','S18','S19','S20'};
for s = 19:19
    filename = [path_dir '/' original_files(s).name];
    fprintf(' loading %s',original_files(s).name)
    load(filename)
    for k=1:20 % loop over connectivities
        m = 1;
        fn = fieldnames(ETI_GC(s).ETI_StressorArm);
        p = 1;
        for i = 1:size(fn,1) % loop over days
            d = cell2str(fn(i));
            if strcmp(fn(p),d)==1  % check if the field (day1,...) exists, for some sub it doesn't exists
                p = p+1;
                fnn = fieldnames(ETI_GC(s).ETI_StressorArm.(d));
                q = 1;  % counter
                for j =1:size(fnn)  % loop over trials
                    t = sprintf('Tr_%d',j);
                    if strcmp(fnn(q),t)==1 % check if the trial exists
                        q = q+1;
                        fprintf('%s %s \n',d,t);
                        GC_indiv(k,m)= ETI_GC(s).ETI_StressorArm.(d).(t).AllCohort(:,[k],:);
                        label_indiv(m) = ETI_GC(s).ETI_StressorArm.(d).(t).AllCohort(:,end);
                        m = m+1;
                    end
                end
            end
        end
    end
    X_trial = 1:size(GC_indiv,2);
    for k=1:20 % loop over connectivities
        close all;
        GC_red = GC_indiv(k,:);
        GC_red(label_indiv>0) = NaN; 
        GC_green = GC_indiv(k,:);
        GC_green(label_indiv<0) = NaN;
%            d1 = movmean(GC_indiv(k,1:10),5); % Moving mean of first day
        d1 = mean(GC_indiv(k,1:exC(1)));
        d2 = mean(GC_indiv(k,exC(1)+1:exC(2)));
        d3 = mean(GC_indiv(k,exC(2)+1:exC(3)));
        if numel(exC)>3
            d4 = mean(GC_indiv(k,exC(3)+1:exC(4)));
            d5 = mean(GC_indiv(k,exC(4)+1:exC(5)));
        end
        f = figure();
        plot(X_trial,GC_red,'r','linewidth',2); hold on
        plot(X_trial,GC_green,'g','linewidth',2); hold on;
        
        plot([1,exC(1)],[d1,d1],'--black', 'linewidth',2.25); hold on;
        plot([exC(1),exC(2)],[d2,d2],'--black', 'linewidth',2.25); hold on;
        plot([exC(2),exC(3)],[d3,d3],'--black', 'linewidth',2.25); hold on;
        if numel(exC) >3
            plot([exC(3),exC(4)],[d4,d4],'--black', 'linewidth',2.25); hold on;
            plot([exC(4),exC(5)],[d5,d5],'--black', 'linewidth',2.25); hold on;
            vline2 ([exC(1), exC(2), exC(3), exC(4)])
            text([5, exC(1)+3, exC(2)+3, exC(3), exC(4)],2*[1, 1, 1, 1, 1], {'Day1','Day2','Day3','Ret','Tra'}) 
        else        
            vline2 ([exC(1), exC(2)])
            text([5, exC(1)+3, exC(2)+3],2*[1, 1, 1], {'Day1','Day2','Day3'}) % NO RETEN and TRANS
        end
        xlabel('Trial')
        ylabel('Granger Causality')
        ylim([-0.01 3])
        legend('Unsuccessful','Successful','mean', 'location','northwest')
        ttl = sprintf('Stressor ETI subject %s :%s ',original_files(s).name(1:3), connection_GC{k});
        title(ttl);
        if numel(exC) < 4
            f.Position = [200 200 650 450];    % sets plot size     
        else
            f.Position = [200 200 1050 600];    % sets plot size     
        end
            
        baseFileName = sprintf('Stressor_%s_conn%d.png',original_files(s).name(1:3), k);
        fullFileName = fullfile(fpath, baseFileName);
        saveas(f, fullFileName)
    end
end
%% 
% %%  ############ Box plot #############
% fpath = '..\Results\\GC_indiviual_subjects\Plots_ETI_stressor';
% % plot for the successful subjects
% for i = 1:size(fieldnames(ETI_StressorArm),1)
%     d = sprintf('Day_%d',i);
%     for j =1:size(fieldnames(ETI_StressorArm.(d)))
%         t = sprintf('Tr_%d',j);
%         close all;
%         f = figure('visible','off');
%         boxplot(ETI_StressorArm.(d).(t).Succ(:,1:end-1),'Labels',connection_GC, 'Colors','b','jitter',0, 'symbol', ''); % 'whisker', inf,  'PlotStyle','Compact'
%         set(gca,'FontSize',10,'XTickLabelRotation',45)
%         xlabel('Connection')
%         ylabel('Granger Causality')
%         ylim([-0.01 3.5])
%         ttl = sprintf('Stressor arm successful: Day-%d, Trial-%d',i,j);
%         title(ttl);
%         baseFileName = sprintf('stressorETI_succ_D%d_T%d.png',i,j);
%         fullFileName = fullfile(fpath, baseFileName);
%         %saveas(f,fullFileName)
%     end
% end
% % plot for the successful subjects
% for i = 1:size(fieldnames(ETI_StressorArm),1)
%     d = sprintf('Day_%d',i);
%     for j =1:size(fieldnames(ETI_StressorArm.(d)))
%         t = sprintf('Tr_%d',j);
%         close all;
%         f = figure('visible','off');
%         boxplot(ETI_StressorArm.(d).(t).Unsucc(:,1:end-1),'Labels',connection_GC, 'Colors','b','jitter',0, 'symbol', ''); % 'whisker', inf,  'PlotStyle','Compact'
%         set(gca,'FontSize',10,'XTickLabelRotation',45)
%         xlabel('Connection')
%         ylabel('Granger Causality')
%         ylim([-0.01 3.5])
%         ttl = sprintf('Stressor arm unsuccessful: Day-%d, Trial-%d',i,j);
%         title(ttl);
%         baseFileName = sprintf('stressorETI_unsucc_D%d_T%d.png',i,j);
%         fullFileName = fullfile(fpath, baseFileName);
%         %saveas(f,fullFileName)
%     end
% end
% %% ############# Plot the means of all the trails ##############
% % plot the boxplot successful
% q = 1; % grand mean of all the trials
% for con= 1:numel(connection_GC)  % loop over connectivity
%     close all
%     p =1; % counter of total number of trials in entire experiment
%     for i=1:size(fieldnames(ETI_StressorArm),1)
%         d = sprintf('Day_%d',i);
%         for j=1:size(fieldnames(ETI_StressorArm.(d)),1)
%             t = sprintf('Tr_%d',j);
%             M(p) = mean(ETI_StressorArm.(d).(t).Succ(:,con));
%             STD(p) = std(ETI_StressorArm.(d).(t).Succ(:,con));
%             GM(q) = M(p);
%             p = p+1;
%             q = q+1;
%         end
%     end
%     f = figure('visible','off');
%     errorbar(M,STD,"-s","MarkerSize",10,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90])
%     ylim([-0.15 2.5])
%     title('Connectivity over trials: successful')
%     legend(connection_GC(con))
%     xlabel('Trial')
%     ylabel('Connectivity')
%     vline2([10, 20])
%     text([5 15 25], 2*[1 1 1], {'Day 1', 'Day 2', 'Day 3'})
%     baseFileName = sprintf('Conn%d_successful.png',con);
%     fullFileName = fullfile(fpath, baseFileName);
%     %saveas(f,fullFileName)
% end
% % plot the boxplot unsuccessful
% q = 1; % grand mean of all the trials
% for con= 1:numel(connection_GC)  % loop over connectivity
%     close all
%     p =1; % counter of total number of trials in entire experiment
%     for i=1:size(fieldnames(ETI_StressorArm),1)
%         d = sprintf('Day_%d',i);
%         for j=1:size(fieldnames(ETI_StressorArm.(d)),1)
%             t = sprintf('Tr_%d',j);
%             M(p) = mean(ETI_StressorArm.(d).(t).Unsucc(:,con));
%             STD(p) = std(ETI_StressorArm.(d).(t).Unsucc(:,con));
%             GM(q) = M(p);
%             p = p+1;
%             q = q+1;
%         end
%     end
%     f = figure('visible','off');
%     errorbar(M,STD,"-s","MarkerSize",10,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90])
%     ylim([-0.15 2.5])
%     title('Connectivity over trials: unsuccessful')
%     legend(connection_GC(con))
%     xlabel('Trial')
%     ylabel('Connectivity')
%     vline2([10, 20])
%     text([5 15 25], 2*[1 1 1], {'Day 1', 'Day 2', 'Day 3'})
%     baseFileName = sprintf('Conn%d_unsuccessful.png',con);
%     fullFileName = fullfile(fpath, baseFileName);
%     %saveas(f,fullFileName)
% end
% % plot histogram of all the connectivities of all the trials
% %figure(6)
% f = figure(); %('visible','off');
% histogram(GM,30)
% title('Distribution: connectivity')
% xlabel('Granger causality')
% ylabel('samples')
% 
% %%
% for d = 1:3  % loop over days
%     D = sprintf('Day_%d',d);
%     for tr =1:10  % loop ove trials
%         Tr = sprintf('Tr_%d',tr);
%         fileName_ = sprintf('C:\\Users\\_Kamat_\\Desktop\\RPI\\ResearchWork\\Papers_\\Effective_Connectivity\\Stressor_arm\\Results\\connectivities\\ETI_stressor\\fNIRS_GC_D%dT%d.csv',d,tr);
%         %temp = [B_Exp.(Tr).Exp(:,:,i) ones(size(B_Exp.(Tr).Exp(:,:,i),1),1); B_Nov.(Tr).Nov(:,:,i) -ones(size(B_Nov.(Tr).Nov(:,:,i),1),1)];
%         %T = array2table(temp);
%         %temp_WCOH = [WCOH_Exp.(Tr).Exp(:,:,i) ;
%         GC_data = [ETI_StressorArm.(D).(Tr).Succ ; ETI_StressorArm.(D).(Tr).Unsucc];
%         T = array2table(GC_data);
%         T.Properties.VariableNames(1:21) = [connection_GC 'Class'];
%         %writetable(T,fileName_B_exp)
%         %writetable(T,fileName_)
%     end
% end
% fprintf('writing to connectivites done!')
% 
% %%  writing connectivity from T1-T5 
% for d = 1:3  % loop over days
%     GC_data1 = [];
%     GC_data2 = [];
%     GC_data = [];
%     D = sprintf('Day_%d',d);
%     fileName_ = sprintf('C:\\Users\\_Kamat_\\Desktop\\RPI\\ResearchWork\\Papers_\\Effective_Connectivity\\Stressor_arm\\Results\\connectivities\\ETI_stressor\\fNIRS_GC_D%dT1_5_succ_vs_unsucc.csv',d);
%     for tr =1:5  % loop ove trials
%         Tr = sprintf('Tr_%d',tr);
%         GC_data1 = [GC_data1 ; ETI_StressorArm.(D).(Tr).Succ];
%         GC_data2 = [GC_data2 ; ETI_StressorArm.(D).(Tr).Unsucc];
%     end
%     GC_data = [GC_data1; GC_data2];
%     T = array2table(GC_data);
%     T.Properties.VariableNames(1:21) = [connection_GC 'Class'];
%     %writetable(T,fileName_B_exp)
%     %writetable(T,fileName_)
% end
% fprintf('writing to connectivites done!')
