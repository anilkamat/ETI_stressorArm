% CUSUM plot of the ETI stressor and non-stressor group. 
%% Installation of BEAST algorithm 
% beastPath = 'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Codes\BEAST\';
% eval( webread('http://b.link/rbeast') ) 

%% ########## % CUSUM plot of the ETI non-stressor group ##########
clc; clear all; close all;

addpath 'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Codes\BEAST\';
fpath = '..\Results\GC_indiviual_subjects\CUSUM_plots_nonStressor';   % CPD plotsaving directory
fpath2 = '..\Results\GC_indiviual_subjects\CUSUM_non_stressor';       % CUSUM plot saving directory
subjID = {'LC01','LC02','LC03','LC04','LC06','LC07','LC08',...
    'LC11','LC12','LC13','LC17','LC18','LC19','LC20','LC21'};
connection_GC = {'LPFC-->RPFC','LPFC-->LPMC','LPFC-->RPMC','LPFC-->SMA',...
    'RPFC-->LPFC','RPFC-->LPMC','RPFC-->RPMC','RPFC-->SMA',...
    'LPMC-->LPFC','LPMC-->RPFC','LPMC-->RPMC','LPMC-->SMA',...
    'RPMC-->LPFC','RPMC-->RPFC','RPMC-->LPMC','RPMC-->SMA',...
    'SMA-->LPFC','SMA-->RPFC','SMA-->LPMC','SMA-->RPMC'};
connection_func = {'LPFC_RPFC','LPFC_LPMC','LPFC_RPMC','LPFC_SMA',...
    'RPFC_LPMC','RPFC_RPMC','RPFC_SMA','LPMC_RPMC','LPMC_SMA','RPMC_SMA'};
for s = 16:19     % Loop over subjects
    subject_ID = subjID{s};
    filename = sprintf('%sETI_normal.mat',subject_ID);
    fullFN = sprintf('C:\\Users\\_Kamat_\\Desktop\\RPI\\ResearchWork\\Papers_\\Effective_Connectivity\\Stressor_arm\\Results\\GC_indiviual_subjects\\Connectivity_non_stressor\\%s',filename);  % Full file name
    load(fullFN)
    for k=1:20 % loop over connectivities
        m = 1;
        for i = 1:size(fieldnames(ETI_normal),1) 
            d = sprintf('Day_%d',i);
            for j =1:size(fieldnames(ETI_normal.(d)))
                t = sprintf('Tr_%d',j);
                GC_indiv(k,m)= ETI_normal.(d).(t).AllCohort(:,[k],:);
                label_indiv(m) = ETI_normal.(d).(t).AllCohort(:,end);
                m = m+1;
            end
        end
        % ## change point detection (based on deviation of mean)
%         close all;
%         gca = figure(2);
%         set(gca,'DefaultLineLineWidth',1.5)
%         findchangepts(GC_indiv(k,:))
%         xlabel('Trials')
%         ylabel('GC')
%         ylim([0,2])
%         ttl = sprintf('Sub: CPD %s conn: %s',subject_ID, connection_GC{k});
%         title(ttl)
%         %sgtitle(ttl)
%         baseFileName = sprintf('CPD2_nonStressor_%s_conn%d.png',subject_ID,  k);
%         fullFileName = fullfile(fpath, baseFileName);
%         saveas(figure(2), fullFileName)
    end
    
    % ###### CUSUM plot  #######
    for k = 1:20 % loop over 20 connectivity
        close all;
        gca = figure(3);
        set(gca,'DefaultLineLineWidth',1.5)
        cusum(GC_indiv(k,:),3,'all')
        xlabel('Trials')
        ylabel('CUSUM score')
        ylim([-6,6])
        ttl = sprintf('Sub:%s_conn:%s',subject_ID, connection_GC{k});
        title(ttl)
        baseFileName = sprintf('CUSUM_nonStressor_%s_conn%d.png',subject_ID,  k);
        fullFileName = fullfile(fpath2, baseFileName);
        saveas(figure(3), fullFileName)
    end
end
%% #########  All subject together  ###########
clc; clear all; close all;
connection_GC = {'LPFC-->RPFC','LPFC-->LPMC','LPFC-->RPMC','LPFC-->SMA',...
    'RPFC-->LPFC','RPFC-->LPMC','RPFC-->RPMC','RPFC-->SMA',...
    'LPMC-->LPFC','LPMC-->RPFC','LPMC-->RPMC','LPMC-->SMA',...
    'RPMC-->LPFC','RPMC-->RPFC','RPMC-->LPMC','RPMC-->SMA',...
    'SMA-->LPFC','SMA-->RPFC','SMA-->LPMC','SMA-->RPMC'};

path_dir = sprintf('C:\\Users\\_Kamat_\\Desktop\\RPI\\ResearchWork\\Papers_\\Effective_Connectivity\\Stressor_arm\\Results\\GC_indiviual_subjects\\Connectivity_non_stressor');
original_files = dir([path_dir, '/*.mat']);
for k=1:length(original_files)
    filename = [path_dir '/' original_files(k).name];
    ETI_GC(k) = load(filename);
end
%%
sub_ID = {original_files([1,2,3,4,7,8,9,11,12,13,20]).name};
for k=1:20 % loop over connectivities
    GC_indiv = zeros(23, 30);
    label_indiv = zeros(23, 30);
    for s = 1:23 % loop over subjects
        if s==5 || s==10 || s==14 || s==15 ||  s==16 || s==17 || s==21 || s==22 || s==23
            continue;
        else
            m = 1;
            for i = 1:size(fieldnames(ETI_GC(s).ETI_normal),1) 
                d = sprintf('Day_%d',i);
                for j =1:size(fieldnames(ETI_GC(s).ETI_normal.(d)))
                    t = sprintf('Tr_%d',j);
                    GC_indiv(s,m)= ETI_GC(s).ETI_normal.(d).(t).AllCohort(:,[k],:);
                    label_indiv(s,m) = ETI_GC(s).ETI_normal.(d).(t).AllCohort(:,end);
                    m = m+1;
                end
            end
        end
    end    
    GC_indiv([5,6,10,14,15,16,17,18,19,21,22,23],:) = [];  % subjects with some issue
    label_indiv([5,6,10,14,15,16,17,18,19,21,22,23],:) = [];    % subjects with some issue
    % GC_indiv(GC_indiv ==0) = [];
    % label_indiv(label_indiv ==0) = [];
    fpath = '..\Results\GC_indiviual_subjects\CUSUM_plots_nonStressor';   % Figure saving directory
    iupper = struct(); ilower = struct(); uppersum = struct(); lowersum = struct();
    subj_ID = struct();
    close all;
    for s = 1:11    % number of subjects after removing problematic ones
        gca = figure(2);
        x0=10;
        y0=10;
        width=950;
        height=600;
        set(gcf,'position',[x0,y0,width,height])
        set(gca,'DefaultLineLineWidth',1.5)
        cusum(GC_indiv(s,:),3,'all'); hold on;
        temp = char(sub_ID(s));
        temp = temp(1:4);
        subj_ID(s).SID = temp;
        %subj_ID(s) = char(sub_ID(s));
        xlabel('Trials')
        ylabel('CUSUM score')
        ylim([-6 6])
    end
    ttl = sprintf(' conn :%s',connection_GC{k});
    title(ttl)
    legend(struct2cell(subj_ID),'Location', 'eastoutside') % 'Orientation','horizontal' )
    baseFileName = sprintf('CUSUM_nonStressor_conn%d.png', k); 
    fullFileName = fullfile(fpath, baseFileName);
    saveas(figure(2), fullFileName)
end
%%
% fpath = '..\Results\GC_indiviual_subjects\CUSUM_plots_nonStressor';   % Figure saving directory
% sub_ID = {original_files([1,2,3,4,7,8,9,11,12,13,20]).name};
% iupper = struct(); ilower = struct(); uppersum = struct(); lowersum = struct();
% subj_ID = struct();
% x = 1:30;
% for k = 1:1 % loop over 20 connectivity
%     close all;
%     for s = 1:11
%         [iupper(s).IU, ilower(s).IL, uppersum(s).US, lowersum(s).LS ] = cusum(GC_indiv(s,:),3,'all');
%         U_sum(s,1:30) = uppersum(s).US;
%         L_sum(s,1:30) = lowersum(s).LS;
%         temp = char(sub_ID(s));
%         temp = temp(1:4);
%         subj_ID(s).SID = temp;
%         %subj_ID(s) = char(sub_ID(s));
%     end
%     gca = figure(2);
%     set(gca,'DefaultLineLineWidth',1.5)
%     plot(U_sum'); hold on;
%     %p = plot(x,y,'o-','MarkerFaceColor','red','MarkerEdgeColor','red','MarkerIndices',10)
%     yline(3,'--')
%     xlabel('Trials')
%     ylabel('CUSUM score')
%     yline(1.5,'--', 'linewidth',1.5)
%     %yline([2 -2],'--',{'UI','LI'})
%     hold on;
%     ttl = sprintf(' conn :%s',connection_GC{k});
%     title(ttl)
%     legend(struct2cell(subj_ID),'Location', 'eastoutside') % 'Orientation','horizontal' )
%     baseFileName = sprintf('nonStressor_conn%d.png', k);
%     fullFileName = fullfile(fpath, baseFileName);
%     %saveas(figure(2), fullFileName)
% end
%% ############## CUSUM plot of the ETI stressor group. ###########
clc; clear all; close all;
subject_ID = 'S12';  %         if s==1 || s==10 || s==11 || s==12 || s==13 || s==14 (problematic subjects)
filename = sprintf('%sETI_StressorArm.mat',subject_ID);
fullFN = sprintf('C:\\Users\\_Kamat_\\Desktop\\RPI\\ResearchWork\\Papers_\\Effective_Connectivity\\Stressor_arm\\Results\\GC_indiviual_subjects\\Connectivity_stressor\\%s',filename);  % Full file name
connection_GC = {'LPFC-->RPFC','LPFC-->LPMC','LPFC-->RPMC','LPFC-->SMA',...
    'RPFC-->LPFC','RPFC-->LPMC','RPFC-->RPMC','RPFC-->SMA',...
    'LPMC-->LPFC','LPMC-->RPFC','LPMC-->RPMC','LPMC-->SMA',...
    'RPMC-->LPFC','RPMC-->RPFC','RPMC-->LPMC','RPMC-->SMA',...
    'SMA-->LPFC','SMA-->RPFC','SMA-->LPMC','SMA-->RPMC'};
connection_func = {'LPFC_RPFC','LPFC_LPMC','LPFC_RPMC','LPFC_SMA',...
    'RPFC_LPMC','RPFC_RPMC','RPFC_SMA','LPMC_RPMC','LPMC_SMA','RPMC_SMA'};
load(fullFN)
for k=1:20 % loop over connectivities
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
    m = 1;
    for i = 1:size(fieldnames(ETI_StressorArm),1) 
        d = sprintf('Day_%d',i);
        for j =1:size(fieldnames(ETI_StressorArm.(d)))
            t = sprintf('Tr_%d',j);
            GC_indiv(s,m)= ETI_ETI_StressorArm.(d).(t).AllCohort(:,[k],:);
            label_indiv(s,m) = ETI_ETI_StressorArm.(d).(t).AllCohort(:,end);
            m = m+1;
        end
    end
end
fpath = '..\Results\GC_indiviual_subjects\CUSUM_plots_Stressor';   % Figure saving directory
for k = 1:20 % loop over 20 connectivity
    close all;
    gca = figure(2);
    set(gca,'DefaultLineLineWidth',1.5)
    cusum(GC_indiv(k,:),3,'all')
    xlabel('Trials')
    ylabel('CUSUM score')
    ylim([-6,6])
    ttl = sprintf('Sub:%s_conn:%s',subject_ID, connection_GC{k});
    title(ttl)
    baseFileName = sprintf('Stressor_%s_conn%d.png',subject_ID,  k);
    fullFileName = fullfile(fpath, baseFileName);
    saveas(figure(2), fullFileName)
end

% for k = 1:20 % loop over 20 connectivity
%     close all;
%     if k==9 || k==9 % some issue with BEAST algo at the forth connectivity 
%         continue;
%     else
%         gca = figure(2);
%         set(gca,'DefaultLineLineWidth',1.5)
%         out = beast(GC_indiv(k,:));
%         plotbeast(out)
%         xlabel('Trials')
%     %     ylim([-6,6])
%         ttl = sprintf('Sub:%s_conn:%s',subject_ID, connection_GC{k});
%         %title(ttl)
%         sgtitle(ttl)
%         baseFileName = sprintf('CPD_Stressor_%s_conn%d.png',subject_ID,  k);
%         fullFileName = fullfile(fpath, baseFileName);
%         saveas(figure(2), fullFileName)
%     end
% end

%% stressor group: all subjects together
clc; clear all; close all;
connection_GC = {'LPFC-->RPFC','LPFC-->LPMC','LPFC-->RPMC','LPFC-->SMA',...
    'RPFC-->LPFC','RPFC-->LPMC','RPFC-->RPMC','RPFC-->SMA',...
    'LPMC-->LPFC','LPMC-->RPFC','LPMC-->RPMC','LPMC-->SMA',...
    'RPMC-->LPFC','RPMC-->RPFC','RPMC-->LPMC','RPMC-->SMA',...
    'SMA-->LPFC','SMA-->RPFC','SMA-->LPMC','SMA-->RPMC'};

path_dir = sprintf('C:\\Users\\_Kamat_\\Desktop\\RPI\\ResearchWork\\Papers_\\Effective_Connectivity\\Stressor_arm\\Results\\GC_indiviual_subjects\\Connectivity_stressor');
original_files = dir([path_dir, '/*.mat']);
for k=1:length(original_files)
    filename = [path_dir '/' original_files(k).name];
    ETI_GC(k) = load(filename);
end

sub_ID = {original_files([1,2,3,4,7,8,9,11,12,13,14,15,16,17,18]).name};
for k=1:20 % loop over connectivities
    GC_indiv = zeros(18, 30);
    label_indiv = zeros(18, 30);
    for s = 1:18 % loop over subjects
        if s==1 || s==10 || s==11 || s==12 || s==13 || s==14
            continue;
        else
            m = 1;
            for i = 1:size(fieldnames(ETI_GC(s).ETI_StressorArm),1) 
                d = sprintf('Day_%d',i);
                for j =1:size(fieldnames(ETI_GC(s).ETI_StressorArm.(d)))
                    t = sprintf('Tr_%d',j);
                    GC_indiv(s,m)= ETI_GC(s).ETI_StressorArm.(d).(t).AllCohort(:,[k],:);
                    label_indiv(s,m) = ETI_GC(s).ETI_StressorArm.(d).(t).AllCohort(:,end);
                    m = m+1;
                end
            end
        end
    end    
    GC_indiv([1,10,11,12,13,14],:) = [];    % subjects with some issue
    label_indiv([1,10,11,12,13,14],:) = []; % subjects with some issue
    % GC_indiv(GC_indiv ==0) = []; 
    % label_indiv(label_indiv ==0) = [];
    fpath = '..\Results\GC_indiviual_subjects\CUSUM_plots_Stressor';   % Figure saving directory
    iupper = struct(); ilower = struct(); uppersum = struct(); lowersum = struct();
    subj_ID = struct();
    close all;
    for s = 1:12   % number of subjects after removing problematic ones
        gca = figure(2);
        x0=10;
        y0=10;
        width=950;
        height=600;
        set(gcf,'position',[x0,y0,width,height])
        set(gca,'DefaultLineLineWidth',1.5)
        cusum(GC_indiv(s,:),3,'all'); hold on;
        temp = char(sub_ID(s));
        temp = temp(1:4);
        subj_ID(s).SID = temp;
        %subj_ID(s) = char(sub_ID(s));
        xlabel('Trials')
        ylabel('CUSUM score')
        ylim([-6 6])
    end
    ttl = sprintf('Stressor conn :%s',connection_GC{k});
    title(ttl)
    legend(struct2cell(subj_ID),'Location', 'eastoutside') % 'Orientation','horizontal' )
    baseFileName = sprintf('allSubj_CUSUM_Stressor_conn%d.png', k); 
    fullFileName = fullfile(fpath, baseFileName);
    saveas(figure(2), fullFileName)
end