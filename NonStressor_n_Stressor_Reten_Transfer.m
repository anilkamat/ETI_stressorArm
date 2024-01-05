clc; clear; warning off; 
connection_GC = {'LPFC-->RPFC','LPFC-->LPMC','LPFC-->RPMC','LPFC-->SMA',...
    'RPFC-->LPFC','RPFC-->LPMC','RPFC-->RPMC','RPFC-->SMA',...
    'LPMC-->LPFC','LPMC-->RPFC','LPMC-->RPMC','LPMC-->SMA',...
    'RPMC-->LPFC','RPMC-->RPFC','RPMC-->LPMC','RPMC-->SMA',...
    'SMA-->LPFC','SMA-->RPFC','SMA-->LPMC','SMA-->RPMC'};
%Load the computed connectivities of Mannequin(Non-stressor) and stressor arm group
load('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Results\ETI_reten_transf_MA.mat')
load('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Results\ETI_reten_transf_SA.mat')

% NOTE: The first three trials of "ETI_reten_transf_MA.mat" and
% "ETI_reten_transf_SA.mat" are for retention task and the last three
% trials are the transfer tasks.

%% Determine the window size for moving mean analysis
q = 1;
for con= 1:1 % numel(connection_GC)  % loop over connectivity
    close all
    p =1; % counter of total number of trials in entire experiment
    for i=1:1 %size(fieldnames(ETI_normal),1)
        d = sprintf('Day_%d',4); % forth day is the retention day on which retention and trasfer both were present.
        for j=1:size(fieldnames(ETI_reten_transf_MA.(d)),1)
            t = sprintf('Tr_%d',j);
            %Mean of successful ETI normal 
            M_sMA(p) = mean(ETI_reten_transf_MA.(d).(t).Succ(:,con));
            STD_sMA(p) = std(ETI_reten_transf_MA.(d).(t).Succ(:,con));
            GM_sMA(q) = M_sMA(p);

            %Mean of successful ETI stressor arm
            M_sSA(p) = mean(ETI_reten_transf_SA.(d).(t).Succ(:,con));
            STD_sSA(p) = std(ETI_reten_transf_SA.(d).(t).Succ(:,con));
            GM_sSA(q) = M_sSA(p);
            
            p = p+1;
            q = q+1;
        end
    end
    for i=1:1:20
        Mmean_sMA = movmean(M_sMA,i+1);
        SAD_sMA(i) = sum(abs(Mmean_sMA-M_sMA)); % sum of absolute difference
        Mmean_sSA = movmean(M_sSA,i+1);
        SAD_sSA(i) = sum(abs(Mmean_sSA-M_sSA)); % sum of absolute difference
    end
    figure(25)
    plot(SAD_sMA, 'color','green', 'LineWidth',2); hold on;
%     plot(SAD_uMA, 'color','blue', 'LineWidth',2); hold on;
    plot(SAD_sSA, 'color','magenta', 'LineWidth',2); hold on;
%     plot(SAD_uSA, 'color','red', 'LineWidth',2); hold on;
    legend('sMA','sSA')
    ylabel('Sum absolute difference between raw and movMean')
    xlabel('Window size')
    title('Effect of window size moving mean')
end

%% Moving average of both the groups
fpath = 'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Results\Plots_Reten_Transf_MA_SA';
q = 1;
for con= 1:20 %20%numel(connection_GC)  % loop over connectivity
    close all;
    p =1; % counter of total number of trials in entire experiment
    for i=1:size(fieldnames(ETI_reten_transf_MA),1)  %loop over days
        d = sprintf('Day_%d',4);
        for j=1:size(fieldnames(ETI_reten_transf_MA.(d)),1)  %loop over trials
            t = sprintf('Tr_%d',j);
            %Mean of successful ETI normal 
            M_sMA(p) = mean(ETI_reten_transf_MA.(d).(t).Succ(:,con));
            num_sMA(p) = size(ETI_reten_transf_MA.(d).(t).Succ(:,con),1);
            STD_sMA(p) = std(ETI_reten_transf_MA.(d).(t).Succ(:,con));
            GM_sMA(q) = M_sMA(p);
            
            %Mean of successful ETI stressor arm
            M_sSA(p) = mean(ETI_reten_transf_SA.(d).(t).Succ(:,con));
            num_sSA(p) = size(ETI_reten_transf_SA.(d).(t).Succ(:,con),1);
            STD_sSA(p) = std(ETI_reten_transf_SA.(d).(t).Succ(:,con));
            GM_sSA(q) = M_sSA(p);
            
            p = p+1;
            q = q+1;
        end
    end
    Mmean_sMA = movmean(M_sMA,2);
    Mmean_sMA_allconn(con,:) = Mmean_sMA;
    Mmean_sSA = movmean(M_sSA,2);
    Mmean_sSA_allconn(con,:) = Mmean_sSA;
    f= figure(1);%f = figure('visible','off');
    errorbar(M_sMA,STD_sMA,"-s",'color','blue',"MarkerSize",10,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90]); hold on;
    errorbar(M_sSA,STD_sSA,"-s","MarkerSize",10,"MarkerEdgeColor","magenta","MarkerFaceColor",[0.9961 0.8 0.9961]); hold on;
    plot(Mmean_sMA,'LineWidth',2, 'color','blue')
    plot(Mmean_sSA,'LineWidth',2, 'color','magenta')
    ylim([-0.15 3.5])
    xlim([0.5 6.5])
    title(['Connectivity over trials:successful MA vs SA' connection_GC(con)])
    legend('Non-stressor','Stressor')
    xlabel('Trial')
    ylabel('Connectivity')
    vline2(3,'color','black')
    text([1 4], 2*[1.25 1.25], {'Retention', 'Transfer'})
    baseFileName = sprintf('Conn%d_successful_MA_vs_SA.png',con);
    fullFileName = fullfile(fpath, baseFileName);
    saveas(f,fullFileName)
end