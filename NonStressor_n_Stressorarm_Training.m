%% This codes plots the strength of GC of stressor arm experimnets across all the trials. It also computs moving average and overlays on the GC strength. 
clc; clear; warning off; 
connection_GC = {'LPFC-->RPFC','LPFC-->LPMC','LPFC-->RPMC','LPFC-->SMA',...
    'RPFC-->LPFC','RPFC-->LPMC','RPFC-->RPMC','RPFC-->SMA',...
    'LPMC-->LPFC','LPMC-->RPFC','LPMC-->RPMC','LPMC-->SMA',...
    'RPMC-->LPFC','RPMC-->RPFC','RPMC-->LPMC','RPMC-->SMA',...
    'SMA-->LPFC','SMA-->RPFC','SMA-->LPMC','SMA-->RPMC'};
%Load the computed connectivities of Mannequin and stressor arm group
load('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Results\ETI_normal.mat')
load('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Results\ETI_StressorArm.mat')
load('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Results\ETI_retentation_MA.mat')
load('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Results\ETI_reten_SA.mat')
load('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Results\ETI_transfer_MA.mat')
load('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Results\ETI_transfer_SA.mat')

%% Determine the window size for moving mean analysis
q = 1;
for con= 1:numel(connection_GC)  % loop over connectivity
    close all
    p =1; % counter of total number of trials in entire experiment
    for i=1:size(fieldnames(ETI_normal),1)
        d = sprintf('Day_%d',i);
        for j=1:size(fieldnames(ETI_normal.(d)),1)
            t = sprintf('Tr_%d',j);
            %Mean of successful ETI normal 
            M_sMA(p) = mean(ETI_normal.(d).(t).Succ(:,con));
            STD_sMA(p) = std(ETI_normal.(d).(t).Succ(:,con));
            GM_sMA(q) = M_sMA(p);
            %Mean of unsuccessful ETI normal 
            M_uMA(p) = mean(ETI_normal.(d).(t).Unsucc(:,con));
            STD_uMA(p) = std(ETI_normal.(d).(t).Unsucc(:,con));
            GM_uMA(q) = M_uMA(p);

            %Mean of successful ETI stressor arm
            M_sSA(p) = mean(ETI_StressorArm.(d).(t).Succ(:,con));
            STD_sSA(p) = std(ETI_StressorArm.(d).(t).Succ(:,con));
            GM_sSA(q) = M_sSA(p);
            %Mean of ETI stressor arm
            M_uSA(p) = mean(ETI_StressorArm.(d).(t).Unsucc(:,con));
            STD_uSA(p) = std(ETI_StressorArm.(d).(t).Unsucc(:,con));
            GM_uSA(q) = M_uSA(p);
            
            p = p+1;
            q = q+1;
        end
    end
    for i=1:1:20
        Mmean_sMA = movmean(M_sMA,i+1);
        SAD_sMA(i) = sum(abs(Mmean_sMA-M_sMA)); % sum of absolute difference
        Mmean_uMA = movmean(M_uMA,i+1);
        SAD_uMA(i) = sum(abs(Mmean_uMA-M_uMA)); % sum of absolute difference
        Mmean_sSA = movmean(M_sSA,i+1);
        SAD_sSA(i) = sum(abs(Mmean_sSA-M_sSA)); % sum of absolute difference
        Mmean_uSA = movmean(M_uSA,i+1);
        SAD_uSA(i) = sum(abs(Mmean_uSA-M_uSA)); % sum of absolute difference
    end
    figure(25)
    plot(SAD_sMA, 'color','green', 'LineWidth',2); hold on;
    plot(SAD_uMA, 'color','blue', 'LineWidth',2); hold on;
    plot(SAD_sSA, 'color','magenta', 'LineWidth',2); hold on;
    plot(SAD_uSA, 'color','red', 'LineWidth',2); hold on;
    legend('sMA','uMA','sSA','uSA')
    ylabel('Sum absolute difference between raw and movMean')
    xlabel('Window size')
    title('Effect of window size moving mean')
end

%% Moving average of both the groups
fpath = 'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Results\Plots';
q = 1;
for con= 1:20 %20%numel(connection_GC)  % loop over connectivity
    close all;
    p =1; % counter of total number of trials in entire experiment
    for i=1:size(fieldnames(ETI_normal),1)  %loop over days
        d = sprintf('Day_%d',i);
        for j=1:size(fieldnames(ETI_normal.(d)),1)  %loop over trials
            t = sprintf('Tr_%d',j);
            %Mean of successful ETI normal 
            M_sMA(p) = mean(ETI_normal.(d).(t).Succ(:,con));
            num_sMA(p) = size(ETI_normal.(d).(t).Succ(:,con),1);
            STD_sMA(p) = std(ETI_normal.(d).(t).Succ(:,con));
            GM_sMA(q) = M_sMA(p);
            %Mean of unsuccessful ETI normal 
            M_uMA(p) = mean(ETI_normal.(d).(t).Unsucc(:,con));
            num_uMA(p) = size(ETI_normal.(d).(t).Unsucc(:,con),1);
            STD_uMA(p) = std(ETI_normal.(d).(t).Unsucc(:,con));
            GM_uMA(q) = M_uMA(p);
            
            %Mean of successful ETI stressor arm
            M_sSA(p) = mean(ETI_StressorArm.(d).(t).Succ(:,con));
            num_sSA(p) = size(ETI_StressorArm.(d).(t).Succ(:,con),1);
            STD_sSA(p) = std(ETI_StressorArm.(d).(t).Succ(:,con));
            GM_sSA(q) = M_sSA(p);
            %Mean of ETI stressor arm
            M_uSA(p) = mean(ETI_StressorArm.(d).(t).Unsucc(:,con));
            num_uSA(p) = size(ETI_StressorArm.(d).(t).Unsucc(:,con),1);
            STD_uSA(p) = std(ETI_StressorArm.(d).(t).Unsucc(:,con));
            GM_uSA(q) = M_uSA(p);
            
            p = p+1;
            q = q+1;
        end
    end
    Mmean_sMA = movmean(M_sMA,3);
    Mmean_sMA_allconn(con,:) = Mmean_sMA;
    Mmean_uMA = movmean(M_uMA,3);
    Mmean_uMA_allconn(con,:) = Mmean_uMA;
    Mmean_sSA = movmean(M_sSA,3);
    Mmean_sSA_allconn(con,:) = Mmean_sSA;
    Mmean_uSA = movmean(M_uSA,3);
    Mmean_uSA_allconn(con,:) = Mmean_uSA;
    %plot for successful MA vs SA
    f= figure(1);%f = figure('visible','off');
    errorbar(M_sMA,STD_sMA,"-s",'color','blue',"MarkerSize",10,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90]); hold on;
    errorbar(M_sSA,STD_sSA,"-s","MarkerSize",10,"MarkerEdgeColor","magenta","MarkerFaceColor",[0.9961 0.8 0.9961]); hold on;
    plot(Mmean_sMA,'LineWidth',2, 'color','blue')
    plot(Mmean_sSA,'LineWidth',2, 'color','magenta')
    ylim([-0.15 3])
    title(['Connectivity over trials:successful MA vs SA' connection_GC(con)])
    legend('sMA','sSA')
    xlabel('Trial')
    ylabel('Connectivity')
    vline2([10, 20],'color',[0.6350 0.0780 0.1840])
    text([5 15 25], 2*[1 1 1], {'Day 1', 'Day 2', 'Day 3'})
    baseFileName = sprintf('Conn%d_successful_MA_vs_SA.png',con);
    fullFileName = fullfile(fpath, baseFileName);
    %saveas(f,fullFileName)

    %plot for unsuccessful MA vs SA
    f = figure(2);%f = figure('visible','off');
    errorbar(M_uMA,STD_uMA,"-s","MarkerSize",10,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90]); hold on;
    errorbar(M_uSA,STD_uSA,"-s","MarkerSize",10,"MarkerEdgeColor","magenta","MarkerFaceColor",[0.9961 0.8 0.9961]); hold on;
    plot(Mmean_uMA,'LineWidth',2, 'color','blue')
    plot(Mmean_uSA,'LineWidth',2, 'color','magenta')
    ylim([-0.15 3])
    title(['Connectivity over trials:unsuccessful MA vs SA' connection_GC(con)])
    legend('uMA','uSA')
    xlabel('Trial')
    ylabel('Connectivity')
    vline2([10, 20],'color',[0.6350 0.0780 0.1840])
    text([5 15 25], 2*[1 1 1], {'Day 1', 'Day 2', 'Day 3'})
    baseFileName = sprintf('Conn%d_unsuccessful_MA_vs_SA.png',con);
    fullFileName = fullfile(fpath, baseFileName);
    %saveas(f,fullFileName)

    %plot for MA successful vs unsuccessful
    f = figure(3);%f = figure('visible','off');
    errorbar(M_sMA,STD_sMA,"-s","MarkerSize",10,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90]); hold on;
    errorbar(M_uMA,STD_uMA,"-s","MarkerSize",10,"MarkerEdgeColor","magenta","MarkerFaceColor",[0.9961 0.8 0.9961]); hold on;
    plot(Mmean_sMA,'LineWidth',2, 'color','blue')
    plot(Mmean_uMA,'LineWidth',2, 'color','magenta')
    ylim([-0.15 3])
    title(['Connectivity over trials:MA successful vs unsuccessful' connection_GC(con)])
    legend('sMA','uMA')
    xlabel('Trial')
    ylabel('Connectivity')
    vline2([10, 20], 'color',[0.6350 0.0780 0.1840])
    text([5 15 25], 2*[1 1 1], {'Day 1', 'Day 2', 'Day 3'})
    baseFileName = sprintf('Conn%d_MA_succ_vs_unsucc.png',con);
    fullFileName = fullfile(fpath, baseFileName);
    %saveas(f,fullFileName)

    %plot for SA successful vs unsuccessful
    f = figure(4);%f = figure('visible','off');
    errorbar(M_sSA,STD_sSA,"-s","MarkerSize",10,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90]); hold on;
    errorbar(M_uSA,STD_uSA,"-s","MarkerSize",10,"MarkerEdgeColor","magenta","MarkerFaceColor",[0.9961 0.8 0.9961]); hold on;
    plot(Mmean_sSA,'LineWidth',2, 'color','blue')
    plot(Mmean_uSA,'LineWidth',2, 'color','magenta')
    ylim([-0.15 3])
    title(['Connectivity over trials:SA successful vs unsuccessful' connection_GC(con)])
    legend('sSA','uSA')
    xlabel('Trial')
    ylabel('Connectivity')
    vline2([10, 20],'color',[0.6350 0.0780 0.1840])
    text([5 15 25], 2*[1 1 1], {'Day 1', 'Day 2', 'Day 3'})
    baseFileName = sprintf('Conn%d_SA_succ_vs_unsucc.png',con);
    fullFileName = fullfile(fpath, baseFileName);
    %saveas(f,fullFileName)    
%     num_sMA
%     num_uMA
%     num_sSA
%     num_uSA
%     f = figure(5);
%     plot(num_sMA, 'linewidth',2, 'color','blue'); hold on;
%     plot(num_uMA, 'linewidth',2, 'color','green'); hold on;
%     plot(num_sSA, 'linewidth',2, 'color','magenta'); hold on;
%     plot(num_uSA, 'linewidth',2, 'color','black');
%     xlabel('Trials')
%     ylabel('Number of Subjects')
%     legend('sMA','uMA','sSA','uSA')
%     
%     f = figure(6);
%     plot(sum([num_sMA; num_uMA],1), 'linewidth',2, 'color','blue'); hold on;
%     plot(sum([num_sSA; num_uSA],1), 'linewidth',2, 'color','magenta'); hold on;
%     xlabel('Trials')
%     ylabel('Number of Subjects')
%     legend('MA','SA')
end

%% Moving average of both the groups: Day-1, Day-2, Day-3, Retention, and Transfer
fpath2 = 'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Results\Plot_D123_reten_transfer_';
ETI_normal.Day_4 = ETI_retention_MA.Day_4; 
ETI_StressorArm.Day_4 = ETI_reten_SA.Day_4;
ETI_normal.Day_5 = ETI_transfer_MA.Day_4; 
ETI_StressorArm.Day_5 = ETI_transfer_SA.Day_4;
q = 1;
for con= 1:20 % 20%numel(connection_GC)  % loop over connectivity
    close all;
    p =1; % counter of total number of trials in entire experiment
    for i=1:size(fieldnames(ETI_normal),1)  %loop over days
        d = sprintf('Day_%d',i);
        for j=1:size(fieldnames(ETI_normal.(d)),1)  %loop over trials
            t = sprintf('Tr_%d',j);
            %Mean of successful ETI normal 
            M_sMA(p) = mean(ETI_normal.(d).(t).Succ(:,con));
            num_sMA(p) = size(ETI_normal.(d).(t).Succ(:,con),1);
            STD_sMA(p) = std(ETI_normal.(d).(t).Succ(:,con));
            GM_sMA(q) = M_sMA(p);
%             %Mean of unsuccessful ETI normal 
%             M_uMA(p) = mean(ETI_normal.(d).(t).Unsucc(:,con));
%             num_uMA(p) = size(ETI_normal.(d).(t).Unsucc(:,con),1);
%             STD_uMA(p) = std(ETI_normal.(d).(t).Unsucc(:,con));
%             GM_uMA(q) = M_uMA(p);
            
            %Mean of successful ETI stressor arm
            M_sSA(p) = mean(ETI_StressorArm.(d).(t).Succ(:,con));
            num_sSA(p) = size(ETI_StressorArm.(d).(t).Succ(:,con),1);
            STD_sSA(p) = std(ETI_StressorArm.(d).(t).Succ(:,con));
            GM_sSA(q) = M_sSA(p);
            %Mean of ETI stressor arm
%             M_uSA(p) = mean(ETI_StressorArm.(d).(t).Unsucc(:,con));
%             num_uSA(p) = size(ETI_StressorArm.(d).(t).Unsucc(:,con),1);
%             STD_uSA(p) = std(ETI_StressorArm.(d).(t).Unsucc(:,con));
%             GM_uSA(q) = M_uSA(p);
            p = p+1;
            q = q+1;
        end
    end
    Mmean_sMA = movmean(M_sMA,3);
    Mmean_sMA_allconn(con,:) = Mmean_sMA;
%     Mmean_uMA = movmean(M_uMA,3);
%     Mmean_uMA_allconn(con,:) = Mmean_uMA;

    Mmean_sSA = movmean(M_sSA,3);
    Mmean_sSA_allconn(con,:) = Mmean_sSA;
    
%     Mmean_uSA = movmean(M_uSA,3);
%     Mmean_uSA_allconn(con,:) = Mmean_uSA;
    %plot for successful MA vs SA
    %f= figure(1);%f = figure('visible','off');
    f = figure('Renderer', 'painters', 'Position', [10 10 900 400]);
    %errorbar(M_sMA,STD_sMA,"-s",'color','blue',"MarkerSize",10,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90]); hold on;
    errorbar(M_sSA,STD_sSA,"-s","MarkerSize",10,"MarkerEdgeColor","magenta","MarkerFaceColor",[0.9961 0.8 0.9961]); hold on;
    %plot(Mmean_sMA,'LineWidth',2, 'color','blue')
    plot(Mmean_sSA,'LineWidth',2, 'color','magenta')
    ylim([-0.15 3])
    xlim([0,37])
    title(['Connectivity over trials:succ Stressor ' connection_GC(con)])
    legend('Stressor', 'Location','northwest') % 'nStressor',
    xlabel('Trial')
    ylabel('Connectivity')
    vline2([10, 20, 30 33],'color',[0.6350 0.0780 0.1840])
    text([5 15 25 31 34], 2*[1 1 1 1 1], {'Day 1', 'Day 2', 'Day 3', 'Reten', 'Trans'})
    baseFileName = sprintf('Conn%d_successful_SA.png',con);
    fullFileName = fullfile(fpath2, baseFileName);
    saveas(f,fullFileName)

%     %plot for unsuccessful MA vs SA
%     f = figure(2);%f = figure('visible','off');
%     errorbar(M_uMA,STD_uMA,"-s","MarkerSize",10,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90]); hold on;
%     errorbar(M_uSA,STD_uSA,"-s","MarkerSize",10,"MarkerEdgeColor","magenta","MarkerFaceColor",[0.9961 0.8 0.9961]); hold on;
%     plot(Mmean_uMA,'LineWidth',2, 'color','blue')
%     plot(Mmean_uSA,'LineWidth',2, 'color','magenta')
%     ylim([-0.15 3])
%     title(['Connectivity over trials:unsuccessful MA vs SA' connection_GC(con)])
%     legend('uMA','uSA')
%     xlabel('Trial')
%     ylabel('Connectivity')
%     vline2([10, 20],'color',[0.6350 0.0780 0.1840])
%     text([5 15 25], 2*[1 1 1], {'Day 1', 'Day 2', 'Day 3'})
%     baseFileName = sprintf('Conn%d_unsuccessful_MA_vs_SA.png',con);
%     fullFileName = fullfile(fpath2, baseFileName);
%     %saveas(f,fullFileName)
% 
%     %plot for MA successful vs unsuccessful
%     f = figure(3);%f = figure('visible','off');
%     errorbar(M_sMA,STD_sMA,"-s","MarkerSize",10,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90]); hold on;
%     errorbar(M_uMA,STD_uMA,"-s","MarkerSize",10,"MarkerEdgeColor","magenta","MarkerFaceColor",[0.9961 0.8 0.9961]); hold on;
%     plot(Mmean_sMA,'LineWidth',2, 'color','blue')
%     plot(Mmean_uMA,'LineWidth',2, 'color','magenta')
%     ylim([-0.15 3])
%     title(['Connectivity over trials:MA successful vs unsuccessful' connection_GC(con)])
%     legend('sMA','uMA')
%     xlabel('Trial')
%     ylabel('Connectivity')
%     vline2([10, 20], 'color',[0.6350 0.0780 0.1840])
%     text([5 15 25], 2*[1 1 1], {'Day 1', 'Day 2', 'Day 3'})
%     baseFileName = sprintf('Conn%d_MA_succ_vs_unsucc.png',con);
%     fullFileName = fullfile(fpath2, baseFileName);
%     %saveas(f,fullFileName)
% 
%     %plot for SA successful vs unsuccessful
%     f = figure(4);%f = figure('visible','off');
%     errorbar(M_sSA,STD_sSA,"-s","MarkerSize",10,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90]); hold on;
%     errorbar(M_uSA,STD_uSA,"-s","MarkerSize",10,"MarkerEdgeColor","magenta","MarkerFaceColor",[0.9961 0.8 0.9961]); hold on;
%     plot(Mmean_sSA,'LineWidth',2, 'color','blue')
%     plot(Mmean_uSA,'LineWidth',2, 'color','magenta')
%     ylim([-0.15 3])
%     title(['Connectivity over trials:SA successful vs unsuccessful' connection_GC(con)])
%     legend('sSA','uSA')
%     xlabel('Trial')
%     ylabel('Connectivity')
%     vline2([10, 20],'color',[0.6350 0.0780 0.1840])
%     text([5 15 25], 2*[1 1 1], {'Day 1', 'Day 2', 'Day 3'})
%     baseFileName = sprintf('Conn%d_SA_succ_vs_unsucc.png',con);
%     fullFileName = fullfile(fpath2, baseFileName);
%     %saveas(f,fullFileName)    

end

%% 
% % plot for the successful subjects
% q = 1; % grand mean of all the trials
% for con= 1:numel(connection_GC)  % loop over connectivity
%     close all
%     p =1; % counter of total number of trials in entire experiment
%     for i=1:size(fieldnames(ETI_normal),1)
%         d = sprintf('Day_%d',i);
%         for j=1:size(fieldnames(ETI_normal.(d)),1)
%             t = sprintf('Tr_%d',j);
%             M(p) = mean(ETI_normal.(d).(t).Succ(:,con));
%             STD(p) = std(ETI_normal.(d).(t).Succ(:,con));
%             GM(q) = M(p);
%             p = p+1;
%             q = q+1;
%         end
%     end
%     f = figure('visible','off');
%     errorbar(M,STD,"-s","MarkerSize",10,"MarkerEdgeColor","magenta","MarkerFaceColor",[0.65 0.85 0.90])
%     ylim([-0.15 3])
%     title('Connectivity over trials:successful')
%     legend(connection_GC(con))
%     xlabel('Trial')
%     ylabel('Connectivity')
%     vline2([10, 20])
%     text([5 15 25], 2*[1 1 1], {'Day 1', 'Day 2', 'Day 3'})
%     baseFileName = sprintf('Conn%d_successful.png',con);
%     fullFileName = fullfile(fpath, baseFileName);
%     saveas(f,fullFileName)
% end
% % plot for the unsuccessful subjects
% q = 1; % grand mean of all the trials
% for con= 1:numel(connection_GC)  % loop over connectivity
%     close all
%     p =1; % counter of total number of trials in entire experiment
%     for i=1:size(fieldnames(ETI_normal),1)
%         d = sprintf('Day_%d',i);
%         for j=1:size(fieldnames(ETI_normal.(d)),1)
%             t = sprintf('Tr_%d',j);
%             M(p) = mean(ETI_normal.(d).(t).Unsucc(:,con));
%             STD(p) = std(ETI_normal.(d).(t).Unsucc(:,con));
%             GM(q) = M(p);
%             p = p+1;
%             q = q+1;
%         end
%     end
%     f = figure('visible','off');
%     errorbar(M,STD,"-s","MarkerSize",10,"MarkerEdgeColor","magenta","MarkerFaceColor",[0.65 0.85 0.90])
%     ylim([-0.15 3])
%     title('Connectivity over trials: unsuccessful')
%     legend(connection_GC(con))
%     xlabel('Trial')
%     ylabel('Connectivity')
%     vline2([10, 20])
%     text([5 15 25], 2*[1 1 1], {'Day 1', 'Day 2', 'Day 3'})
%     baseFileName = sprintf('Conn%d_unsuccessful.png',con);
%     fullFileName = fullfile(fpath, baseFileName);
%     saveas(f,fullFileName)
% end