function [Data_Reg] = get_LS_channelData(subj_data)
Data = subj_data.data;
channels = subj_data.textdata; % list of channels in the current loaded data
% req_HbO_channels = {'S1_D1_HbO','S1_D15_HbO','S2_D2_HbO','S2_D3_HbO','S3_D1_HbO','S3_D3_HbO','S3_D5_HbO','S4_D1_HbO','S4_D3_HbO','S4_D15_HbO','S5_D2_HbO','S5_D3_HbO','S5_D4_HbO','S5_D5_HbO','S6_D4_HbO','S6_D5_HbO','S6_D6_HbO','S7_D2_HbO','S7_D4_HbO','S7_D7_HbO','S8_D4_HbO','S8_D6_HbO','S8_D7_HbO','S9_D8_HbO','S9_D15_HbO','S10_D8_HbO','S10_D9_HbO','S10_D13_HbO','S11_D9_HbO','S11_D10_HbO','S12_D8_HbO','S12_D9_HbO','S12_D15_HbO','S13_D9_HbO','S13_D10_HbO','S13_D12_HbO','S13_D13_HbO','S14_D12_HbO','S14_D13_HbO','S14_D14_HbO','S15_D10_HbO','S15_D11_HbO','S15_D12_HbO','S16_D11_HbO','S16_D12_HbO','S16_D14_HbO'};   %set of required channels
% [~, Seq] = CStrAinBP(req_HbO_channels, channels); % compares two set of channels names and Seq store the index of req.indecs
% Data = Data(Seq,:);
%LPFC = Data([ 1 2 5 7 8 9 10 11 12],:); % list from suturing montage
Lpfc_channels = {'S1_D1_HbO','S1_D15_HbO','S2_D3_HbO','S3_D1_HbO','S3_D3_HbO'...
    ,'S3_D5_HbO','S4_D1_HbO','S4_D3_HbO','S2_D15_HbO'};
[~, Seq] = CStrAinBP(Lpfc_channels, channels);
LPFC = Data(:,Seq);
LPFC = nanmean(LPFC,2);

%RPFC = Data([28 29 31 32 33 35 37 38 40],:); list from suturing montage
Rpfc_channels = {'S9_D8_HbO','S9_D15_HbO','S10_D8_HbO','S10_D9_HbO','S10_D13_HbO'...
    ,'S11_D9_HbO','S12_D8_HbO','S12_D9_HbO','S13_D9_HbO'};
[~, Seq] = CStrAinBP(Rpfc_channels, channels);
RPFC = Data(:, Seq);
RPFC = nanmean(RPFC,2);      

% LSMA  = Data([4 13 14 16 21],:);
lsma_channels = {'S2_D2_HbO','S5_D2_HbO','S5_D3_HbO','S5_D5_HbO','S7_D2_HbO'};
[~, Seq] = CStrAinBP(lsma_channels, channels);
LSMA = Data(:, Seq);
LSMA = nanmean(LSMA,2);


%RSMA  = Data([36 41 43 48],:); list from suturing montage
rsma_channels = {'S11_D10_HbO','S13_D10_HbO','S13_D13_HbO','S15_D10_HbO'};
[~, Seq] = CStrAinBP(rsma_channels, channels);
RSMA = Data(:, Seq);
RSMA = nanmean(RSMA,2);
SMA = [LSMA RSMA];
SMA = nanmean(SMA,2);

%LPMC = Data([15 17 18 19 22 23 24 25 26],:); list from suturing montage
lpmc_channels = {'S5_D4_HbO','S6_D4_HbO','S6_D5_HbO','S6_D6_HbO','S7_D4_HbO',...
    'S7_D7_HbO','S8_D4_HbO','S8_D6_HbO','S8_D7_HbO'};
[~, Seq] = CStrAinBP(lpmc_channels, channels);
LPMC = Data(:, Seq);
LPMC = nanmean(LPMC,2);

% RPMC = Data([42 45 46 47 49 50 51 52 53],:); list from suturing montage
rpmc_channels = {'S13_D12_HbO','S14_D12_HbO','S14_D13_HbO','S14_D14_HbO',...
    'S15_D11_HbO','S15_D12_HbO','S16_D11_HbO','S16_D12_HbO','S16_D14_HbO'};
[~, Seq] = CStrAinBP(rpmc_channels, channels);
RPMC = Data(:, Seq);
RPMC = nanmean(RPMC,2);
Data_Reg = [LPFC RPFC LPMC RPMC SMA]'; % CARE for the order of each region
%valid_Ch_names = Cnames;
NaN_rows = find(all(isnan(Data_Reg),1));
%valid_Ch_names(:,NaN_rows) = [];

end