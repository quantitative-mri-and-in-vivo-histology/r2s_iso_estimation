function  Binarised_Data = AnisotropicBinning_InSilico_fromExVivo(varargin)

%31.08.2020: Small code who creates 4 different figures with the "loaded"
% matlab file. If it is OLS or MLE, it should apply the same principle. But
% it is required that the entire dataset is calculated.
% OLS: (Max TE data (1,3 and 5),(DispNoMye,NoDispNoMye,DispMye,NoDispMye),
% 67 kappa values,45 angles, all models)

load(varargin{1});

if nargin == 2
    params = varargin{2};
    clearvars -except 'OLS_Fitting' 'params'
else
    error('Params variable is needed for proper plotting');
end

% Ex vivo study: These are defined as default but change it accordingly:
% The way that it is estimated goes as follows: 
% 1.- Obtain the \kappa distribution per anisotropic bin in ex vivo data:
% h = histogram(ExtraValues.kappa_aniso_bin{bin_indx,kappa_range},[beg:0.1:end);
% 2.- Obtain the cumulative bin (sum across all the bins).
% 3.- Replace the number of bins with an array of equal length to this
% number with variable = to the half of the bin.
% full_ar = cat(2,full_ar,beg + 0.05*bin_indx*ones(1,bins(bins_indx)));
% 4.- To check if beta distribution is suitable, normalise the array to be
% between 0 and 1.
% h = histfit((full_ar-beg)/(end-beg),bin_indx,'beta');
% 5.- Obtain the beta parameters (or from any other distribution if
% required).
% fitdist((full_ar'- beg)/(end-beg),'beta')

% \kappa distribution for \kappa < 1.0: beta_function(3.145,1.234)
w_values_k1 = betapdf(0:0.1:1.0,3.145,1.234);
w_values_k1 = w_values_k1/sum(w_values_k1);
% Obtained manually from fitted study
std_w_k1 = [0.0016, 0.0021, 0.0054, 0.0104, 0.0158, 0.0108, 0.0135, 0.0108, 0.0177, 0.0166];

% \kappa distribution for 1.0 < \kappa < 2.5: beta_function(1.34,4.33).
pd_k1_k2p5 = makedist('HalfNormal','mu',0,'sigma',0.4498);
w_values_k1_k2p5 = pdf(pd_k1_k2p5,((1.0:0.1:2.5)-1.0)/(2.5-1.0));
w_values_k1_k2p5 = w_values_k1_k2p5/sum(w_values_k1_k2p5);
% Obtained manually from fitted study
std_w_k1_k2p5 = [0.0082, 0.0060, 0.0064, 0.0073, 0.0080, 0.0068,...
                 0.0063, 0.0054, 0.0062, 0.0049, 0.0031, 0.0054,...
                 0.0040, 0.0039, 0.0042, 0.0042];

% \kappa distribution for 2.5 < \kappa < 6: beta_function(1.34,4.33)
pd_k2p5 = makedist('Exponential','mu',0.2241);
w_values_k2p5 = pdf(pd_k2p5,((2.5:0.1:6.0)-2.5)/(4.5-2.5));
w_values_k2p5 = w_values_k2p5/sum(w_values_k2p5);
% Obtained manually from fitted study
std_w_k2p5 = [0.0219, 0.0153, 0.0228, 0.0163, 0.0218, 0.0177,...
              0.0214, 0.0165, 0.0108, 0.0068, 0.0103, 0.0071,...
              0.0093, 0.0065, 0.0035, 0.0068, 0, 0.0044, 0.0035,...
              zeros(1,17)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET THE AVERAGE ACROSS SAMPLES.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kappa_indx = 1:size(OLS_Fitting,2)-1
    for angle_indx = 1:size(OLS_Fitting,3)
        for temax_indx = 1:size(OLS_Fitting,1)
            for model_indx = 1:size(OLS_Fitting,4)
                DataPlot_Mean(kappa_indx,angle_indx,temax_indx,model_indx,:) = mean(OLS_Fitting{temax_indx,kappa_indx,angle_indx,model_indx}.full_b.*1000.^([0:2]'),2);
                DataPlot_Std(kappa_indx,angle_indx,temax_indx,model_indx,:) = std(OLS_Fitting{temax_indx,kappa_indx,angle_indx,model_indx}.full_b.*1000.^([0:2]'),[],2);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET THE AVERAGE ACROSS K-RANGES based on the ex vivo data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the weighted average in ex vivo scenario for all TE's and data_type
% sum_weight = sum(data * distribution)/sum(distribution);
% std (error_expanded) = 1/sum(distribution) *
% sqrt(sum((std_data*distribution)^2));
wvalues = [w_values_k1(2:end)';w_values_k1_k2p5';w_values_k2p5'];
kappa_indexes_lK1 = 1:numel(w_values_k1(2:end));
kappa_indexes_bK1_K2p5 = (kappa_indexes_lK1(end)+1):numel([w_values_k1(2:end)';w_values_k1_k2p5']);
kappa_indexes_gK2p5 = (kappa_indexes_bK1_K2p5(end)+1):numel(wvalues);
wvalues_mean = repmat(wvalues,[1,size(DataPlot_Mean,2:5)]);

sd_wvalues = [std_w_k1';std_w_k1_k2p5';std_w_k2p5'];
sd_wvalues_estimation = repmat(sd_wvalues,[1,size(DataPlot_Mean,2:5)]);

DataPlot_ExVivo_wKappa = DataPlot_Mean.*wvalues_mean;
DataPlot_ExVivo_MeanKappa(:,1,:,:,:) = squeeze(sum(DataPlot_ExVivo_wKappa(kappa_indexes_lK1,:,:,:,:)));
DataPlot_ExVivo_MeanKappa(:,2,:,:,:) = squeeze(sum(DataPlot_ExVivo_wKappa(kappa_indexes_bK1_K2p5,:,:,:,:)));
DataPlot_ExVivo_MeanKappa(:,3,:,:,:) = squeeze(sum(DataPlot_ExVivo_wKappa(kappa_indexes_gK2p5,:,:,:,:)));

DataPlot_ExVivo_wSdKappa = (DataPlot_Mean.*sd_wvalues_estimation).^2 + (DataPlot_Std.*wvalues_mean).^2;
DataPlot_ExVivo_SdKappa(:,1,:,:,:) = squeeze(sqrt(sum(DataPlot_ExVivo_wSdKappa(kappa_indexes_lK1,:,:,:,:))));
DataPlot_ExVivo_SdKappa(:,2,:,:,:) = squeeze(sqrt(sum(DataPlot_ExVivo_wSdKappa(kappa_indexes_bK1_K2p5,:,:,:,:))));
DataPlot_ExVivo_SdKappa(:,3,:,:,:) = squeeze(sqrt(sum(DataPlot_ExVivo_wSdKappa(kappa_indexes_gK2p5,:,:,:,:))));


% for model_indx = 1:size(DataPlot_Mean,3)
%     for beta_indx = 1:size(DataPlot_Mean,1)
%         for temax_indx = 1:size(DataPlot_Mean,2)
%             % Mean per range
%             DataPlot_ExVivo_K1_Mean(:,temax_indx,model_indx,beta_indx) = sum(squeeze(DataPlot_Mean(beta_indx,temax_indx,model_indx,1:10,:))'.*w_values_k1(2:end),2);
%             DataPlot_ExVivo_K1_K2p5_Mean(:,temax_indx,model_indx,beta_indx) = sum(squeeze(DataPlot_Mean(beta_indx,temax_indx,model_indx,10:25,:))'.*w_values_k1_k2p5,2);
%             DataPlot_ExVivo_K2p5_Mean(:,temax_indx,model_indx,beta_indx) = sum(squeeze(DataPlot_Mean(beta_indx,temax_indx,model_indx,25:60,:))'.*w_values_k2p5,2);
%           
%             % Std due to weight-std
%             Weight_K1_Std = sum((squeeze(DataPlot_Mean(beta_indx,temax_indx,model_indx,1:10,:))'.*std_w_k1).^2,2);
%             Weight_K1_K2p5_Std = sum((squeeze(DataPlot_Mean(beta_indx,temax_indx,model_indx,10:25,:))'.*std_w_k1_k2p5).^2,2);            
%             Weight_K2p5_Std = sum((squeeze(DataPlot_Mean(beta_indx,temax_indx,model_indx,25:60,:))'.*std_w_k2p5).^2,2);
%             
%             % Std due to beta-std
%             Beta_K1_Std = sum((squeeze(DataPlot_Std(beta_indx,temax_indx,model_indx,1:10,:))'.*w_values_k1(2:end)).^2,2);
%             Beta_K1_K2p5_Std = sum((squeeze(DataPlot_Std(beta_indx,temax_indx,model_indx,10:25,:))'.*w_values_k1_k2p5).^2,2);            
%             Beta_K2p5_Std = sum((squeeze(DataPlot_Std(beta_indx,temax_indx,model_indx,25:60,:))'.*w_values_k2p5).^2,2);
%             
%             % Final std
%             DataPlot_ExVivo_K1_Std(:,temax_indx,model_indx,beta_indx) = sqrt(Weight_K1_Std + Beta_K1_Std);
%             DataPlot_ExVivo_K1_K2p5_Std(:,temax_indx,model_indx,beta_indx) = sqrt(Weight_K1_K2p5_Std + Beta_K1_K2p5_Std);            
%             DataPlot_ExVivo_K2p5_Std(:,temax_indx,model_indx,beta_indx) = sqrt(Weight_K2p5_Std + Beta_K2p5_Std);
%         end
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET THE AVERAGE ACROSS IRREGULAR ANGULAR ORIENTATIONS
% based on the ex vivo data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FJv21(08.02): For all the bins, except the first, the probability
% distribution is uniform, therefore DataPlot_ExVivo.. is weighted, again,
% by the corresponding probability distribution per bin. For the first bin,
% unfortunately, it is a linear increment. Therefore, the probability
% distribution is defined as 2*x/(x_fin^2 - x_in^2) with x the angular
% "position".

if ~exist('AngleDensityExVivo_AnisotropicBin.mat','file')
    error('Do not run this section of the code if you do not have the file: AngleDensityExVivo_AnisotropicBin');
else
    load('AngleDensityExVivo_AnisotropicBin.mat');
end

AngMidValue_Indexes = [squeeze(floor(MatrixForInSilico(:,1,2:4)/2)) squeeze(ceil(MatrixForInSilico(:,2,2:4)/2))];
AngMidValue_Indexes(AngMidValue_Indexes(:,1:3) <= 0) = 1;
AngMidValue_Indexes(AngMidValue_Indexes(:,4:6) > 45) = 45;

SumIndexes = (AngMidValue_Indexes(:,4:6).*(AngMidValue_Indexes(:,4:6)-1))/2 - (AngMidValue_Indexes(:,1:3).*(AngMidValue_Indexes(:,1:3)-1))/2;

XAnisoBin_Mean = AngMidValue_Indexes(:,4:6) + AngMidValue_Indexes(:,1:3);
XAnisoBin_Sd = XAnisoBin_Mean - 2*AngMidValue_Indexes(:,1:3);

for bin_indx = 1:20 % 20 anisotropic bins
    % If the .mat file is loaded, then variable name is applied properly
    Weight_Coefficients_lK1 = AngMidValue_Indexes(bin_indx,1):AngMidValue_Indexes(bin_indx,4);
    Weight_Coefficients_bK1_K2p5 = AngMidValue_Indexes(bin_indx,2):AngMidValue_Indexes(bin_indx,5);
    Weight_Coefficients_gK2p5 = AngMidValue_Indexes(bin_indx,3):AngMidValue_Indexes(bin_indx,6);

    switch bin_indx == 1
        case 1
            std_weight_K1 = Weight_Coefficients_lK1'./SumIndexes(bin_indx,1);
            std_weight_K1_K2p5 = Weight_Coefficients_bK1_K2p5'./SumIndexes(bin_indx,2);
            std_weight_K2p5 = Weight_Coefficients_gK2p5'./SumIndexes(bin_indx,3);
            % XBin** are the mean and sd re-estimated irregular angular
            % binning.
            XAnisoBin_Mean(bin_indx,:) = 2*[sum(Weight_Coefficients_lK1.*std_weight_K1),...
                                            sum(Weight_Coefficients_bK1_K2p5.*std_weight_K1_K2p5),...
                                            sum(Weight_Coefficients_gK2p5.*std_weight_K2p5)];
            XAnisoBin_Sd(bin_indx,:) = sqrt([sum((2*Weight_Coefficients_lK1 - XAnisoBin_Mean(bin_indx,1)).^2.*std_weight_K1),...
                                             sum((2*Weight_Coefficients_bK1_K2p5 - XAnisoBin_Mean(bin_indx,2)).^2.*std_weight_K1_K2p5),...
                                             sum((2*Weight_Coefficients_gK2p5 - XAnisoBin_Mean(bin_indx,3)).^2.*std_weight_K2p5)]);

        case 0
            std_weight_K1 = ones(numel(Weight_Coefficients_lK1),1)./numel(Weight_Coefficients_lK1);
            std_weight_K1_K2p5 = ones(numel(Weight_Coefficients_bK1_K2p5'),1)./numel(Weight_Coefficients_bK1_K2p5);
            std_weight_K2p5 = ones(numel(Weight_Coefficients_gK2p5'),1)./numel(Weight_Coefficients_gK2p5);
    end

    AnisoBin_DataPlot_Mean(bin_indx,1,:,:,:) = squeeze(sum(DataPlot_ExVivo_MeanKappa(Weight_Coefficients_lK1,:,:,:,:).*std_weight_K1));
    AnisoBin_DataPlot_Mean(bin_indx,2,:,:,:) = squeeze(sum(DataPlot_ExVivo_MeanKappa(Weight_Coefficients_bK1_K2p5,:,:,:,:).*std_weight_K1_K2p5));
    AnisoBin_DataPlot_Mean(bin_indx,3,:,:,:) = squeeze(sum(DataPlot_ExVivo_MeanKappa(Weight_Coefficients_gK2p5,:,:,:,:).*std_weight_K2p5));
    
    AnisoBin_DataPlot_Sd(bin_indx,1,:,:,:) = squeeze(sqrt(sum((DataPlot_ExVivo_SdKappa(Weight_Coefficients_lK1,:,:,:,:).*std_weight_K1).^2)));
    AnisoBin_DataPlot_Sd(bin_indx,2,:,:,:) = squeeze(sqrt(sum((DataPlot_ExVivo_SdKappa(Weight_Coefficients_bK1_K2p5,:,:,:,:).*std_weight_K1_K2p5).^2)));
    AnisoBin_DataPlot_Sd(bin_indx,3,:,:,:) = squeeze(sqrt(sum((DataPlot_ExVivo_SdKappa(Weight_Coefficients_gK2p5,:,:,:,:).*std_weight_K2p5).^2)));    
end

%Table creation here:
tabVarNames = ["Gratio" "DataType" "T2ValueReference" "Dispersion_Range" "TEmax" "Mean_bin_angle" "Sd_bin_angle"...
               "Alpha0_mean" "Alpha1_mean" "Alpha0_sd" "Alpha1_sd" "Beta0_mean" "Beta1_mean" "Beta2_mean" "Beta0_sd",...
               "Beta1_sd" "Beta2_sd"];
MeanDataTable = reshape(AnisoBin_DataPlot_Mean,[prod(size(AnisoBin_DataPlot_Mean,1:3)),prod(size(AnisoBin_DataPlot_Mean,4:5))]);
SdDataTable = reshape(AnisoBin_DataPlot_Sd,[prod(size(AnisoBin_DataPlot_Sd,1:3)),prod(size(AnisoBin_DataPlot_Sd,4:5))]);
GratioTable = repmat(params.gratio,size(MeanDataTable,1),1);
DataTypeTable = repmat(params.datatype,size(MeanDataTable,1),1);
T2ValueReferenceTable = repmat(params.T2Reference,size(MeanDataTable,1),1);

Dispersion_perTEmax = repmat(["Highly" "Mildly" "Negligible"],size(AnisoBin_DataPlot_Mean,1),1);
DispersionTable = repmat(Dispersion_perTEmax(:),size(AnisoBin_DataPlot_Mean,3),1);

TEmaxTable = repmat(params.TEmaxTable,size(MeanDataTable,1),1);

CellTable = {GratioTable,DataTypeTable,T2ValueReferenceTable,DispersionTable,TEmaxTable,...
             XAnisoBin_Mean(:),XAnisoBin_Sd(:),MeanDataTable(:,1:2),SdDataTable(:,1:2),...
             MeanDataTable(:,4:6),SdDataTable(:,4:6)};

Binarised_Data = cell2table(CellTable,tabVarNames);
Binarised_Data.Properties.VariableUnits = ["","","","","ms","degrees","degrees",...
                                          "a.u.","s^{-1}","a.u.","s^{-1}","a.u.","s^{-1}","s^{2}",...
                                          "a.u.","s^{-1}","s^{2}"];
% FJv21(08.02): Here will be added a small test plot. In theory, the
% anisotropy bins should coincide with the main curve.
%     markers = {{'ob--', 'b', 'b'},{'sy:', 'y', 'y'},{'xr-.', 'r', 'r'}};
%     sgtitle('\beta variation in silico analysis for all models at max TE different ex vivo \kappa.');
%     for model_indx = 1:4
%         for beta_indx = 1:3
%             figure(100);
%             subplot(3,4,model_indx + 4*(beta_indx-1)), errorbar(2:2:90,squeeze(DataPlot_ExVivo_K1_Mean(:,end,model_indx,beta_indx)),...
%                                                                 squeeze(DataPlot_ExVivo_K1_Std(:,end,model_indx,beta_indx)),'o',...
%                                                                 'DisplayName','\kappa =< 1.0'); grid on; box off; xlim([-5,95]); hold on; 
%             subplot(3,4,model_indx + 4*(beta_indx-1)), errorbarxy(XBinK1(:,1)',squeeze(AnisoBinK1_Mean_InSilico(:,end,model_indx,beta_indx)),...
%                                                                   XBinK1(:,2)',squeeze(AnisoBinK1_Std_InSilico(:,end,model_indx,beta_indx)),markers{1}); ...
%                                                                   grid on; box off; xlim([-5,95]); hold on;                                                
%             figure(101);
%             subplot(3,4,model_indx + 4*(beta_indx-1)), errorbar(2:2:90,squeeze(DataPlot_ExVivo_K1_K2p5_Mean(:,end,model_indx,beta_indx)),...
%                                                                 squeeze(DataPlot_ExVivo_K1_K2p5_Std(:,end,model_indx,beta_indx)),'o',...
%                                                                 'DisplayName','\kappa =< 1.0'); grid on; box off; xlim([-5,95]); hold on; 
%             subplot(3,4,model_indx + 4*(beta_indx-1)), errorbarxy(XBinK1_K2p5(:,1)',squeeze(AnisoBinK1_K2p5_Mean_InSilico(:,end,model_indx,beta_indx)),...
%                                                                   XBinK1_K2p5(:,2)',squeeze(AnisoBinK1_K2p5_Std_InSilico(:,end,model_indx,beta_indx)),markers{2}); ...
%                                                                   grid on; box off; xlim([-5,95]); hold on;                                                
%             figure(102);
%             subplot(3,4,model_indx + 4*(beta_indx-1)), errorbar(2:2:90,squeeze(DataPlot_ExVivo_K2p5_Mean(:,end,model_indx,beta_indx)),...
%                                                                squeeze(DataPlot_ExVivo_K2p5_Std(:,end,model_indx,beta_indx)),'x',...
%                                                                'DisplayName','\kappa >= 2.5'); grid on; box off; xlim([-5,95]); hold on;     
%             subplot(3,4,model_indx + 4*(beta_indx-1)), errorbarxy(XBinK2p5(:,1)',squeeze(AnisoBinK2p5_Mean_InSilico(:,end,model_indx,beta_indx))',...
%                                                                   XBinK2p5(:,2)',squeeze(AnisoBinK2p5_Std_InSilico(:,end,model_indx,beta_indx))',markers{3});...
%                                                                   grid on; box off; xlim([-5,95]); hold on;                                                 
%             subplot(3,4,model_indx), title(['Model ' num2str(model_indx) ' analysis']);
%         end
%         %subplot(3,4,1+4*(beta_indx-1)), ylabel(['\beta_' num2str(beta_indx-1) ' parameter ' label_units{beta_indx}]);
%     end
%     subplot(3,4,11), xlabel('Angular orientation (\theta_\mu) [degrees]'); 

Binarised_Data{1} = AnisoBinK1_Mean_InSilico;
Binarised_Data{2} = AnisoBinK1_Std_InSilico;
Binarised_Data{3} = AnisoBinK1_K2p5_Mean_InSilico;
Binarised_Data{4} = AnisoBinK1_K2p5_Std_InSilico;
Binarised_Data{5} = AnisoBinK2p5_Mean_InSilico;
Binarised_Data{6} = AnisoBinK2p5_Std_InSilico;

end