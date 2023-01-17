function [ExVivoTable,HistInfo] = ...
         AnisotropyBinning_Submission(Angle_Meas,Kappa_Meas,ICVF_Meas,BetaParams,Mask,params)
%ANISOTROPYBINNING Here the beta values are binned anisotropically to
%ensure equal number of "voxels" as a function of angular orientation.
%Several filters are performed.
%   This function uses the angular map obtained across all the measurements
%   (16), the dispersion map,the ICVF map and mask to filter the obtained 
%   \beta parameters in the OC. This allows to bin the data anisotropically
%    with a number of bins defined by params.nbins variable. The size of
%    each input is given below. Four different voxel filters were used in 
%    this code: (1) all angles with ICVF > 0.8, (2) all angles with ICVF >  
%    0.8 and 0 < \kappa <= 1.0, (3) all angles with ICVF > 0.8 and 1.0 < 
%    \kappa < 2.5, and (4) all angles with ICVF > 0.8 and 2.5 <= \kappa.

% Input:
% Angle_Meas: Matrix with size (dim1,dim2,dim3,total_meas)
% This variable contains the voxel-wise angular orientation per angular
% measurement. This is filtered here to bound it between 0 and 90°.
% Kappa_Meas: Matrix with size (dim1,dim2,dim3)
% This variable contains the voxel-wise \kappa value.
% ICVF_Meas: Matrix with size (dim1,dim2,dim3) 
% This variable contains the voxel-wise ICVF value.
% Beta_Params: Matrix with size
% (dim1,dim2,dim3,total_meas,max_TEs,models,betas)
% This variable contains the voxel-wise \beta parameter values per each
% angular orientation (dim4), max. TE used for analysis (dim5), models
% applied (dim6) and for each beta value (3). This format is obtained from
% FittingEquations function and OC_AnalysisCode.m. If this is not used,
% please accomodate any input in this format. It is important that dim4
% here is equal in size than dim4 of Angle_Meas variable.
% Mask: Matrix with size (dim1,dim2,dim3)
% Binary mask to reduce unwanted voxels.
% Params: Structure with fields:
%   Params.nbins = number of desired bins.
%
% Output:
% XBin_mean: Matrix with size (nbins,filters)
% This variable contains the mean value of the bin (with a total of nbins
% defined in params.nbins) for each filter.
% YBin_mean: Matrix with size (nbins,max_TEs,models,betas,filters)
% This variable contains the mean value of the data in each defined bin 
% (with a total of nbins defined in params.nbins) and each filter. This is
% also applied for each max TE, model and beta values (the same size that
% the input data).
% XBin_std: Matrix with size (nbins,filters)
% This variable contains the std value of the bin (with a total of nbins
% defined in params.nbins) for each filter.
% YBin_std: Matrix with size (nbins,max_TEs,models,betas,filters)
% This variable contains the std value of the data in each defined bin 
% (with a total of nbins defined in params.nbins) and each filter. This is
% also applied for each max TE, model and beta values (the same size that
% the input data).
% Extra_Values: Structure parameter with fields:
%   Extra_Values.ang_aniso_bin: Cell variable which contains all the
%   angular values used for each specific bin and filter.
%   Extra_Values.kappa_aniso_bin: Cell variable which contains all the
%   \kappa values used for each specific bin and filter.

%% Start of the code
    % Masking all the input parameters:
    % Mask the OC not only with the manual masking but also for unrealistic
    % values in \kappa and ICVF.
    FullOC_indx = (Mask == 1 & Kappa_Meas > 0 & ICVF_Meas > 0 & ICVF_Meas < 1);
    
    % Reshape the angular maps to leave it in n-total voxels x m-measurements.
    theta_Z = reshape(Angle_Meas,[prod(size(Angle_Meas,1:3)),size(Angle_Meas,4)]);
    theta_Z_ROI = theta_Z(FullOC_indx(:),:);
    clear 'theta_Z'
    
    % Check the angular values and modify accordingly: negative values are
    % changed, angles with > pi/2 are reflected (pi - angle value) and
    % unrealistic values (angle < 0.5° or pi/360 rads) are converted to 0.
    theta_Z_ROI(theta_Z_ROI < 0) = abs(theta_Z_ROI(theta_Z_ROI < 0));
    theta_Z_ROI(theta_Z_ROI > pi/2) = pi - theta_Z_ROI(theta_Z_ROI > pi/2);
    theta_Z_ROI(theta_Z_ROI < 0) = abs(theta_Z_ROI(theta_Z_ROI < 0));
    theta_Z_ROI(theta_Z_ROI < pi/360) = 0;
    % If an specific voxel has the same angular value across ALL the
    % angular measurements, then it is expected that these voxels are
    % irrelevant for this study. Those voxels are also used to filter the
    % other maps.
    remove_indx = all(theta_Z_ROI(:,:) == theta_Z_ROI(:,1),2);
    theta_Z_ROI(remove_indx,:) = [];   

    Kvalues = Kappa_Meas(FullOC_indx(:));
    Kvalues(remove_indx) = [];

    Ficvfvalues = ICVF_Meas(FullOC_indx(:));
    Ficvfvalues(remove_indx) = [];

    beta_space = reshape(BetaParams,[prod(size(BetaParams,1:3)),size(BetaParams,4:7)]);
    % The beta_space is converted to i-total voxels (x-y-z axis), j-angular
    % measurements (16), k-max. TE's (5), l-models (4) and m-betas (3).
    beta_space = beta_space(FullOC_indx,:,:,:,:);
    beta_space(remove_indx,:,:,:,:) = []; 

    % To generate the cumulative sampling, it is necessary to repeat the
    % ICVF and \kappa maps equal to the number of angular measurements.
    ficvf_fullsampling_Z = repmat(Ficvfvalues,16,1);
    k_values_Z = repmat(Kvalues,16,1);
    ang_points_Z = theta_Z_ROI(:)*180/pi;

    % Anisotropic binning for equal number of voxels:
    % Five different binnings are performed: 
    % (1) All angles, no \kappa and not ICVF filtering. This is defined but
    % not used, so it can be removed.
    %h = histogram(ang_points_Z,0:0.1:90);
    %cumsum_h = cumsum(h.Values)./sum(h.Values);    
    % (2) All angles, ICVF > 0.8 and no \kappa filtering
    hficvf_0p8 = histogram(ang_points_Z(ficvf_fullsampling_Z > 0.8),0:0.1:90);
    cumsum_hficvf_0p8 = cumsum(hficvf_0p8.Values)./sum(hficvf_0p8.Values);
    % (3) All angles, ICVF > 0.8 and 0 < \kappa <= 1 filtering
    hficvf_0p8_k1 = histogram(ang_points_Z(ficvf_fullsampling_Z > 0.8 & k_values_Z < 1),0:0.1:90);
    cumsum_hficvf_0p8_k1 = cumsum(hficvf_0p8_k1.Values)./sum(hficvf_0p8_k1.Values);
    % (4) All angles, ICVF > 0.8 and 1 < \kappa < 2.5 filtering
    hficvf_0p8_k1_k2p5 = histogram(ang_points_Z(ficvf_fullsampling_Z > 0.8 & k_values_Z >= 1 & k_values_Z < 2.5),0:0.1:90);
    cumsum_hficvf_0p8_k1_k2p5 = cumsum(hficvf_0p8_k1_k2p5.Values)./sum(hficvf_0p8_k1_k2p5.Values);
    % (5) All angles, ICVF > 0.8 and 2.5 <= \kappa filtering (due to
    % distribution, max(\kappa) ~ 6.
    hficvf_0p8_k2p5 = histogram(ang_points_Z(ficvf_fullsampling_Z > 0.8 & k_values_Z >= 2.5),0:0.1:90);
    cumsum_hficvf_0p8_k2p5 = cumsum(hficvf_0p8_k2p5.Values)./sum(hficvf_0p8_k2p5.Values);
    
    % This is only for the plotting of Figure 5
    % (3) All angles, ICVF > 0.8 and 0 < \kappa < 1 filtering
    hist_k = histogram(ang_points_Z(ficvf_fullsampling_Z > 0.8 & k_values_Z < 1),0:2:90);
    HistogramRegularBin(1,:) = hist_k.BinEdges;
    HistogramRegularCount(1,:) = hist_k.BinCounts;
    % (4) All angles, ICVF > 0.8 and 1 <= \kappa < 2.5 filtering
    hist_k = histogram(ang_points_Z(ficvf_fullsampling_Z > 0.8 & k_values_Z >= 1 & k_values_Z < 2.5),0:2:90);
    HistogramRegularBin(2,:) = hist_k.BinEdges;
    HistogramRegularCount(2,:) = hist_k.BinCounts;
    % (5) All angles, ICVF > 0.8 and 2.5 <= \kappa filtering (due to
    % distribution, max(\kappa) ~ 6.
    hist_k = histogram(ang_points_Z(ficvf_fullsampling_Z > 0.8 & k_values_Z >= 2.5),0:2:90);
    HistogramRegularBin(3,:) = hist_k.BinEdges;
    HistogramRegularCount(3,:) = hist_k.BinCounts;
    
    HistInfo = {HistogramRegularBin; HistogramRegularCount};
    
    % Here the binning is made. The number of binnings are defined in the
    % structure params.nbins. A max index is searched and used to define
    % the lower edge value. The adding 0.05 is to define the "middle" value
    % of the bin. The use of 0.1 and 0.05 are arbitrary because the
    % histograms were sampled between 0 to 90° in steps of 0.1°.
    for anbin_indx = 1:params.nbins
        ang_NoRegularBin_ICVF0p8(anbin_indx) = max(find(cumsum_hficvf_0p8 <= 1/params.nbins*anbin_indx))*0.1 + 0.05;
        ang_NoRegularBin_ICVF0p8_K1(anbin_indx) = max(find(cumsum_hficvf_0p8_k1 <= 1/params.nbins*anbin_indx))*0.1 + 0.05;
        ang_NoRegularBin_ICVF0p8_K1_K2p5(anbin_indx) = max(find(cumsum_hficvf_0p8_k1_k2p5 <= 1/params.nbins*anbin_indx))*0.1 + 0.05;    
        ang_NoRegularBin_ICVF0p8_K2p5(anbin_indx) = max(find(cumsum_hficvf_0p8_k2p5 <= 1/params.nbins*anbin_indx))*0.1 + 0.05;
    end
    
    % The lowest edge value is added (0).
    ang_NoRegularBin_edge_ICVF0p8 = [0,ang_NoRegularBin_ICVF0p8];
    ang_NoRegularBin_edge_ICVF0p8_K1 = [0,ang_NoRegularBin_ICVF0p8_K1];
    ang_NoRegularBin_edge_ICVF0p8_K1_K2p5 = [0,ang_NoRegularBin_ICVF0p8_K1_K2p5];
    ang_NoRegularBin_edge_ICVF0p8_K2p5 = [0,ang_NoRegularBin_ICVF0p8_K2p5];
    
    % And here all the voxels are searched under the different filtering
    % approaches (as defined before).
    for ang_bin_indx = 1:length(ang_NoRegularBin_edge_ICVF0p8)-1
        indx_bin_ICVF0p8 = ang_points_Z >= ang_NoRegularBin_edge_ICVF0p8(ang_bin_indx) &...
                           ang_points_Z < ang_NoRegularBin_edge_ICVF0p8(ang_bin_indx+1) &...
                           ficvf_fullsampling_Z > 0.8;
        indx_bin_ICVF0p8_K1 = ang_points_Z >= ang_NoRegularBin_edge_ICVF0p8_K1(ang_bin_indx) &...
                              ang_points_Z < ang_NoRegularBin_edge_ICVF0p8_K1(ang_bin_indx+1) &...
                              ficvf_fullsampling_Z > 0.8 & k_values_Z < 1;
        indx_bin_ICVF0p8_K1_K2p5 = ang_points_Z >= ang_NoRegularBin_edge_ICVF0p8_K1_K2p5(ang_bin_indx) &...
                              ang_points_Z < ang_NoRegularBin_edge_ICVF0p8_K1_K2p5(ang_bin_indx+1) &...
                              ficvf_fullsampling_Z > 0.8 & k_values_Z >= 1 & k_values_Z < 2.5;
        indx_bin_ICVF0p8_K2p5 = ang_points_Z >= ang_NoRegularBin_edge_ICVF0p8_K2p5(ang_bin_indx) &...
                                ang_points_Z < ang_NoRegularBin_edge_ICVF0p8_K2p5(ang_bin_indx+1) & ...
                                ficvf_fullsampling_Z > 0.8 & k_values_Z >= 2.5;
        
        % Define mean value of each bin in the x-axis
        XBin_mean(ang_bin_indx,1) = mean(ang_points_Z(indx_bin_ICVF0p8));
        XBin_mean(ang_bin_indx,2) = mean(ang_points_Z(indx_bin_ICVF0p8_K1));
        XBin_mean(ang_bin_indx,3) = mean(ang_points_Z(indx_bin_ICVF0p8_K1_K2p5));
        XBin_mean(ang_bin_indx,4) = mean(ang_points_Z(indx_bin_ICVF0p8_K2p5));
        % Define std value of each bin in the x-axis
        XBin_std(ang_bin_indx,1) = std(ang_points_Z(indx_bin_ICVF0p8));
        XBin_std(ang_bin_indx,2) = std(ang_points_Z(indx_bin_ICVF0p8_K1));
        XBin_std(ang_bin_indx,3) = std(ang_points_Z(indx_bin_ICVF0p8_K1_K2p5));
        XBin_std(ang_bin_indx,4) = std(ang_points_Z(indx_bin_ICVF0p8_K2p5));
        
        % Define extra elements like the total amount of angles and \kappa
        % values selected per bin.
        Num_elements_ang(ang_bin_indx,1) = numel(ang_points_Z(indx_bin_ICVF0p8));
        Num_elements_ang(ang_bin_indx,2) = numel(ang_points_Z(indx_bin_ICVF0p8_K1));
        Num_elements_ang(ang_bin_indx,3) = numel(ang_points_Z(indx_bin_ICVF0p8_K1_K2p5));
        Num_elements_ang(ang_bin_indx,4) = numel(ang_points_Z(indx_bin_ICVF0p8_K2p5));

        %Extra_Values.kappa_aniso_bin{ang_bin_indx,1} = k_values_Z(indx_bin_ICVF0p8);
        %Extra_Values.kappa_aniso_bin{ang_bin_indx,2} = k_values_Z(indx_bin_ICVF0p8_K1);
        %Extra_Values.kappa_aniso_bin{ang_bin_indx,3} = k_values_Z(indx_bin_ICVF0p8_K1_K2p5);
        %Extra_Values.kappa_aniso_bin{ang_bin_indx,4} = k_values_Z(indx_bin_ICVF0p8_K2p5);
        
        % Reshape data to be used in the anisotropy binning.
        data_plot = reshape(beta_space, [prod(size(beta_space,1:2)),size(beta_space,3:5)]);
        % Define mean value of each bin in the y-axis
        YBin_mean(ang_bin_indx,:,:,:,1) = mean(data_plot(indx_bin_ICVF0p8,:,:,:));
        YBin_mean(ang_bin_indx,:,:,:,2) = mean(data_plot(indx_bin_ICVF0p8_K1,:,:,:));
        YBin_mean(ang_bin_indx,:,:,:,3) = mean(data_plot(indx_bin_ICVF0p8_K1_K2p5,:,:,:));
        YBin_mean(ang_bin_indx,:,:,:,4) = mean(data_plot(indx_bin_ICVF0p8_K2p5,:,:,:));
        % Define std value of each bin in the y-axis
        YBin_std(ang_bin_indx,:,:,:,1) = std(data_plot(indx_bin_ICVF0p8,:,:,:));
        YBin_std(ang_bin_indx,:,:,:,2) = std(data_plot(indx_bin_ICVF0p8_K1,:,:,:));
        YBin_std(ang_bin_indx,:,:,:,3) = std(data_plot(indx_bin_ICVF0p8_K1_K2p5,:,:,:));
        YBin_std(ang_bin_indx,:,:,:,4) = std(data_plot(indx_bin_ICVF0p8_K2p5,:,:,:));
    end  
        
    % Table creation here
    tabVarNames = ["Dispersion Range" "TEmax" "MinEdge_AngBin" "MaxEdge_AngBin",...
                   "Frequency_AngBin" "Mean_AngBin" "Sd_AngBin",...
                   "Alpha0_mean" "Alpha1_mean" "Alpha0_sd" "Alpha1_sd",...
                   "Beta0_mean" "Beta1_mean" "Beta2_mean" "Beta0_sd" "Beta1_sd" "Beta2_sd",...
                   "nRMSD_Alpha0" "nRMSD_Alpha1" "nRMSD_Beta0" "nRMSD_Beta1" "nRMSD_Beta2",...
                   "MWF_from_meanBeta1" "sdMWF_from_Beta1" "Gratio_from_meanBeta1" "sd_Gratio_from_Beta1"];
               
    DispersionTable = repmat(["Highly" "Mildly" "Negligible"],prod(size(YBin_mean,1:2)),1);
    TEmax_perDispersion = repmat(18:9:54, size(YBin_mean,1),1);
    TEmaxTable = repmat(TEmax_perDispersion(:),size(YBin_mean,5)-1,1);
    ang_NoRegularBin_edge_ICVF0p8_K1 = [0,ang_NoRegularBin_ICVF0p8_K1];
    ang_NoRegularBin_edge_ICVF0p8_K1_K2p5 = [0,ang_NoRegularBin_ICVF0p8_K1_K2p5];
    ang_NoRegularBin_edge_ICVF0p8_K2p5 = [0,ang_NoRegularBin_ICVF0p8_K2p5];
    
    MinEdge_Bin = repmat([ang_NoRegularBin_edge_ICVF0p8_K1(1:end-1);...
                          ang_NoRegularBin_edge_ICVF0p8_K1_K2p5(1:end-1);...
                          ang_NoRegularBin_edge_ICVF0p8_K2p5(1:end-1)]',size(YBin_mean,2),1);
    MaxEdge_Bin = repmat([ang_NoRegularBin_edge_ICVF0p8_K1(2:end);...
                          ang_NoRegularBin_edge_ICVF0p8_K1_K2p5(2:end);...
                          ang_NoRegularBin_edge_ICVF0p8_K2p5(2:end)]',size(YBin_mean,2),1);
    Frequency_Bin = repmat(Num_elements_ang(:,2:4),size(YBin_mean,2),1);
    
    MeanAngleTable = repmat(XBin_mean(:,2:4),size(YBin_mean,2),1);
    SdAngleTable = repmat(XBin_std(:,2:4),size(YBin_mean,2),1);
    
    ParamsMeanTable = reshape(permute(YBin_mean(:,:,:,:,2:4),[1,2,5,3,4]),[prod(size(YBin_mean,[1,2]))*3,size(YBin_mean,[3,4])]);
    ParamsSdTable = reshape(permute(YBin_std(:,:,:,:,2:4),[1,2,5,3,4]),[prod(size(YBin_mean,[1,2]))*3,size(YBin_mean,[3,4])]);
    
    % nRMSD estimation from the irregular binnings
    nRMSD_ExVivo = squeeze(sqrt(sum((YBin_mean(2:end,:,:,:,:) - YBin_mean(1,:,:,:,:)).^2)./(size(YBin_mean,1)-1))./YBin_mean(1,:,:,:,:));
    AuxnRMSDTable(1,:,:,:,:) = nRMSD_ExVivo(:,:,:,2:4);
    nRMSDTable = reshape(repmat(permute(AuxnRMSDTable,[1,2,5,3,4]),[size(YBin_mean,1),1,1,1,1]),[prod(size(YBin_mean,[1,2]))*3,size(YBin_mean,[3,4])]);
    
    % MWF and g-ratio estimations from the irregular binning
    MWF_MeanTable = (ParamsMeanTable(:,2,2) - 18.5)./(75.4 - 18.5);
    MWF_sdTable = 1./(75.4 - 18.5) .* ParamsSdTable(:,2,2);
    
    Gratio_MeanTable = sqrt(1 - (1./0.7).*(MWF_MeanTable)./(0.5.*(MWF_MeanTable.*(1./0.7-1)+1)));
    Gratio_sdTable = abs(1./Gratio_MeanTable).*(140/9./(MWF_MeanTable + 7/3).^2).^2.*MWF_sdTable;
    
    ExVivoTable = table(DispersionTable(:),TEmaxTable(:),MinEdge_Bin(:),MaxEdge_Bin(:),Frequency_Bin(:),MeanAngleTable(:),SdAngleTable(:),...
                        ParamsMeanTable(:,1,1),ParamsMeanTable(:,1,2),ParamsSdTable(:,1,1),ParamsSdTable(:,1,2),...
                        ParamsMeanTable(:,2,1),ParamsMeanTable(:,2,2),ParamsMeanTable(:,2,3),ParamsSdTable(:,2,1),...
                        ParamsSdTable(:,2,2),ParamsSdTable(:,2,3),nRMSDTable(:,1,1),nRMSDTable(:,1,2),...
                        nRMSDTable(:,2,1),nRMSDTable(:,2,2),nRMSDTable(:,2,3),MWF_MeanTable,MWF_sdTable,...
                        Gratio_MeanTable,Gratio_sdTable,'VariableNames',tabVarNames);
end