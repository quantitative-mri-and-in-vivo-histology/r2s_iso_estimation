%% Main code for Figures.
% This code creates some of the Figures entirely (e.g. Figures 4 to 12) or part of them,
% which will require further "edition" (e.g. Figures 3B-C). Even though
% some figures are obtained directly from Matlab, it required further
% manual arrangement (e.g. moving legends or figure position arrangement). 
clear all; clc; close all;
addpath(genpath('Colormaps_ColorBlind'));

%Load data for all figures
load(fullfile(fileparts(pwd),'InSilico','Directions1500_Cartesian.mat')); % For Figure 4

% The following maps for Figure 5
%KappaMap = niftiread(fullfile(fileparts(pwd),'ExVivo','Datasets','example_kappa_maskedx10_x_manual_adapted.nii'));
%ICVFMap = niftiread(fullfile(fileparts(pwd),'ExVivo','Datasets','example_ficvf_masked_x_manual_adapted.nii'));
%AngleMaps = load(fullfile(fileparts(pwd),'ExVivo','Datasets','OC_AnglesLoaded.mat'));
%MaskOC = double(niftiread(fullfile(fileparts(pwd),'ExVivo','Datasets','Mask_OC.nii')));

% This dataset is used for all Figures except 4.
date_exvivo = '230104'; %Change this date accordingly.
load(fullfile(fileparts(pwd),'ExVivo','Results',[date_exvivo '_M1M2_SingleMeas'],'ExVivoResults_IrregularBinning_NotNormalised.mat'));

% This dataset is used for Figure 6 onwards.
date_insilico = '221229'; %Change this date accordingly.
load(fullfile(fileparts(pwd),'InSilico',['Myelin_NoMyelin_AnisotropicBinarised_InSilicoMR_' date_insilico '.mat']));

%% Figure 4 (Cylinders space and cylinders weighted by the Watson distribution).
% Note: It needs to "re-arrange" vertically.
figure(1);
[phi_axons, theta_axons, ~] = cart2sph(DAxons1_Norm(:,1),DAxons1_Norm(:,2),DAxons1_Norm(:,3));
theta_axons = pi/2 - theta_axons; % Matlab defines pi/2 in the +z axis.
kappa_values = [0.001,1,2.5,4.25,6];
Watson = @(k, mean_theta, angle_theta, mean_phi, angle_phi) (hypergeom(1/2,3/2,k)).^(-1) .* exp(k .* (sin(mean_theta).*sin(angle_theta).*cos(mean_phi - angle_phi) + cos(mean_theta).*cos(angle_theta)).^2);
P1 = [0,0,0];
P2 = [0,0,1.5];

subplot(2,5,[1,2,6,7]), scatter3(DAxons1_Norm(:,1),DAxons1_Norm(:,2),DAxons1_Norm(:,3),...
                           40*ones(size(DAxons1_Norm,1),1),'filled');hold on;
subplot(2,5,[1,2,6,7]), quiver3(0,0,0,1/sqrt(3),-1/sqrt(3),1/sqrt(3),'Color',[0.494,0.184,0.556],'LineWidth',5);
subplot(2,5,[1,2,6,7]), quiver3(0,0,0,cos(20/180*pi),sin(20/180*pi),1.5,'Color',[0.466,0.674,0.188],'LineWidth',5);
subplot(2,5,[1,2,6,7]), zlim([-1,1]);
subplot(2,5,[1,2,6,7]), quiver3(-0.9,0.9,0,0,0,1,'Color',[0.929,0.694,0.125],'LineWidth',5);

indx_subplot = [3:5,8,9];

for kappa_index = 1:length(kappa_values)
    Watson_weight = Watson(kappa_values(kappa_index),0,theta_axons,0,phi_axons);
    subplot(2,5,indx_subplot(kappa_index)),scatter3(DAxons1_Norm(:,1),DAxons1_Norm(:,2),DAxons1_Norm(:,3),...
                                      200*ones(size(DAxons1_Norm,1),1),Watson_weight/max(Watson_weight),'filled');hold on;
    caxis([0,1]);
    subplot(2,5,indx_subplot(kappa_index)), quiver3(0,0,1,0,0,1,'Color',[0,0,0],'LineWidth',5);
    axis off
    subplot(2,5,indx_subplot(kappa_index)),title(['\kappa: ' num2str(kappa_values(kappa_index))],'FontSize',24);
    colormap('viridis');
end

%% Figure 5 - Irregular binning and kappa mapping.
figure('DefaultAxesFontSize',18,'DefaultAxesFontWeight','Bold');
% colormap("gray");
% subplot(2,4,1), imagesc(KappaMap(:,:,62).*MaskOC(:,:,62),[0,5]);
% subplot(2,4,1), axis off;
% 
% DispersionMask = MaskOC .* (ICVFMap >= 0.8 & ICVFMap < 1).*...
%                             (logical(KappaMap < 1 & KappaMap > 0) + ...
%                             2*logical(KappaMap >= 1 & KappaMap < 2.5) + ...
%                             3*logical(KappaMap >= 2.5)) - (MaskOC==0);
% 
% ax1 = subplot(2,4,5);
% im1 = imagesc(KappaMap(:,:,62)'.*MaskOC(:,:,62)');
% im1.AlphaData = 1;
% hold all;
% colormap(ax1,"gray");
% 
% ax2 = subplot(2,4,5);
% im2 = imagesc(DispersionMask(:,:,62)');
% im2.AlphaData = 0.1;
% hold all;
% ax2.Visible = 'off';
% ax2.XTick = [];
% ax2.YTick = [];
% Background = [1 1 1];
% RedHighDispersed = [0.635 0.078 0.184];
% GreenMildlyDisp = [0.466 0.674 0.188];
% BlueNegliDisp = [0 0.447 0.741];
% colormap(ax2,[Background;RedHighDispersed;GreenMildlyDisp;BlueNegliDisp]);
% set(ax2,'CLim',[-1,3]);
% ax2.UserData = linkprop([ax1,ax2],{'Position','InnerPosition','DataAspectRatio',...
%                                    'xtick','ytick','ydir','xdir'});
for k_indx = 1:3
    switch k_indx 
        case 1
            kappa_str = 'Highly';
        case 2
            kappa_str = 'Mildly';
        case 3
            kappa_str = 'Negligible';
    end

    TableExVivoConditions = ExVivoTableSummary.("Dispersion Range") == kappa_str & ...
                            ExVivoTableSummary.("TEmax") == 54;

    subplot(2,4,[1+k_indx,5+k_indx]), histogram('BinEdges',HistogramFigure5{1}(k_indx,:),'BinCounts',HistogramFigure5{2}(k_indx,:)); hold on;
    subplot(2,4,[1+k_indx,5+k_indx]), histogram('BinEdges',[ExVivoTableSummary.MinEdge_AngBin(TableExVivoConditions);90]','BinCounts',[ExVivoTableSummary.Frequency_AngBin(TableExVivoConditions)]'); grid minor; box off;
end
subplot(2,4,[2,6]), ylabel({'Number of voxels across all measurements','(cumulated data) [frequency]'}); title('\kappa < 1');
subplot(2,4,[3,7]), xlabel('Angular orientation (\theta_\mu) [degrees]'); title('1 \leq \kappa < 2.5');
subplot(2,4,[4,8]), title('2.5 \leq \kappa');


%% Figures 6 to 11 - Plots
% Figure 6
figure('DefaultAxesFontSize',20,'DefaultAxesFontWeight','Bold');
for k_indx = 1:3
    switch k_indx
        case 1
            kappa_str = 'Highly';
        case 2
            kappa_str = 'Mildly';
        case 3
            kappa_str = 'Negligible';
    end

    TableExVivoConditions = ExVivoTableSummary.("Dispersion Range") == kappa_str & ...
                            ExVivoTableSummary.("TEmax") == 54;

    for g_indx = 1:3
        switch g_indx
            case 1
                gratio_val = 0.66;
            case 2
                gratio_val = 0.73;
            case 3
                gratio_val = 0.8;
        end

        TableConditions = Myelin_InSilicoResults.DataType == 'Mye' & ...
                          Myelin_InSilicoResults.Dispersion_Range == kappa_str & ...
                          Myelin_InSilicoResults.Gratio == gratio_val & ...
                          Myelin_InSilicoResults.TEmax == 54 & ...
                          Myelin_InSilicoResults.T2ValueReference == 'T2ExVivo';
        
        InSilicoDataType_mAlpha1 = Myelin_InSilicoResults.Alpha1_mean(TableConditions,:);
        InSilicoDataType_sdAlpha1 = Myelin_InSilicoResults.Alpha1_sd(TableConditions,:);
        InSilicoDataType_mBeta1 = Myelin_InSilicoResults.Beta1_mean(TableConditions,:);
        InSilicoDataType_sdBeta1 = Myelin_InSilicoResults.Beta1_sd(TableConditions,:);

        subplot(2,3,k_indx),errorbar(ExVivoTableSummary.Mean_AngBin(TableExVivoConditions),InSilicoDataType_mAlpha1,...
                                     InSilicoDataType_sdAlpha1,'Color',ColorRGB_ForPlot(g_indx),'LineWidth',2.5,'Marker','x'); ...
                                     xlim([-5,95]); ylim([10,60]); hold on;
        subplot(2,3,3 + k_indx),errorbar(ExVivoTableSummary.Mean_AngBin(TableExVivoConditions),InSilicoDataType_mBeta1,...
                                        InSilicoDataType_sdBeta1,'Color',ColorRGB_ForPlot(g_indx),'LineWidth',2.5,'Marker','x'); ...
                                        xlim([-5,95]); ylim([10,60]); hold on;
    end
    
    subplot(2,3,k_indx),errorbar(ExVivoTableSummary.Mean_AngBin(TableExVivoConditions),ExVivoTableSummary.Alpha1_mean(TableExVivoConditions),...
                                 ExVivoTableSummary.Alpha1_sd(TableExVivoConditions),'Color',ColorRGB_ForPlot(6),'LineWidth',2.5,'Marker','x'); hold on;...
                                 grid minor; box off; xlim([-5,95]); ylim([10,60]);
    subplot(2,3,3 + k_indx),errorbar(ExVivoTableSummary.Mean_AngBin(TableExVivoConditions),ExVivoTableSummary.Beta1_mean(TableExVivoConditions),...
                                     ExVivoTableSummary.Beta1_sd(TableExVivoConditions),'Color',ColorRGB_ForPlot(6),'LineWidth',2.5,'Marker','x'); hold on;...
                                     grid minor; box off; xlim([-5,95]); ylim([10,60]); 
     
end
subplot(2,3,1), title({'Highly dispersed fibres','(\kappa < 1)'});
subplot(2,3,2), title({'Mildly dispersed fibres','(1 \leq \kappa < 2.5)'});
subplot(2,3,3), title({'Negligibly dispersed fibres','(2.5 \leq \kappa)'});
subplot(2,3,1), ylabel({'\alpha_1 parameter [1/s]'});
subplot(2,3,4), ylabel('\beta_1 parameter [1/s]');
subplot(2,3,5), xlabel('Angular orientation (\theta_\mu) [degrees]');

%% Figure 7
figure('DefaultAxesFontSize',20,'DefaultAxesFontWeight','Bold');
for k_indx = 1:3
    switch k_indx
        case 1
            kappa_str = 'Highly';
        case 2
            kappa_str = 'Mildly';
        case 3
            kappa_str = 'Negligible';
    end

    TableExVivoConditions = ExVivoTableSummary.("Dispersion Range") == kappa_str & ...
                            ExVivoTableSummary.("TEmax") == 54;
    
    %-FJv23(16.01): Since it is the same value, taking the mean is the same
    %that taking the first index.
    nRMSD_ExVivo_Alpha1(k_indx) = mean(ExVivoTableSummary.nRMSD_Alpha1(TableExVivoConditions,:));
    nRMSD_ExVivo_Beta1(k_indx) = mean(ExVivoTableSummary.nRMSD_Beta1(TableExVivoConditions,:));

    for g_indx = 1:3
        switch g_indx
            case 1
                gratio_val = 0.66;
            case 2
                gratio_val = 0.73;
            case 3
                gratio_val = 0.8;
        end

        TableConditions = Myelin_InSilicoResults.DataType == 'Mye' & ...
                          Myelin_InSilicoResults.Dispersion_Range == kappa_str & ...
                          Myelin_InSilicoResults.Gratio == gratio_val & ...
                          Myelin_InSilicoResults.TEmax == 54 & ...
                          Myelin_InSilicoResults.T2ValueReference == 'T2ExVivo';
        
        InSilicoDataType_mAlpha1 = Myelin_InSilicoResults.Alpha1_mean(TableConditions,:);
        %InSilicoDataType_sdAlpha1 = Myelin_InSilicoResults.Alpha1_sd(TableConditions,:);
        InSilicoDataType_mBeta1 = Myelin_InSilicoResults.Beta1_mean(TableConditions,:);
        %InSilicoDataType_sdBeta1 = Myelin_InSilicoResults.Beta1_sd(TableConditions,:);

        nRMSD_Alpha1(g_indx,k_indx) = sqrt(sum((InSilicoDataType_mAlpha1(2:end)-InSilicoDataType_mAlpha1(1)).^2)/19)./InSilicoDataType_mAlpha1(1);
        nRMSD_Beta1(g_indx,k_indx) = sqrt(sum((InSilicoDataType_mBeta1(2:end)-InSilicoDataType_mBeta1(1)).^2)/19)./InSilicoDataType_mBeta1(1);
    end


end
delta_nRMSD_InSilico = nRMSD_Beta1 - nRMSD_Alpha1;

delta_nRMSD_ExVivo = nRMSD_ExVivo_Beta1 - nRMSD_ExVivo_Alpha1;

plot1 = bar(subplot(1,3,1),[nRMSD_Alpha1;nRMSD_ExVivo_Alpha1]'*100,'grouped','FaceColor','flat');
plot2 = bar(subplot(1,3,2),[nRMSD_Beta1;nRMSD_ExVivo_Beta1]'*100,'grouped','FaceColor','flat');
plot3 = bar(subplot(1,3,3),[delta_nRMSD_InSilico;squeeze(delta_nRMSD_ExVivo)]'*100,'grouped','FaceColor','flat');

for group_indx = 1:4
    Color_bar = logical(group_indx < 4)*ColorRGB_ForPlot(group_indx) + logical(group_indx == 4)*ColorRGB_ForPlot(6);
    for bar_indx = 1:3
        plot1(group_indx).CData(bar_indx,:) = Color_bar;
        plot2(group_indx).CData(bar_indx,:) = Color_bar;
        plot3(group_indx).CData(bar_indx,:) = Color_bar;
    end
end

subplot(1,3,1), ylabel('nRMSD \alpha_1 [%]'); 
subplot(1,3,2), ylabel('nRMSD \beta_1 [%]'); 
subplot(1,3,3), ylabel('\Delta nRMSD [%-points]'); 

for plot_indx = 1:3
    subplot(1,3,plot_indx), xlabel('Dispersed fibres'); xticklabels({'Highly','Mildly','Negligibly'}); xtickangle(45);...
                            ylim(logical(plot_indx < 3)*[0 55] + logical(plot_indx == 3)*[-55 0]); grid minor; box off; hold on;
end

%% Figure 8
gratio_truth = [30.8, 28.5, 26.1];
%Small calculation outside, it is easy to plot afterwards and get some
%values ;)
TableConditions_Rows = Myelin_InSilicoResults.DataType == 'Mye' & ...
                       Myelin_InSilicoResults.TEmax == 54 & ...
                       Myelin_InSilicoResults.T2ValueReference == 'T2ExVivo';
% The resulting array is ordered first across the different dispersion ranges,
% then across g-ratios. That results in 20 x 3 x 3 = 180 values.
InSilicoDataType_mBeta1 = reshape(Myelin_InSilicoResults.Beta1_mean(TableConditions_Rows,:),20,3,3);
InSilicoDataType_sdBeta1 = reshape(Myelin_InSilicoResults.Beta1_sd(TableConditions_Rows,:),20,3,3);

% The "single" \beta1 truth (\beta1,nm) is 18.5.
epsilon_nm = (1 - squeeze(InSilicoDataType_mBeta1/18.5))*100;
sd_eps_nm = 1/18.5*squeeze(InSilicoDataType_sdBeta1)*100;

% The "variable" \beta1 truth (\beta1,m) as a function of g-ratio are 30.8,
% 28.5 and 26.1.
beta1m = reshape(repmat([30.8, 28.5, 26.1],60,1),20,3,3);
epsilon_m = (1 - InSilicoDataType_mBeta1./beta1m)*100;
sd_eps_m = 1./beta1m.*InSilicoDataType_sdBeta1*100;    

figure('DefaultAxesFontSize',20,'DefaultAxesFontWeight','Bold');
for k_indx = 1:3
    switch k_indx
        case 1
            kappa_str = 'Highly';
        case 2
            kappa_str = 'Mildly';
        case 3
            kappa_str = 'Negligible';
    end

    TableExVivoConditions = ExVivoTableSummary.("Dispersion Range") == kappa_str & ...
                            ExVivoTableSummary.("TEmax") == 54;
    for g_indx = 1:3
        subplot(2,4,k_indx),plot(ExVivoTableSummary.Mean_AngBin(TableExVivoConditions),epsilon_nm(:,k_indx,g_indx),...
                                 'Color',ColorRGB_ForPlot(g_indx),'LineWidth',2.5,'Marker','o'); ...
                                 xlim([-5,95]); ylim([-100,22.5]); hold on; grid minor; box off;
        subplot(2,4,4 + k_indx),plot(ExVivoTableSummary.Mean_AngBin(TableExVivoConditions),epsilon_m(:,k_indx,g_indx),...
                                     'Color',ColorRGB_ForPlot(g_indx),'LineWidth',2.5,'Marker','o'); ...
                                     xlim([-5,95]); ylim([-100,22.5]); hold on; grid minor; box off;
    end  
end
subplot(2,4,1), title({'Highly dispersed fibres','(\kappa < 1)'});
subplot(2,4,2), title({'Mildly dispersed fibres','(1 \leq \kappa < 2.5)'});
subplot(2,4,3), title({'Negligibly dispersed fibres','(2.5 \leq \kappa)'});
subplot(2,4,1), ylabel({'\epsilon_{nm} estimation [%]'});
subplot(2,4,5), ylabel('\epsilon_{m} estimation [%]');
subplot(2,4,6), xlabel('Angular orientation (\theta_\mu) [degrees]');

%For error bar plot
mean_epsilon_nm = squeeze(mean(epsilon_nm,1));
sd_epsilon_nm = squeeze(std(epsilon_nm,[],1));

mean_epsilon_m = squeeze(mean(epsilon_m,1));
sd_epsilon_m = squeeze(std(epsilon_m,[],1));

model_series = [mean_epsilon_nm mean_epsilon_m]; 
b = bar(subplot(2,4,[4,8]),model_series, 'grouped','FaceColor','flat');

for group_indx = 1:6
    Color_bar = logical(group_indx < 4)*ColorRGB_ForPlot(group_indx) + logical(group_indx >= 4)*ColorRGB_ForPlot(group_indx+3);
    for bar_indx = 1:3
        b(group_indx).CData(bar_indx,:) = Color_bar;
    end
end

hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
%In this point, resize manually this plot to make the x-axis label to
%appear
subplot(2,4,[4,8]), errorbar(x',model_series,[sd_epsilon_nm sd_epsilon_m],'k','linestyle','none','linewidth',2); grid minor;
subplot(2,4,[4,8]), xticklabels({'Highly','Mildly','Negligibly'}); xtickangle(45); xlabel('Dispersed fibres');
subplot(2,4,[4,8]), ylabel('Mean (<\epsilon>) and sd(\epsilon) estimation [%]');
hold off

%% Figure 9
%MWF_fun = @(x) (x-18.5)./(75.4-18.5);
%MWF_exvivo = MWF_fun(squeeze(AnisoBinY_mean(:,end,2,2,2:4)));

figure('DefaultAxesFontSize',20,'DefaultAxesFontWeight','Bold');
for k_indx = 1:3
    switch k_indx
        case 1
            kappa_str = 'Highly';
        case 2
            kappa_str = 'Mildly';
        case 3
            kappa_str = 'Negligible';
    end

    TableExVivoConditions = ExVivoTableSummary.("Dispersion Range") == kappa_str & ...
                            ExVivoTableSummary.("TEmax") == 54;
    
    subplot(2,3,[1,2,4,5]),plot(ExVivoTableSummary.Mean_AngBin(TableExVivoConditions,:),...
                                ExVivoTableSummary.MWF_from_meanBeta1(TableExVivoConditions,:),'LineWidth',2,'Marker','o','Color',ColorRGB_ForPlot(k_indx+3)); hold on; grid minor;
    MWF_exvivo(:,k_indx) = ExVivoTableSummary.MWF_from_meanBeta1(TableExVivoConditions,:);
end


subplot(2,3,[1,2,4,5]), ylabel('Estimated myelin water fraction (MWF) [n.u.]');
subplot(2,3,[1,2,4,5]), xlabel('Angular orientation (\theta_\mu) [degrees]'); xlim([5,95]);
subplot(2,3,[3,6]), boxplot(MWF_exvivo), ylabel('Median and sd. of the estimated MWF [n.u.]'); xlabel('Dispersed fibres'); ...
                    xticklabels({'Highly','Mildly','Negligibly'}); xtickangle(45); hold on; grid minor;

%% Figure 10 - TE var for \beta_1 M1 and \beta_1 M2
TableConditions_Rows = Myelin_InSilicoResults.DataType == 'Mye' & ...
                       Myelin_InSilicoResults.Gratio == 0.8 & ...
                       Myelin_InSilicoResults.T2ValueReference == 'T2ExVivo' & ...
                       Myelin_InSilicoResults.Dispersion_Range == 'Negligible';

% The resulting array is ordered from the lowest to highest TEmax. Since it
% is 5 TEmax, then the total data is 20 x 5 = 100.
InSilicoDataType_mAlpha1 = reshape(Myelin_InSilicoResults.Alpha1_mean(TableConditions_Rows,:),20,5);
InSilicoDataType_sdAlpha1 = reshape(Myelin_InSilicoResults.Alpha1_sd(TableConditions_Rows,:),20,5);

InSilicoDataType_mBeta1 = reshape(Myelin_InSilicoResults.Beta1_mean(TableConditions_Rows,:),20,5);
InSilicoDataType_sdBeta1 = reshape(Myelin_InSilicoResults.Beta1_sd(TableConditions_Rows,:),20,5);

figure('DefaultAxesFontSize',20,'DefaultAxesFontWeight','Bold');
for te_indx = 1:3
    TableExVivoConditions = ExVivoTableSummary.("Dispersion Range") == 'Negligible' & ...
                            ExVivoTableSummary.("TEmax") == te_indx*18;

    subplot(2,3,4 - te_indx),errorbar(ExVivoTableSummary.Mean_AngBin(TableExVivoConditions,:),squeeze(InSilicoDataType_mAlpha1(:,2*te_indx-1)),...
                                    squeeze(InSilicoDataType_sdAlpha1(:,2*te_indx-1)),'Color',ColorRGB_ForPlot(3),'LineWidth',2.5,'Marker','x'); hold on;...
                                    xlim([-5,95]); ylim([20,45]);
    subplot(2,3,4 - te_indx),errorbar(ExVivoTableSummary.Mean_AngBin(TableExVivoConditions,:),ExVivoTableSummary.Alpha1_mean(TableExVivoConditions),...
                                    ExVivoTableSummary.Alpha1_sd(TableExVivoConditions),'Color',ColorRGB_ForPlot(6),'LineWidth',2.5,'Marker','x'); hold on;...
                                    grid minor; box off; xlim([-5,95]); ylim([20,45]);
    subplot(2,3,7 - te_indx),errorbar(ExVivoTableSummary.Mean_AngBin(TableExVivoConditions,:),squeeze(InSilicoDataType_mBeta1(:,2*te_indx-1)),...
                                    squeeze(InSilicoDataType_sdBeta1(:,2*te_indx-1)),'Color',ColorRGB_ForPlot(3),'LineWidth',2.5,'Marker','x'); hold on;...
                                    xlim([-5,95]); ylim([20,45]);        
    subplot(2,3,7 - te_indx),errorbar(ExVivoTableSummary.Mean_AngBin(TableExVivoConditions,:),ExVivoTableSummary.Beta1_mean(TableExVivoConditions),...
                                    ExVivoTableSummary.Beta1_sd(TableExVivoConditions),'Color',ColorRGB_ForPlot(6),'LineWidth',2.5,'Marker','x'); hold on;...
                                    grid minor; box off; xlim([-5,95]); ylim([20,45]);        
end
subplot(2,3,1), title('TE_{max}: 54 ms');
subplot(2,3,2), title('TE_{max}: 36 ms');
subplot(2,3,3), title('TE_{max}: 18 ms');
subplot(2,3,3), xline(55,'--','Linewidth',2.5,'Color', [0.635 0.078 0.184]);
subplot(2,3,6), xline(55,'--','Linewidth',2.5,'Color', [0.635 0.078 0.184]);
subplot(2,3,1), ylabel({'\alpha_1 parameter [1/s]'});
subplot(2,3,4), ylabel('\beta_1 parameter [1/s]');
subplot(2,3,5), xlabel('Angular orientation (\theta_\mu) [degrees]');

%% Figure 11
figure('DefaultAxesFontSize',20,'DefaultAxesFontWeight','Bold');

for te_indx = 1:3
    TEvalue = 54 - (te_indx - 1)*18;
    
    TableExVivoConditions = ExVivoTableSummary.("Dispersion Range") == 'Negligible' & ...
                            ExVivoTableSummary.("TEmax") == TEvalue;
    
    %-FJv23(16.01): Since it is the same value, taking the mean is the same
    %that taking the first index.
    nRMSD_ExVivo_TEmax_Alpha1(1,te_indx) = mean(ExVivoTableSummary.nRMSD_Alpha1(TableExVivoConditions,:));
    nRMSD_ExVivo_TEmax_Beta1(1,te_indx) = mean(ExVivoTableSummary.nRMSD_Beta1(TableExVivoConditions,:));

    TableConditions = Myelin_InSilicoResults.DataType == 'Mye' & ...
                      Myelin_InSilicoResults.Dispersion_Range == 'Negligible' & ...
                      Myelin_InSilicoResults.Gratio == 0.8 & ...
                      Myelin_InSilicoResults.TEmax == TEvalue & ...
                      Myelin_InSilicoResults.T2ValueReference == 'T2ExVivo';
        
    InSilicoDataType_mAlpha1 = Myelin_InSilicoResults.Alpha1_mean(TableConditions,:);
    InSilicoDataType_sdAlpha1 = Myelin_InSilicoResults.Alpha1_sd(TableConditions,:);
    InSilicoDataType_mBeta1 = Myelin_InSilicoResults.Beta1_mean(TableConditions,:);
    InSilicoDataType_sdBeta1 = Myelin_InSilicoResults.Beta1_sd(TableConditions,:);

    nRMSD_TEmax_Alpha1(1,te_indx) = sqrt(sum((InSilicoDataType_mAlpha1(2:end)-InSilicoDataType_mAlpha1(1)).^2)/19)./InSilicoDataType_mAlpha1(1);
    nRMSD_TEmax_Beta1(1,te_indx) = sqrt(sum((InSilicoDataType_mBeta1(2:end)-InSilicoDataType_mBeta1(1)).^2)/19)./InSilicoDataType_mBeta1(1);
end

delta_nRMSD_InSilico = nRMSD_TEmax_Beta1 - nRMSD_TEmax_Alpha1;
delta_nRMSD_ExVivo = nRMSD_ExVivo_TEmax_Beta1 - nRMSD_ExVivo_TEmax_Alpha1;

plot1 = bar(subplot(1,3,1),[nRMSD_TEmax_Alpha1;nRMSD_ExVivo_TEmax_Alpha1]'*100,'grouped','FaceColor','flat');
plot2 = bar(subplot(1,3,2),[nRMSD_TEmax_Beta1;nRMSD_ExVivo_TEmax_Beta1]'*100,'grouped','FaceColor','flat');
plot3 = bar(subplot(1,3,3),[delta_nRMSD_InSilico;squeeze(delta_nRMSD_ExVivo)]'*100,'grouped','FaceColor','flat');

for group_indx = 1:2
    Color_bar = ColorRGB_ForPlot(3*group_indx);
    for bar_indx = 1:3
        plot1(group_indx).CData(bar_indx,:) = Color_bar;
        plot2(group_indx).CData(bar_indx,:) = Color_bar;
        plot3(group_indx).CData(bar_indx,:) = Color_bar;
    end
end

subplot(1,3,1), ylabel('nRMSD \alpha_1 [%]'); 
subplot(1,3,2), ylabel('nRMSD \beta_1 [%]'); 
subplot(1,3,3), ylabel('\Delta nRMSD [%-points]'); 

for plot_indx = 1:3
    subplot(1,3,plot_indx), xlabel('TE_{max}'); xticklabels({'54 ms','36 ms','18 ms'}); xtickangle(45);...
                            ylim(logical(plot_indx < 3)*[0 50] + logical(plot_indx == 3)*[-40 40]); grid minor; box off; hold on;
end

%% Figure 12
% The mean wAICc for ex vivo has to be done here sadly. So it requires a
% little bit of masking since wAICc was determined per voxel/angular
% orientation/max TE.
figure('DefaultAxesFontSize',20,'DefaultAxesFontWeight','Bold');

for te_indx = 1:3
    TEvalue = 54 - (te_indx - 1)*18;
    
    TableConditions = Myelin_InSilicoResults.DataType == 'Mye' & ...
                      Myelin_InSilicoResults.Dispersion_Range == 'Negligible' & ...
                      Myelin_InSilicoResults.Gratio == 0.8 & ...
                      Myelin_InSilicoResults.TEmax == TEvalue & ...
                      Myelin_InSilicoResults.T2ValueReference == 'T2ExVivo';
        
    InSilicoDataType_wAICc(:,te_indx) = Myelin_InSilicoResults.wAICc_M2_mean_uniform(TableConditions,:);
    InSilicoDataType_sdwAICc(:,te_indx) = Myelin_InSilicoResults.wAICc_M2_sd_uniform(TableConditions,:);

    TableExVivoConditions = ExVivoTableSummary.("Dispersion Range") == 'Negligible' & ...
                            ExVivoTableSummary.("TEmax") == TEvalue;
    
    %-FJv23(16.01): Since it is the same value, taking the mean is the same
    %that taking the first index.
    ExVivo_mean_wAICc(1,te_indx) = mean(ExVivoTableSummary.wAICc_mean_M2(TableExVivoConditions,:));
    ExVivo_sd_wAICc(1,te_indx) = mean(ExVivoTableSummary.wAICc_sd_M2(TableExVivoConditions,:));

end

plot1 = bar([InSilicoDataType_wAICc(1,:);ExVivo_mean_wAICc]','grouped','FaceColor','flat');

for group_indx = 1:2
    Color_bar = ColorRGB_ForPlot(3*group_indx);
    for bar_indx = 1:3
        plot1(group_indx).CData(bar_indx,:) = Color_bar;
    end
end

hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size([InSilicoDataType_wAICc(1,:);ExVivo_mean_wAICc]');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = plot1(i).XEndPoints;
end
% Plot the errorbars
%In this point, resize manually this plot to make the x-axis label to
%appear
errorbar(x',[InSilicoDataType_wAICc(1,:);ExVivo_mean_wAICc]',[InSilicoDataType_sdwAICc(1,:);ExVivo_sd_wAICc]','k','linestyle','none','linewidth',2); grid minor;
yline(0.5,'--','Linewidth',2.5,'Color', [0.85 0.325 0.098]);
yline(0.75,'--','Linewidth',2.5,'Color', [0.494 0.184 0.556]);
xlabel('Maximum echo time (TE_{max})'); xticklabels({'54 ms','36 ms','18 ms'}); xtickangle(45);
ylabel('Average of the wAICc for M2 [n.u.]'); ylim([-0.2 1.2]);

%% Supplementary Figure 1
% This is ONLY executed when data is available.
% meas_indx_pos = [1,11,12,10,8,7,6,13];
% t = tiledlayout(3,8,'TileSpacing','Compact','Padding','Compact');
% 
% for ang_indx = 1:8
%     meas_indx = meas_indx_pos(ang_indx);
%     AngleMap = load_untouch_nii(fullfile(fileparts(pwd),'Dataset','AngularOrientation',...
%                          ['ThetatoB0Meas' num2str(meas_indx) '_DiffSpace_Zdir_0toPiHalf_x_manual_adapted.nii']));  
%     nexttile
%     imagesc(imrotate(squeeze(AngleMap.img(35:85,25:75,62)*180/pi .* (FICV_map(35:85,25:75,62) > 0.8)),90),[0,90]); colormap('parula'); axis off;
% end
% for ang_indx = 1:8
%     meas_indx = meas_indx_pos(ang_indx);
%     BetaImage = load_untouch_nii(fullfile(results_path,['Model_' num2str(1)],'NonNormalisedFitting',...
%                          ['NoReg_Meas' num2str(meas_indx) '_Beta' num2str(1) '_TEmax' num2str(5) '.nii']));  
%     nexttile
%     imagesc(imrotate(squeeze(BetaImage.img(35:85,25:75,62) .* (FICV_map(35:85,25:75,62) > 0.8)),90),[10,40]); colormap('gray'); axis off;
% end
% for ang_indx = 1:8
%     meas_indx = meas_indx_pos(ang_indx);
%     BetaImage2 = load_untouch_nii(fullfile(results_path,['Model_' num2str(2)],'NonNormalisedFitting',...
%                          ['NoReg_Meas' num2str(meas_indx) '_Beta' num2str(1) '_TEmax' num2str(5) '.nii']));  
%     nexttile
%     imagesc(imrotate(squeeze(BetaImage2.img(35:85,25:75,62) .* (FICV_map(35:85,25:75,62) > 0.8)),90),[10,30]); colormap('gray'); axis off;
% end
