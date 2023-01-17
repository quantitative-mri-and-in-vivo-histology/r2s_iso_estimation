%% OC_AnalysisCode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code analyses the optic chiasm (OC) data 4D_X27P.nii R2* weighted
% data. This dataset has 16 angular measurements and 16 echo times per each
% angular measurement. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;

% Add required paths and load the required volumes and variables for analysis
% All those folders are located in the same directory this code is
% executed.
addpath(fullfile(pwd,'spm12','spm12'));
addpath(genpath(fullfile(pwd,'spm12','spm12','toolbox','ACID')));
addpath(genpath(fullfile(pwd,'Nifti_tool')));
addpath(fullfile(pwd,'Functions'));

dataset_path = fullfile(fileparts(pwd),'ExVivo','Dataset');
if ~exist(dataset_path,'dir')
    error('Folder Dataset does not exist, but must exist and contain all the required ex vivo information');
end

% Define paths for results 
%date_string = datestr(today,'yymmdd');
date_string = '230104';
results_path = fullfile(fileparts(pwd),'ExVivo','Results',[date_string '_M1M2_SingleMeas']);
if ~exist(results_path,'dir')
    mkdir(results_path);
end

% Loading *ALREADY COREGISTERED* diffusion maps (in R2* space): \kappa
% (dispersion) and ICVF maps. Also, a BaseFile is used to open .nii files 
% with the ACID toolbox and a manual segmented mask.
BaseFile = fullfile(dataset_path,'sPMB_Chlange-0003-00001-000001-01-001.nii');
Kappa_map = ACID_read_vols(spm_vol(fullfile(dataset_path,'example_kappa_maskedx10_x_manual_adapted.nii')),spm_vol(BaseFile),1);
FICV_map = ACID_read_vols(spm_vol(fullfile(dataset_path,'example_ficvf_masked_x_manual_adapted.nii')),spm_vol(BaseFile),1);
Mask = logical(ACID_read_vols(spm_vol(fullfile(dataset_path,'Mask_OC.nii')),spm_vol(BaseFile),1));

% Ex vivo data results are saved not only in .mat file but also in
% NII-format.
BaseNii = load_untouch_nii(BaseFile);
BaseNii.hdr.dime.datatype = 64;

% Some variables used in the code and can be modified:
exp_time = 3.4:3.34:53.5; % The 16 echo times in OC - in ms.
params.lambda = zeros(1,4); % Regularisation parameters (comment properly)

%SNR calculation
Background = logical(ACID_read_vols(spm_vol(fullfile(dataset_path,'Background_OC.nii')),spm_vol(BaseFile),1));
Im_OC = ACID_read_vols(spm_vol(BaseFile),spm_vol(BaseFile),1);
VoxelsInBack = Im_OC(Background);
SNR_full = Im_OC/std(VoxelsInBack(VoxelsInBack < 100));

%Reported value in manuscript
SNR_in_WM = mean(SNR_full(Mask(:) == 1 & FICV_map(:) > 0.8));

% (Data Loading) Loading data used for analysis
% Loading OC data from volumes or from already .mat file

if exist(fullfile(dataset_path,'Chiasm1_allecho'),'dir') &&...
  ~exist(fullfile(dataset_path,'OC_VolsLoaded.mat'),'file') 
    disp('Mat-file of the R2s OC data does not exist, so it can be created because volumes are available');
    code_path = pwd; % saves the code path to return at it afterwards.
    for ang_meas_indx = 1:16 % The 16 angular measurements.
        path = fullfile(dataset_path,'Chiasm1_allecho',num2str(ang_meas_indx));
        cd(path)

        if ang_meas_indx == 1 
            % First measurement is not corregistered, so it has this
            % prefix.
            files = dir('sPMB*-001.nii');
            BaseFile = fullfile(path,files(1).name);
        else
            % Corregistered measurements have this double rr- prefix added.
            files = dir('rrsPMB*-001.nii');
        end

        for echo_time_indx = 1:size(files,1)
            % Going through all the echo times read in the section above.
            DataImages_OC(:,:,:,ang_meas_indx,echo_time_indx) = ...
            double(ACID_read_vols(spm_vol(files(echo_time_indx).name),spm_vol(BaseFile),1));
        end
    end

    save(fullfile(dataset_path,'OC_VolsLoaded.mat'), 'DataImages_OC', '-v7.3');
    cd(code_path); % Here it returns to the code's path.
else
    disp('Mat-file of the R2s OC data does exist, so it is being loaded');
    load(fullfile(dataset_path,'OC_VolsLoaded.mat'));
end

% Loading angular maps data from volumes or from already .mat file. These
% angular maps were already estimated. See code "".
if exist(fullfile(dataset_path,'AngularOrientation'),'dir') && ~exist(fullfile(dataset_path,'OC_AnglesLoaded.mat'),'file')
    disp('Mat-file of the angular orientations of the OC data measurements does not exist, so it can be created because data is available');
    for ang_meas_indx = 1:16 % The 16 angular measurements.
        theta_map_path = fullfile(dataset_path,'AngularOrientation',...
                         ['ThetatoB0Meas' num2str(ang_meas_indx) '_DiffSpace_Zdir_0toPiHalf_x_manual_adapted.nii']);
        theta_map = ACID_read_vols(spm_vol(theta_map_path),spm_vol(BaseFile),1);
        
        % Some corrections: All negative angles becomes positive (due to
        % cos(angle) symmetry, and given also the symmetry in the
        % measurement, all angles higher than pi/2 are corrected.
        theta_map(theta_map < 0) = abs(theta_map(theta_map < 0)); % Due to symmetry, negative angle become positive. 
        theta_map(theta_map > pi/2) = pi - theta_map(theta_map > pi/2); % All angles over pi/2 gets "corrected".
        theta_map(theta_map < 0) = abs(theta_map(theta_map < 0)); % Due to symmetry, negative angle become positive.
        Theta_OC(:,:,:,ang_meas_indx) = theta_map;
    end
    save(fullfile(dataset_path,'OC_AnglesLoaded.mat'), 'Theta_OC', '-v7.3');
elseif exist(fullfile(dataset_path,'AngularOrientation'),'dir') && exist(fullfile(dataset_path,'OC_AnglesLoaded.mat'),'file')
    disp('Mat-file of the angular measurement OC data does exist, so it is being loaded');
    load(fullfile(dataset_path,'OC_AnglesLoaded.mat'));
else
    disp('Mat-file and .nii of the angular orientations of the OC data measurements do not exist, so new ones are calculated.');
    if ~exist(fullfile(dataset_path,'AngularOrientation_DiffSpace'),'dir')
        mkdir(fullfile(dataset_path,'AngularOrientation_DiffSpace'));
    end
    ExVivo_Setup_FromDiffusion(dataset_path);
    error('In this point, manual transformation from diffusion to R2* in AngularOrientation_DiffSpace must be done outside Matlab');
end

%% (Data fitting) OC data analysis
clc;
h = waitbar(0,'Initializing waitbar...');

reg_string = 'NoReg';
params.kappa = Kappa_map(Mask(:));
        
for indx_data = 1 % Non-normalised (1) and normalised fitting (2)
    for indx_angle = 1:16 % All angular orientations
        % Update the angular orientation per measurement.
        loaded_theta_map = squeeze(Theta_OC(:,:,:,indx_angle));
        params.theta = loaded_theta_map(Mask(:));
        
        DataA = reshape(squeeze(DataImages_OC(:,:,:,indx_angle,:)), prod(size(DataImages_OC,1:3)),size(DataImages_OC,5));
        DataA = DataA(Mask(:),:);
        
        if indx_data == 1 % Non-normalised analysis.
            DataA = real(log(DataA));
            % These parameters are used into "FittingEquations.m" to
            % analyse this data in a non-normalised approach.
            params.norm = 0;
            params.ref = 0;
            params.remove = 0;
            string_file = 'NonNormalisedFitting';
        else % Normalised analysis with respect the first echo.
            DataA = real(log(DataA./DataA(:,1)));
            % These parameters are used into "FittingEquations.m" to
            % analyse this data in a normalised approach.
            params.norm = 1;
            params.ref = 0;
            params.remove = 1;
            string_file = 'NormalisedFitting';
        end
       
        for te_indx = 1:5 % As In-silico analysis: max TE of 18, 27, 36, 45 and 54 ms are used for analysis.
            % "Filter" the time points until a maximum of the chosen max
            % TE.
            time = exp_time(exp_time < 9*(1 + te_indx))';
            
            for model_indx = 1:2 % All models, from M1 to M4 - FJv20(23.11): M3 is divided in M3a and M3b, as well as M4a and M4b.
                params.model = model_indx - 1; 
                [results, gof] = FittingEquations_Submission(time, DataA(:, 1:length(time))',[],params);
                % From the resulting data, filter the NaN's, INF's, -INF's
                % and Imag values.
                gof.full_res(isnan(results.full_b) | ...
                    isinf(real(results.full_b)) | imag(results.full_b) ~= 0 | ...
                    isinf(abs(real(results.full_b)))) = 0;
                results.full_b(isnan(results.full_b) | ...
                    isinf(real(results.full_b)) | imag(results.full_b) ~= 0 | ...
                    isinf(abs(real(results.full_b)))) = 0;
                
                % Here the "statistics" are saved properly: BIC, R^2 and
                % R^2 predicted:
                Auxmap = zeros(112,112,112);
                Auxmap(Mask) = length(time)*log(gof.full_res/length(time)) + (model_indx+1)*log(length(time));
                BIC_fit(:,:,:,indx_angle,te_indx,model_indx) = Auxmap;
                
                Auxmap = zeros(112,112,112);
                Auxmap(Mask) = 1 - gof.full_res./gof.sst;
                SquaredR_fit(:,:,:,indx_angle,te_indx,model_indx) = Auxmap;
                
                Auxmap = zeros(112,112,112);
                Auxmap(Mask) = 1 - length(time)*gof.cross_validation_value./gof.sst;
                SquaredRPred_fit(:,:,:,indx_angle,te_indx,model_indx) = Auxmap;
                
                Auxmap = zeros(112,112,112);
                Auxmap(Mask) = 2*(model_indx+1) + length(time)*log(gof.full_res/length(time)) + (2*(model_indx+1)^2+2*(model_indx+1))/(length(time)-(model_indx+1)-1);
                AICc_fit(:,:,:,indx_angle,te_indx,model_indx) = Auxmap;
                
                Auxmap = zeros(112,112,112);
                Auxmap(Mask) = gof.full_res;
                SSE_fit(:,:,:,indx_angle,te_indx,model_indx) = Auxmap;
                                
                for indx_stats = 4   
                    switch indx_stats
                        case 1 % BIC
                            BaseNii.img = Auxmap;
                            BaseNii.hdr.dime.glmax = 10000;
                            BaseNii.hdr.dime.glmin = -10000;
                            save_untouch_nii(BaseNii,fullfile(results_path,['BIC_Model_' num2str(model_indx) '_' string_file '_' reg_string '_Meas' num2str(indx_angle) '_TEmax' num2str(te_indx)]));
                        case 2 % R^2
                            BaseNii.img = Auxmap;
                            BaseNii.hdr.dime.glmax = 1;
                            BaseNii.hdr.dime.glmin = 0;
                            save_untouch_nii(BaseNii,fullfile(results_path,['SquaredR_Model_' num2str(model_indx) '_' string_file '_' reg_string '_Meas' num2str(indx_angle) '_TEmax' num2str(te_indx)]));
                        case 3 % Predicted R^2
                            BaseNii.img = Auxmap;
                            BaseNii.hdr.dime.glmax = 1;
                            BaseNii.hdr.dime.glmin = 0;
                            save_untouch_nii(BaseNii,fullfile(results_path,['SquaredRPred_Model_' num2str(model_indx) '_' string_file '_' reg_string '_Meas' num2str(indx_angle) '_TEmax' num2str(te_indx)]));
                        case 4 % AICc
                            BaseNii.img = Auxmap;
                            BaseNii.hdr.dime.glmax = 1;
                            BaseNii.hdr.dime.glmin = 0;
                            save_untouch_nii(BaseNii,fullfile(results_path,['AICc_Model_' num2str(model_indx) '_' string_file '_' reg_string '_Meas' num2str(indx_angle) '_TEmax' num2str(te_indx)]));
                        case 5 % SSE
                            BaseNii.img = Auxmap;
                            BaseNii.hdr.dime.glmax = 1;
                            BaseNii.hdr.dime.glmin = 0;
                            save_untouch_nii(BaseNii,fullfile(results_path,['SSE_Model_' num2str(model_indx) '_' string_file '_' reg_string '_Meas' num2str(indx_angle) '_TEmax' num2str(te_indx)]));
                    end
                end
                
                % Here the beta parameters are saved in the .mat file and
                % .nii accordingly.
                for beta_indx = 1:3
                    Auxmap = zeros(112,112,112);
                    % Since the beta estimations were performed in ms,
                    % \beta_1 and \beta_2 must be multiplied by 1000 and
                    % 1000000 (that leaves them in 1/s and
                    % 1/s^2 units) respectively.
                    Auxmap(Mask) = squeeze(real(results.full_b(beta_indx,:).*1000^(beta_indx-1)));
                    
                    if indx_data == 1
                        beta_nn(:,:,:,indx_angle,te_indx,model_indx,beta_indx) = Auxmap;
                    else
                        beta_n(:,:,:,indx_angle,te_indx,model_indx,beta_indx) = Auxmap;
                    end
                    
                    if beta_indx == 1 
                        %\beta_0 is defined analytically as ~ log(S_0).
                        BaseNii.img = exp(Auxmap);
                    elseif beta_indx == 2
                        Auxmap(Auxmap < 0 | Auxmap > 100) = 0;
                        BaseNii.img = Auxmap;
                        BaseNii.hdr.dime.glmax = 100;
                        BaseNii.hdr.dime.glmin = 0;
                    elseif beta_indx == 3
                        Auxmap(Auxmap < -2000 | Auxmap > 2000) = 0;
                        BaseNii.img = Auxmap;
                        BaseNii.hdr.dime.glmax = 2000;
                        BaseNii.hdr.dime.glmin = -2000;
                    end
                    
                    if indx_data == 1 
                        if ~exist(fullfile(results_path,['Model_' num2str(model_indx)],'NonNormalisedFitting'),'dir')
                            mkdir(fullfile(results_path,['Model_' num2str(model_indx)],'NonNormalisedFitting'));
                        end
                        save_untouch_nii(BaseNii,fullfile(results_path,['Model_' num2str(model_indx)],'NonNormalisedFitting',[reg_string '_Meas' num2str(indx_angle) '_Beta' num2str(beta_indx-1) '_TEmax' num2str(te_indx)]));
                    else
                        if ~exist(fullfile(results_path,['Model_' num2str(model_indx)],'NormalisedFitting'),'dir')
                            mkdir(fullfile(results_path,['Model_' num2str(model_indx)],'NormalisedFitting'));
                        end
                        save_untouch_nii(BaseNii,fullfile(results_path,['Model_' num2str(model_indx)],'NormalisedFitting',[reg_string '_Meas' num2str(indx_angle) '_Beta' num2str(beta_indx-1) '_TEmax' num2str(te_indx)]));  
                    end
                end
            end
              
            waitbar((te_indx + 5*(indx_angle-1) + 5*16*(indx_data-1))/(2*5*16),h,...
                ['Index #' num2str(indx_data) ', angle index #' num2str(indx_angle) ' with TEmax #' num2str(te_indx)]);   
        end
    end
    if indx_data == 1
        save(fullfile(results_path,[reg_string '_Betas_M1M2_NonNormalised.mat']),...
        'beta_nn', 'BIC_fit', 'SquaredR_fit', 'SquaredRPred_fit', 'AICc_fit', 'SSE_fit', '-v7.3');
        clear 'beta_nn' 'BIC_fit' 'SquaredR_fit' 'SquaredRPred_fit' 'AICc_fit' 'SSE_fit'
    else
        save(fullfile(results_path,[reg_string '_Betas_M1M2_NormalisedFirstEchoRemoved.mat']),...
        'beta_n', 'BIC_fit', 'SquaredR_fit', 'SquaredRPred_fit', 'AICc_fit', 'SSE_fit', '-v7.3');
        clear 'beta_n' 'BIC_fit' 'SquaredR_fit' 'SquaredRPred_fit' 'AICc_fit' 'SSE_fit'
    end
end

close(h);

%% (Data post-processing and Figures) - In the entire OC.
% Anisotropic binning is performed here *AND* the obtained values are used
% to create the figures.
clc;

% load here the required data (change manually if required):
load(fullfile(results_path,'NoReg_Betas_M1M2_NonNormalised.mat'));
clear 'res_nn' 'BIC_fit' 'SquaredR_fit' 'SquaredRPred_fit' 'SSE_fit'; % IF Normalised is used, change this to res_n.

% Anisotropic binning function has the following inputs: angle maps, \kappa
% map, ICVF map, beta values from the section before and mask. More
% information IN the function.
if ~ exist('Theta_OC','var') % Theta_OC data is needed.
    load(fullfile(fileparts(pwd),'Dataset','OC_AnglesLoaded.mat'));
end % Otherwise change accordingly in the next function if another name is used.
params_anisotropy.nbins = 20;

% Given the loaded .mat file, beta_nn is used. If Normalised .mat file is
% used, change beta_nn for beta_n. And for any other name, change
% accordingly.
[ExVivoTableSummary,HistogramFigure5] = ...
    AnisotropyBinning_Submission(Theta_OC,Kappa_map,FICV_map,beta_nn,Mask,params_anisotropy);

[wAICc_mean_estimation,wAICc_sd_estimation] = ...
    wAICcEstimationExVivo_Submission(AICc_fit,Theta_OC,Kappa_map,FICV_map,Mask,params_anisotropy);

ExVivoTableSummary.wAICc_mean_M2 = reshape(permute(wAICc_mean_estimation(:,:,:,2),[1,3,2,4]),[prod(size(wAICc_mean_estimation,1:3)),1]);
ExVivoTableSummary.wAICc_sd_M2 = reshape(permute(wAICc_sd_estimation(:,:,:,2),[1,3,2,4]),[prod(size(wAICc_mean_estimation,1:3)),1]);

save(fullfile(results_path,'ExVivoResults_IrregularBinning_NotNormalised.mat'),...
     'ExVivoTableSummary','HistogramFigure5','-v7.3');
