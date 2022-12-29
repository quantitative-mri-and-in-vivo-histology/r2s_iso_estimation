%% General cylinders study

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was created to obtain R2* synthetic data for any cylinder
% configuration and be fitted afterwards by R2* signal models.
% The code operates as follows: A set of
% directions (defined by a file or created on fly) are created, 
% each one representing a cylinder. Each one contributes to the
% total signal by its intra-axonal and myelin (if not neglected)
% compartments, and its effect to the extra-axonal compartment.
% All variables are recommended to be defined at the beginning of
% the code. When the signal is created, it is fitted immediately 
% by each model (classic and quadratics).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Path, folder and variable names (and checking if they are normalised):
path_for_data = fileparts(pwd);

% If data has not be created or want to run the code entirely at once,
% executed this line
date_for_folder_filenames = datestr(today,'yymmdd');
% Otherwise, it has been created previously or want to run some sections of
% the code, change the date_for_folder_filenames accordingly:
% e.g. date_for_folder_filenames = '210811';
folder_signal_decay = fullfile(path_for_data,'InSilico',['SignalDecay_' date_for_folder_filenames]);
folder_fitted_signal = fullfile(path_for_data,'InSilico',['FittedSignal_' date_for_folder_filenames]);

DispNoMye_Sampling_points = [];
NoDispNoMye_Sampling_points = [];
DispMye_Sampling_points = [];
NoDispMye_Sampling_points = [];
datatypes = {'DispNoMye','NoDispNoMye','DispMye','NoDispMye'};
%For fitting
reg_string = 'NoReg';

% Cylinders
Data = load('Directions1500_Cartesian.mat'); % If defined by text file
DAxons1 = Data.DAxons1;
DAxons1_Norm = [DAxons1(:,1), DAxons1(:,2), DAxons1(:,3)]./sqrt(sum(DAxons1.^2,2));

% Parameter simulation space.
angle_values = pi/90:pi/90:pi/2; % In radians. 
kappa_values = [0.1:0.1:6.0,7.5,10,20]; % dispersion factor                                            
gratio_values = [0.66,0.73,0.8];
sim_time_range = (0:0.00025:0.06)*1000; % Oversampling time in msec.
samples = 5000; %number of samples created with added random noise

% Adaptable variables for study:
datatype_indexes = [1,3]; % Change HERE the indexes used for data analysis.
% The order, once more, goes as follows:
% datatype_indexes = 1 -> No myelin with dispersion
% datatype_indexes = 2 -> No myelin with no dispersion
% datatype_indexes = 3 -> Myelin with dispersion
% datatype_indexes = 4 -> Myelin with no dispersion

R2a_test = [1/53.96,1/36];
R2e_test = [1/53.96,1/36];
R2m_test = [1/13.26,1/8];
% Extra parameters useful for renaming the file name accordingly based on
% the previous listed g_samples_for_study and R2x_test's.
suffix_string = {'_T2ExVivo_GratioExVivo';'_T2ExVivo_GratioMiddle';...
                 '_T2ExVivo_GratioWB';'_T2WB_GratioExVivo';...
                 '_T2WB_GratioMiddle';'_T2WB_GratioWB'};

exp_time_range = 3.4:3.34:53.5; % Experimental undersampled time from ex vivo.
% TEmax study:
max_time = [18, 27, 36, 45, 54];

% FJv21(16.08): List here the SNRs to be used for study. As before, the
% noise is added on-fly based on the SNR to be studied.
SNR_study = 112; %OC WM showed to have a mean SNR of 112. Value changed accordingly.

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Creating synthetic R2* signal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FJv21(28.01): Ex vivo setup (Note: find the references for using these
% numbers. 
% FJv21(11.02): Ok, the "ex vivo setup" was performed with the following
% parameters: T2m ~ 13.26 ms and T2a = T2e ~ 53.96 ms from (https://onlinelibrary.
%wiley.com/doi/epdf/10.1002/mrm.22267) which, doing a monoexponential study
%on it, the "fitted" mean T2 is similar of what Evgeniya reported (T2_mean ~
%45 ms, from https://advances.sciencemag.org/content/6/41/eaaz9281), FVF and G-
%ratio were from Laurin's presentation, resulting in FVF ~ (0.213 + 0.276) = 0.49
% and G-ratio ~ sqrt(1 - 0.276/0.49) ~ 0.66. 
% FJv21(16.08): Due to space and consistency, only the R2* signal decay at
% SNR = infty will be saved. Noise addition and such is "added on fly" in
% data fitting.

for R2_indx = 1:numel(R2a_test)
    params.R2m = R2m_test(R2_indx);
    params.R2a = R2a_test(R2_indx);
    params.R2e = R2e_test(R2_indx);
    h = waitbar(0,'Initializing waitbar...');
    switch R2_indx
        case 1
            T2Reference = 'T_2 ex vivo values';
        case 2
            T2Reference = 'T_2 values from Wharton and Bowtell 2013';
    end
    for g_indx = 1:numel(gratio_values)     
        for k = 1:numel(kappa_values)  
            for j = 1:numel(angle_values)
                params.theta = angle_values(j);
                params.kappa = kappa_values(k);
                params.FiberLimit = 0.5; %This is the upper limit for "total FVF" after adding all cylinders.
                % Getting equations for datapoints
                Signal_time = SignalModelR2_Submission(DAxons1_Norm,...
                                                       0.5,... %FVF per cylinder = 0.5
                                                       gratio_values(g_indx),...
                                                       params);

                % Signal at SNR = infty is in complex, so complex noise is
                % added afterwards
                DispNoMye_Sampling_points(k,j,1:numel(sim_time_range)) = Signal_time{1}(sim_time_range);
                NoDispNoMye_Sampling_points(k,j,1:numel(sim_time_range)) = Signal_time{3}(sim_time_range);
                DispMye_Sampling_points(k,j,1:numel(sim_time_range)) = Signal_time{2}(sim_time_range);
                NoDispMye_Sampling_points(k,j,1:numel(sim_time_range)) = Signal_time{4}(sim_time_range);
                
                waitbar((j + numel(angle_values)*(k-1) + numel(angle_values)*numel(kappa_values)*(g_indx-1))...
                    /(numel(gratio_values)*numel(kappa_values)*numel(angle_values)),h,...
                    {'Creating the in silico data with the following microstructural properties: ',...
                    ['FVF-Gratio: 0.5 - ' num2str(gratio_values(g_indx)) ', with \kappa = ' num2str(kappa_values(k)) ', \theta_\mu ' num2str(angle_values(j)*180/pi) ' using ' T2Reference]});
            end
        end
        
        if ~exist(folder_signal_decay,'dir')
            mkdir(folder_signal_decay);
        end
        
        % data saving per signal equation and SNR:
        save(fullfile(folder_signal_decay,['DispNoMye_data_Gratio' num2str(gratio_values(g_indx)*100) '_FVF5' suffix_string{g_indx + (R2_indx-1)*numel(gratio_values)} '.mat']), 'DispNoMye_Sampling_points','-v7.3');
        save(fullfile(folder_signal_decay,['NoDispNoMye_data_Gratio' num2str(gratio_values(g_indx)*100) '_FVF5' suffix_string{g_indx + (R2_indx-1)*numel(gratio_values)} '.mat']), 'NoDispNoMye_Sampling_points','-v7.3');   
        save(fullfile(folder_signal_decay,['DispMye_data_Gratio' num2str(gratio_values(g_indx)*100) '_FVF5' suffix_string{g_indx + (R2_indx-1)*numel(gratio_values)} '.mat']), 'DispMye_Sampling_points','-v7.3');   
        save(fullfile(folder_signal_decay,['NoDispMye_data_Gratio' num2str(gratio_values(g_indx)*100) '_FVF5' suffix_string{g_indx + (R2_indx-1)*numel(gratio_values)} '.mat']), 'NoDispMye_Sampling_points','-v7.3');    
        
        clear 'DispNoMye_Sampling_points' 'NoDispNoMye_Sampling_points' 'DispMye_Sampling_points' 'NoDispMye_Sampling_points'               
    end
close(h);
end
%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fitting synthetic R2* signal (if created, loading file. Otherwise,
% it continues with the signals created above.) per sample and average!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;

% This section (the for loops) must be adapted AS EQUAL as the previous
% section.
for SNR_indx = 1:numel(SNR_study)
    
    SNR_to_study = SNR_study(SNR_indx);
    h = waitbar(0,'Initializing waitbar...');

    for R2_indx = 1:numel(R2a_test)
        switch R2_indx
            case 1
                T2Reference = 'T_2 ex vivo values';
            case 2
                T2Reference = 'T_2 values from Wharton and Bowtell 2013';
        end
        for g_indx = 1:numel(gratio_values)       
            for dt_indx = 1:numel(datatype_indexes) % Go through all desired datatypes.
                load(fullfile(folder_signal_decay,[datatypes{datatype_indexes(dt_indx)} '_data_Gratio' num2str(gratio_values(g_indx)*100) '_FVF5' suffix_string{g_indx + (R2_indx-1)*numel(gratio_values)} '.mat']));

                if exist('DispNoMye_Sampling_points','var') && ~isempty(DispNoMye_Sampling_points)
                    SNR_Sampling_points = DispNoMye_Sampling_points;
                    DispNoMye_Sampling_points = [];
                elseif exist('NoDispNoMye_Sampling_points','var') && ~isempty(NoDispNoMye_Sampling_points)
                    SNR_Sampling_points = NoDispNoMye_Sampling_points;
                    NoDispNoMye_Sampling_points = [];
                elseif exist('DispMye_Sampling_points','var') && ~isempty(DispMye_Sampling_points)
                    SNR_Sampling_points = DispMye_Sampling_points;
                    DispMye_Sampling_points = [];
                elseif exist('NoDispMye_Sampling_points','var') && ~isempty(NoDispMye_Sampling_points)
                    SNR_Sampling_points = NoDispMye_Sampling_points;
                    NoDispMye_Sampling_points = [];
                end
                                
                for norm_indx = 1:2 % Without/with normalisation with first echo        
                    for max_indx = 1:numel(max_time) % Maximum time index
                        time_end = max(find(exp_time_range < max_time(max_indx)));
                        time_indx_exp = zeros(1,time_end);

                        for ti = 1:time_end
                            time_indx_exp(ti) = max(find(time_range <= exp_time_range(ti)));
                        end

                        for k = 1:numel(kappa_values)
                            for j = 1:numel(angle_values)
                                Noise_addition = NoiseAddition(squeeze(SNR_Sampling_points(k,j,1)),SNR_to_study);
                                Sampling_points = squeeze(abs(squeeze(SNR_Sampling_points(k,j,time_indx_exp))' + squeeze(Noise_addition(samples,length(time_indx_exp)))));
                                
                                data_avg = squeeze(mean(Sampling_points));

                                if norm_indx == 1
                                    data_points = log(Sampling_points)';
                                    data_points_avg = log(squeeze(mean(Sampling_points)))';
                                    params_fit.norm = 0;
                                    params_fit.ref = 0;
                                    params_fit.remove = 0;
                                else                                    
                                    data_points = log(Sampling_points./Sampling_points(:,1))';
                                    data_points_avg = log(data_avg/data_avg(1))';
                                    params_fit.norm = 1; % If data is normalised
                                    params_fit.ref = 0; % If takes the first (0) or last (1) echo for normalisation
                                    params_fit.remove = 1; % If the echo used for normalisation is removed (1 if yes)
                                end

                                params_fit.kappa = kappa_values(k);
                                params_fit.theta = angle_values(j);

                                time_exp = time_range(time_indx_exp)';

                                for model_indx = 1:2
                                    params_fit.model = model_indx-1;

                                    % Full samples
                                    [OLS_Fitting{max_indx,k,j,model_indx}, OLS_gof{max_indx,k,j,model_indx}] = ...
                                        FittingEquations_Submission(time_exp, data_points, [], params_fit);

                                    % Average signal from samples
                                    [OLS_Fitting_avg{max_indx,k,j,model_indx}, OLS_gof_avg{max_indx,k,j,model_indx}] = ...
                                        FittingEquations_Submission(time_exp, squeeze(data_points_avg), [], params_fit);
                                end

                            waitbar((j + (k-1)*numel(angle_values) + (max_indx-1)*numel(angle_values)*numel(kappa_values))/(numel(angle_values)*numel(kappa_values)*5),...
                                h,{['Fitting in silico data with ' T2Reference ' and with the following microstructural characteristics: '],...
                                ['FVF-Gratio: 0.5 - ' num2str(gratio_values(g_indx)*100) ', Fitting SNR ' num2str(SNR_to_study), ', Max TE ' num2str(max_time(max_indx)) ', (\kappa,\theta_\mu): (' num2str(kappa_values(k)) ',' num2str(angle_values(j)*180/pi) ' deg)']})
                            end
                        end
                    end

                    if ~exist(folder_fitted_signal,'dir')
                        mkdir(folder_fitted_signal);
                    end

                    if norm_indx == 1
                        save(fullfile(folder_fitted_signal,['FittedSignal_' datatypes{datatype_indexes(dt_indx)} '_1500cylinders_TimeES_' reg_string '_NNormalised_SNR' num2str(SNR_to_study) '_FVF5' suffix_string{g_indx + (R2_indx-1)*numel(gratio_values)} '.mat']), 'OLS_Fitting', 'OLS_Fitting_avg', 'OLS_gof', 'OLS_gof_avg', '-v7.3');
                    else
                        save(fullfile(folder_fitted_signal,['FittedSignal_' datatypes{datatype_indexes(dt_indx)} '_1500cylinders_TimeES_' reg_string '_NormalisedRemoved1stEcho_SNR' num2str(SNR_to_study) '_FVF5' suffix_string{g_indx + (R2_indx-1)*numel(gratio_values)} '.mat']), 'OLS_Fitting', 'OLS_Fitting_avg', 'OLS_gof', 'OLS_gof_avg', '-v7.3');
                    end
                    clear 'OLS_Fitting' 'OLS_Fitting_avg' 'OLS_gof' 'OLS_gof_avg' 'time_end' 'time_indx_exp' 
                end
            end
            clear 'Sampling_points'
        end
    end
    close(h)
end

%% Figure for papers:
clc;

%This is used ONLY to load the data
g_samples = repmat(gratio_values,1,2);
params_fig.onlyAniso = 1;
params_fig.TEmaxTable = max_time;

%This is saved in a table.
for i = 1:6
    params_fig.gratio = g_samples(i);
    params_fig.datatype = 'NoMye';
    switch i < 4
        case 1
            params_fig.T2Reference = 'T2ExVivo';
        case 0
            params_fig.T2Reference = 'T2WB';
    end
    NoMyelin_InSilicoResults{i} = AnisotropicBinning_InSilico_fromExVivo(fullfile(folder_fitted_signal,['FittedSignal_' datatypes{datatype_indexes(1)} '_1500cylinders_TimeES_' reg_string '_NNormalised_SNR' num2str(SNR_study) '_Gratio' num2str(g_samples(i)*100) '_FVF' num2str(FVF_samples(1)*10) suffix_string{i} '.mat']),params_fig);
    params_fig.datatype = 'Mye';
    Myelin_InSilicoResults{i} = AnisotropicBinning_InSilico_fromExVivo(fullfile(folder_fitted_signal,['FittedSignal_' datatypes{datatype_indexes(2)} '_1500cylinders_TimeES_' reg_string '_NNormalised_SNR' num2str(SNR_study) '_Gratio' num2str(g_samples(i)*100) '_FVF' num2str(FVF_samples(1)*10) suffix_string{i} '.mat']),params_fig);
end

save('Myelin_NoMyelin_AnisotropicBinarised_InSilicoMR_NewWeighting_20092021.mat','NoMyelin_InSilicoResults','Myelin_InSilicoResults');

load('Myelin_NoMyelin_AnisotropicBinarised_InSilicoMR_NewWeighting_20092021.mat');
r11082021_FiguresForPaper; % Main Figures 6 and 9.
params.DynRange = 'Aniso';
% Main results, Figures 7B, 8 and 10.
[RMSD_ratio_B1, rCoV_B2, DynRange_Data] = Stats_calculation(NoMyelin_InSilicoResults,Myelin_InSilicoResults,'ExVivoAnisotropyBinResults.mat',params);

% Supplementary Figures:
SupplementaryFigures_Code;
