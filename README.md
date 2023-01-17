# r2s_iso_estimation
Source code and data used for the paper "Deciphering the fibre-orientation independent component of R2* (R2,iso*) in the human brain with a single multi-echo gradient-recalled-echo measurement under varying microstructural conditions" from Fritz et al 202X."
The two main codes are located in "Code/InSilicoStudy_Submission.m" and "Code/OC_AnalysisCode_Submission.m". The remaining functions are contained in the folder "Code/Functions".

In simple words, the codes function in the following way:
OC_AnalysisCode_Submission expects a dataset (meGRE and diffusion) in > ExVivo/Datasets folder. After analysis, everything is saved in ExVivo/Results/[date]_M1M2_SingleMeas folder. The [date] is automatically defined by the code (which uses "today" as default). Then a table is created which summarises all the results from the ex vivo analysis (i.e. model parameters, irregular binning, wAICc, nRMSD, etc). This table is uploaded accordingly (ExVivo/[date]_M1M2_SingleMeas) and it is loaded in Code/FiguresForPaper.m function. The Dataset folder will not be uploaded with data except when it is correspondingly required.

InSilicoStudy_Submission expects the summary table from the ex vivo analysis and the 1500 directions. After that, it creates the meGRE signal without noise and saves it in InSilico/SignalDecay_[date]. After that, the created meGRE signal per simulated parameter (i.e. g-ratio, dispersion, mean angular orientation and T2 values) is fitted but with added noise (the noise addition is performed on fly). Then the fitted model parameters (plus other metrics) are saved in FittedSignal_[date]. Then the resulted data is irregularly-binarised to become similar to the ex vivo data and the summary table is created. This table contains similar information as the ex vivo summary table (plus extra simulation parameter dimensions like g-ratio). This table is uploaded accordingly (in InSilico/) and loaded accordingly in Code/FiguresForPaper.m

FigureForPaper creates (most) of the figures based on the loaded summary tables described before.

It is important to:
1. Execute any of the main codes in their corresponding path working directory. Otherwise the code will not find the corresponding ExVivo and/or InSilico paths accordingly.
2. Download the spm12 folder into the Code folder.
