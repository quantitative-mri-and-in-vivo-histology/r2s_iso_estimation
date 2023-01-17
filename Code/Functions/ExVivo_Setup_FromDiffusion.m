function ExVivo_Setup_FromDiffusion(path)
%% Changing path names from Sebastian's data:
clc;
load(fullfile(path,'Chiasm1_allecho','Image_lists_R2s_1stO.mat'));

%figure;
for indx = 1:16
    [pthname, fname, ext] = fileparts(Porg{indx,1});
    new_pthname = fullfile(path,'Chiasm1_allecho',num2str(indx));
    fname_base = strsplit(fname, '_R');
    Porg_Frank{indx,1} = fullfile(new_pthname,[fname_base{1,1} ext]);
    [pthname, fname, ext] = fileparts(Prot{indx,1});
    new_pthname = fullfile(path,'Chiasm1_allecho',num2str(indx));
    fname_base = strsplit(fname, '_R');
    Prot_Frank{indx,1} = fullfile(new_pthname,[fname_base{1,1} ext]);
end

% Loading vVector in diffusion space. The out is the "angular orientation"
% per angle. The estimations are done "inside".
Estimatetheta_Voxelwise(Porg_Frank,Prot_Frank,...
                        fullfile(path,'DiffusionResults_NODDI','example_all_vector.nii'),...
                        fullfile(path,'DiffusionResults_NODDI','MaskToBeUsed_DiffusionOC.nii'));

end
