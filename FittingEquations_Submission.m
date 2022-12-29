function [params_fitting, gof_fitting] = ...
         FittingEquations_Submission(data_x, data_y, ~, params)

% Since the equation to be fitted is quite extensive and complicated, it is
% then necessary to perform in a different code to avoid overwritting in
% MatlabCode_for_Figures.m

% UPDATES:
% 15.07.2022: Reduced version created (FittingEquations_Submission.m),
% where only M1 and M2 are fitted + no regularization is performed.
% 22.02.2020: Opt_norm added. This tells to the fitting procedure if the
% data is normalised or not with respect the first time value.
% 26.05.2020: All the extra parameters were saved in "params" variables. It
% should contain:
% params.kappa: kappa value,
% params.theta: angular orientation in rads.
% params.model: model number (1: linear, 2: quadratic, 3: quadratic*sin^4 and
%  4: quadratic with dispersion.
% params.norm: boolean input if data and fitting is normalised or not.
% params.ref: boolean input if echo used for normalisation is the first (0) or
%  the last (1).
% params.remove: boolean input if the echo used for normalisation is removed or not.
% 09.09.2020: Code cleaning for in silico and ex vivo analysis. Also, the
% params.lambda was also incorporated in case of regularisation study is needed.

% Variables initialisation:
lambda_used = zeros(1,size(data_y,2));
resnorm = Inf*ones(1,size(data_y,2));
x_final = []; %zeros(3,size(data_y,2));
x_new = []; %zeros(3,size(data_y,2));
resnorm_new = zeros(1,size(data_y,2));

% Data initialisation:
if params.norm
    if params.ref % Last echo normalisation
        data_x_1st = data_x - data_x(end);
        data_x_2nd = data_x.^2 - data_x(end).^2;
        
        if params.remove
            data_x_1st = data_x_1st(1:end-1);
            data_x_2nd = data_x_2nd(1:end-1);
            data_y = data_y(1:end-1,:);
        end
    else % first echo normalisation
        data_x_1st = data_x - data_x(1);
        data_x_2nd = data_x.^2 - data_x(1).^2;
        
        if params.remove
            data_x_1st = data_x_1st(2:end);
            data_x_2nd = data_x_2nd(2:end);
            data_y = data_y(2:end,:);
        end
    end
    A = [zeros(length(data_x_1st),1), -data_x_1st, -data_x_2nd];
else
    data_x_1st = data_x;
    data_x_2nd = data_x.^2;
    A = [ones(length(data_x_1st),1), -data_x_1st, -data_x_2nd];
end

switch params.model
    case 0 % Model 1 or linear model
        if length(params.theta) == 1 % same design matrix for all voxels
            A(:,3) = 1.*A(:,3);
        else % design matrix in which the last component changes per voxel
            A_extended = A(:,3)*ones(1,size(data_y,2)); % n times points x m voxels.
        end
        
        if params.norm
            A = A(:,2);
        else
            A = A(:,1:2);
        end
        
    case 1 % Model 2 or quadratic model
        if length(params.theta) == 1 % same design matrix for all voxels
            A(:,3) = 1.*A(:,3);
        else % design matrix in which the last component changes per voxel
            A_extended = A(:,3)*ones(1,size(data_y,2)); % n times points x m voxels.
        end
        
        if params.norm
            A = A(:,2:3);
        end
end

aux_matrix = zeros(size(A,2));
aux_matrix(end,end) = data_x_1st(end).^(numel(A(1,:))-1);
% H matrix for calculating the cross validation function
diag_h_hat = diag(A*inv(A'*A)*A');

if length(params.theta) == 1 % one angle for "all voxels" -> mostly in silico data case
    x_new = A\data_y;
    error_fit = A*x_new - data_y;
    resnorm_new = sum(error_fit.^2);
    crossV_fit = 1/size(data_y,1)*sum((error_fit./(1-diag_h_hat)).^2);
    total_sum_squares = sum((data_y - mean(data_y)).^2);
else % mostly ex vivo case
    if params.model < 2 % No regularisation and same design matrix for all voxels
        x_new = A\data_y;
        error_fit = A*x_new - data_y;
        resnorm_new = sum(error_fit.^2);
        crossV_fit = 1/size(data_y,1)*sum((error_fit./(1-diag_h_hat)).^2);
        total_sum_squares = sum((data_y - mean(data_y)).^2);
    else
        for vox_indx = 1:length(params.theta) % iteration per voxel
            x_new(:,vox_indx) = A\data_y(:,vox_indx);
            error_fit = A*x_new(:,vox_indx) - data_y(:,vox_indx);
            resnorm_new(vox_indx) = sum(error_fit.^2);
            crossV_fit(vox_indx) = 1/size(data_y,1)*sum((error_fit./(1-diag_h_hat)).^2);
            total_sum_squares(vox_indx) = sum((data_y(:,vox_indx) - mean(data_y(:,vox_indx))).^2);
        end
    end
end
        
res_indx = find(resnorm > resnorm_new);
resnorm(res_indx) = resnorm_new(res_indx);
x_final(:,res_indx) = x_new(:,res_indx);

if ~isempty(x_new)
    if params.norm
        params_fitting.b0 = zeros(size(x_final(1,:)));
        params_fitting.b1 = x_final(1,:);
        if (size(x_final(:,1),1) == 2)
            params_fitting.b2 = x_final(2,:);
        else
            params_fitting.b2 = zeros(size(x_final(1,:)));
        end
    else
        params_fitting.b0 = x_final(1,:);
        params_fitting.b1 = x_final(2,:);
        if (size(x_final(:,1),1) == 3)
            params_fitting.b2 = x_final(3,:);
        else
            params_fitting.b2 = zeros(size(x_final(1,:)));
        end
    end
end
                
params_fitting.full_b = [params_fitting.b0;...
                         params_fitting.b1;...
                         params_fitting.b2];        
    
gof_fitting.res = resnorm;
gof_fitting.full_res = resnorm_new;
gof_fitting.cross_validation_value = crossV_fit;
gof_fitting.sst = total_sum_squares;
end