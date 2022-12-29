function [noise_function] = NoiseAddition(model_at_t0,SNR)
%NOISEADDITION This is a small function which generates the background
%noise in the in silico data with a specific SNR at time = 0 ms.
%   Detailed explanation goes here
% Input:
% model_at_t0 - double
% This variable is the value of the model/function at time = 0. This
% value can be in complex or real number.
% SNR - double:
% This is a double variable that specifies the required SNR.
%
% Output:
% noise_function - function handle: noise_values @(n,m) 
% This is a function handle that generates a (complex) random gaussian
% noise with n-rows and m-columns.
% model_with_noise - function handle: model_with_noise @(t,n)
% This is a function handle that contains the model with added (complex)
% noise. Ideally the input values in model must be in columns, because the
% total number of elements are used to evaluate the m-elements in the
% noise_function.

% Standard deviation of the noise in complex space.
noise_std_complex = abs(model_at_t0)/SNR;% * sqrt(2 /(4 - pi));
noise_function = @(n,m) normrnd(zeros([n,size(model_at_t0),m]),permute(repmat(noise_std_complex,1,1,n,m),[3,1,2,4]))...
                        + 1i*normrnd(zeros([n,size(model_at_t0),m]),permute(repmat(noise_std_complex,1,1,n,m),[3,1,2,4]));
end

