% --------------------------------------------------------------------
% example code for Causality-driven Trustworthy Multi-View maPping approach (Cad-TMVP)
%------------------------------------------
% Author: jinzhang@mail.nwpu.edu.cn
% Date created: 06-07-2025
% @Northwestern Ploytechnical University.
% ------------------------------------------------------------------------------------

% This study used data from the Alzheimer's Disease Neuroimaging Initiative (ADNI) database (http://adni.loni.usc.edu/). 
% ADNI was established to investigate the use of MRI, PET, and other measures to track the progression of mild cognitive impairment (MCI) and Alzheimer's disease (AD). 
% Researchers can access data by registering and agreeing to the ADNI data use agreement on the website. 
% The other cancer omics can be obtained from the referred paper in our manuscript.

close all; 
clearvars -except round;
clc;

% Load data
load('Data.mat');  % From http://adni.loni.usc.edu/
% Normalization
Y{1} = getNormalization(Y{1}, 'normalize');
Y{2} = getNormalization(Y{2}, 'normalize');
Y{3} = getNormalization(Y{3}, 'normalize');
DX = getNormalization(DX, 'normalize');

X{1} = Y{1};
X{2} = Y{2};
X{3} = Y{3};

% time_start = cputime;
% Tuned parameters
% set parameters

opts.CadTMVP.lambda1 = 1;
opts.CadTMVP.lambda2 = 0.1;

% Kfold Cross validation
n = size(X{1}, 1);
k_fold = 2;
indices = crossvalind('Kfold', n, k_fold);

fprintf('===================================\n');
for k = 1 : k_fold
    fprintf('Current fold: %d\n', k);
    % Split training data and test data
    test = (indices == k);
    train = ~test;
    for i = 1 : numel(X)
        trainData.X{i} = normalize(X{i}(train, :), 'norm');
        trainData.DX = normalize(DX(train, :), 'norm');
        testData.X{i} = normalize(X{i}(test, :), 'norm');
        testData.DX = normalize(DX(test, :), 'norm');
    end
 % Training step
 % tic;
    [W.CadTMVP{k}, u.CadTMVP(:, k), v.CadTMVP(:, k), w.CadTMVP(:, k)] = FastCadTMVP(trainData, opts.CadTMVP);
    if k ~= k_fold
        fprintf('\n');
    end
end

fprintf('===================================\n');

u.CadTMVP_mean = mean(u.CadTMVP, 2);
v.CadTMVP_mean = mean(v.CadTMVP, 2);
w.CadTMVP_mean = mean(w.CadTMVP, 2);
idx_selected1 = f_sel_top_K_features(abs(u.CadTMVP_mean), top_k_selected);
idx_selected2 = f_sel_top_K_features(abs(v.CadTMVP_mean), top_k_selected);
idx_selected3 = f_sel_top_K_features(abs(w.CadTMVP_mean), top_k_selected);



