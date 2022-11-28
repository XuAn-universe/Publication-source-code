%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Jan 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%%
corr_pawl.rho = [corrY_pawlF.rho; corrY_pawlP.rho];
corr_pawr.rho = [corrY_pawrF.rho; corrY_pawrP.rho];
corr_pawl.p = [corrY_pawlF.p; corrY_pawlP.p];
corr_pawr.p = [corrY_pawrF.p; corrY_pawrP.p];

%%
corr_pawl.rho = [corrX_pawlF.rho; corrX_pawlP.rho];
corr_pawr.rho = [corrX_pawrF.rho; corrX_pawrP.rho];
corr_pawl.p = [corrX_pawlF.p; corrX_pawlP.p];
corr_pawr.p = [corrX_pawrF.p; corrX_pawrP.p];

%%
corr_pawl.rho = [corrZ_pawlF.rho; corrZ_pawlP.rho];
corr_pawr.rho = [corrZ_pawrF.rho; corrZ_pawrP.rho];
corr_pawl.p = [corrZ_pawlF.p; corrZ_pawlP.p];
corr_pawr.p = [corrZ_pawrF.p; corrZ_pawrP.p];

%%
clc;
figure;
subplot(1, 2, 1);
disp('Rho R NI');
paired_plot(corr_pawr.rho(:, 1), corr_pawl.rho(:, 1), 'Corr coef', [], 'dot', [1 2]);
disp('Rho L NI');
paired_plot(corr_pawl.rho(:, 2), corr_pawr.rho(:, 2), 'Corr coef', [], 'dot', [3 4]);
xlim([0.4 4.6]);
subplot(1, 2, 2);
disp('Rho R I');
paired_plot(corr_pawr.rho(:, 3), corr_pawl.rho(:, 3), 'Corr coef', [], 'dot', [1 2]);
disp('Rho L I');
paired_plot(corr_pawl.rho(:, 4), corr_pawr.rho(:, 4), 'Corr coef', [], 'dot', [3 4]);
xlim([0.4 4.6]);

figure;
subplot(1, 2, 1);
disp('Rho R NI');
paired_plot(corr_pawr.rho(corr_pawr.p(:, 1) < 0.05, 1), corr_pawl.rho(corr_pawr.p(:, 1) < 0.05, 1), 'Corr coef', [], 'dot', [1 2]);
disp('Rho L NI');
paired_plot(corr_pawl.rho(corr_pawl.p(:, 2) < 0.05, 2), corr_pawr.rho(corr_pawl.p(:, 2) < 0.05, 2), 'Corr coef', [], 'dot', [3 4]);
xlim([0.4 4.6]);
subplot(1, 2, 2);
disp('Rho R I');
paired_plot(corr_pawr.rho(corr_pawr.p(:, 3) < 0.05, 3), corr_pawl.rho(corr_pawr.p(:, 3) < 0.05, 3), 'Corr coef', [], 'dot', [1 2]);
disp('Rho L I');
paired_plot(corr_pawl.rho(corr_pawl.p(:, 4) < 0.05, 4), corr_pawr.rho(corr_pawl.p(:, 4) < 0.05, 4), 'Corr coef', [], 'dot', [3 4]);
xlim([0.4 4.6]);