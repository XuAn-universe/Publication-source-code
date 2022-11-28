%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, April 2020
% xan@cshl.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%%
clear;
clc;
[filename, pathname] = uigetfile('*.mat', 'Pick all of the individual data set', 'MultiSelect', 'on');
if isequal(filename, 0)
    return;
end
names1 = {'zlight', 'zspeedlight', 'xyzspeedlight', 'zaccelerationlight', 'xyzaccelerationlight'};
names2 = {'dzlight', 'dxyzlight', 'orientationxylight'};
N = numel(filename);
load([pathname filename{1}]);
names = fieldnames(result);
ninhibition = numel(eval(['result.' names{1}]));
switch numel(names)
    case 22
        nvar = numel(names1);
        for i = 1:nvar
            eval([names1{i} '_mean_NI = zeros(N, ninhibition);']);
            eval([names1{i} '_mean_I = zeros(N, ninhibition);']);
            eval([names1{i} '_var_NI = zeros(N, ninhibition);']);
            eval([names1{i} '_var_I = zeros(N, ninhibition);']);
        end
    case 12
        nvar = numel(names2);
        for i = 1:nvar
            eval([names2{i} '_mean_NI = zeros(N, ninhibition);']);
            eval([names2{i} '_mean_I = zeros(N, ninhibition);']);
            eval([names2{i} '_var_NI = zeros(N, ninhibition);']);
            eval([names2{i} '_var_I = zeros(N, ninhibition);']);
        end
end
for i = 1:N
    load([pathname filename{i}]);
    switch numel(names)
        case 22
            for j = 1:nvar
                eval([names1{j} '_mean_NI(i, :) = result.' names1{j} '_mean_NI;']);
                eval([names1{j} '_mean_I(i, :) = result.' names1{j} '_mean_I;']);
                eval([names1{j} '_var_NI(i, :) = result.' names1{j} '_var_NI;']);
                eval([names1{j} '_var_I(i, :) = result.' names1{j} '_var_I;']);
            end
        case 12
            for j = 1:nvar
                eval([names2{j} '_mean_NI(i, :) = result.' names2{j} '_mean_NI;']);
                eval([names2{j} '_mean_I(i, :) = result.' names2{j} '_mean_I;']);
                eval([names2{j} '_var_NI(i, :) = result.' names2{j} '_var_NI;']);
                eval([names2{j} '_var_I(i, :) = result.' names2{j} '_var_I;']);
            end
    end
end

for i = 1:ninhibition
    figure('Name', ['Inhibition #' num2str(i)]);
    switch numel(names)
        case 22
            for j = 1:nvar
                subplot(2, nvar, j);
                disp([names1{j} ' mean']);
                paired_plot(eval([names1{j} '_mean_NI(:, i)']), eval([names1{j} '_mean_I(:, i)']), [names1{j} ' mean'], filename, 'bar');
                subplot(2, nvar, nvar+j);
                disp([names1{j} ' var']);
                paired_plot(eval([names1{j} '_var_NI(:, i)']), eval([names1{j} '_var_I(:, i)']), [names1{j} ' var'], filename, 'bar');
            end
        case 12
            for j = 1:nvar
                subplot(2, nvar, j);
                disp([names2{j} ' mean']);
                paired_plot(eval([names2{j} '_mean_NI(:, i)']), eval([names2{j} '_mean_I(:, i)']), [names2{j} ' mean'], filename, 'bar');
                subplot(2, nvar, nvar+j);
                disp([names2{j} ' var']);
                paired_plot(eval([names2{j} '_var_NI(:, i)']), eval([names2{j} '_var_I(:, i)']), [names2{j} ' var'], filename, 'bar');
            end
    end
end