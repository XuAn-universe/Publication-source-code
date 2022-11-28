function [risingindex, fallingindex, rise2fall, fall2rise] = ZeroCrossingDetection(data)
data(data == 0) = eps;
data_shift = circshift(data, -1);
zcindex = find(data.*data_shift < 0);
zcindex(zcindex == numel(data)) = [];
risingindex = zcindex(data(zcindex+1) > 0);
fallingindex = zcindex(data(zcindex+1) < 0);
try
    if risingindex(1) < fallingindex(1)
        rise2fall = zeros(numel(fallingindex), 2);
        rise2fall(:, 1) = risingindex(1:numel(fallingindex));
        rise2fall(:, 2) = fallingindex;
        
        fall2rise = zeros(numel(risingindex)-1, 2);
        fall2rise(:, 1) = fallingindex(1:numel(risingindex)-1);
        fall2rise(:, 2) = risingindex(2:end);
    else
        rise2fall = zeros(numel(fallingindex)-1, 2);
        rise2fall(:, 1) = risingindex(1:numel(fallingindex)-1);
        rise2fall(:, 2) = fallingindex(2:end);
        
        fall2rise = zeros(numel(risingindex), 2);
        fall2rise(:, 1) = fallingindex(1:numel(risingindex));
        fall2rise(:, 2) = risingindex;
    end
catch
    rise2fall = [];
    fall2rise = [];
end
