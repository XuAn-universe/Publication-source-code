function chew2handleID = compute_boutstartHMM(raise_interval, distancedata)
window_smooth = 5;
chew2handleID = raise_interval(3, :);
distancedata = mean(distancedata);
distancedata = smoothdata(distancedata, 'gaussian', window_smooth);
distancespeed = [NaN diff(distancedata)];
distancespeed = fillmissing(distancespeed, 'spline');
[~, ~, rise2fall, ~] = ZeroCrossingDetection(distancespeed);
for i = 1:size(raise_interval, 2)
    rise2fall_temp = rise2fall(rise2fall(:, 1) >= raise_interval(1, i) & rise2fall(:, 2) <= raise_interval(2, i), :);
    if ~isempty(rise2fall_temp)
        riseheight = distancedata(rise2fall_temp(:, 2))-distancedata(rise2fall_temp(:, 1));
        [~, ID] = max(riseheight);
        chew2handleID(i) = rise2fall_temp(ID(end), 1);
    end
end