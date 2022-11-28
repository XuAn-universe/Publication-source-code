function timestamp = get_laser_timestamp_sound(data, state, n, pmethod)
trials = length(data);
switch state
    case 'start'
        for i = 1:trials
            if size(data(i).laser_timestamps, 1) == 1
                timestamp(i, :) = data(i).laser_timestamps(1);
            else
                timestamp(i, :) = data(i).laser_timestamps(1, :);
            end
        end
    case 'stop'
        for i = 1:trials
            if size(data(i).laser_timestamps, 1) == 1
                timestamp(i, :) = data(i).laser_timestamps(2);
            else
                timestamp(i, :) = data(i).laser_timestamps(2, :);
            end
        end
end

timestamp = eval([pmethod '(timestamp);']);
timestamp = timestamp(n);
