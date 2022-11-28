function [t_aligned, signal_aligned, signal_aligned_random] = AlignSignal2Event(t, signal, t_event, pre, post)
nchannel = size(signal, 2);
SampleRate = 1/(mean(diff(t)));
% SampleRate = 19.6; %
pre = ceil(pre*SampleRate);
post = ceil(post*SampleRate);
t_aligned = (-pre:post)'/SampleRate;
if ~isempty(t_event) && ~isempty(signal)
    signal_aligned = nan(pre+post+1, numel(t_event), nchannel);
    signal_aligned_random = nan(size(signal_aligned));
    for i = 1:numel(t_event)
        timestamp = t_event(i);
        [~, id] = min(abs(t-timestamp));
        if id-pre > 0 && id+post <= numel(t)
            signal_aligned(:, i, :) = reshape(signal(id-pre:id+post, :), [pre+post+1 1, nchannel]);
        elseif id+post <= numel(t)
            signal_aligned(end-id-post+1:end, i, :) = reshape(signal(1:id+post, :), [id+post 1, nchannel]);
        elseif id-pre > 0
            signal_aligned(1:numel(t)+pre-id+1, i, :) = reshape(signal(id-pre:numel(t), :), [numel(t)+pre-id+1 1, nchannel]);
        else
            signal_aligned(end-id-post+1:numel(t)+pre-id+1, i, :) = reshape(signal(1:numel(t), :), [numel(t) 1, nchannel]);
        end
        
        rand_index = randi(numel(t), 1);
        if rand_index-pre > 0 && rand_index+post <= numel(t)
            signal_aligned_random(:, i, :) = reshape(signal(rand_index-pre:rand_index+post, :), [pre+post+1 1, nchannel]);
        elseif rand_index+post <= numel(t)
            signal_aligned_random(end-rand_index-post+1:end, i, :) = reshape(signal(1:rand_index+post, :), [rand_index+post 1, nchannel]);
        elseif rand_index-pre > 0
            signal_aligned_random(1:numel(t)+pre-rand_index+1, i, :) = reshape(signal(rand_index-pre:numel(t), :), [numel(t)+pre-rand_index+1 1, nchannel]);
        else
            signal_aligned_random(end-rand_index-post+1:numel(t)+pre-rand_index+1, i, :) = reshape(signal(1:numel(t), :), [numel(t) 1, nchannel]);
        end
    end

%     n = 0;
%     for i = 1:numel(t_event)
%         timestamp = t_event(i);
%         [~, id] = min(abs(t-timestamp));
%         if id-pre > 0 && id+post <= numel(t)
%             n = n+1;
%             signal_aligned(:, n, :) = reshape(signal(id-pre:id+post, :), [pre+post+1 1, nchannel]);
%             rand_index = randi(numel(t)-pre-post, 1);
%             signal_aligned_random(:, n, :) = reshape(signal(rand_index:rand_index+pre+post, :), [pre+post+1 1, nchannel]);
%         end
%     end
end