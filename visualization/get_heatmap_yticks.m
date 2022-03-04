function [y_ticks, ylabels] = get_heatmap_yticks(percent_replicated, varargin)

    percent_replicated = percent_replicated .* 100;

    if isempty(varargin)
        ylabels = 10:10:100;
    else
        ylabels = varargin{1};
    end
    
    y_ticks = NaN(length(ylabels), 1);
    for y = 1:length(ylabels) - 1
        is_in_bin = percent_replicated > ylabels(y) & ...
            percent_replicated < ylabels(y + 1);
        if sum(is_in_bin) > 0
            y_ticks(y, 1) = find(is_in_bin, 1);
        end
    end

    keep_label = y_ticks > 1 & y_ticks < length(percent_replicated);
    y_ticks = y_ticks(keep_label);
    ylabels = ylabels(keep_label);
end
