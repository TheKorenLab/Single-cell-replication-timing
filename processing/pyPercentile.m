function percentile = pyPercentile(vector, p)

    if isempty(vector)
        percentile = NaN;
        return
    end

    vector = sort(vector,'ascend');
    vector = vector(~isnan(vector));
    position = (p/100) * (length(vector)-1) + 1;
    if floor(position) < length(vector)
        percentile = interp1(vector, position);
    else
        percentile = vector(position);
    end
    
end
