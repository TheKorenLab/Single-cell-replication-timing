function mu = fit_two_component_poisson_mixture(seq)

    mixed_poisson = @(seq, proportion1, logmu1, logmu2) ...
        proportion1 * poisspdf(seq, exp(logmu1)) + (1 - proportion1) * poisspdf(seq, exp(logmu2));

    initial_guess = log([nanmedian(seq) 2*nanmedian(seq)]);
    
    if any(~isfinite(initial_guess)) || any(initial_guess == 0) % median count after aggregation = 0
        mu = [];
        return
    end
    
    p_hat = mle(seq, 'pdf', mixed_poisson, 'start', [0.5 initial_guess], 'LowerBound', [0 0 0], ...
        'UpperBound', [1 Inf Inf], 'Options', statset('MaxIter', 1e3, 'MaxFunEvals', 1e3, ...
        'FunValCheck', 'off'));
    
    mu = sort(exp(p_hat(2:3)), 'ascend');
    
end
