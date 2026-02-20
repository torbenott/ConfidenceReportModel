function c = MakeBeliefFunction(p, s, x, Px, Correct)
% MakeBeliefFunction  Map percepts to confidence via P(correct | percept).
%
%   c = make_belief_function(p, s)
%   c = make_belief_function(p, s, x, Px, Correct)
%
%   Computes the belief function c(p, sigma) = P(correct | percept, choice),
%   which links the decision maker's internal percept to their confidence.
%
%   Analytic form (uniform DV distribution):
%     Under Gaussian noise and a uniform stimulus distribution, P(correct)
%     has a closed-form solution as a half-Gaussian:
%       c = 0.5 * (1 + erf(|p| / (sigma * sqrt(2))))
%
%   Monte Carlo form (arbitrary DV distribution):
%     When the stimulus distribution is non-uniform, the analytic form is
%     no longer valid. Instead, P(correct | percept) is estimated empirically
%     by binning simulated percepts and averaging correct outcomes per bin.
%
%   Inputs:
%     p       - percept values (internal evidence), scalar or vector
%     s       - perceptual sigma (sensitivity), used in analytic form
%     x       - percept bin centers from simulation (unused placeholder, for API consistency)
%     Px      - percept distribution from simulation; if empty, use analytic form
%     Correct - binary correct/error outcomes for each simulated percept in p
%
%   Output:
%     c       - confidence values in [0.5, 1], same size as p

% Maximum confidence (can be < 1 for probe/catch trial variants)
pmax = 1;

if nargin < 3 || isempty(Px)
    %% --- Analytic form (uniform DV distribution) ------------------------
    % Valid when stimulus distribution is symmetric and uniform.
    % Ref: Lak et al. (2014)
    c = pmax * 0.5 * (1 + erf(abs(p) / (s * sqrt(2))));

else
    %% --- Monte Carlo form (arbitrary DV distribution) -------------------
    % Bin simulated percepts and compute mean accuracy per bin as a proxy
    % for P(correct | percept). Assumes symmetric percept distribution.
    p = abs(p);
    edges = linspace(-max(abs(p)) - 1000*eps, max(abs(p)) + 1000*eps, 501);
    p_i   = discretize(p, edges);
    bins  = unique(p_i);

    c = nan(size(p));
    for k = 1:length(bins)
        inBin      = (p_i == bins(k));
        c(inBin)   = mean(Correct(inBin));
    end
end

end
