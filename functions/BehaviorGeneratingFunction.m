function [Sample, Choices, Percept, Correct, Confidences, BeliefF] = BehaviorGeneratingFunction(x, bias, sigma, n, xpx, Px)
% BehaviorGeneratingFunction  Simulate decisions and confidence for a noisy decision agent.
%
%   [Sample, Choices, Percept, Correct, Confidences, BeliefF] = ...
%       BehaviorGeneratingFunction(x, bias, sigma, n, xpx, Px)
%
%   Simulates a decision agent that observes stimulus evidence x with
%   Gaussian perceptual noise, makes binary choices based on the sign of
%   the percept, and computes confidence as P(correct | percept) via the
%   belief function (see MakeBeliefFunction).
%
%   Inputs:
%     x     - stimulus evidence values, one per trial (vector)
%     bias  - perceptual bias (shifts the Gaussian mean)
%     sigma - perceptual noise (Gaussian SD)
%     n     - total number of simulated trials (across all stimulus values)
%     xpx   - percept bin centers for empirical belief function (optional)
%     Px    - percept distribution for empirical belief function (optional)
%             if xpx/Px are omitted, the analytic half-Gaussian belief
%             function is used (valid for uniform stimulus distributions)
%
%   Outputs:
%     Sample      - stimulus value for each simulated trial
%     Choices     - binary choice: +1 (percept > 0) or -1 (percept <= 0)
%     Percept     - sampled percept (stimulus + Gaussian noise) per trial
%     Correct     - logical, true if sign(Percept) == sign(Sample)
%     Confidences - P(correct | percept) for each trial, in [0.5, 1]
%     BeliefF     - function handle: confidence as a function of |percept|
%                   via interpolation over the simulated distribution

%% --- Simulate percepts --------------------------------------------------
x       = x(:);
nStim   = length(x);
repsPerStim = ceil(n / nStim);  % repetitions per stimulus value

% For each stimulus, draw repsPerStim percepts from N(x + bias, sigma)
Sample  = repmat(x, 1, repsPerStim);
Percept = normrnd(repmat(x + bias, 1, repsPerStim), sigma);

% Flatten to trial vectors
Sample  = Sample(:);
Percept = Percept(:);

%% --- Decisions and accuracy ---------------------------------------------
% Choice is determined by the sign of the percept
Choices = sign(Percept);           % +1 or -1
Choices(Choices == 0) = 1;         % break ties arbitrarily

% Correct if percept and stimulus agree in sign
Correct = sign(Sample) == sign(Percept);

% For zero-evidence stimuli, accuracy is chance (randomly assign correct)
Correct(Sample == 0) = rand(sum(Sample == 0), 1) < 0.5;

%% --- Confidence via belief function -------------------------------------
if nargin < 5 || isempty(xpx)
    % Analytic form: assumes uniform stimulus distribution
    Confidences = MakeBeliefFunction(Percept, sigma);
else
    % Empirical form: uses simulated percept distribution Px at xpx
    Confidences = MakeBeliefFunction(Percept, sigma, xpx, Px, Correct);
end

% Return belief function as an interpolant over the simulated distribution
BeliefF = @(query) interp1(abs(Percept), Confidences, query);

end
