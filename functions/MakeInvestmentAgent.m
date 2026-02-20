% MakeInvestmentAgent
% Simulates a decision-making agent with confidence-guided time investment.
%
% The agent:
%   (1) receives sensory evidence x drawn from a uniform distribution
%   (2) forms an internal percept x_hat with Gaussian noise (s.d. sigma)
%   (3) makes a binary decision d by thresholding x_hat at bias b
%   (4) computes confidence c = P(correct | x_hat) via the analytic
%       half-Gaussian belief function 
%   (5) invests time ti = m(c) via the optimal mapping function,
%       with sub-optimal noise controlled by efficacy parameter alpha
%
% Output is saved to data/InvestmentAgentData.mat as structs Data and Meta.
%
% Ott et al. 2026
% torben.ott@bccn-berlin.de
%

SCRIPTPATH = fileparts(mfilename('fullpath'));
DATAPATH   = fullfile(SCRIPTPATH, '..', 'data');
DATAFILE   = 'InvestmentAgentData.mat';

rng(1);  % fix random seed for reproducibility

%% --- Agent parameters ---------------------------------------------------
sigma = 0.4;  % perceptual noise SD of Gaussian likelihood P(x_hat | x)
b     = 0;    % decision threshold (b=0: unbiased)
alpha = 0.7;  % confidence efficacy: fraction of trials using optimal TI
              % (alpha=1: fully optimal; alpha=0: fully random)

%% --- Simulate perceptual decisions --------------------------------------
N     = 100000;                      % number of trials
x     = rand(1, N) * 2 - 1;         % stimulus evidence, uniform on [-1, 1]
x_hat = x + randn(1, N) * sigma;    % noisy percept
d     = x_hat > b;                  % binary decision (true = positive evidence)
correct = (d*2 - 1) == sign(x);     % correct if decision matches stimulus sign

%% --- Confidence and optimal time investment -----------------------------
% Analytic belief function: P(correct | x_hat) for Gaussian noise (Lak et al. 2014)
c = 0.5 * (1 + erf(abs(x_hat) / (sigma * sqrt(2))));

% Optimal mapping function ti = m(c) (Lak et al. 2014)
tau   = 1.5;   % time scaling parameter
kappa = 0.08;  % reward rate parameter
offset = 0.6;  % offset of truncated reward delay distribution
ti_optimal = tau .* log((c*0.9 ./ (1 - c*0.9)) .* ((1 - kappa*tau) / (kappa*tau))) + offset;

%% --- Confidence efficacy noise ------------------------------------------
% On a fraction (1-alpha) of trials, time investment is replaced by a
% random draw from the same distribution (randreplace model)
ti          = ti_optimal;
replaceIdx  = rand(size(ti)) > alpha;
randomDraw  = randi(length(ti), sum(replaceIdx), 1);
ti(replaceIdx) = ti_optimal(randomDraw);

%% --- Behavioral task structure: reward interruption --------------------
% In correct trials, reward is delivered after a random delay drawn from a
% truncated exponential. If the reward arrives before the agent stops
% waiting (ti), the trial is rewarded. Probe trials (omissions) are never
% rewarded regardless of outcome.
reward_delay_dist = truncate(makedist('Exponential', 'mu', 1.5), 0.6, 8);
reward_delay      = random(reward_delay_dist, 1, N);
omission          = rand(1, N) < 0.1;              % 10% probe (omission) trials
unrewarded        = (~correct) | omission | (reward_delay > ti);

%% --- Save output --------------------------------------------------------
Data.x          = x;
Data.x_hat      = x_hat;
Data.c          = c;
Data.d          = d;
Data.correct    = correct;
Data.ti_optimal = ti_optimal;
Data.ti         = ti;
Data.omission   = omission;
Data.unrewarded = unrewarded;

Meta.sigma = sigma;
Meta.b     = b;
Meta.alpha = alpha;
Meta.N     = N;
Meta.tau   = tau;
Meta.kappa = kappa;

save(fullfile(DATAPATH, DATAFILE), 'Data', 'Meta');
fprintf('Saved agent data to %s\n', fullfile(DATAPATH, DATAFILE));
