function Data = ConfidenceReportModel(Data, Params)
% ConfidenceReportModel  Fit and evaluate the confidence-guided time investment model.
%
%   Data = ConfidenceReportModel(Data, Params)
%
%   This function implements a two-step procedure:
%
%   Step 1 — Confidence mapping: estimates the mapping function ti = m(c)
%     that converts statistically optimal confidence c (P(correct|percept))
%     to time investment ti, by matching the model confidence distribution
%     to the empirical waiting time distribution (assuming monotonicity).
%
%   Step 2 — Confidence efficacy (alpha): fits or applies a single
%     sub-optimality parameter alpha that controls how faithfully the
%     decision maker uses their confidence to guide time investment.
%     alpha = 1: optimal use of confidence; alpha = 0: random time investment.
%
%   Required Data fields:
%     .DV      - decision variable (evidence) per trial
%     .Conf    - behavioral confidence report (e.g. waiting time) per trial
%     .Outcome - binary accuracy (1 = correct, 0 = error)
%     .Fit     - struct with psychometric fit results (.Sigma, .Bias)
%
%   Optional Data fields:
%     .Omission - logical flag for reward-omission trials (default: all true)
%
%   Key Params fields:
%     .plotmapping        - plot the mapping function m(c) (default: true)
%     .Resample           - resample omission trials (default: false)
%     .MinConf / .MaxConf - confidence range (default: 0 / 15)
%     .Hmapping           - axes handle for mapping plot
%     .Halpha             - axes handle for alpha/GOF plot
%     .Model.N            - number of simulated trials (default: 100000)
%     .Model.nhist        - bins for KS density estimation (default: 30)
%     .Model.alpha        - fixed alpha if not fitting (default: 1)
%     .Model.fit_alpha    - fit alpha by MLE (default: false)
%     .Model.Nalpha       - alpha grid size for fitting (default: 10)
%     .Model.Nfitboot     - bootstrap iterations for alpha CI (default: 199)
%     .Model.efficacy_type- noise model for alpha ('randreplace')
%     .Model.Shuffle      - shuffle model confidence as control (default: false)
%     .Model.GOFplot      - plot log-likelihood curve (default: true)
%     .Model.UniformDV    - assume uniform DV distribution (default: true)
%
% Ott et al. 2026
% torben.ott@bccn-berlin.de
%

%% --- Default general parameters -----------------------------------------
if ~isfield(Params, 'plotmapping'),  Params.plotmapping  = true;  end
if ~isfield(Params, 'Resample'),     Params.Resample     = false; end
if ~isfield(Params, 'UniqueDV'),     Params.UniqueDV     = false; end
if ~isfield(Params, 'MinConf'),      Params.MinConf      = 0;     end
if ~isfield(Params, 'MaxConf'),      Params.MaxConf      = 100;   end
if ~isfield(Params, 'CrossVal'),     Params.CrossVal     = false; end
if ~isfield(Params, 'WriteTI'),      Params.WriteTI      = [];    end

% Determine whether to create new figures or plot into provided axes
if isfield(Params, 'Hmapping') && ~isempty(Params.Hmapping)
    axes(Params.Hmapping(1));
    Params.create_match_fig = false;
else
    Params.create_match_fig = true;
end
if isfield(Params, 'Halpha') && ~isempty(Params.Halpha)
    axes(Params.Halpha(1));
    Params.create_alpha_fig = false;
else
    Params.create_alpha_fig = true;
end

%% --- Default model parameters -------------------------------------------
if ~isfield(Params.Model, 'N'),              Params.Model.N              = 100000;       end
if ~isfield(Params.Model, 'nhist'),          Params.Model.nhist          = 30;           end
if ~isfield(Params.Model, 'alpha'),          Params.Model.alpha          = 1;            end
if ~isfield(Params.Model, 'fit_alpha'),      Params.Model.fit_alpha      = false;        end
if ~isfield(Params.Model, 'Nalpha'),         Params.Model.Nalpha         = 10;           end
if ~isfield(Params.Model, 'Nfitboot'),       Params.Model.Nfitboot       = 199;          end
if ~isfield(Params.Model, 'efficacy_type'),  Params.Model.efficacy_type  = 'randreplace';end
if ~isfield(Params.Model, 'Shuffle'),        Params.Model.Shuffle        = false;        end
if ~isfield(Params.Model, 'GOFplot'),        Params.Model.GOFplot        = true;         end
if ~isfield(Params.Model, 'Label'),          Params.Model.Label          = true;         end
if ~isfield(Params.Model, 'Lapse'),          Params.Model.Lapse          = false;        end
if ~isfield(Params.Model, 'UniformDV'),      Params.Model.UniformDV      = true;         end
if ~isfield(Params.Model, 'discrete')
    Params.Model.discrete = length(unique(Data.Conf(~isnan(Data.Conf)))) < 20;
end

%% --- Data validation ----------------------------------------------------
if ~isfield(Data, 'Omission')
    Data.Omission = true(length(Data.DV), 1);
end
if ~isfield(Data, 'Fit')
    error('ConfidenceReportModel: no psychometric fit found in Data.Fit.');
end

%% -----------------------------------------------------------------------
%  STEP 1a: Simulate percepts and model confidence
%  Run BehaviorGeneratingFunction to generate N simulated trials, producing
%  percepts drawn from N(DV + bias, sigma), binary choices, accuracy, and
%  confidence via the belief function P(correct | percept).
%% -----------------------------------------------------------------------
sigma = Data.Fit.Sigma;
bias  = Data.Fit.Bias;
n     = Params.Model.N;
x     = Data.DV;

% Build percept distribution for belief function estimation.
% If UniformDV=true, the analytic half-Gaussian belief function is used.
% Otherwise, the empirical DV distribution is passed to MakeBeliefFunction.
if Params.Model.UniformDV
    xpx = [];
    Px  = [];
else
    dvBinEdges = linspace(-max(abs(Data.DV)), max(abs(Data.DV)), 50);
    xpx = (dvBinEdges(1:end-1) + dvBinEdges(2:end)) / 2;  % bin centers
    Px  = histcounts(Data.DV, dvBinEdges, 'Normalization', 'probability');
end

fprintf('Generating model behavior...');
tic
[SamplesModel, ChoicesModel, PerceptModel, CorrectModel, ConfidenceModel, BeliefF] = ...
    BehaviorGeneratingFunction(x, bias, sigma, n, xpx, Px);
fprintf(' Done. (%.2f s)\n', toc);

% Shuffle control: randomly permute model confidence (breaks DV-confidence link)
if Params.Model.Shuffle
    shuffleIdx     = randperm(length(ConfidenceModel));
    ConfidenceModel = ConfidenceModel(shuffleIdx);
    PerceptModel    = PerceptModel(shuffleIdx);
end

% Store the zero-alpha (optimal) model as a reference baseline
Data.Model.Zero.DV      = SamplesModel;
Data.Model.Zero.ChoiceA = ChoicesModel;
Data.Model.Zero.Percept = PerceptModel;
Data.Model.Zero.Correct = CorrectModel;
Data.Model.Zero.Conf    = ConfidenceModel;
Data.Model.Zero.BeliefF = BeliefF;

%% -----------------------------------------------------------------------
%  STEP 1b: Estimate confidence mapping function ti = m(c)
%  Calibrates model confidence to time investment units by matching the
%  model confidence CDF to the empirical waiting time distribution (assuming monotonicity).
%% -----------------------------------------------------------------------
[ConfmatchOptimal, Data] = calibrate_confidence(abs(PerceptModel), Data, Params);

%% -----------------------------------------------------------------------
%  STEP 2: Confidence efficacy (alpha)
%  Fits or applies the sub-optimality parameter alpha, which controls how
%  much random noise is mixed into the confidence-based time investment.
%% -----------------------------------------------------------------------
if Params.Model.fit_alpha
    fprintf('Fitting confidence efficacy (alpha)...');
    tic
    alphaGrid = linspace(0, 1, Params.Model.Nalpha);
    Data = fit_confidence_efficacy(ConfmatchOptimal, CorrectModel, SamplesModel, ...
                                   ConfidenceModel, Data, Params, alphaGrid);
    alpha = Data.Model.alpha_best;
    fprintf(' Done. (%.2f s)\n', toc);
else
    % Use fixed alpha supplied by caller (alpha=1: optimal, alpha=0: random)
    alpha = Params.Model.alpha;
end

% Re-run mapping with efficacy noise applied at the chosen alpha
ConfmatchAlpha = confidence_efficacy(ConfmatchOptimal, alpha, Data, Params);

%% --- Assemble model output ----------------------------------------------
Data.Model.Conf    = ConfmatchAlpha;
Data.Model.DV      = SamplesModel;
Data.Model.ChoiceA = ChoicesModel == 1;
Data.Model.Outcome = CorrectModel == 1;

end % main function


%% =========================================================================
function [ModelConfCal, Data] = calibrate_confidence(ModelConf, Data, Params)
% calibrate_confidence  Estimate the mapping function ti = m(c).
%
%   Maps model confidence (P(correct|percept), in [0.5,1]) to time investment
%   units by matching quantiles of the model confidence distribution to the
%   empirical waiting time (confidence) distribution on omission trials.
%   This enforces a monotonic mapping that preserves the rank-order of confidence.
%
%   Inputs:
%     ModelConf - absolute model percept values (used as proxy for confidence)
%     Data      - data struct with .Conf and .Omission fields
%     Params    - parameter struct
%
%   Outputs:
%     ModelConfCal - model confidence values remapped to time investment units
%     Data         - updated with Data.Stats.ModelFit.ConfM and Data.Confp

% Empirical waiting time distribution on omission trials
DataConf  = Data.Conf(Data.Omission);
edgehist  = linspace(min(DataConf), max(DataConf), Params.Model.nhist + 1);
centerhist = (edgehist(1:end-1) + edgehist(2:end)) / 2;

% KS density estimate of empirical waiting time distribution
[ConfDensity, ~] = ksdensity(DataConf, centerhist);
ConfHistNorm     = histcounts(DataConf, edgehist, 'Normalization', 'probability');

% Record peak of waiting time distribution
[~, peakIdx] = max(ConfDensity);
Data.Stats.ModelFit.ConfM = centerhist(peakIdx);

% Compute CDF of empirical waiting time, floored to avoid interpolation issues
ConfCDF = cumsum(ConfDensity);
ConfCDF = ConfCDF / ConfCDF(end);
ConfCDF = floor(ConfCDF * 100000) / 100000;  % floor to 5 decimal places

% Map model confidence quantiles to waiting time bin edges
q = quantile(ModelConf, [0, ConfCDF]);

% Ensure interpolation points are unique (duplicate quantiles can arise with
% discrete or very peaked distributions)
if length(unique(q)) ~= length(q)
    [q, uniqueIdx]  = unique(q);
    edgehist        = edgehist(uniqueIdx);
    centerhist      = (edgehist(1:end-1) + edgehist(2:end)) / 2;
    ConfHistNorm    = ConfHistNorm(uniqueIdx(1:end-1));
end

% Forward mapping: model confidence → time investment
ModelConfCal     = interp1(q, edgehist, ModelConf);

% Inverse mapping: time investment → model confidence (for Data.Confp)
Data.Confp = interp1(edgehist, q, Data.Conf);

%% --- Plot mapping function m(c) -----------------------------------------
if Params.plotmapping
    if Params.create_match_fig
        figure; hold on;
        ax = axes;
    else
        ax = Params.Hmapping(1);
    end
    axes(ax); hold on;

    % Evaluate mapping over a grid of confidence values via belief function
    [xx, vv] = unique(Data.Model.Zero.BeliefF(q));
    mapping_f = @(c) interp1(xx, edgehist(vv), c);
    c_grid  = linspace(0.5, 1, 21);
    ti_grid = mapping_f(c_grid);
    plot(c_grid, ti_grid, 'Color', [0, 0, 0], 'LineWidth', 2);

    xlim([0.5, 1]);
    set(gca, 'XTick', [0.5, 0.75, 1]);
    if Params.Model.Label
        ylabel('Time investment (s)');
        xlabel('Confidence');
        title('Mapping function m(c)');
    end
end

end % calibrate_confidence


%% =========================================================================
function Data = fit_confidence_efficacy(ConfmatchOptimal, CorrectModel, SamplesModel, ...
                                         ConfidenceModel, Data, Params, alphaGrid)
% fit_confidence_efficacy  Fit the confidence efficacy parameter alpha by MLE.
%
%   Finds the alpha in [0,1] that maximizes the log-likelihood of the observed
%   confidence distribution given the model, using a polynomial fit + fmincon
%   to find a smooth optimum. Bootstrap resampling provides 95% CIs.
%
%   alpha = 1: subject uses confidence optimally (low noise)
%   alpha = 0: subject time investment is random (high noise / no confidence use)
%
%   Results stored in Data.Model:
%     .alpha_best             - MLE estimate of alpha
%     .alpha_best_ci_low/high - 95% bootstrap CI bounds
%     .p_alpha_0              - bootstrap p-value for alpha > 0
%     .logL                   - full log-likelihood matrix (boot x alpha)

% Fixed binning parameters for the likelihood computation
NBinDV        = 6;
NBinConf      = 10;
MinCondPoints = 10;

% Remove trials where calibrated confidence is undefined (NaN after mapping)
validModel      = ~isnan(ConfmatchOptimal);
ConfidenceModel = ConfmatchOptimal(validModel);
CorrectModel    = CorrectModel(validModel);
SamplesModel    = SamplesModel(validModel);

% Valid data trials
validData  = ~isnan(Data.DV) & ~isnan(Data.Conf) & ~isnan(Data.Outcome);
DataConf   = Data.Conf(validData);
DataOutcome= Data.Outcome(validData);
DataDV     = abs(Data.DV(validData));

% Bin edges for confidence and DV
ConfBin = linspace(min(DataConf) - 0.0001, max(DataConf) + 0.0001, NBinConf + 1);
DVBin   = linspace(0, max(abs(Data.DV)) + 0.0001, NBinDV + 1);

% Discretize data and model into bins
DataConfIndex  = discretize(DataConf,   ConfBin);
DataDVIndex    = discretize(DataDV,     DVBin);
ModelConfIndex = discretize(ConfidenceModel, ConfBin);
ModelDVIndex   = discretize(SamplesModel,    DVBin);

% Marginal confidence distribution of the model (for reference / random baseline)
pModelMarginal = histcounts(ModelConfIndex, 0.5:1:NBinConf+1);
pModelMarginal(pModelMarginal == 0) = 1;         % avoid log(0)
pModelMarginal = pModelMarginal / sum(pModelMarginal);

%% --- Bootstrap MLE over alpha grid --------------------------------------
nboot  = Params.Model.Nfitboot;
nalpha = length(alphaGrid);
logL   = nan(nboot + 1, nalpha);   % data log-likelihood
logLM  = nan(nboot + 1, nalpha);   % model-to-model log-likelihood (for recovery CI)

parfor n = 1:nboot + 1
    % First iteration (n=1) uses full data; subsequent iterations resample
    if n > 1
        % Bootstrap resample data and model trials independently
        ridxData  = randi(length(DataDVIndex),  length(DataDVIndex),  1);
        ridxModel = randi(length(ModelDVIndex), length(DataDVIndex),  1);
        DataDVIndex_boot   = DataDVIndex(ridxData);
        DataConfIndex_boot = DataConfIndex(ridxData);
        DataOutcome_boot   = DataOutcome(ridxData);
        ModelDVIndex_boot  = ModelDVIndex(ridxModel);
        ModelConfIndex_boot= ModelConfIndex(ridxModel);
        CorrectModel_boot  = CorrectModel(ridxModel);
    else
        DataDVIndex_boot   = DataDVIndex;
        DataConfIndex_boot = DataConfIndex;
        DataOutcome_boot   = DataOutcome;
        ModelDVIndex_boot  = ModelDVIndex;
        ModelConfIndex_boot= ModelConfIndex;
        CorrectModel_boot  = CorrectModel;
    end

    for a = 1:nalpha
        alpha = alphaGrid(a);

        % Apply efficacy noise at this alpha level and re-bin
        ConfidenceModela  = confidence_efficacy(ConfidenceModel, alpha, Data, Params);
        ModelConfIndexa   = discretize(ConfidenceModela, ConfBin);

        % Log-likelihood of data given model at this alpha
        logL(n, a) = monte_carlo_model_likelihood( ...
            DataDVIndex_boot, DataConfIndex_boot, DataOutcome_boot, ...
            ModelDVIndex, ModelConfIndexa, CorrectModel, ...
            DVBin, ConfBin, MinCondPoints);

        % Model-to-model log-likelihood (for recovery interval)
        logLM(n, a) = monte_carlo_model_likelihood( ...
            ModelDVIndex_boot, ModelConfIndex_boot, CorrectModel_boot, ...
            ModelDVIndex, ModelConfIndexa, CorrectModel, ...
            DVBin, ConfBin, MinCondPoints);
    end
end

%% --- Find best alpha via polynomial fit of the LL curve -----------------
% Fit a degree-4 polynomial to smooth the LL curve, then minimize its negative
polyCoeffs = polyfit(alphaGrid, logL(1,:), 4);
negPolyLL  = @(a) -(polyCoeffs(1).*a.^4 + polyCoeffs(2).*a.^3 + ...
                    polyCoeffs(3).*a.^2 + polyCoeffs(4).*a    + polyCoeffs(5));
alphaBest  = fmincon(negPolyLL, 0.99, [], [], [], [], 0, 1);
LLbest     = -negPolyLL(alphaBest);

Data.Model.alpha_best = alphaBest;
Data.Model.logL       = logL;

%% --- Bootstrap confidence intervals for alpha ---------------------------
[~, bootBestIdx] = max(logL, [], 2);
AlphaBoot = sort(alphaGrid(bootBestIdx));

Data.Model.alpha_best_ci_low  = AlphaBoot(ceil((nboot+1) * 0.025));
Data.Model.alpha_best_ci_high = AlphaBoot(ceil((nboot+1) * 0.975));
Data.Model.p_alpha_0          = find(AlphaBoot > 0, 1, 'first') / (nboot + 1);

% Model recovery interval (model-to-model bootstrap)
[~, bootBestIdxM] = max(logLM, [], 2);
AlphaBootM = sort(alphaGrid(bootBestIdxM));
Data.Model.alpha_best_model_ci_low  = AlphaBootM(ceil((nboot+1) * 0.025));
Data.Model.alpha_best_model_ci_high = AlphaBootM(ceil((nboot+1) * 0.975));

%% --- Plot log-likelihood curve ------------------------------------------
if Params.Model.GOFplot
    % Bootstrap CI band over alpha values
    logLLow  = nan(1, nalpha);
    logLHigh = nan(1, nalpha);
    for a = 1:nalpha
        bootSorted   = sort(logL(:, a));
        logLLow(a)   = bootSorted(ceil((nboot+1) * 0.025));
        logLHigh(a)  = bootSorted(ceil((nboot+1) * 0.975));
    end

    if Params.create_alpha_fig
        figure('Color', [1,1,1]);
        ax = axes;
    else
        ax = Params.Halpha(1);
        axes(ax);
    end
    hold on

    % Shaded CI band for bootstrap LL
    hCI = fill([alphaGrid, alphaGrid(end:-1:1)], [logLLow, logLHigh(end:-1:1)], [.8,.8,.8]);
    hCI.EdgeColor = 'none';

    % Y-axis limits (sign-aware padding)
    sig    = [sign(min(logLLow)), sign(max(logLHigh))];
    yLimits = [min(logLLow) * (1 - sig(1)*0.05), max(logLHigh) * (1 + sig(2)*0.05)];

    % Shaded model recovery interval
    hRec = fill([Data.Model.alpha_best_model_ci_low,  Data.Model.alpha_best_model_ci_high, ...
                 Data.Model.alpha_best_model_ci_high,  Data.Model.alpha_best_model_ci_low], ...
                [yLimits(1), yLimits(1), yLimits(2), yLimits(2)], [.8,.8,.8]);
    hRec.EdgeColor = 'none';

    % Polynomial LL curve
    alphaFine = linspace(min(alphaGrid), max(alphaGrid), 100);
    plot(alphaFine, -negPolyLL(alphaFine), '-k');

    % Best alpha marker and CI
    plot(alphaBest, LLbest, '.r', 'MarkerSize', 12);
    plot([alphaBest, alphaBest], [min(get(gca,'YLim')), LLbest], '--r');
    plot([Data.Model.alpha_best_ci_low, Data.Model.alpha_best_ci_high], ...
         [LLbest, LLbest], '-r', 'LineWidth', 2);

    if Params.Model.Label
        xlabel('Time investment efficiency (\alpha)');
        ylabel('Sum(log L) / N');
    end
    title(['Log-likelihood (\alpha_{best} = ', num2str(round(alphaBest*100)/100), ')']);
    leg = legend(ax.Children(3:5), {'\alpha_{best}', 'LL polyfit', 'LL 95% CI'});
    leg.Location = 'northwest';
end

end % fit_confidence_efficacy


%% =========================================================================
function L = monte_carlo_model_likelihood(DataDVIndex, DataConfIndex, DataOutcome, ...
                                           ModelDVIndex, ModelConfIndex, ModelOutcome, ...
                                           DVBin, ConfBin, MinCondPoints)
% monte_carlo_model_likelihood  Compute normalized log-likelihood of data given model.
%
%   Estimates L = sum(log P(conf | DV, decision)) / N over all data trials,
%   where P(conf | DV, decision) is estimated from the model's simulated
%   confidence distribution within each (DV bin, correct/error) cell.
%
%   A small probability floor (1/N_cell) is applied to avoid log(0) for
%   confidence bins not visited by the model.
%
%   Inputs:
%     DataDVIndex   - DV bin index for each data trial
%     DataConfIndex - confidence bin index for each data trial
%     DataOutcome   - binary outcome (1=correct, 0=error) for each data trial
%     ModelDVIndex  - DV bin index for each model trial
%     ModelConfIndex- confidence bin index for each model trial
%     ModelOutcome  - binary outcome for each model trial
%     DVBin         - DV bin edges (used only for size)
%     ConfBin       - confidence bin edges (used only for size)
%     MinCondPoints - minimum model trials per cell; fewer → skip with warning

NBinDV   = length(DVBin)   - 1;
NBinConf = length(ConfBin) - 1;

% Reindex to integer bins for histogram counting
DVBinIdx   = 0.5:(0.5 + NBinDV);
ConfBinIdx = 0.5:(0.5 + NBinConf);

% Initialize likelihood and count arrays
Lcond = nan(2, NBinDV, NBinConf);   % P(conf | DV, decision) from model
Ncond = zeros(2, NBinDV, NBinConf); % observed counts in data

for dv = 1:NBinDV
    for decIdx = 1:2
        outcome = 2 - decIdx;  % decIdx=1 → correct (outcome=1), decIdx=2 → error (outcome=0)

        % Model trials in this (DV, decision) cell
        modelCell = ModelDVIndex == dv & ModelOutcome == outcome;

        if sum(modelCell) > MinCondPoints
            % Estimate P(conf | DV, decision) from model histogram
            pModel = histcounts(ModelConfIndex(modelCell), ConfBinIdx, 'Normalization', 'probability');
            pModel(pModel == 0) = 1 / sum(modelCell);  % floor: avoid log(0)
            pModel = pModel / sum(pModel);              % re-normalize after flooring
            Lcond(decIdx, dv, :) = pModel;
        else
            outcomeLabels = {'correct', 'error'};
            warning('ConfidenceReportModel: model underspecified for DV bin %i, %s trials.', ...
                dv, outcomeLabels{decIdx});
        end

        % Count data trials in this cell per confidence bin
        for c = 1:NBinConf
            Ncond(decIdx, dv, c) = sum(DataDVIndex == dv & DataOutcome == outcome & DataConfIndex == c);
        end
    end
end
% Final log-likelihood 
L = nansum(log(Lcond(:)) .* Ncond(:)) / length(DataDVIndex);
end % monte_carlo_model_likelihood


%% =========================================================================
function confOut = confidence_efficacy(conf_0, alpha, Data, Params)
% confidence_efficacy  Apply sub-optimality noise to model confidence.
%
%   Mixes optimal confidence conf_0 with a random sample from the same
%   distribution, controlled by alpha:
%     alpha = 1 → conf unchanged (fully optimal)
%     alpha = 0 → conf fully randomized
%
%   Current supported efficacy type:
%     'randreplace' : each trial is replaced by a random draw with
%                     probability (1 - alpha)

switch Params.Model.efficacy_type
    case 'randreplace'
        % For each trial, with probability (1-alpha), replace its confidence
        % with a randomly sampled confidence from the full distribution
        replaceIdx       = rand(size(conf_0)) > alpha;
        randomDrawIdx    = randi(length(conf_0), sum(replaceIdx), 1);
        confOut          = conf_0;
        confOut(replaceIdx) = conf_0(randomDrawIdx);

    otherwise
        error('ConfidenceReportModel: unknown efficacy type ''%s''.', Params.Model.efficacy_type);
end

end % confidence_efficacy
