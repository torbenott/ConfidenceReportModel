function Data = PsychometricModel(Data, Params)
% PsychometricModel  Fit and optionally plot a psychometric curve.
%
%   Data = PsychometricModel(Data, Params)
%
%   Fits a psychometric function to choice data (Data.DV vs Data.ChoiceA)
%   using FitPsychometric, stores fit parameters in Data.Fit, and
%   optionally plots the binned data and fitted curve.
%
%   Required Data fields:
%     .DV      - decision variable (stimulus evidence)
%     .ChoiceA - binary choice toward option A (0/1)
%
%   Fit results stored in Data.Fit:
%     .Sigma / .Sigma95  - sensitivity and its 95% CI half-width
%     .Bias  / .Bias95   - bias and its 95% CI half-width
%     .Lapse, .Guess     - lapse/guess rates (model-dependent)
%
%   Params fields (all optional):
%     .Function    - fit type passed to FitPsychometric (default: 'erf')
%     .Bin         - number of bins for plotting (default: 7)
%     .UniqueDV    - treat DV as discrete values (default: false)
%     .fitc        - run the fit (default: true)
%     .plotc       - plot the result (default: true)
%     .H           - axes handle to plot into (default: new figure)
%     .Label       - add axis labels (default: true)
%     .Error       - error bar type: 'sem','ci','sd','none' (default: 'sem')
%     .DataMarker  - marker style (default: 'o')
%     .MarkerSize  - marker size (default: 3)
%     .DataColor   - data marker/line color (default: dark grey)
%     .LineColor   - fit line color (default: light grey)
%     .LineWidth   - line width for error bars (default: 1)

%% --- Default parameters -------------------------------------------------
if ~isfield(Params, 'plotc'),       Params.plotc       = true;            end
if ~isfield(Params, 'fitc'),        Params.fitc        = true;            end
if ~isfield(Params, 'Label'),       Params.Label       = true;            end
if ~isfield(Params, 'Function'),    Params.Function    = 'erf';           end
if ~isfield(Params, 'Bin'),         Params.Bin         = 7;               end
if ~isfield(Params, 'UniqueDV'),    Params.UniqueDV    = false;           end
if ~isfield(Params, 'Error'),       Params.Error       = 'sem';           end
if ~isfield(Params, 'DataMarker'),  Params.DataMarker  = 'o';             end
if ~isfield(Params, 'MarkerSize'),  Params.MarkerSize  = 3;               end
if ~isfield(Params, 'DataColor'),   Params.DataColor   = [0.2, 0.2, 0.2]; end
if ~isfield(Params, 'LineColor'),   Params.LineColor   = [0.8, 0.8, 0.8]; end
if ~isfield(Params, 'LineWidth'),   Params.LineWidth   = 1;               end

% Determine whether to use an existing axes or create a new figure
if isfield(Params, 'H') && ~isempty(Params.H)
    axes(Params.H(1));
    create_fig = false;
else
    create_fig = true;
end

xData = Data.DV;
yData = Data.ChoiceA;

%% --- Fit psychometric function ------------------------------------------
if Params.fitc
    fprintf('Fitting psychometric function...');
    tic

    excl  = isnan(xData(:)) | isnan(yData(:));
    fitp  = FitPsychometric(xData(~excl), yData(~excl), Params.Function);

    Data.Fit.Sigma = fitp.sigma;

    % Bias parameter (not present in all model variants)
    try
        Data.Fit.Bias = fitp.m;
    catch
        Data.Fit.Bias = 0;
    end

    % Extract lapse/guess rates depending on fit type
    switch Params.Function
        case 'log'
            Data.Fit.Lapse = fitp.lapse;
            Data.Fit.Guess = fitp.guess;
        case 'erf-lapse'
            Data.Fit.Lapse = fitp.lapse;
            Data.Fit.Guess = fitp.lapse;  % symmetric lapse acts as both
        case 'erf-guess'
            Data.Fit.Lapse = 0;
            Data.Fit.Guess = fitp.guess;
        otherwise
            % Compute 95% CI half-widths for bias and sigma from the cfit object
            if isa(fitp, 'cfit')
                CI = confint(fitp);
                try
                    Data.Fit.Bias95  = (CI(1,2) - CI(1,1)) / 2;
                    sigmaRow = 2;
                catch
                    Data.Fit.Bias95 = NaN;
                    sigmaRow = 1;
                    CI = CI';
                end
                Data.Fit.Sigma95 = (CI(sigmaRow,2) - CI(sigmaRow,1)) / 2;
            end
    end

    fprintf(' Done. (%.2f s)\n', toc);
end

%% --- Plot psychometric curve --------------------------------------------
if Params.plotc
    if create_fig
        figure; hold on;
    end

    % Bin the choice data over the DV range
    edges = linspace(-max(abs(xData)) - 10*eps, max(abs(xData)) + 10*eps, Params.Bin + 1);
    opt = [];
    if Params.UniqueDV
        opt.unique = true;
    end
    [x, y, e, ~, v, ci] = BinData(xData, yData, edges, opt);

    % Suppress error bars when no markers are shown
    if strcmp(Params.DataMarker, 'none')
        e  = nan(size(e));
        ci = nan(size(ci));
        v  = nan(size(v));
    end

    % Select error bar type
    switch Params.Error
        case 'sem',  plote = e;
        case 'ci',   plote = ci;
        case 'sd',   plote = v;
        otherwise,   plote = nan(size(e));
    end

    if Params.Label
        xlabel('Evidence (delta-click / sum-click)');
        ylabel('Choice A');
    end

    % Overlay fitted curve
    if Params.fitc
        xf = linspace(min(x), max(x), 100);
        if isa(fitp, 'cfit')
            yf = feval(fitp, xf);
        else
            yf = feval(fitp.f, xf, Data.Fit.Sigma, Data.Fit.Bias);
        end
        plot(xf, yf, '-', 'Color', Params.LineColor, 'LineWidth', 2);
    end

    % Plot binned data with error bars
    errorbar(x, y, plote, ...
        'LineStyle',       'none', ...
        'Marker',          Params.DataMarker, ...
        'MarkerSize',      Params.MarkerSize, ...
        'MarkerFaceColor', Params.DataColor, ...
        'MarkerEdgeColor', Params.DataColor, ...
        'Color',           Params.DataColor, ...
        'LineWidth',       Params.LineWidth, ...
        'CapSize',         0);

    ylim([-0.05, 1.05]);
    xlim([min(x) - 0.01, max(x) + 0.01]);
end

end
