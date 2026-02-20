function [Data, Params] = PlotConfidenceSignatures(Data, Params)
% PlotConfidenceSignatures  Plot the three canonical confidence signatures.
%
%   [Data, Params] = PlotConfidenceSignatures(Data, Params)
%
%   Produces three plots characterizing the relationship between behavioral
%   confidence (waiting time / time investment) and task performance,
%   as defined in Ott et al. (2018) and Masset, Ott et al. (2020):
%
%     1. Calibration curve  : P(correct) as a function of confidence,
%                             computed on reward-omission trials.
%     2. Vevaiometric curve : Mean confidence as a function of evidence (DV),
%                             split by correct and error trials.
%     3. Conditioned psychometric : P(correct) vs DV, split by high/low
%                             confidence (median split on omission trials).
%
%   Required Data fields:
%     .DV      - decision variable (evidence); positive = toward choice A
%     .ChoiceA - binary choice (1 = A, 0 = B)
%     .Conf    - confidence report (e.g. waiting time in seconds)
%     .Outcome - binary accuracy (1 = correct, 0 = error)
%
%   Optional Data fields:
%     .Omission - logical flag for reward-omission trials (default: all true)
%
%   Params fields (all optional):
%     .Figure       - figure handle to plot into (default: new figure)
%     .Axes         - 1x3 array of axes handles (default: new subplots)
%     .NBinDV       - number of DV bins (default: 8)
%     .NBinConf     - number of confidence bins (default: 6)
%     .EvidenceBin  - explicit DV bin edges (overrides NBinDV)
%     .ConfidenceBin- explicit confidence bin edges (overrides NBinConf)
%     .AbsDV        - use absolute DV, split by outcome (default: false)
%     .UniqueDV     - treat DV as discrete values (default: false)
%     .ErrorBar     - error bar type: 'sem', 'ci', or 'none' (default: 'none')
%     .MinPoints    - min trials per bin; fewer → NaN (default: 0)
%     .ConfMed      - fixed median for high/low confidence split (default: [])
%     .Label        - add axis labels and titles (default: true)
%     .LineStyle    - line style for curves (default: '-')
%     .LineWidth    - line width; 0 = no line (default: 2)
%     .MarkerSize   - marker size (default: 5)
%     .ColorsC      - calibration curve color (default: grey)
%     .ColorsV      - vevaiometric colors {correct, error} (default: green/red)
%     .ColorsP      - conditioned psychometric colors {high, low} (default: dark/light)

%% --- Default parameters -------------------------------------------------
if ~isfield(Params, 'Figure'),    Params.Figure    = false;          end
if ~isfield(Params, 'Axes'),      Params.Axes      = nan(1, 3);      end
if ~isfield(Params, 'NBinDV'),    Params.NBinDV    = 8;              end
if ~isfield(Params, 'NBinConf'),  Params.NBinConf  = 6;              end
if ~isfield(Params, 'AbsDV'),     Params.AbsDV     = false;          end
if ~isfield(Params, 'ErrorBar'),  Params.ErrorBar  = 'none';         end
if ~isfield(Params, 'MarkerSize'),Params.MarkerSize= 5;              end
if ~isfield(Params, 'LineWidth'), Params.LineWidth  = 2;             end
if ~isfield(Params, 'MinPoints'), Params.MinPoints  = 0;             end
if ~isfield(Params, 'ConfMed'),   Params.ConfMed   = [];             end
if ~isfield(Params, 'UniqueDV'),  Params.UniqueDV  = false;          end
if ~isfield(Params, 'Label'),     Params.Label     = true;           end
if ~isfield(Params, 'LineStyle'), Params.LineStyle  = '-';           end
if ~isfield(Params, 'ColorsZ'),   Params.ColorsZ   = {[230,140,0]/255}; end

%% --- Figure and style setup ---------------------------------------------
if ~isa(Params.Figure, 'matlab.ui.Figure')
    Params.Figure = figure('Units', 'normalized', 'Position', [0.26, 0.37, 0.45, 0.41]);
else
    figure(Params.Figure);
end

% Colors and markers depend on whether error bars (data mode) or lines (model mode)
noErrorBar = strcmp(Params.ErrorBar, 'none') || ...
             (islogical(Params.ErrorBar) && ~Params.ErrorBar);
if noErrorBar
    % Line/model style: lighter colors, no markers
    if ~isfield(Params, 'ColorsC'), Params.ColorsC = {[.7, .7, .7]};                      end
    if ~isfield(Params, 'ColorsV'), Params.ColorsV = {[.5,.9,.5], [.9,.5,.5], [.8,.7,.1]}; end
    if ~isfield(Params, 'ColorsP'), Params.ColorsP = {[.4,.4,.4], [.85,.85,.85]};          end
    marker = 'none';
else
    % Data style: darker colors, circle markers
    if ~isfield(Params, 'ColorsC'), Params.ColorsC = {[.2, .2, .2]};                      end
    if ~isfield(Params, 'ColorsV'), Params.ColorsV = {[.1,.9,.1], [.9,.1,.1], [.8,.7,.1]}; end
    if ~isfield(Params, 'ColorsP'), Params.ColorsP = {[.2,.2,.2], [.7,.7,.7], [.8,.8,.1]}; end
    marker = 'o';
end

% Line style: LineWidth=0 suppresses lines entirely
if Params.LineWidth == 0
    linestyle = 'none';
    linewidth  = 1;
else
    linestyle = Params.LineStyle;
    linewidth  = Params.LineWidth;
end

% If Omission field is missing or mismatched in size, treat all trials as omissions
if ~isfield(Data, 'Omission') || any(size(Data.DV) ~= size(Data.Omission))
    Data.Omission = true(numel(Data.DV), 1);
end

%% --- DV preprocessing ---------------------------------------------------
% AbsDV=true: fold DV to positive axis and split by correct/error (outcome)
% AbsDV=false: keep signed DV and split by choice direction
if Params.AbsDV
    DV     = abs(Data.DV);
    Choice = Data.Outcome;
    xlim_dv = [0, max(DV) + 0.05];
else
    DV     = Data.DV;
    Choice = Data.ChoiceA;
    xlim_dv = [-max(abs(DV)) - 0.05, max(abs(DV)) + 0.05];
end

% BinData option for discrete DV values
opt.unique = Params.UniqueDV;

%% --- Evidence bin edges -------------------------------------------------
if isfield(Params, 'EvidenceBin') && ~isempty(Params.EvidenceBin)
    EvidenceBin = Params.EvidenceBin;
else
    uniqueDV = unique(DV);
    if length(uniqueDV) < 10
        % Discrete DV: place bin edges midway between unique values
        dd  = diff(uniqueDV);
        idx = find(dd > 0.001);
        EvidenceBin      = zeros(1, length(idx) + 2);
        EvidenceBin(1)   = uniqueDV(1) - 10*eps;
        EvidenceBin(end) = uniqueDV(end) + 10*eps;
        for k = 2:length(idx) + 1
            EvidenceBin(k) = uniqueDV(idx(k-1)) + dd(idx(k-1)) / 2;
        end
    else
        % Continuous DV: uniform bins
        EvidenceBin = linspace(min(DV) - 10*eps, max(DV) + 10*eps, Params.NBinDV);
    end
end
Params.EvidenceBin = EvidenceBin;

%% --- Confidence bin edges -----------------------------------------------
if isfield(Params, 'ConfidenceBin') && ~isempty(Params.ConfidenceBin)
    ConfidenceBin = Params.ConfidenceBin;
else
    ConfidenceBin = linspace(min(Data.Conf) - 10*eps, max(Data.Conf) + 10*eps, Params.NBinConf);
    Params.ConfidenceBin = ConfidenceBin;
end

%% --- Median confidence split for conditioned psychometric ---------------
% Split omission trials into high/low confidence halves using median.
% If Params.ConfMed is set, use that as a fixed threshold instead.
if isempty(Params.ConfMed)
    ConfLow  = nanmedian(Data.Conf(Data.Omission));
    ConfHigh = quantile(Data.Conf(Data.Omission), 0.5);  % same as median
else
    ConfLow  = Params.ConfMed;
    ConfHigh = Params.ConfMed;
end

%% -----------------------------------------------------------------------
%  PLOT 1: Vevaiometric (mean confidence vs evidence, correct vs error)
%  Ref: Ott et al. 2018, Masset & Ott et al. 2020
%% -----------------------------------------------------------------------
if ~isa(Params.Axes(2), 'matlab.graphics.axis.Axes')
    Params.Axes(2) = subplot(1, 3, 2);
else
    axes(Params.Axes(2));
end
hold on

% Bin mean confidence separately for correct and error trials
[xc, yc, ec, nc, ~, cic] = BinData(DV(Data.Outcome==1), Data.Conf(Data.Outcome==1), EvidenceBin, opt);
[xe, ye, ee, ne, ~, cie] = BinData(DV(Data.Outcome==0), Data.Conf(Data.Outcome==0), EvidenceBin, opt);

% Mask bins with too few trials
xc(nc<Params.MinPoints) = NaN;  yc(nc<Params.MinPoints) = NaN;
ec(nc<Params.MinPoints) = NaN;  cic(nc<Params.MinPoints) = NaN;
xe(ne<Params.MinPoints) = NaN;  ye(ne<Params.MinPoints) = NaN;
ee(ne<Params.MinPoints) = NaN;  cie(ne<Params.MinPoints) = NaN;

% Select error bar values
if strcmp(Params.ErrorBar, 'ci')
    ec = cic;  ee = cie;
elseif noErrorBar
    ec = nan(1, numel(xc));  ee = nan(1, numel(xe));
end

vis_c = ~isnan(yc);  vis_e = ~isnan(ye);
Data.Axes.Vevaio.Correct = errorbar(xc(vis_c), yc(vis_c), ec(vis_c), ...
    'LineStyle', linestyle, 'Marker', marker, 'MarkerSize', Params.MarkerSize, ...
    'MarkerFaceColor', Params.ColorsV{1}, 'MarkerEdgeColor', Params.ColorsV{1}, ...
    'Color', Params.ColorsV{1}, 'LineWidth', linewidth, 'CapSize', 0);
Data.Axes.Vevaio.Error = errorbar(xe(vis_e), ye(vis_e), ee(vis_e), ...
    'LineStyle', linestyle, 'Marker', marker, 'MarkerSize', Params.MarkerSize, ...
    'MarkerFaceColor', Params.ColorsV{2}, 'MarkerEdgeColor', Params.ColorsV{2}, ...
    'Color', Params.ColorsV{2}, 'LineWidth', linewidth, 'CapSize', 0);

if Params.Label
    xlabel('Evidence (delta-click / sum-click)');
    ylabel('Time investment (s)');
    title('Vevaiometric');
    l=legend({'Correct','Error'});
    l.Location='northwest';
end
xlim(xlim_dv);

%% -----------------------------------------------------------------------
%  PLOT 2: Calibration curve (P(correct) vs confidence, omission trials)
%% -----------------------------------------------------------------------
if ~isa(Params.Axes(1), 'matlab.graphics.axis.Axes')
    Params.Axes(1) = subplot(1, 3, 1);
else
    axes(Params.Axes(1));
end
hold on

[xCal, yCal, eCal, nCal, ~, ciCal] = BinData( ...
    Data.Conf(Data.Omission), Data.Outcome(Data.Omission), ConfidenceBin);
xCal(nCal < Params.MinPoints) = NaN;

if strcmp(Params.ErrorBar, 'ci')
    eCal = ciCal;
elseif noErrorBar
    eCal = nan(1, numel(xCal));
end

vis_cal = ~isnan(yCal);
Data.Axes.Calibration = errorbar(xCal(vis_cal), yCal(vis_cal), eCal(vis_cal), ...
    'LineStyle', linestyle, 'Marker', marker, 'MarkerSize', Params.MarkerSize, ...
    'MarkerFaceColor', Params.ColorsC{1}, 'MarkerEdgeColor', Params.ColorsC{1}, ...
    'Color', Params.ColorsC{1}, 'LineWidth', linewidth, 'CapSize', 0);

if Params.Label
    xlabel('Time investment (s)');
    ylabel('P(correct)');
    title('Calibration curve');
end
ylim([0.4, 1]);
set(gca, 'YTick', 0.5:0.25:1, 'YTickLabel', 0.5:0.25:1);

%% -----------------------------------------------------------------------
%  PLOT 3: Conditioned psychometric (P(correct) vs evidence, high vs low conf)
%% -----------------------------------------------------------------------
if ~isa(Params.Axes(3), 'matlab.graphics.axis.Axes')
    Params.Axes(3) = subplot(1, 3, 3);
else
    axes(Params.Axes(3));
end
hold on

% High confidence trials: above median; Low: at or below median
highConf = Data.Omission & Data.Conf >  ConfHigh;
lowConf  = Data.Omission & Data.Conf <= ConfLow;

[xh, yh, eh, nh, ~, cih] = BinData(DV(highConf), Choice(highConf), EvidenceBin, opt);
[xl, yl, el, nl, ~, cil] = BinData(DV(lowConf),  Choice(lowConf),  EvidenceBin, opt);

xh(nh<Params.MinPoints) = NaN;  yh(nh<Params.MinPoints) = NaN;  eh(nh<Params.MinPoints) = NaN;
xl(nl<Params.MinPoints) = NaN;  yl(nl<Params.MinPoints) = NaN;  el(nl<Params.MinPoints) = NaN;

if strcmp(Params.ErrorBar, 'ci')
    eh = cih;  el = cil;
elseif noErrorBar
    eh = nan(1, numel(xh));  el = nan(1, numel(xl));
end

vis_h = ~isnan(yh);  vis_l = ~isnan(yl);
Data.Axes.Cond.High = errorbar(xh(vis_h), yh(vis_h), eh(vis_h), ...
    'LineStyle', linestyle, 'Marker', marker, 'MarkerSize', Params.MarkerSize, ...
    'MarkerFaceColor', Params.ColorsP{1}, 'MarkerEdgeColor', Params.ColorsP{1}, ...
    'Color', Params.ColorsP{1}, 'LineWidth', linewidth, 'CapSize', 0);
Data.Axes.Cond.Low = errorbar(xl(vis_l), yl(vis_l), el(vis_l), ...
    'LineStyle', linestyle, 'Marker', marker, 'MarkerSize', Params.MarkerSize, ...
    'MarkerFaceColor', Params.ColorsP{2}, 'MarkerEdgeColor', Params.ColorsP{2}, ...
    'Color', Params.ColorsP{2}, 'LineWidth', linewidth, 'CapSize', 0);

if Params.Label
    xlabel('Evidence (delta-click / sum-click)');
    ylabel('P(correct)');
    title('Conditioned psychometric');
    l=legend({['TI > ', num2str(round(ConfHigh*100)/100), ' s'], ...
            ['TI \leq ', num2str(round(ConfLow*100)/100), ' s']}, ...
           'Location', 'SouthEast');
    l.Box='on';
end
xlim(xlim_dv);
ylim([0.5, 1]);
set(gca, 'YTick', 0.5:0.25:1, 'YTickLabel', 0.5:0.25:1);

RedoTicks(gcf);

end
