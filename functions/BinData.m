function [x, y, e, n, v, ci] = BinData(datax, datay, edges, opt)
% BinData  Bin datay by datax and return group statistics.
%
%   [x, y, e, n, v, ci] = BinData(datax, datay, edges, opt)
%
%   Inputs:
%     datax  - predictor variable (used for binning)
%     datay  - outcome variable to summarize within each bin
%     edges  - bin edges (passed to discretize), or unique values if opt.unique
%     opt    - options struct (all fields optional):
%                .unique    - if true, treat datax as discrete values (default: false)
%                .bin_edge  - if true, clamp out-of-range values to edge bins (default: false)
%
%   Outputs:
%     x   - bin centers (or unique datax values)
%     y   - per-bin mean
%     e   - per-bin SEM
%     n   - per-bin count
%     v   - per-bin SD (e * sqrt(n))
%     ci  - 95% confidence interval half-width (t-based); NaN if not requested

%% --- Default options ----------------------------------------------------
if nargin < 4
    opt = struct;
end
if ~isfield(opt, 'unique'),   opt.unique   = false; end
if ~isfield(opt, 'bin_edge'), opt.bin_edge = false; end

%% --- Assign each datax value to a bin index (idx) ----------------------
if opt.unique
    % Discrete mode: round to 3 decimal places and use unique values as bins
    idx = round(datax * 1000) / 1000;
    x   = unique(idx);
    x   = x(~isnan(x));
else
    % Continuous mode: bin by edges
    idx = discretize(datax, edges);
    x   = mean([edges(1:end-1); edges(2:end)]);  % bin centers

    if opt.bin_edge
        % Clamp values outside the edge range into the first/last bin
        idx(idx == 1)              = 2;
        idx(idx == length(edges)+1) = length(edges);
        idx = idx - 1;
    end
end

%% --- Compute per-bin statistics ----------------------------------------
[y_temp, n_temp, gnames] = grpstats(datay, idx, {'mean', 'numel', 'gname'});
e_temp = grpstats(datay, idx, 'sem');

% grpstats omits empty bins; pad e_temp with NaN if shorter than y_temp
if length(e_temp) < length(y_temp)
    e_temp = nan(size(y_temp));
end

%% --- Expand to full bin count, filling missing bins with NaN -----------
if ~isempty(gnames) && ~opt.unique
    gnames = cellfun(@str2double, gnames);

    % Identify bins with no data
    missing = setxor(1:length(edges)-1,gnames);
    missingdata = nan(length(missing),1);
    
    %Fill in nan data
    namesfull = [gnames;missing(:)];
    y = [y_temp;missingdata];
    y = sortrows([y,namesfull],2); %sort
    y=y(:,1);
    
    e = [e_temp;missingdata];
    e = sortrows([e,namesfull],2); %sort
    e=e(:,1);
    
    n=[n_temp;zeros(length(missing),1)];
    n = sortrows([n,namesfull],2); %sort
    n=n(:,1);
    n=n';

    v = e.*sqrt(n(:));

elseif opt.unique
    y = y_temp;
    e = e_temp;
    n = length(y);
    v = e .* sqrt(n);

else  % no data at all
    y = nan(length(edges)-1, 1);
    e = nan(length(edges)-1, 1);
    n = zeros(length(edges)-1, 1);
    v = nan(length(edges)-1, 1);
end

%% --- Confidence interval (only computed if caller requests it) ---------
if nargout > 5
    fa = tinv(0.975, n);
    ci = e(:) .* fa(:);
else
    ci = nan(size(v));
end

%% --- Ensure column vectors ---------------------------------------------
x = x(:); y = y(:); e = e(:);

end
