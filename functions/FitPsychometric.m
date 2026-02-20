function fitresult = FitPsychometric(x, y, fittype_str)
% FitPsychometric  Fit a psychometric function to binned choice data.
%
%   fitresult = FitPsychometric(x, y, fittype_str)
%
%   Inputs:
%     x           - stimulus values (e.g. DV or contrast)
%     y           - proportion correct / choice probability per bin
%     fittype_str - model to fit (string):
%                     'log'       : logistic with guess and lapse
%                     'erf-full'  : cumulative Gaussian (bias + sigma only)
%                     'erf-guess' : cumulative Gaussian with lower asymptote
%                     'erf-lapse' : cumulative Gaussian with symmetric lapse
%
%   Output:
%     fitresult - MATLAB cfit object with fitted parameters

% Fallback dummy data used when too few points to fit reliably
MIN_POINTS = 3;

switch fittype_str

    case 'log'
        % Logistic with guess rate and lapse rate
        % P(x) = guess + (1 - guess - lapse) / (1 + exp(-(x - m) / sigma))
        ft   = fittype('guess+(1-guess-lapse)./(1+exp(-(x-m)/sigma))', ...
                       'independent', 'x', 'dependent', 'y');
        opts = fitoptions(ft);
        opts.Lower      = [0,    0,    -1,   0  ];
        opts.StartPoint = [0.01, 0.01,  0.01, 0.5];
        opts.Upper      = [1,    1,     1,   10 ];
        nParams = 4;

    case 'erf-full'
        % Cumulative Gaussian: bias (m) and sensitivity (sigma)
        % P(x) = normcdf(x, m, sigma)
        ft   = fittype('normcdf(x, m, sigma)', ...
                       'independent', 'x', 'dependent', 'y');
        opts = fitoptions(ft);
        opts.Lower      = [-100,  0.001];
        opts.StartPoint = [   0,  1    ];
        opts.Upper      = [ 100, 100   ];
        nParams = 4;  % dummy data size kept consistent with other branches

    case 'erf-guess'
        % Cumulative Gaussian with lower asymptote (guess rate)
        % P(x) = guess + (1 - guess) * normcdf(x, m, sigma)
        ft   = fittype('guess + (1 - guess)*normcdf(x, m, sigma)', ...
                       'independent', 'x', 'dependent', 'y');
        opts = fitoptions(ft);
        opts.Lower      = [0, -20,  10];
        opts.StartPoint = [0, -10,  10];
        opts.Upper      = [1,  20,  30];
        nParams = 3;

    case 'erf-lapse'
        % Cumulative Gaussian with symmetric lapse rate
        % P(x) = (1 - 2*lapse) * normcdf(x, m, sigma) + lapse
        ft   = fittype('(1-2*lapse)*normcdf(x, m, sigma)+lapse', ...
                       'independent', 'x', 'dependent', 'y');
        opts = fitoptions(ft);
        opts.Lower      = [0, -10,  0.001];
        opts.StartPoint = [0,   0,  1    ];
        opts.Upper      = [1,  10, 10    ];
        nParams = 4;

    otherwise
        error('FitPsychometric: unknown fit type ''%s''.', fittype_str);
end

opts.Display = 'Off';

% Fit model, or substitute dummy data if too few observations
if length(x) > MIN_POINTS
    [fitresult, ~] = fit(x(:), double(y(:)), ft, opts);
else
    dummy = 0.5 * ones(nParams, 1);
    [fitresult, ~] = fit(dummy, dummy, ft, opts);
end

end
