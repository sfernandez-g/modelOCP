%GUESSMSMR Estimate initial MSMR parameters from OCV data.
%
%   [U0, X, w] = guessMSMR(ocvrel, TdegC)
%   [U0, X, w, diag] = guessMSMR(ocvrel, TdegC, ...)
%
%   [U0, X, w] = guessMSMR(ocvrel, J, TdegC)
%   [U0, X, w, diag] = guessMSMR(ocvrel, J, TdegC, ...)
%
%   Estimates an initial set of MSMR parameters from the average of the
%   charge and discharge OCV curves contained in ocvrel.
%
%   If J is provided, the smoothing resolution is increased until the
%   number of detected dZ/dU peaks is less than or equal to J.
%
%   If J is not provided, the number of detected peaks at the selected
%   smoothing resolution is used as the number of MSMR galleries.
%
%   Inputs
%   ------
%   ocvrel : struct
%       OCV structure containing charge and discharge estimates.
%   J : integer, optional
%       Target number of MSMR galleries.
%
%   TdegC : scalar
%       Temperature [degC].
%
%   Optional name-value arguments
%   -----------------------------
%   'dvInitial' : scalar
%       Initial voltage resolution [V]. Default: 1e-3 V.
%   'dvIncrement' : scalar
%       Increment in voltage resolution [V]. Default: 0.5e-3 V.
%   'dvFinal' : scalar
%       Maximum voltage resolution [V]. Default: 0.1 V.
%   'vmin' : scalar
%       Minimum voltage used for the initialization. Default: [].
%   'vmax' : scalar
%       Maximum voltage used for the initialization. Default: [].
%
%   Outputs
%   -------
%   U0 : vector
%       Initial MSMR potential parameters [V].
%   X : vector
%       Initial gallery fractions. sum(X) = 1.
%   w : vector
%       Initial MSMR width parameters [V].
%   diag : struct
%       Diagnostic information used to construct the initialization.

function [U0, X, w, diag] = guessMSMR(ocvrel, varargin)
    if isempty(varargin)
        error('Temperature TdegC must be provided.');
    end
    if isscalar(varargin{1}) && isnumeric(varargin{1}) && ...
            varargin{1} == round(varargin{1}) && varargin{1} >= 1 && ...
            numel(varargin) >= 2 % J explicitly provided.
        J = varargin{1};
        TdegC = varargin{2};
        nameValueStart = 3;
    else % J omitted.
        J = [];
        TdegC = varargin{1};
        nameValueStart = 2;
    end
    if ~isscalar(TdegC) || ~isnumeric(TdegC)
        error('TdegC must be a numeric scalar.');
    end

    %% Parse optional name-value arguments
    parser = inputParser;
    addParameter(parser, 'dvInitial', 1e-3, ...
        @(x) isscalar(x) && x > 0);
    addParameter(parser, 'dvIncrement', 0.5e-3, ...
        @(x) isscalar(x) && x > 0);
    addParameter(parser, 'dvFinal', 0.1, ...
        @(x) isscalar(x) && x > 0);
    addParameter(parser, 'vmin', [], ...
        @(x) isempty(x) || (isscalar(x) && isnumeric(x)));
    addParameter(parser, 'vmax', [], ...
        @(x) isempty(x) || (isscalar(x) && isnumeric(x)));
    parse(parser, varargin{nameValueStart:end});
    arg = parser.Results;

    %% Constants
    R = 8.3144598;      % Molar gas constant [J/mol/K]
    F = 96485.3329;     % Faraday constant [C/mol]

    %% Extract charge and discharge OCV curves
    names = string({ocvrel.estimates.name});
    idxChg = find(names == "charge", 1);
    idxDch = find(names == "discharge", 1);
    if isempty(idxChg)
        error('ocvrel does not contain a "charge" OCV curve.');
    end
    if isempty(idxDch)
        error('ocvrel does not contain a "discharge" OCV curve.');
    end
    Zchg = ocvrel.estimates(idxChg).Z(:);
    Uch = ocvrel.estimates(idxChg).U(:);
    Zdch = ocvrel.estimates(idxDch).Z(:);
    Udch = ocvrel.estimates(idxDch).U(:);

    %% Common composition vector
    Zmin = max([min(Zchg), min(Zdch)]);
    Zmax = min([max(Zchg), max(Zdch)]);
    if Zmin >= Zmax
        error(['Charge and discharge OCV curves do not have an ' ...
               'overlapping composition range.']);
    end
    Z = linspace(Zmin, Zmax, 1000)';
    Uch = linearinterp(Zchg, Uch, Z);
    Udch = linearinterp(Zdch, Udch, Z);

    % Average charge/discharge OCV.
    U = (Uch + Udch) / 2;

    %% Restrict voltage range if requested
    valid = true(size(U));
    if ~isempty(arg.vmin)
        valid = valid & U > arg.vmin;
    end
    if ~isempty(arg.vmax)
        valid = valid & U < arg.vmax;
    end
    Z = Z(valid);
    U = U(valid);
    if numel(Z) < 3
        error('Not enough OCV data points remain after applying the voltage limits.');
    end

    %% Estimate dZ/dU
    found = false;
    for dv = arg.dvInitial:arg.dvIncrement:arg.dvFinal
        [Z0, U0curve, dZdU] = smoothdiff(Z, U, dv);
        % Peaks in dZ/dU versus voltage.
        [~, idxUpeaks] = findpeaks(dZdU);
        % If J is specified, increase smoothing until the number of
        % detected peaks is compatible with the requested number.
        if ~isempty(J)
            if numel(idxUpeaks) <= J
                found = true;
                break;
            end
        else
            % Without a target J, use the first detected set of peaks.
            found = true;
            break;
        end
    end
    if ~found
        warning(['The requested number of galleries was not reached. ' ...
                 'Using the final smoothing resolution.']);
    end

    %% Determine J automatically if it was not provided
    nPeaks = numel(idxUpeaks);
    if isempty(J)
        J = nPeaks;
        if J < 1
            error('No dZ/dU peaks were detected.');
        end
    elseif nPeaks < J
        warning(['Only %d dZ/dU peaks were detected, but J = %d was ' ...
                 'requested. The MSMR initialization is incomplete.'], ...
                 nPeaks, J);
        U0 = [];
        X = [];
        w = [];
    elseif nPeaks > J
        error(['The number of detected peaks (%d) is larger than J (%d). ' ...
               'This should not occur after the smoothing loop.'], ...
               nPeaks, J);
    end

    %% Determine U0 from dZ/dU peaks
    U0 = U0curve(idxUpeaks);
    dZdUpeaks = dZdU(idxUpeaks);
    % -------------------------------- Sort U0 from low to high voltage.
    [U0, order] = sort(U0, 'ascend');
    idxUpeaks = idxUpeaks(order);
    dZdUpeaks = dZdUpeaks(order);

    %% Determine X from valleys between consecutive peaks
    idxBoundaries = zeros(J-1, 1);
    for j = 1:J-1
        idx1 = idxUpeaks(j);
        idx2 = idxUpeaks(j+1);
        if idx1 > idx2
            tmp = idx1;
            idx1 = idx2;
            idx2 = tmp;
        end
        [~, localIdx] = min(dZdU(idx1:idx2)); % Search for the minimum dZ/dU between consecutive peaks.
        idxBoundaries(j) = idx1 + localIdx - 1;
    end
    Zvalleys = Z0(idxBoundaries);
    % --------------------------------- Sort boundaries by composition.
    Zvalleys = sort(Zvalleys);
    % ---------------------------------- Gallery fractions.
    X = diff([min(Z0); Zvalleys; max(Z0)]);
    % ---------------------------------- Normalize to ensure sum(X)=1.
    X = X / sum(X);

    %% Estimate w
    f = F / R / (TdegC + 273.15);
    w = f .* X ./ dZdUpeaks / 4;
    w = abs(w);

    %% Diagnostic output
    diag = struct();
    diag.Z = Z;
    diag.U = U;
    diag.dZdU = dZdU;
    diag.Ucurve = U0curve;
    diag.Upeaks = U0;
    diag.dZdUpeaks = dZdUpeaks;
    diag.Zvalleys = Zvalleys;
    diag.dv = dv;
    diag.found = found;
    diag.nPeaks = nPeaks;
    diag.J = J;

end