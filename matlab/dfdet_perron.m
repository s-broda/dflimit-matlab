function D = dfdet_perron(model,lam)
%DFDET_PERRON  Perron's trend-break models as designs of the general engine.
%
%   D = DFDET_PERRON(model,lam): model 'A' (1, r, 1{r > lambda}), a break in the
%   level; 'B' (1, r, (r - lambda)+), a break in the slope; 'C' (1, r,
%   1{r > lambda}, (r - lambda)+), both; the breaks at the known fractions
%   lam (a scalar or a vector 0 < lambda_1 < ... < lambda_q < 1, each with its
%   own dummies).  For a break in the slope alone the statistic that detrends
%   first (DFDET_DETREND_FIRST) is the one free of the trend coefficients.
%
%   Example: the 5% critical value with a break in level and slope at 0.3
%      dfdet_icdf(0.05, dfdet_perron('C', 0.3))
%
%   See also DFDET_TREND, DFDET_DETREND_FIRST, DFDET_CDF.
if ~(ischar(model) || isstring(model)) || ~any(strcmpi(char(model),{'A','B','C'})), error('dfdet:args','dfdet_perron: model A, B or C.'); end
model = upper(char(model));
if ~(isnumeric(lam) && isreal(lam) && ~isempty(lam)), error('dfdet:args','dfdet_perron: 0 < lambda_1 < ... < lambda_q < 1.'); end
lam = sort(double(lam(:).'));
if ~(all(lam > 0 & lam < 1) && all(diff(lam) > 0)), error('dfdet:args','dfdet_perron: 0 < lambda_1 < ... < lambda_q < 1.'); end
q = numel(lam);
if model == 'C', K = 2 + 2*q; else, K = 2 + q; end
C = cell(1,q + 1);
for j = 1:q+1                                          % piece j: after the first j - 1 breaks
    Cj = zeros(2,K); Cj(1,1) = 1; Cj(2,2) = 1;
    col = 3;
    for i = 1:q
        on = i <= j - 1;
        if model == 'A' || model == 'C'
            if on, Cj(1,col) = 1; end
            col = col + 1;
        end
        if model == 'B' || model == 'C'
            if on, Cj(1,col) = -lam(i); Cj(2,col) = 1; end
            col = col + 1;
        end
    end
    C{j} = Cj;
end
D = dfdet_determ([0 lam 1],C,sprintf('Perron model %s, lambda = %s',model,strjoin(arrayfun(@(v) sprintf('%g',v),lam,'UniformOutput',false),', ')));
end
