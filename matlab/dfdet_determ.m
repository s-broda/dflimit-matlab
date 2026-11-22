function D = dfdet_determ(br,C,name,loc,detrend)
%DFDET_DETERM  Piecewise polynomial deterministic regressors on [0,1] (a design of the general engine).
%
%   D = DFDET_DETERM(br,C,name,loc,detrend): br the piece boundaries
%   0 = br(1) < ... < br(end) = 1, C a cell array whose j-th element is the
%   (deg+1) x K matrix of the coefficients of r^k (k = 0..deg, rows) of the K
%   regressors (columns) on piece j.  A regressor is left-continuous at a
%   break: piece j covers br(j) < r <= br(j+1) (r = 0 belongs to piece 1), so
%   that a level-shift dummy 1{r > lambda} is 0 at r = lambda, as 1{t > T_B}
%   in the regression.  loc (false by default): the coefficients are those of
%   the local variable r - br(j) on piece j (the form DFDET_FDESIGN builds, well
%   conditioned at high degree on short pieces).  detrend: 'onestep' (the
%   default), the t-ratio of y_{t-1} in the regression of y_t on d(t/n) and
%   y_{t-1}, or 'first', the statistic that detrends first
%   (DFDET_DETREND_FIRST).  Pieces with fewer rows are padded with zeros.
%
%   The design is a structure with the fields br, C, name, loc, detrend; the
%   constructors DFDET_TREND, DFDET_PERRON, DFDET_FDESIGN, DFDET_LTREND,
%   DFDET_DESIGN_FROM_MATRIX build the common ones.
%
%   Example: the regressors (1, 1{r > 0.4}), a level shift at 0.4
%      D = dfdet_determ([0 0.4 1], {[1 0], [1 1]}, 'level shift at 0.4');
%      dfdet_cdf(-3, D)
%
%   See also DFDET_TREND, DFDET_PERRON, DFDET_FDESIGN, DFDET_CDF.
if nargin < 3 || isempty(name), name = 'design'; end
if nargin < 4 || isempty(loc), loc = false; end
if nargin < 5 || isempty(detrend), detrend = 'onestep'; end
if ~(isnumeric(br) && isreal(br) && isvector(br) && numel(br) >= 2), error('dfdet:args','dfdet_determ: br must be a real vector 0 = br(1) < ... < br(end) = 1.'); end
br = double(br(:).');
if br(1) ~= 0 || br(end) ~= 1 || any(diff(br) <= 0), error('dfdet:args','dfdet_determ: br must satisfy 0 = br(1) < ... < br(end) = 1.'); end
if isnumeric(C), C = {C}; end
if ~iscell(C) || numel(C) ~= numel(br) - 1, error('dfdet:args','dfdet_determ: C must be a cell array with one coefficient matrix per piece.'); end
K = size(C{1},2); R = 0;
for j = 1:numel(C)
    if ~(isnumeric(C{j}) && isreal(C{j}) && size(C{j},2) == K), error('dfdet:args','dfdet_determ: every coefficient matrix must be real with the same number K of columns.'); end
    R = max(R,size(C{j},1));
end
R = max(R,1);
for j = 1:numel(C)
    Cj = double(C{j});
    if size(Cj,1) < R, Cj = [Cj; zeros(R - size(Cj,1),K)]; end
    C{j} = Cj;
end
detrend = lower(char(detrend));
if ~any(strcmp(detrend,{'onestep','first'})), error('dfdet:args','dfdet_determ: detrend must be ''onestep'' or ''first''.'); end
D = struct('br',br,'C',{reshape(C,1,[])},'name',char(name),'loc',logical(loc),'detrend',detrend);
end
