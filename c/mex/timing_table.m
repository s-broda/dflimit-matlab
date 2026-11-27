% TIMING_TABLE  The first table of c/mex/README.md: seconds per value of the distribution function by one evaluation, the MATLAB
% programs (one thread; single = the line integral, double = the conditional integral) against the C MEX (1 and 10
% threads); the smallest of 7 warm calls (3 from 1 s on).  Then the second table (c/mex/timing_fixed.m, a fixed start).
% Needs the MEX file in this folder and MATLAB; run from anywhere: run('c/mex/timing_table.m').
here = fileparts(fileparts(fileparts(mfilename('fullpath')))); cd(here); mexdir = fullfile(here,'c','mex'); addpath(fullfile(here,'matlab'));
pts = {'DK',0,0,-1.94; 'DU',0,0,-2.86; 'DT',0,0,-3.41; 'DU',-5,0,-2.86; 'DU',10,0,20; 'DU',0,1,-1; 'DT',-2.5,-5,-1.1; ...
       'DK',2,2,0; 'DU',5,-5,126.85; 'SU',-5,0,-2; 'ST',-10,0,-3};
fprintf('MATLAB %s, %s, maxNumCompThreads before: %d\n', version, datestr(now), maxNumCompThreads);
fprintf('%-18s %10s %10s %10s %10s %10s %10s   %18s %18s %18s %18s\n', 'case c c0 x', 'MAT sgl', 'MAT dbl', 'MEX1 sgl', 'MEX1 dbl', ...
        'MEX10 sgl', 'MEX10 dbl', 'F MAT sgl', 'F MAT dbl', 'F MEX sgl', 'F MEX dbl');
for k = 1:size(pts,1)
    cas = pts{k,1}; c = pts{k,2}; c0 = pts{k,3}; x = pts{k,4}; o = struct('c0',c0);
    rmpath(mexdir); addpath(fullfile(here,'matlab')); nt = maxNumCompThreads(1);
    ts = tmin(@() ltratio_cdf(x, cas, c, setfield(o,'method','single')));
    td = tmin(@() ltratio_cdf(x, cas, c, setfield(o,'method','double')));
    Fs = ltratio_cdf(x, cas, c, setfield(o,'method','single')); Fd = ltratio_cdf(x, cas, c, setfield(o,'method','double'));
    maxNumCompThreads(nt);
    addpath(mexdir);
    ltratio('threads', 1);
    s1 = tmin(@() ltratio('cdf', x, cas, c, setfield(o,'method','single')));
    d1 = tmin(@() ltratio('cdf', x, cas, c, setfield(o,'method','double')));
    ltratio('threads', 10);
    s10 = tmin(@() ltratio('cdf', x, cas, c, setfield(o,'method','single')));
    d10 = tmin(@() ltratio('cdf', x, cas, c, setfield(o,'method','double')));
    Gs = ltratio('cdf', x, cas, c, setfield(o,'method','single')); Gd = ltratio('cdf', x, cas, c, setfield(o,'method','double'));
    fprintf('%-18s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f   %18.15f %18.15f %18.15f %18.15f\n', ...
            sprintf('%s %g %g %g', cas, c, c0, x), ts, td, s1, d1, s10, d10, Fs, Fd, Gs, Gd);
end
rmpath(mexdir);
fprintf('\n--- the second table (timing_fixed.m)\n');
run(fullfile(mexdir, 'timing_fixed.m'));

function t = tmin(f)
% the smallest of 7 timed calls after a warm-up call (of 3 where a call takes 1 s or more)
f(); t0 = tic; f(); t = toc(t0); n = 7; if t >= 1, n = 3; end
for k = 2:n, t0 = tic; f(); t = min(t, toc(t0)); end
end
