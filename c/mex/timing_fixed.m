% TIMING_FIXED  Seconds per value of the distribution function with a fixed start: the MATLAB programs (the conditional
% integral, and the contour integral in cases DU, DT; 1 thread) against the C MEX (the line integral with the
% noncentral kernel, the conditional integral and the contour integral), 1 and 10 threads; timeit medians (the
% fixed-start rows of the table in README.md and the cost paragraph of the Supplementary Material), and the values.
%   Needs the MEX file in this folder and MATLAB; run from anywhere.
here = fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(fullfile(here,'matlab')); mexdir = fullfile(here,'c','mex');
cd(here);   % not this folder: its engine replacements would shadow the MATLAB programs from the current directory
pts = {'DK',-5,0,1,-2; 'DK',2,0,1,1; 'DU',-5,0,1,-2.86; 'DU',-5,0,1,0.5; 'DU',2,0,1,-3; 'DT',-5,0,2,-3; 'DT',2,0,1,-2; 'DU',-5,-5,1,-2; 'DT',-2.5,-5,1,-1.1};
fprintf('%-18s %10s %10s %10s %10s %10s %10s %10s %10s   %18s %18s %18s %18s %18s\n', 'case c c0 g0 x', 'MATLAB dbl', 'MATLAB ctr', ...
        'MEX1 sgl', 'MEX1 dbl', 'MEX1 ctr', 'MEX10 sgl', 'MEX10 dbl', 'MEX10 ctr', 'F single', 'F double', 'F contour', 'MATLAB double', 'MATLAB contour');
for k = 1:size(pts,1)
    cas = pts{k,1}; c = pts{k,2}; c0 = pts{k,3}; g0 = pts{k,4}; x = pts{k,5};
    o = struct('c0',c0,'g0',g0); ctr = any(strcmp(cas,{'DU','DT'}));
    % the MATLAB programs: matlab/ ahead of c/mex on the path, one thread
    rmpath(mexdir); addpath(fullfile(here,'matlab')); maxNumCompThreads(1);
    tm = timeit(@() ltratio_cdf(x, cas, c, setfield(o,'method','double')));
    Fm = ltratio_cdf(x, cas, c, setfield(o,'method','double'));
    tmc = NaN; Fmc = NaN;
    if ctr
        tmc = timeit(@() ltratio_cdf(x, cas, c, setfield(o,'method','contour')));
        Fmc = ltratio_cdf(x, cas, c, setfield(o,'method','contour'));
    end
    % the MEX: c/mex ahead
    addpath(mexdir);
    ltratio('threads', 1);
    ts1 = timeit(@() ltratio('cdf', x, cas, c, setfield(o,'method','single')));
    td1 = timeit(@() ltratio('cdf', x, cas, c, setfield(o,'method','double')));
    tc1 = NaN; if ctr, tc1 = timeit(@() ltratio('cdf', x, cas, c, setfield(o,'method','contour'))); end
    ltratio('threads', 10);
    ts10 = timeit(@() ltratio('cdf', x, cas, c, setfield(o,'method','single')));
    td10 = timeit(@() ltratio('cdf', x, cas, c, setfield(o,'method','double')));
    tc10 = NaN; if ctr, tc10 = timeit(@() ltratio('cdf', x, cas, c, setfield(o,'method','contour'))); end
    Fs = ltratio('cdf', x, cas, c, setfield(o,'method','single')); Fd = ltratio('cdf', x, cas, c, setfield(o,'method','double'));
    Fc = NaN; if ctr, Fc = ltratio('cdf', x, cas, c, setfield(o,'method','contour')); end
    fprintf('%-18s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f   %18.15f %18.15f %18.15f %18.15f %18.15f\n', ...
            sprintf('%s %g %g %g %g', cas, c, c0, g0, x), tm, tmc, ts1, td1, tc1, ts10, td10, tc10, Fs, Fd, Fc, Fm, Fmc);
end
