function dfdet_make_paths(srcdir,outfile)
%DFDET_MAKE_PATHS  The path tables of the general engine (dfdet_paths.mat) from those of the Julia package.
%
%   DFDET_MAKE_PATHS(srcdir) reads the outer paths of cases DK, DU and DT at the
%   unit root with a zero start that the Julia package ships (data/paths of
%   DFLimit.jl: DK_c0_zero.bin, DU_c0_zero.bin, DT_c0_zero.bin, written by
%   write_path_table of DFDet.jl: a magic string, the number of nodes K, the
%   number of roots r and of positive roots rp as Int64, the scalars c, tau,
%   t_*, a, b, h, phi, minpos, maxneg, then the nodes t, the derivatives t',
%   the roots (r x K), log P and the weights, little-endian Float64) and saves
%   them as the structures DK, DU, DT in dfdet_paths.mat beside this file (or
%   outfile), from which DFDET_CPATH loads the paths of these cases in place of
%   computing them.  DFDET_TEST checks the tables against the computed paths.
%
%   Example: dfdet_make_paths('/path/to/DFLimit.jl/data/paths')
%
%   See also DFDET_CPATH.
if nargin < 2 || isempty(outfile), outfile = fullfile(fileparts(mfilename('fullpath')),'dfdet_paths.mat'); end
S = struct;
for nm = {'DK','DU','DT'}
    f = fullfile(srcdir,[nm{1} '_c0_zero.bin']);
    fid = fopen(f,'r','ieee-le');
    if fid < 0, error('dfdet:table','dfdet_make_paths: cannot open %s.',f); end
    cl = onCleanup(@() fclose(fid));
    magic = fread(fid,[1 10],'*char');
    if ~strcmp(magic,'DFDetPath1'), error('dfdet:table','dfdet_make_paths: %s is not a path table.',f); end
    n = fread(fid,3,'int64'); K = n(1); r = n(2);
    sc = fread(fid,9,'double');
    cx = @(m) complex(m(1:2:end),m(2:2:end));
    t = cx(fread(fid,2*K,'double')); dt = cx(fread(fid,2*K,'double'));
    R = reshape(cx(fread(fid,2*r*K,'double')),r,K); lP = cx(fread(fid,2*K,'double')); wt = fread(fid,K,'double');
    clear cl
    S.(nm{1}) = struct('c',sc(1),'tau',sc(2),'tstar',sc(3),'a',sc(4),'b',sc(5),'h',sc(6),'phi',sc(7),'minpos',sc(8),'maxneg',sc(9), ...
        't',t,'dt',dt,'R',R,'lP',lP,'wt',wt,'rp',n(3));
end
save(outfile,'-struct','S','-v7');
fprintf('wrote %s\n',outfile);
end
