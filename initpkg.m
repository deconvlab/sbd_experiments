function initpkg(addrm)
%INITPKG  Initialize the SBD-iPALM package
%
%  initpkg()  initializes the SBD-iPALM package by adding the appropriate
%  files to the MATLAB path.
%
%  initpkg(ADDRM)  removes the same files from the MATLAB path if ADDRM == FALSE.
%
%  Directories to be managed can be removed by changing the DIRS variable
%  at the beginning of the funciton.
%

% Directories + subdirectories to add to path
dirs = {
    '.',        false;
    'solvers',  true;
    'utils',    true;
};

warning('OFF', 'MATLAB:addpath:DirNotFound');
warning('OFF', 'MATLAB:rmpath:DirNotFound');

if nargin < 1 || isempty(addrm);  addrm = true;  end

pkgpath = mfilename('fullpath');
d = strsplit(pkgpath, {'\','/'});
pkgpath = pkgpath(1:end-numel(d{end}));

for d = 1:size(dirs,1)
    pathop([pkgpath dirs{d,1}], addrm, dirs{d,2});
end

warning('ON', 'MATLAB:addpath:DirNotFound');
warning('ON', 'MATLAB:rmpath:DirNotFound');

end

function pathop(path, addrm, recursive)
    if recursive
        path = genpath(path);
    end
    
    if addrm
        addpath(path);
    else
        rmpath(path);
    end
end