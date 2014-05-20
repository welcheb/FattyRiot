%% setup FattyRiot
%
% add parent folder of FattyRiot.m to MATLAB path
%
% setup_FattyRiot.m should remain in the same folder as FattyRiot.m
%

%% parse location of setup_FattyRiot.m
[pathstr,name,ext] = fileparts(mfilename('fullpath'));

%% Store present working directory and change directory temporarily to matlabroot
pwd_path = pwd;
cd(matlabroot);

if exist('FattyRiot')~=2,
    addpath(pathstr);
    savepath;
    if exist('FattyRiot')~=2,
        disp( sprintf('FattyRiot.m appears to NOT be setup. Check that ''%s'' is added to MATLAB path.', pathstr) );
    else
        disp( sprintf('FattyRiot.m is now accessible at ''%s''', which('FattyRiot') ) );
    end
else
    disp( sprintf('FattyRiot.m is already accessible at ''%s''', which('FattyRiot') ) );
end

%% change back to the original directory
cd(pwd_path);
