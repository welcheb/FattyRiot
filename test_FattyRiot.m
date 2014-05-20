%% Test FattyRiot
%clear all; close all; clc;

%% start clock
test_FattyRiot_start_time = tic;

%% detect run location details
[mfile_pathstr, mfile_name, mfile_ext] = fileparts(mfilename('fullpath'));

%% check for test_results folder
test_results_folder = sprintf('%s/test_results', mfile_pathstr);
if exist(test_results_folder)~=7,
    mkdir(test_results_folder);
end

%% start diary
diary_file_name = sprintf('%s/test_FattyRiot_diary.txt', test_results_folder);
warning off; delete(diary_file_name); warning on;
diary(diary_file_name);
diary on;

%% folder locations
case_folder = sprintf('%s/test_cases', mfile_pathstr);

%% loop through cases
%for c=4, % Quick test
%for c=1:10, % Phase I cases
%for c=11:17, % Phase II cases
for c=1:17, % All cases
    
    clear imDataParams;
    case_matfilename = sprintf('%s/%02d.mat', case_folder, c);
    load(case_matfilename);
    [nx ny nz ncoils nTE] = size(imDataParams.images);
        
    [FW,INFO] = FattyRiot(imDataParams);
    
    FW = FW - min(FW(:));
    FW = FW / max(FW(:));
    F = FW(:,:,[1:nz]);
    W = FW(:,:,[1:nz]+nz);

    figure(1000+c);
    imagesc([F(:,:) ; W(:,:)]);
    axis image;
    colormap(gray);
    title(sprintf('CASE %02d',c));
    drawnow;
    
    FWpngfilename = sprintf('%s/%02d_FattyRiot.png', test_results_folder, c);
    imwrite([F(:,:) ; W(:,:)], FWpngfilename, 'PNG', 'BitDepth', 8, 'Author', 'FattyRiot', 'Description', sprintf('2012 ISMRM Fat Water Challenge : Phase I : Case %02d : %s', c, datestr(now)) );
    
    disp( sprintf('COMPLETED CASE %02d',c) );

end

%% completion time
disp(sprintf('Completed test_FattyRiot in %0.2f seconds', toc(test_FattyRiot_start_time) ));

%% stop diary
diary off;
