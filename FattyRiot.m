%% FATTYRIOT Top-level function for 2012 ISMRM Fat-Water Challenge
%
%  FW = FattyRiot(imDataParams);
%
% Input: structure imDataParams
%   - imDataParams.images                : acquired images, array of size [nx, ny, nz, ncoils, nTE]
%   - imDataParams.TE                    : echo times (in seconds), vector of length nTE
%   - imDataParams.FieldStrength         : (in Tesla), (default 1.5)
%   - imDataParams.PrecessionIsClockwise : ==1 is clockwise, ~=1 is counterclockwise, (default 1) 
%   - imDataParams.mask                  : logical [nx, ny, nz], (default true([nx, ny, nz]) )
%
% Output: FW (fat and water magnitude images) size [nx, ny, 2*nz]
%
% Fat magnitude images are stored in the first nz slices of FW
% Water magnitude images are stored in the second nz slices of FW
%
%  [FW,INFO] = FattyRiot(imDataParams);
%
% Returns an INFO structure containing the following fields
%   - INFO.imDataParams           : originally passed imDataParams structure
%   - INFO.imDataParams_FattyRiot : FattyRiot modified imDataParams structure
%   - INFO.F                      : complex F images [nx, ny, nz]
%   - INFO.W                      : complex W images [nx, ny, nz]
%   - INFO.FM_HZ                  : fieldmap in Hz [nx, ny, nz]
%   - INFO.R2                     : R2* map in s^-1 [nx, ny, nz]
%   - INFO.ERROR                  : residual error map [nx, ny, nz]
%   - INFO.nTEeven                : number of evenly spaced echo times
%   - INFO.idx_TEeven             : indices of evenly spaced echo times
%   - INFO.dTEeven                : delta TE for evenly space echo times (in seconds)
%   - INFO.FattyRiot_version      : version number of FattyRiot.m
%
% Example:
%
% load('10.mat'); % challenge case 10
% FW = FattyRiot(imDataParams);
% [nx ny nz ncoils nTE] = size(imDataParams.images);
% F = FW(:,:,[1:nz]);
% W = FW(:,:,[1:nz]+nz);
% figure;
% subplot(2,1,1); imagesc(F(:,:)); axis image; colormap(gray); title('FAT');
% subplot(2,1,2); imagesc(W(:,:)); axis image; colormap(gray); title('WATER');
%

%
% last modified: 2013.03.29 by E. Brian Welch <brian.welch@vanderbilt.edu>
%
function [FW,INFO] = FattyRiot(imDataParams)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% overall algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1. initializations
%   a. start clock
%   b. FattyRiot version
%   d. detect path to this m-file
%   d. FW
%   e. INFO
%
% 2. validate input structure
%   a. imDataParams
%   b. imDataParams.images
%   c. imDataParams.TE
%   d. imDataParams.FieldStrength
%   e. imDataParams.PrecessionIsClockwise
%   f. imDataParams.mask
%   g. unhandled nTE values 
%
% 3. sanitize input
%   a. combine coils
%   b. account for precession direction
%   c. sort TE times
%   d. average redundant echo times
%   e. customize mask
%   f. detect crop indices
%
% 4. detect indices and number of evenly spaced echoes, idx_TEeven and nTEeven
% 
% 5. setup common algorithm settings (fat signal model)
%
% 6. R2* Calibration
%    * if (nTEeven>=4), use QPBO_2D (should always be satisfied for ISMRM 2012 Challenge)
%    * elseif (nTE>=4), use GC
%    * else   GLOBAL_maxR2 = 0; (should not happen for ISMRM 2012 Challenge)
%
% 7. if (nTEeven>=3) apply Berglund QPBO slicewise
%    * FM_HZ_QPBO_2D yielded as long as nTEeven>=3 (guaranteed for ISMRM 2012 Challenge) 
%    * R2_QPBO_2D yielded as long as nTEeven>=4 (not guaranteed)
%
% 8. if (nTEeven<=3), apply Hernando GraphCut slicewise with all echoes
%    * R2_GC yielded (needed because no R2 available from QPBO)
%
% 9. apply decomposeGivenFieldMapAndDampings using all echoes
%    * if (nTEeven>=4), use FM_HZ_QPBO_2D and R2_QPBO_2D
%    * else
%    *   if (nTEeven==3), use FM_HZ_QPBO_2D and R2_GC
%    *   else use FM_HZ_GC and R2_GC (nTEeven==2, should not happen in 2012 ISMRM Challenge)
%
% 10. wrap-up
%    a. build rest of INFO structure if necessary
%    b. print execution time
%
% 11. internal functions
%    a. FattyRiot_elapsed_time
%    b. print success summary
%    c. FattyRiot_time_remaining
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. initializations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1a. start clock
FattyRiot_start_time_original = tic;

% this is the time cushion required after R2* calibration
% it is not imperative that all slices be processed to estimate GLOBAL_maxR2
% the time is needed afterward for full QPBO_2D, GC, DECOMPOSE and DECOMPOSE SELECTION
% will spend at most 300 seconds (5 minutes) on R2* calibration
% one slice will always be done
time_cushion_very_long = 1500;

% this is the time cushion required to safely complete DECOMPOSE and DECOMPOSE SELECTION
time_cushion_long = 200;

% this is the time chusion to safely wrap-up and exit the function
time_cushion_short = 10;

out_of_time_cushion_very_long = 0;
out_of_time_cushion_long = 0;
out_of_time_cushion_short = 0;

% increase slice_runtime estimate by a little bit for timing calculations
slice_runtime_amplificaton_factor = 1.1;

% 1b. FattyRiot version
FattyRiot_version = 'v1.2_20130329';
disp( sprintf('FattyRiot: Version %s', FattyRiot_version) );

% 1c. detect path to this m-file
[FattyRiot_pathstr, name, ext] = fileparts(mfilename('fullpath'));
disp( sprintf('FattyRiot: Running from %s', FattyRiot_pathstr) );

% 1d. FW
FW = [];

% 1e. INFO
INFO = {};
INFO.FattyRiot_version = FattyRiot_version;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. validate input structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2a. imDataParams
if nargin~=1,
    error('FattyRiot: detected %d inputs (expected 1 input)', nargin);
else
    % initialize INFO
    INFO.imDataParams = imDataParams;
end

% 2b. imDataParams.images
if ~isfield(imDataParams,'images'),
    error('FattyRiot: input must contain field ''images'' with size [nx ny nz ncoils nTE]');
else
    if ndims(imDataParams.images)~=5,
        error('FattyRiot: imDataParams.images should have dimensions [nx ny nz ncoils nTE]'); % also produces error for nTE=1
    else
        [nx ny nz ncoils nTE] = size(imDataParams.images);
        
        % initialize FW
        [nx ny nz ncoils nTE] = size(imDataParams.images);
        FW = zeros([nx ny 2*nz]);
        FW(:,:,[1:nz]+nz) = ones([nx ny nz]); % dumbest answer possible, i.e. everything is 100% water
        FW_dumb_flag = ones(nz,1);
        success_GLOBAL_maxR2_calibration = zeros(nz,1);
        success_QPBO_2D = zeros(nz,1);
        success_QPBO_3D = zeros(nz,1);
        success_GC = zeros(nz,1);
        success_DECOMPOSITION_SELECTION = zeros(nz,1);
        attempted_GLOBAL_maxR2_calibration = 0;
        attempted_QPBO_2D = 0;
        attempted_QPBO_3D = 0;
        attempted_GC = 0;
        attempted_DECOMPOSE = 0;
    end
end

% 2c. imDataParams.TE
if ~isfield(imDataParams,'TE'),
    error('FattyRiot: input must contain field ''TE'' with number of elements nTE');
end
if length(imDataParams.TE(:))~=nTE,
    error('FattyRiot; imDataParams.TE has length %d (expected %d)', length(imDataParams.TE(:)), nTE);
end

% 2d. imDataParams.FieldStrength
if ~isfield(imDataParams,'FieldStrength'),
    imDataParams.FieldStrength = 1.5;
    warning('FattyRiot: defaulting imDataParams.FieldStrength to 1.5');
end

% 2e. imDataParams.PrecessionIsClockwise
if ~isfield(imDataParams,'PrecessionIsClockwise'),
    imDataParams.PrecessionIsClockwise = 1;
    warning('FattyRiot: defaulting imDataParams.PrecessionIsClockwise to 1');
end

% 2f. imDataParams.mask
if ~isfield(imDataParams,'mask'),
    imDataParams.mask = true([nx, ny, nz]);
    warning( sprintf('FattyRiot: defaulting imDataParams.mask to true([%d, %d, %d])', nx, ny, nz) );
else
    imDataParams.mask = logical(imDataParams.mask);
end

% 2g. unhandled nTE values
if nTE<3,
    warning( sprintf('FattyRiot: only %d echoes provided, but ISMRM 2012 Challenge organizers promised at least 3 evenly spaced echoes', nTE) );
    return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. sanitize input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3a. combine coils
if (ncoils>1),
    imDataParams.images = FattyRiot_coilCombine3D(imDataParams.images);
    ncoils = 1;
end

% 3b. account for precession direction
if imDataParams.PrecessionIsClockwise ~= 1, 
    imDataParams.images = conj(imDataParams.images);
    imDataParams.PrecessionIsClockwise = 1;
end

% 3c. sort TE times (just to be sure)
imDataParams.TE = imDataParams.TE(:);
[dummy,sorted_idx] = sort(imDataParams.TE);
imDataParams.TE = imDataParams.TE(sorted_idx);
imDataParams.images = imDataParams.images(:,:,:,1,sorted_idx);

% 3d. average redundant echo times
TE_usec_unique = unique( round(1e6*(imDataParams.TE(:))) );
nTE_usec_unique = length(TE_usec_unique(:));
if nTE_usec_unique~=nTE,
    disp( sprintf('FattyRiot: Found %d unique echo times out of %d provided echo times : Redundant echo times will be averaged', nTE_usec_unique, nTE) );
    images_tmp = zeros([nx ny nz 1 nTE_usec_unique]);
    TE_tmp = zeros(nTE_usec_unique,1);
    for kTE=1:nTE_usec_unique,
        idx_TE = find(TE_usec_unique(kTE)==imDataParams.TE(:));
        images_tmp(:,:,:,1,kTE) = sum(imDataParams.images(:,:,:,1,idx_TE),5)/length(idx_TE(:));
        TE_tmp(kTE) = imDataParams.TE(idx_TE(1));
    end
    imDataParams.images = images_tmp;
    imDataParams.TE = TE_tmp;
    clear images_tmp TE_tmp;
end

% 3e. customize mask
abs_echo_sum = zeros([nx ny nz]);
for kTE = 1:nTE,
    abs_echo_sum = abs_echo_sum + abs( squeeze(imDataParams.images(:,:,:,1,kTE)) );
end
abs_echo_sum = abs_echo_sum / max(abs_echo_sum(:));
mask_otsu = abs_echo_sum > graythresh(abs_echo_sum(:));
imDataParams.mask = imDataParams.mask | mask_otsu;
clear abs_echo_sum mask_otsu;

% 4f. detect crop indices
crop_cushion = 1;
for sl=1:nz,
    [idx_x idx_y] = find( imDataParams.mask(:,:,sl) );
    idx_crop_x_start = max( 1,min(idx_x(:))-crop_cushion);
    idx_crop_x_stop  = min(nx,max(idx_x(:))+crop_cushion);
    idx_crop_y_start = max( 1,min(idx_y(:))-crop_cushion);
    idx_crop_y_stop  = min(ny,max(idx_y(:))+crop_cushion);
    idx_crop_x{sl} = [idx_crop_x_start:idx_crop_x_stop];
    idx_crop_y{sl} = [idx_crop_y_start:idx_crop_y_stop];
    disp( sprintf('FattyRiot: Cropping slice %d to save %.2f%%', sl, 100*(1-length(idx_crop_x{sl})*length(idx_crop_y{sl})/(nx*ny)) ) );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. detect indices and number of evenly spaced echoes, idx_TEeven and nTEeven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dTE_usec = round(1e6*diff(imDataParams.TE(:)));
dTE_candidates = dTE_usec(:);
for idx_dTE = 1:length(dTE_usec(:)),
    dTE_candidates = unique([dTE_candidates(:) ; cumsum(dTE_usec(idx_dTE:end)) ]);
end
dTE_candidates = unique(dTE_candidates(:));

best_echoes = [];
for idx_dTE_candidate = 1:length(dTE_candidates),
    target_dTE_usec = dTE_candidates(idx_dTE_candidate);
    ne = length(imDataParams.TE);

    % initialize dTE_groups
    dTE_groups = {};
    for idx=1:(ne-1),
        dTE_groups{idx}.dTE_usec = dTE_usec(idx);
        dTE_groups{idx}.echoes = [idx idx+1];
    end

    % combine dTE_groups
    for idx=1:(ne-1),
        if (dTE_groups{idx}.dTE_usec == target_dTE_usec),
            % do nothing
        elseif (dTE_groups{idx}.dTE_usec < target_dTE_usec)
            % combine with next group
                if idx<(ne-1),
                    dTE_groups{idx+1}.dTE_usec = dTE_groups{idx}.dTE_usec + dTE_groups{idx+1}.dTE_usec;
                    dTE_groups{idx+1}.echoes = [dTE_groups{idx}.echoes(1:end-1) dTE_groups{idx+1}.echoes(2:end)];
                    dTE_groups{idx}.dTE_usec = NaN;
                    dTE_groups{idx}.echoes = [];
                else
                    dTE_groups{idx}.dTE_usec = NaN;
                    dTE_groups{idx}.echoes = [];
                end
        else
            % echo spacing is too large
            dTE_groups{idx}.dTE_usec = NaN;
            dTE_groups{idx}.echoes = [];
        end
    end

    % find usable TE set
    echoes_to_use = [];
    idx = 1;
    while idx <= (ne-1),
        if ( isnan(dTE_groups{idx}.dTE_usec) && (length(echoes_to_use)>0) ),
            break;
        end
        if ~isnan(dTE_groups{idx}.dTE_usec),
            echoes_to_use = [echoes_to_use dTE_groups{idx}.echoes];
        end
        idx = idx + 1;
    end

    echoes_to_use = sort( unique(echoes_to_use) );

    % test if this is better (longer) echo set
    if length(echoes_to_use) > length(best_echoes),
        best_echoes = echoes_to_use;
    end
end

idx_TEeven = best_echoes;
nTEeven = length(idx_TEeven);

% double check that detect evenly space echoes are truly evenly spaced
dTEeven_usec = round(1e6*diff(imDataParams.TE(idx_TEeven)));
if all(dTEeven_usec==dTEeven_usec(1)),
    dTEeven = imDataParams.TE(idx_TEeven(2)) - imDataParams.TE(idx_TEeven(1));
    disp( sprintf('FattyRiot: Available evenly space echoes confirmed with delta TE = %.6f seconds', dTEeven) );
end

if nTEeven<3,
    warning('FattyRiot: Could not find at least 3 evenly spaced echoes');
else
    disp( sprintf('FattyRiot: Found %d of %d evenly spaced echoes. idx_nTEeven = [ %s]', nTEeven, nTE, sprintf('%d ', idx_TEeven) ) );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. setup common algorithm settings (fat signal model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma_Hz_per_Tesla = 42.577481e6;
%gamma_Hz_per_Tesla = 42.58e6;
algoParams = {};
algoParams.species(1).name = 'water' ; % Water
algoParams.species(1).frequency = [0] ;
algoParams.species(1).relAmps = [1] ;  
algoParams.species(2).name = 'fat' ; % Fat
algoParams.species(2).frequency = 1e6 * [-242.7060, -217.1580, -166.0620, -123.9078, -24.9093, 38.3220] / (1.5 * gamma_Hz_per_Tesla); 
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048] ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. R2* Calibration
%    * if (nTEeven>=4), use QPBO_2D (should always be satisfied for ISMRM 2012 Challenge)
%    * elseif (nTE>=4), use GC
%    * else   GLOBAL_maxR2 = 0; (should not happen for ISMRM 2012 Challenge)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GLOBAL_maxR2 = 0; % will be calibrated in a moment
calibration_maxR2 = 800; % use a big R2* value for calibration
%calibration_maxR2 = 1000; % use a big R2* value for calibration

try
    FattyRiot_start_time_GLOBAL_maxR2_calibration = tic;
    attempted_GLOBAL_maxR2_calibration = 1;
    R2_GLOBAL_maxR2_calibration = zeros([nx ny nz]);
    
    if ( (nTEeven>=4) ),
    %if ( (nTEeven>=4) && (nTEeven==nTE) ),

        disp( sprintf('FattyRiot: Attempting R2* Calibration using QPBO_2D with calibration_maxR2 = %.2f', calibration_maxR2) );
        
        % QPBO algorithm settings
        algoParams_QPBO = algoParams;
        algoParams_QPBO.decoupled_estimation = true; % flag for decoupled R2 estimation
        algoParams_QPBO.Fibonacci_search = true; % flag for Fibonacci search
        algoParams_QPBO.B0_smooth_in_stack_direction = false; % flag for B0 smooth in stack direction
        algoParams_QPBO.multigrid = true; % flag for multi-level resolution pyramid
        algoParams_QPBO.estimate_R2 = true; % flag to estimate R2star
        %algoParams_QPBO.verbose = true; % flag for verbose status messages
        algoParams_QPBO.verbose = false; % flag for verbose status messages
        algoParams_QPBO.ICM_iterations = 2; % ICM iterations
        algoParams_QPBO.num_B0_labels = 100; % number of discretized B0 values
        algoParams_QPBO.mu = 10; % regularization parameter
        %algoParams_QPBO.R2_stepsize = 1; % R2 stepsize in s^-1
        algoParams_QPBO.R2_stepsize = 2; % R2 stepsize in s^-1
        algoParams_QPBO.max_R2 = calibration_maxR2; % maximum R2 in s^-1
        algoParams_QPBO.max_label_change = 0.1; % 
        %algoParams_QPBO.fine_R2_stepsize = 1.0; % Fine stepsize for discretization of R2(*) [sec-1] (used in decoupled fine-tuning step)
        algoParams_QPBO.fine_R2_stepsize = 2.0; % Fine stepsize for discretization of R2(*) [sec-1] (used in decoupled fine-tuning step)    
        %algoParams_QPBO.coarse_R2_stepsize = 10.0; % Coarse stepsize for discretization of R2(*) [sec-1] (used for joint estimation step, value 0.0 --> Decoupled estimation only
        algoParams_QPBO.coarse_R2_stepsize = 20.0; % Coarse stepsize for discretization of R2(*) [sec-1] (used for joint estimation step, value 0.0 --> Decoupled estimation only

        % QPBO image changes
        imDataParams_QPBO = imDataParams;
        imDataParams_QPBO.images = imDataParams_QPBO.images(:,:,:,1,idx_TEeven);
        imDataParams_QPBO.TE = imDataParams_QPBO.TE(idx_TEeven);        

        imDataParams_QPBO_2D.TE = imDataParams_QPBO.TE;
        imDataParams_QPBO_2D.FieldStrength = imDataParams.FieldStrength;
        imDataParams_QPBO_2D.PrecessionIsClockwise = imDataParams.PrecessionIsClockwise;

        try
            if exist('FattyRiot_fw_i3cm1i_3pluspoint_berglund_QPBO.m')~=2,
                addpath( sprintf('%s/FattyRiot_toolbox/berglund/QPBO', FattyRiot_pathstr) );
            end
        catch
            warning('FattyRiot: Encountered error adding path to FattyRiot_fw_i3cm1i_3pluspoint_berglund_QPBO.m');
            %rethrow(lasterror); 
        end        
        
        try
            for sl=1:nz,

                imDataParams_QPBO_2D.images = imDataParams_QPBO.images(idx_crop_x{sl},idx_crop_y{sl},sl,1,:);
                result = FattyRiot_fw_i3cm1i_3pluspoint_berglund_QPBO(imDataParams_QPBO_2D,algoParams_QPBO);
                R2_GLOBAL_maxR2_calibration(idx_crop_x{sl},idx_crop_y{sl},sl) = result.r2starmap;

                % update success flag
                success_GLOBAL_maxR2_calibration(sl) = 1;

                % update average slice runtime
                average_slice_runtime = FattyRiot_elapsed_time(FattyRiot_start_time_GLOBAL_maxR2_calibration) / sl;
                disp( sprintf('FattyRiot: R2* Calibration with QPBO_2D average slice runtime (%d/%d) = %.2f seconds', sl, nz, average_slice_runtime) );
                
                % check time
                if (FattyRiot_time_remaining(time_cushion_very_long+slice_runtime_amplificaton_factor*average_slice_runtime)<0),
                    out_of_time_cushion_very_long = 1;
                    error('FattyRiot: OUT OF TIME!');
                end

            end
        catch
            if (out_of_time_cushion_very_long==1),
                disp( sprintf('FattyRiot: R2* calibration with QPBO_2D ran out of time at slice (%d/%d)', sl, nz) );
            else
                warning('FattyRiot: Encountered unexpected error performing maximum R2* calibration with QPBO_2D');
                %rethrow(lasterror);
            end
        end    

    elseif (nTE>=4), % Use GC (will not happen for Phase I and Phase II Cases of ISMRM Challenge 2012)
        warning('FattyRiot: GC called to perform GLOBAL_maxR2 calibration');
        
        disp( sprintf('FattyRiot: Attempting R2* Calibration using GC with calibration_maxR2 = %.2f', calibration_maxR2) );
        
        %% Graph Cut algorithm settings
        algoParams_GC = algoParams;
        algoParams_GC.size_clique = 1;           % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
        algoParams_GC.range_r2star = [0 calibration_maxR2];    % Range of R2* values
        algoParams_GC.NUM_R2STARS = 1;           % Number of R2* values for quantization
        algoParams_GC.range_fm = [-400 400];     % Range of field map values
        algoParams_GC.NUM_FMS = 301;             % Number of field map values to discretize
        algoParams_GC.NUM_ITERS = 40;            % Number of graph cut iterations
        algoParams_GC.SUBSAMPLE = 1;             % Spatial subsampling for field map estimation (for speed)
        algoParams_GC.DO_OT = 1;                 % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
        algoParams_GC.LMAP_POWER = 2;            % Spatially-varying regularization (2 gives ~ uniformn resolution)
        algoParams_GC.lambda = 0.05;             % Regularization parameter
        algoParams_GC.LMAP_EXTRA = 0.05;         % More smoothing for low-signal regions
        algoParams_GC.TRY_PERIODIC_RESIDUAL = 0; % Take advantage of periodic residual if uniform TEs (will change range_fm)  

        try
            if exist('FattyRiot_fw_i2cm1i_3pluspoint_hernando_graphcut.m')~=2,
                addpath( sprintf('%s/FattyRiot_toolbox/hernando/graphcut', FattyRiot_pathstr) );
            end
        catch
            warning('FattyRiot: Encountered error adding path to FattyRiot_fw_i2cm1i_3pluspoint_hernando_graphcut.m');
            %rethrow(lasterror); 
        end

        try
            if exist('FattyRiot_fw_i2cm0i_3plusploint_hernando_optimtransfer.m')~=2,
                addpath( sprintf('%s/FattyRiot_toolbox/hernando/descent', FattyRiot_pathstr) );
            end
        catch
            warning('FattyRiot: Encountered error adding path to FattyRiot_fw_i2cm0i_3plusploint_hernando_optimtransfer.m');
            %rethrow(lasterror); 
        end

        try
            if exist('FattyRiot_decomposeGivenFieldMapAndDampings.m')~=2,
                addpath( sprintf('%s/FattyRiot_toolbox/hernando/common', FattyRiot_pathstr) );
            end
        catch
            warning('FattyRiot: Encountered error adding path to FattyRiot_decomposeGivenFieldMapAndDampings.m');
            %rethrow(lasterror); 
        end

        try
            addpath( sprintf('%s/FattyRiot_toolbox/hernando/matlab_bgl', FattyRiot_pathstr) );
        catch
            warning('FattyRiot: Encountered error adding path to matlab_bgl');
            %rethrow(lasterror); 
        end    
    
        try

            imDataParams_GC_slice.FieldStrength = imDataParams.FieldStrength;
            imDataParams_GC_slice.TE = imDataParams.TE;
            imDataParams_GC_slice.PrecessionIsClockwise = imDataParams.PrecessionIsClockwise;

            for sl=1:nz,

                % use all echoes
                imDataParams_GC_slice.images = imDataParams.images(idx_crop_x{sl},idx_crop_y{sl},sl,1,1:nTE);
                imDataParams_GC_slice.mask = imDataParams.mask(idx_crop_x{sl},idx_crop_y{sl},sl);

                warning off;
                result = FattyRiot_fw_i2cm1i_3pluspoint_hernando_graphcut(imDataParams_GC_slice, algoParams_GC);
                warning on;

                R2_GLOBAL_maxR2_calibration(idx_crop_x{sl},idx_crop_y{sl},sl) = result.r2starmap;
                        
                % update success flag
                success_GLOBAL_maxR2_calibration(sl) = 1;

                % update average slice runtime
                average_slice_runtime = FattyRiot_elapsed_time(FattyRiot_start_time_GLOBAL_maxR2_calibration) / sl;
                disp( sprintf('FattyRiot: R2* Calibration with GC average slice runtime (%d/%d) = %.2f seconds', sl, nz, average_slice_runtime) );
                
                % check time
                if (FattyRiot_time_remaining(time_cushion_very_long+slice_runtime_amplificaton_factor*average_slice_runtime)<0),
                    out_of_time_cushion_very_long = 1;
                    error('FattyRiot: OUT OF TIME!');
                end

            end
            
        catch
            if (out_of_time_cushion_very_long==1),
                disp( sprintf('FattyRiot: R2* calibration with GC ran out of time at slice (%d/%d)', sl, nz) );
            else
                warning('FattyRiot: Encountered unexpected error performing maximum R2* calibration with GC');
                %rethrow(lasterror);
            end
        end 
        
    else
        % should not happen for ISMRM 2012 Challenge, GLOBAL_maxR2 will remain 0
        warning('FattyRiot: Not enough echoes to perform GLOBAL_maxR2 calibration');
    end

    try
        if (any(success_GLOBAL_maxR2_calibration)),
            
            mask_R2cal = zeros([nx ny nz]);
            for sl=1:nz,
                if (success_GLOBAL_maxR2_calibration(sl)==1),
                    mask_R2cal(:,:,sl) = imDataParams.mask(:,:,sl);
                end
            end
            
            idx_mask_R2cal = find(mask_R2cal(:)==1);
            numel_idx_mask_R2cal = numel(idx_mask_R2cal);
            num_bins = 1000;
            %num_bins = 10000;
            bins = linspace(0,calibration_maxR2,num_bins);
            [bin_counts , bin_centers] = hist(R2_GLOBAL_maxR2_calibration(idx_mask_R2cal), bins );
            bin_cumsum_percent = (cumsum(bin_counts)/numel_idx_mask_R2cal);
            INFO.bin_cumsum_percent = bin_cumsum_percent;
            cdf_threshold = 0.9975;
            idx_cumsum_above_threshold = find( bin_cumsum_percent >= cdf_threshold);
            
            % check to see if cdf_threshold should be relaxed
            while ( (idx_cumsum_above_threshold(1) == num_bins) && (cdf_threshold>0) ),
                cdf_threshold = cdf_threshold - 0.0150;
                %cdf_threshold = cdf_threshold - 0.0125;
                %cdf_threshold = cdf_threshold - 0.0050;
                idx_cumsum_above_threshold = find( bin_cumsum_percent >= cdf_threshold);
            end
            GLOBAL_maxR2  = ceil( bin_centers(idx_cumsum_above_threshold(1)) );
            
            % default to something close to reasonable if this cdf method completely bombs
            if (idx_cumsum_above_threshold(1) == num_bins),
                if (imDataParams.FieldStrength>1.6),
                    GLOBAL_maxR2 = 350;
                else
                    GLOBAL_maxR2 = 200;
                end
                warning('FattyRiot: R2* Calibration ... could not find a GLOBAL_maxR2 less than calibration_maxR2');
            end
                    
            disp( sprintf('FattyRiot: R2* Calibration ... CDF threshold %.4f satisfied with GLOBAL_maxR2 = %.2f', cdf_threshold, GLOBAL_maxR2) );
        else
            if (imDataParams.FieldStrength>1.6),
                GLOBAL_maxR2 = 350;
            else
                GLOBAL_maxR2 = 200;
            end
            warning('FattyRiot: R2* Calibration ... could not find a GLOBAL_maxR2 less than calibration_maxR2');
        end
        
        disp( sprintf('FattyRiot: GLOBAL_maxR2 calibration completed in %.2f seconds', FattyRiot_elapsed_time(FattyRiot_start_time_GLOBAL_maxR2_calibration) ) );
    
    catch
        warning('FattyRiot: Encountered unexpected error finding R2* calibration CDF threshold');
        %rethrow(lasterror);
    end
    
catch
    warning('FattyRiot: Encountered unexpected error performing maximum R2* calibration');
    %rethrow(lasterror); 
end
INFO.GLOBAL_maxR2 = GLOBAL_maxR2;
disp( sprintf('FattyRiot: GLOBAL_maxR2 = %.2f', GLOBAL_maxR2) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. if (nTEeven>=3) apply Berglund QPBO slicewise
%    * FM_HZ_QPBO_2D yielded as long as nTEeven>=3 (guaranteed for ISMRM 2012 Challenge) 
%    * R2_QPBO_2D yielded as long as nTEeven>=4 (not guaranteed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nTEeven>=3),

    FattyRiot_start_time_QPBO_2D = tic;
    attempted_QPBO_2D = 1;
    disp('FattyRiot: Attempting QPBO_2D');
    
    % QPBO algorithm settings
    algoParams_QPBO = algoParams;
    algoParams_QPBO.decoupled_estimation = true; % flag for decoupled R2 estimation
    algoParams_QPBO.Fibonacci_search = true; % flag for Fibonacci search
    algoParams_QPBO.B0_smooth_in_stack_direction = false; % flag for B0 smooth in stack direction
    algoParams_QPBO.multigrid = true; % flag for multi-level resolution pyramid
    algoParams_QPBO.estimate_R2 = true; % flag to estimate R2star
    %algoParams_QPBO.verbose = true; % flag for verbose status messages
    algoParams_QPBO.verbose = false; % flag for verbose status messages
    algoParams_QPBO.ICM_iterations = 2; % ICM iterations
    algoParams_QPBO.num_B0_labels = 100; % number of discretized B0 values
    algoParams_QPBO.mu = 10; % regularization parameter
    algoParams_QPBO.R2_stepsize = 1; % R2 stepsize in s^-1
    algoParams_QPBO.max_R2 = GLOBAL_maxR2; % maximum R2 in s^-1
    algoParams_QPBO.max_label_change = 0.1; % 
    algoParams_QPBO.fine_R2_stepsize = 1.0; % Fine stepsize for discretization of R2(*) [sec-1] (used in decoupled fine-tuning step)
    algoParams_QPBO.coarse_R2_stepsize = 10.0; % Coarse stepsize for discretization of R2(*) [sec-1] (used for joint estimation step, value 0.0 --> Decoupled estimation only

    % QPBO image changes
    imDataParams_QPBO = imDataParams;
    imDataParams_QPBO.images = imDataParams_QPBO.images(:,:,:,1,idx_TEeven);
    imDataParams_QPBO.TE = imDataParams_QPBO.TE(idx_TEeven);
    
    try
        if exist('FattyRiot_fw_i3cm1i_3pluspoint_berglund_QPBO.m')~=2,
            addpath( sprintf('%s/FattyRiot_toolbox/berglund/QPBO', FattyRiot_pathstr) );
        end
    catch
        warning('FattyRiot: Encountered error adding path to FattyRiot_fw_i3cm1i_3pluspoint_berglund_QPBO.m');
        %rethrow(lasterror); 
    end

    try
        W_QPBO_2D = zeros([nx ny nz]);
        F_QPBO_2D = zeros([nx ny nz]);
        FM_DEGREES_QPBO_2D = zeros([nx ny nz]);
        FM_HZ_QPBO_2D = zeros([nx ny nz]);
        R2_QPBO_2D = zeros([nx ny nz]);
        imDataParams_QPBO_2D.TE = imDataParams_QPBO.TE;
        imDataParams_QPBO_2D.FieldStrength = imDataParams.FieldStrength;
        imDataParams_QPBO_2D.PrecessionIsClockwise = imDataParams.PrecessionIsClockwise;
        for sl=1:nz,
            
            imDataParams_QPBO_2D.images = imDataParams_QPBO.images(idx_crop_x{sl},idx_crop_y{sl},sl,1,:);
            result = FattyRiot_fw_i3cm1i_3pluspoint_berglund_QPBO(imDataParams_QPBO_2D,algoParams_QPBO);
            W_QPBO_2D(idx_crop_x{sl},idx_crop_y{sl},sl) = result.species(1).amps;
            F_QPBO_2D(idx_crop_x{sl},idx_crop_y{sl},sl) = result.species(2).amps;
            FM_DEGREES_QPBO_2D(idx_crop_x{sl},idx_crop_y{sl},sl) = result.fieldmap;
            if (nTEeven>=4),
                R2_QPBO_2D(idx_crop_x{sl},idx_crop_y{sl},sl) = result.r2starmap;
            end
                        
            FM_HZ_QPBO_2D(:,:,sl) = FM_DEGREES_QPBO_2D(:,:,sl) * ( pi / 180.0 ) / (2*pi*dTEeven);
            
            % overwrite dumb (100% water) result if nothing has replaced it yet
            if FW_dumb_flag(sl)==1,
                FW(:,:,sl)    = abs(F_QPBO_2D(:,:,sl));
                FW(:,:,sl+nz) = abs(W_QPBO_2D(:,:,sl));
                FW_dumb_flag(sl) = 0;
            end
            
            % update success flag
            success_QPBO_2D(sl) = 1;
            
            % update average slice runtime
            average_slice_runtime = FattyRiot_elapsed_time(FattyRiot_start_time_QPBO_2D) / sl;
            disp( sprintf('FattyRiot: QPBO_2D average slice runtime (%d/%d) = %.2f seconds', sl, nz, average_slice_runtime) );
                
            % check time
            if (FattyRiot_time_remaining(time_cushion_long++slice_runtime_amplificaton_factor*average_slice_runtime)<0),
                out_of_time_cushion_long = 1;
                error('FattyRiot: OUT OF TIME!');
            end
            
        end
        clear F_QPBO_2D W_QPBO_2D FM_DEGREES_QPBO_2D;
        
        disp( sprintf('FattyRiot: QPBO_2D completed in %.2f seconds', FattyRiot_elapsed_time(FattyRiot_start_time_QPBO_2D) ) );
    catch
        if (out_of_time_cushion_long==1),
            disp( sprintf('FattyRiot: QPBO_2D ran out of time at slice (%d/%d)', sl, nz) );
        else
            warning('FattyRiot: Encountered unexpected problem running QPBO_2D');
            %rethrow(lasterror);
        end
    end

else    
    warning( sprintf('FattyRiot: Could not run QPBO_2D because nTEeven=%d', nTEeven) );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8. if ( (out_of_time_cushion_long==0) ), apply Hernando GraphCut slicewise with all echoes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( (out_of_time_cushion_long==0) ),
    
    FattyRiot_start_time_GC = tic;
    attempted_GC = 1;
    disp('FattyRiot: Attempting GC');
    
    if (nTEeven<=3),
        disp('FattyRiot: Need to run GC because nTEeven<=3');
    end

    if all(success_QPBO_2D==0),
        disp('FattyRiot: Need to run GC because QPBO_2D failed to succeed');
    end
    
    %% Graph Cut algorithm settings
    algoParams_GC = algoParams;
    algoParams_GC.size_clique = 1;           % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
    algoParams_GC.range_r2star = [0 GLOBAL_maxR2];    % Range of R2* values
    algoParams_GC.NUM_R2STARS = 1;           % Number of R2* values for quantization
    algoParams_GC.range_fm = [-400 400];     % Range of field map values
    algoParams_GC.NUM_FMS = 301;             % Number of field map values to discretize
    algoParams_GC.NUM_ITERS = 40;            % Number of graph cut iterations
    algoParams_GC.SUBSAMPLE = 1;             % Spatial subsampling for field map estimation (for speed)
    algoParams_GC.DO_OT = 1;                 % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
    algoParams_GC.LMAP_POWER = 2;            % Spatially-varying regularization (2 gives ~ uniformn resolution)
    algoParams_GC.lambda = 0.05;             % Regularization parameter
    algoParams_GC.LMAP_EXTRA = 0.05;         % More smoothing for low-signal regions
    algoParams_GC.TRY_PERIODIC_RESIDUAL = 0; % Take advantage of periodic residual if uniform TEs (will change range_fm)  

    W_GC = zeros([nx ny nz]);
    F_GC = zeros([nx ny nz]);
    FM_HZ_GC = zeros([nx ny nz]);
    R2_GC = zeros([nx ny nz]);

    try
        if exist('FattyRiot_fw_i2cm1i_3pluspoint_hernando_graphcut.m')~=2,
            addpath( sprintf('%s/FattyRiot_toolbox/hernando/graphcut', FattyRiot_pathstr) );
        end
    catch
        warning('FattyRiot: Encountered error adding path to FattyRiot_fw_i2cm1i_3pluspoint_hernando_graphcut.m');
        %rethrow(lasterror); 
    end

    try
        if exist('FattyRiot_fw_i2cm0i_3plusploint_hernando_optimtransfer.m')~=2,
            addpath( sprintf('%s/FattyRiot_toolbox/hernando/descent', FattyRiot_pathstr) );
        end
    catch
        warning('FattyRiot: Encountered error adding path to FattyRiot_fw_i2cm0i_3plusploint_hernando_optimtransfer.m');
        %rethrow(lasterror); 
    end

    try
        if exist('FattyRiot_decomposeGivenFieldMapAndDampings.m')~=2,
            addpath( sprintf('%s/FattyRiot_toolbox/hernando/common', FattyRiot_pathstr) );
        end
    catch
        warning('FattyRiot: Encountered error adding path to FattyRiot_decomposeGivenFieldMapAndDampings.m');
        %rethrow(lasterror); 
    end

    try
        addpath( sprintf('%s/FattyRiot_toolbox/hernando/matlab_bgl', FattyRiot_pathstr) );
    catch
        warning('FattyRiot: Encountered error adding path to matlab_bgl');
        %rethrow(lasterror); 
    end    
    
    try
        %error('FattyRiot: Forcing skip of GC');

        imDataParams_GC_slice.FieldStrength = imDataParams.FieldStrength;
        imDataParams_GC_slice.TE = imDataParams.TE;
        imDataParams_GC_slice.PrecessionIsClockwise = imDataParams.PrecessionIsClockwise;

        for sl=1:nz,
            
            % use all echoes
            imDataParams_GC_slice.images = imDataParams.images(idx_crop_x{sl},idx_crop_y{sl},sl,1,1:nTE);
            imDataParams_GC_slice.mask = imDataParams.mask(idx_crop_x{sl},idx_crop_y{sl},sl);
            
            warning off;
            result = FattyRiot_fw_i2cm1i_3pluspoint_hernando_graphcut(imDataParams_GC_slice, algoParams_GC);
            warning on;
            
            W_GC(idx_crop_x{sl},idx_crop_y{sl},sl) = result.species(1).amps;
            F_GC(idx_crop_x{sl},idx_crop_y{sl},sl) = result.species(2).amps;
            FM_HZ_GC(idx_crop_x{sl},idx_crop_y{sl},sl) = result.fieldmap;
            R2_GC(idx_crop_x{sl},idx_crop_y{sl},sl) = result.r2starmap;
            
            % overwrite dumb (100% water) result if nothing has replaced it yet
            if FW_dumb_flag(sl)==1,
                FW(:,:,sl)    = abs(F_GC(:,:,sl));
                FW(:,:,sl+nz) = abs(W_GC(:,:,sl));
                FW_dumb_flag(sl) = 0;
            end
            
            % update success flag
            success_GC(sl) = 1;
            
            % update average slice runtime
            average_slice_runtime = FattyRiot_elapsed_time(FattyRiot_start_time_GC) / sl;
            disp( sprintf('FattyRiot: GC average slice runtime (%d/%d) = %.2f seconds', sl, nz, average_slice_runtime) );
            
            % check time
            if (FattyRiot_time_remaining(time_cushion_long+slice_runtime_amplificaton_factor*average_slice_runtime)<0),
                out_of_time_cushion_long = 1;
                error('FattyRiot: OUT OF TIME!');
            end
            
        end
        clear W_GC F_GC;
        
        disp( sprintf('FattyRiot: GC completed in %.2f seconds', FattyRiot_elapsed_time(FattyRiot_start_time_GC) ) );
    catch
        if (out_of_time_cushion_long==1),
            disp( sprintf('FattyRiot: GC ran out of time at slice (%d/%d)', sl, nz) );
        else
            warning('FattyRiot: Encountered unexpected problem running GC');
            %rethrow(lasterror);
        end
    end

else
    disp( sprintf('FattyRiot: No need to run GC because nTEeven=%d (>=4)', nTEeven) );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 9. apply DECOMPOSITIONS using all echoes
%    * FM_HZ_QPBO_2D and R2_QPBO_2D
%    * FM_HZ_GC and R2_GC
%    * FM_HZ_QPBO_2D and R2_GC
%    * FM_HZ_GC and R2_QPBO_2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% begin DECOMPOSITIONS
try
    
    FattyRiot_start_time_DECOMPOSITIONS = tic;
    attempted_DECOMPOSITIONS = 1;
    disp('FattyRiot: Attempting DECOMPOSE');
     
    % intitialize empty FW_options structure
    FW_options = {};

    % PURE DECOMPOSITIONS
    
    % FM_HZ_QPBO_2D and R2_QPBO_2D
    if( any(success_QPBO_2D==1) ),
        try
            DECOMPOSE_struct = FattyRiot_DECOMPOSE(FM_HZ_QPBO_2D,R2_QPBO_2D,success_QPBO_2D);
            FW_options{end+1} = DECOMPOSE_struct;
            FW_options{end}.name = 'FM_HZ_QPBO_2D,R2_QPBO_2D'; 
        catch
            warning('FattyRiot: Encountered unexpected problem running DECOMPOSITIONS with FM_HZ_QPBO_2D and R2_QPBO_2D');
        end
    else
        disp('FattyRiot: No QPBO_2D slices to DECOMPOSE');
    end
 
    % FM_HZ_GC and R2_GC
    %success_GC(end)=0; % test slice skipping
    %success_GC(:)=0;   % test GC failure, e.g. no time
    if( any(success_GC==1) ),
        try
            DECOMPOSE_struct = FattyRiot_DECOMPOSE(FM_HZ_GC,R2_GC,success_GC);
            FW_options{end+1} = DECOMPOSE_struct;
            FW_options{end}.name = 'FM_HZ_GC,R2_GC'; 
        catch
            warning('FattyRiot: Encountered unexpected problem running DECOMPOSITIONS with FM_HZ_GC and R2_GC');
        end
    else
        disp('FattyRiot: No GC slices to DECOMPOSE');
    end

    % MIXED DECOMPOSITIONS   

    % FM_HZ_QPBO_2D and R2_GC
%     try
%         DECOMPOSE_struct = FattyRiot_DECOMPOSE(FM_HZ_QPBO_2D,R2_GC);
%         FW_options{end+1} = DECOMPOSE_struct;
%         FW_options{end}.name = 'FM_HZ_QPBO_2D,R2_GC'; 
%     catch
%         warning('FattyRiot: Encountered problem running DECOMPOSITIONS with FM_HZ_QPBO_2D and R2_GC');
%     end       
    
    % FM_HZ_GC and R2_QPBO_2D
%     try
%         DECOMPOSE_struct = FattyRiot_DECOMPOSE(FM_HZ_GC,R2_QPBO_2D);
%         FW_options{end+1} = DECOMPOSE_struct;
%         FW_options{end}.name = 'FM_HZ_GC,R2_QPBO_2D'; 
%     catch
%         warning('FattyRiot: Encountered problem running DECOMPOSITIONS with FM_HZ_GC and R2_QPBO_2D');
%     end        
      
    disp( sprintf('FattyRiot: DECOMPOSITIONS completed in %.2f seconds', FattyRiot_elapsed_time(FattyRiot_start_time_DECOMPOSITIONS) ) );
catch    
    warning('FattyRiot: Encountered unexpected problem running DECOMPOSITIONS');
    %rethrow(lasterror);
end

% select between decompositions
try
    
    FattyRiot_start_time_DECOMPOSITION_SELECTION = tic;
    attempted_DECOMPOSITION_SELECTION = 1;
    disp('FattyRiot: Attempting DECOMPOSITION_SELECTION');
    
    num_FW_options = length(FW_options);
    INFO.FW_options_REMERROR_DECOMPOSE = zeros([nx ny nz]);
    INFO.FW_options_map = zeros([nx ny nz]);
    INFO.FW_options_map_filtered = zeros([nx ny nz]);
    
    INFO.F = zeros([nx ny nz]);
    INFO.W = zeros([nx ny nz]);
    INFO.FM_HZ = zeros([nx ny nz]);
    INFO.R2 = zeros([nx ny nz]);
    INFO.ERROR = zeros([nx ny nz]);
    
    FW_options_names_str = '';
    for idx_FW_option = 1:num_FW_options,
        INFO.FW_options_names{idx_FW_option} = FW_options{idx_FW_option}.name;
        FW_options_names_str = sprintf('%s %s ', FW_options_names_str, INFO.FW_options_names{idx_FW_option});
    end
    
    if (num_FW_options>1),
    
        disp( sprintf('FattyRiot: DECOMPOSITION_SELECTION choosing among [ %s ]', FW_options_names_str ) );

        for kz=1:nz,
            % cover entire crop boundary box
            idx_x = idx_crop_x{kz};
            idx_y = idx_crop_y{kz};

            for kx=idx_x,
                for ky=idx_y,

                        best_error = Inf;
                        best_FW_option = 0;

                        for idx_FW_option = 1:num_FW_options,

                            if (FW_options{idx_FW_option}.REMERROR_DECOMPOSE(kx,ky,kz) < best_error),
                                best_error = FW_options{idx_FW_option}.REMERROR_DECOMPOSE(kx,ky,kz);
                                best_FW_option = idx_FW_option;
                            end

                        end % end idx_FW_option

                        INFO.FW_options_REMERROR_DECOMPOSE(kx,ky,kz) = best_error;
                        INFO.FW_options_map(kx,ky,kz) = best_FW_option;
                        %FW(kx,ky,kz)    = FW_options{best_FW_option}.FW(kx,ky,kz);
                        %FW(kx,ky,kz+nz) = FW_options{best_FW_option}.FW(kx,ky,kz+nz);

                end % end ky
            end % end kx

            % mode filter the FW_options_map
            neighborhood_radius = 1; % 3x3
            min_idx_x = min(idx_x);
            max_idx_x = max(idx_x);
            min_idx_y = min(idx_y);
            max_idx_y = max(idx_y);

            for kx=idx_x,
                for ky=idx_y,

                    dx_start = max(min_idx_x,kx-neighborhood_radius);
                    dx_stop  = min(max_idx_x,kx+neighborhood_radius);
                    dy_start = max(min_idx_y,ky-neighborhood_radius);
                    dy_stop  = min(max_idx_y,ky+neighborhood_radius);
                    vote_count = zeros(num_FW_options,1);
                    for dx=dx_start:dx_stop,
                        for dy=dy_start:dy_stop,
                            this_FW_option = INFO.FW_options_map(dx,dy,kz);
                            vote_count( this_FW_option ) = vote_count( this_FW_option ) + 1;
                        end
                    end
                    [vote_max, idx_vote_max] = max(vote_count);

                    INFO.FW_options_map_filtered(kx,ky,kz) = idx_vote_max;

                    INFO.F(kx,ky,kz)     = FW_options{idx_vote_max}.F_DECOMPOSE(kx,ky,kz);
                    INFO.W(kx,ky,kz)     = FW_options{idx_vote_max}.W_DECOMPOSE(kx,ky,kz);
                    INFO.FM_HZ(kx,ky,kz) = FW_options{idx_vote_max}.FM_HZ(kx,ky,kz);
                    INFO.R2(kx,ky,kz)    = FW_options{idx_vote_max}.R2(kx,ky,kz);;
                    INFO.ERROR(kx,ky,kz) = FW_options{idx_vote_max}.REMERROR_DECOMPOSE(kx,ky,kz);;

                    FW(kx,ky,kz)    = FW_options{idx_vote_max}.FW(kx,ky,kz);
                    FW(kx,ky,kz+nz) = FW_options{idx_vote_max}.FW(kx,ky,kz+nz);

                end
            end

            % update success flag
            success_DECOMPOSITION_SELECTION(kz) = 1;

            % update average slice runtime
            average_slice_runtime = FattyRiot_elapsed_time(FattyRiot_start_time_DECOMPOSITION_SELECTION) / kz;
            disp( sprintf('FattyRiot: DECOMPOSITION_SELECTION average slice runtime (%d/%d) = %.2f seconds', kz, nz, average_slice_runtime) );

            if (FattyRiot_time_remaining(time_cushion_short+slice_runtime_amplificaton_factor*average_slice_runtime)<0),
                out_of_time_cushion_short = 1;
                error('FattyRiot: OUT OF TIME!');
            end

        end % end kz

    else
        disp('FattyRiot: Only one FW_option available for DECOMPOSITION_SELECTION');
        FW(:,:,[1:nz])    = FW_options{1}.FW(:,:,[1:nz]);
        FW(:,:,[1:nz]+nz) = FW_options{1}.FW(:,:,[1:nz]+nz);
    end
    
    disp( sprintf('FattyRiot: DECOMPOSITION_SELECTION completed in %.2f seconds', FattyRiot_elapsed_time(FattyRiot_start_time_DECOMPOSITION_SELECTION) ) );
catch
    
    if (out_of_time_cushion_short==1),
        disp( sprintf('FattyRiot: DECOMPOSITION_SELECTION ran out of time at slice (%d/%d)', kz, nz) );
    else
        warning('FattyRiot: Encountered unexpected problem running DECOMPOSITION_SELECTION');
        %rethrow(lasterror);
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 10. wrap-up
%    a. build rest of INFO structure if necessary
%    b. print success summary
%    c. print execution time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 10a. build rest of INFO structure if necessary
if nargout>1,
    try INFO.imDataParams_FattyRiot = imDataParams; catch end
    %try INFO.F = F_DECOMPOSE; catch end
    %try INFO.W = W_DECOMPOSE; catch end
    %try INFO.FM_HZ = FM_HZ; catch end
    %try INFO.R2 = R2; catch end
    %try INFO.ERROR = REMERROR_DECOMPOSE; catch end
    try INFO.nTEeven = nTEeven; catch end
    try INFO.idx_TEeven = idx_TEeven; catch end
    try INFO.dTEeven = dTEeven; catch end   
end

% 10b. print success summary
if (attempted_GLOBAL_maxR2_calibration==1),
    disp( sprintf('FattyRiot: Successfully completed R2* Calibration in %d of %d slices', sum(success_GLOBAL_maxR2_calibration), nz) );
end
if (attempted_QPBO_2D==1),
    disp( sprintf('FattyRiot: Successfully completed QPBO_2D in %d of %d slices', sum(success_QPBO_2D), nz) );
end
if (attempted_GC==1),
    disp( sprintf('FattyRiot: Successfully completed GC in %d of %d slices', sum(success_GC), nz) );
end
if (attempted_DECOMPOSITION_SELECTION==1),
    disp( sprintf('FattyRiot: Successfully completed DECOMPOSITION_SELECTION in %d of %d slices', sum(success_DECOMPOSITION_SELECTION), nz) );
end

% 10c. print execution time
disp( sprintf('FattyRiot: Completed in %.2f seconds', FattyRiot_elapsed_time(FattyRiot_start_time_original) ) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 11. internal functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 11a. FattyRiot_elapsed_time    
function elapsed_time_seconds = FattyRiot_elapsed_time(reference_time)
    elapsed_time_seconds = toc(reference_time);
end

% 11b. FattyRiot_time_remaining
function time_remaining_seconds = FattyRiot_time_remaining(time_cushion_seconds)
    allowed_time_seconds = 30 * 60; % 30 minutes per case
    time_remaining_seconds = allowed_time_seconds - time_cushion_seconds - FattyRiot_elapsed_time(FattyRiot_start_time_original);
end

% 11c. FattyRiot_DECOMPOSE
function DECOMPOSE_struct = FattyRiot_DECOMPOSE(FM_HZ,R2,slice_success_array),

    try
        FattyRiot_start_time_DECOMPOSE = tic;
          
        % common DECOMPOSE settings
        algoParams_DECOMPOSE = algoParams;
        algoParams_DECOMPOSE.species(1).frequency = [0.0];

        try 
            ampW = algoParams_DECOMPOSE.species(1).relAmps;
        catch
            ampW = 1.0;
        end

        deltaF = [0 ; gamma_Hz_per_Tesla/1e6*(algoParams_DECOMPOSE.species(2).frequency(:) - algoParams_DECOMPOSE.species(1).frequency(1))*(imDataParams.FieldStrength)];
        relAmps = algoParams_DECOMPOSE.species(2).relAmps;
        t = imDataParams.TE(:); % uses all echoes!
        relAmps = reshape(relAmps,1,[]);  
        
        B1 = zeros(nTE,2);
        B = zeros(nTE,2);
        
        % DECOMPOSE_struct
        DECOMPOSE_struct.FM_HZ = FM_HZ;
        DECOMPOSE_struct.R2 = R2;
        DECOMPOSE_struct.F_DECOMPOSE = zeros([nx ny nz]);
        DECOMPOSE_struct.W_DECOMPOSE = zeros([nx ny nz]);
        DECOMPOSE_struct.FW = zeros([nx ny nz*2]);
        DECOMPOSE_struct.REMERROR_DECOMPOSE = Inf*ones([nx ny nz]);
        DECOMPOSE_struct.success_DECOMPOSE = zeros(nz,1);

        for kz=1:nz,

            if (slice_success_array(kz)==1),
            
                % cover entire crop boundary box
                idx_x = idx_crop_x{kz};
                idx_y = idx_crop_y{kz};

                for kx=idx_x,
                    for ky=idx_y,

                        for n=1:nTE,
                            B1(n,:) = [ampW*exp(1i*2*pi*deltaF(1)*t(n)),sum(relAmps(:).*exp(1i*2*pi*deltaF(2:end)*t(n)))];
                        end
                        s = reshape( imDataParams.images(kx,ky,kz,1,:), [nTE 1]);
                        B(:,1) = B1(:,1).*exp(1i*2*pi*FM_HZ(kx,ky,kz)*t(:) - R2(kx,ky,kz)*t(:));
                        B(:,2) = B1(:,2).*exp(1i*2*pi*FM_HZ(kx,ky,kz)*t(:) - R2(kx,ky,kz)*t(:));
                        amps = B\s;
                        DECOMPOSE_struct.W_DECOMPOSE(kx,ky,kz) = amps(1);
                        DECOMPOSE_struct.F_DECOMPOSE(kx,ky,kz) = amps(2);
                        DECOMPOSE_struct.REMERROR_DECOMPOSE(kx,ky,kz) = norm(s - B*amps,'fro');

                    end % end ky
                end % end kx

                DECOMPOSE_struct.FW(:,:,kz)    = abs(DECOMPOSE_struct.F_DECOMPOSE(:,:,kz));
                DECOMPOSE_struct.FW(:,:,kz+nz) = abs(DECOMPOSE_struct.W_DECOMPOSE(:,:,kz));

                % update success flag
                DECOMPOSE_struct.success_DECOMPOSE(kz) = 1;

                % update average slice runtime
                average_slice_runtime = FattyRiot_elapsed_time(FattyRiot_start_time_DECOMPOSE) / kz;
                disp( sprintf('FattyRiot: DECOMPOSE average slice runtime (%d/%d) = %.2f seconds', kz, nz, average_slice_runtime) );

                if (FattyRiot_time_remaining(time_cushion_short+slice_runtime_amplificaton_factor*average_slice_runtime)<0),
                    out_of_time_cushion_short = 1;
                    error('FattyRiot: OUT OF TIME!');
                end

            else
                disp( sprintf('FattyRiot: DECOMPOSE skipping unsuccessful slice (%d/%d)', kz, nz) );
            end
            
        end % end kz

    catch
        
        if (out_of_time_cushion_short==1),
            disp( sprintf('FattyRiot: DECOMPOSE ran out of time at slice (%d/%d)', kz, nz) );
        else
            warning('FattyRiot: Encountered unexpected problem running DECOMPOSE');
            %rethrow(lasterror);
        end
    
    end

    disp( sprintf('FattyRiot: DECOMPOSE completed in %.2f seconds', FattyRiot_elapsed_time(FattyRiot_start_time_DECOMPOSE) ) );

end % end FattyRiot_DECOMPOSE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % end FattyRiot function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function: FattyRiot_coilCombine
%
% Description: combine multi-coil image sequences
%
% Based on: Walsh DO, Gmitro AF, Marcellin MW. Adaptive reconstruction of
% phased array MR imagery. Magn Reson Med 2000;43:682-690
% 
% Parameters:
% im1: the multi-coil images (size [nx,ny,nz,ncoils,nTE])
%
% Returns: 
% im2: the coil-combined images (size [nx,ny,nz,1,nTE])
%
% Author: Diego Hernando
% Date created: August 13, 2011
% Date last modified: December 8, 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function im2 = FattyRiot_coilCombine( im1 )

% Let's make the coil dimension the fourth one and the TE the third
im1 = permute(im1,[1 2 5 4 3]);

% Get image dimensions and set filter size
[sx,sy,N,C] = size(im1);
filtsize = 7;

% Initialize
im2 = zeros(sx,sy,1,1,N);
Rs = zeros(sx,sy,C,C);

% Get correlation matrices
for kc1=1:C
  for kc2=1:C
    for kn=1:N
      Rs(:,:,kc1,kc2) = Rs(:,:,kc1,kc2) + filter2(ones(filtsize),im1(:,:,kn,kc1).*conj(im1(:,:,kn,kc2)),'same');
    end
  end
end

% Compute and apply filter at each voxel
for kx=1:sx
  for ky=1:sy
% $$$     [U,S] = eig(squeeze(Rs(kx,ky,:,:)));
% $$$     s = diag(S);
% $$$     [maxval,maxind] = max(abs(s));
% $$$     myfilt = U(:,maxind);    
% $$$     im2(kx,ky,:) = myfilt'*reshape(squeeze(im1(kx,ky,:,:)).',[C N]);

    % Change suggested by Mark Bydder
    [U,S] = svd(squeeze(Rs(kx,ky,:,:)));
    myfilt = U(:,1); 
    im2(kx,ky,1,1,:) = myfilt'*reshape(squeeze(im1(kx,ky,:,:)).',[C N]);
  end
end

% In case the input data are single
if strcmp(class(im1),'single')
  im2 = single(im2);
end

end % end FattyRiot_coilCombine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function: FattyRiot_coilCombine3D
%
% Description: combine multi-coil image sequences
%
% Based on: Walsh DO, Gmitro AF, Marcellin MW. Adaptive reconstruction of
% phased array MR imagery. Magn Reson Med 2000;43:682-690
% 
% Parameters:
% im1: the multi-coil images (size [nx,ny,ncoils,nTE])
%
% Returns: 
% im2: the coil-combined images (size [nx,ny,nTE])
%
% Author: Diego Hernando
% Date created: August 13, 2011
% Date last modified: February 7, 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function im2 = FattyRiot_coilCombine3D( im1 )

[sx,sy,sz,C,N] = size(im1);

% Maintain the data type (e.g., single, double) of the input data
ims = zeros([sx,sy,sz,1,N],class(im1));
for kz=1:sz
  im2(:,:,kz,1,:) = FattyRiot_coilCombine(im1(:,:,kz,:,:));
end

end % end FattyRiot_coilCombine3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
