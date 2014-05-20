%% Function name: FattyRiot_fw_i3cm1i_3pluspoint_berglund_QPBO
%%
%% Description: Fat-water separation from three plus complex echoes with uniform
%%              echo time spacing, using a whole image optimization algorithm described in:
%%
%% Berglund J, Kullberg J. Three-dimensional water/fat separation and T2* estimation based on whole-image
%% optimization--application in breathhold liver imaging at 1.5 T. Magn Reson Med. 2012, Jun;67(6):1684-93.
%% 
%% Some properties:
%%   - Image-space
%%   - 2 species (water-fat)
%%   - Complex-fitting
%%   - Multi-peak fat (pre-calibrated)
%%   - R2* (given >3 echoes)
%%   - Independent water/fat phase
%%   - Requires 3 or more uniformly spaced echoes
%%
%% Input: structures imDataParams and algoParams
%%   - imDataParams.images: acquired images, array of size[nx,ny,nz,ncoils,,nTE]
%%   - imDataParams.TE: echo times (in seconds)
%%   - imDataParams.FieldStrength: (in Tesla)
%%   - imDataParams.voxelSize: (mm x mm x mm)
%%
%%   - algoParams.species(ii).name = name of species ii (string)
%%   - algoParams.species(ii).frequency = frequency shift in ppm of each peak within species ii
%%   - algoParams.species(ii).relAmps = relative amplitude of each peak within species ii
%%
%%   Example
%%   - algoParams.species(1).name = 'water';
%%   - algoParams.species(1).frequency = 4.70;
%%   - algoParams.species(1).relAmps = 1;
%%   - algoParams.species(2).name = 'fat';
%%   - algoParams.species(2).frequency = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29];
%%   - algoParams.species(2).relAmps = [88 642 58 62 58 6 39 10 37];
%% 
%%   - algoParams.decoupled_estimation = true; % flag for decoupled R2 estimation
%%   - algoParams.Fibonacci_search = true; % flag for Fibonacci search
%%   - algoParams.B0_smooth_in_stack_direction = false; % flag for B0 smooth in stack direction
%%   - algoParams.multigrid = true; % flag for multi-level resolution pyramid
%%   - algoParams.estimate_R2 = true; % flag to estimate R2star
%%   - algoParams.verbose = false; % flag for verbose status messages

%%   - algoParams.ICM_iterations = 2; % ICM iterations
%%   - algoParams.num_B0_labels = 100; % number of discretized B0 values
%%   - algoParams.mu = 10; % regularization parameter
%%   - algoParams.R2_stepsize = 1; % R2 stepsize in s^-1
%%   - algoParams.max_R2 = 120; % maximum R2 in s^-1
%%   - algoParams.max_label_change = 0.1; % 
%%   - algoParams.fine_R2_stepsize = 1.0; % Fine stepsize for discretization of R2(*) [sec-1] (used in decoupled fine-tuning step)
%%   - algoParams.coarse_R2_stepsize = 10.0; % Coarse stepsize for discretization of R2(*) [sec-1] (used for joint estimation step, value 0.0 --> Decoupled estimation only
%%   - algoParams.water_R2 = 0.0; % Water R2 [sec-1]
%%   - algoParams.fat_R2s = zeros(size(1,9)); % fat peak R2s [sec-1]
%%
%% Output: structure outParams
%%   - outParams.species(ii).name: name of the species (taken from algoParams)
%%   - outParams.species(ii).amps: estimated water/fat images, size [nx,ny,nz] 
%%   - outParams.fieldmap: field map (in degrees, size [nx,ny,nz]), fieldmap is NOT unwrapped
%%   - outParams.r2starmap: R2* map (in s^{-1}, size [nx,ny,nz])
%%   - outParams.status: QPBO status map (size [nx,ny,nz])
%%   - outParams.execution_time_seconds
%%
%% Author: E. Brian Welch
%% Date created: January 23, 2013
%% Date last modified: January 23, 2013

function outParams = FattyRiot_fw_i3cm1i_3pluspoint_berglund_QPBO( imDataParams, algoParams ),
time_start = tic;
outParams = [];

%% check validity of params, and set default algorithm parameters if not provided
[validParams, imDataParams, algoParams] = checkParamsAndSetDefaults_QPBO( imDataParams, algoParams );
if validParams==0
  disp(['Exiting -- data not processed']);
  return;
end

%% Build a signal input structure for DixonQPBO
DixonQPBO_input.images = imDataParams.images;
DixonQPBO_input.TE = imDataParams.TE;
DixonQPBO_input.field_strength = imDataParams.FieldStrength;
DixonQPBO_input.voxelSize = imDataParams.voxelSize;

% put water in species(1) and fat in species(2)
if strcmpi(algoParams.species(2).name,'water'),
    DixonQPBO_input.water_chemical_shift  = algoParams.species(2).frequency(:);
    DixonQPBO_input.fat_chemical_shifts   = algoParams.species(1).frequency(:);
    DixonQPBO_input.fat_amps              = algoParams.species(1).relAmps(:);
else
    DixonQPBO_input.water_chemical_shift  = algoParams.species(1).frequency(:);
    DixonQPBO_input.fat_chemical_shifts   = algoParams.species(2).frequency(:);
    DixonQPBO_input.fat_amps              = algoParams.species(2).relAmps(:);    
end
DixonQPBO_input.fat_R2s = algoParams.fat_R2s;

DixonQPBO_input.verbose = algoParams.verbose;
DixonQPBO_input.decoupled_estimation = algoParams.decoupled_estimation;
DixonQPBO_input.Fibonacci_search = algoParams.Fibonacci_search;
DixonQPBO_input.B0_smooth_in_stack_direction = algoParams.B0_smooth_in_stack_direction;
DixonQPBO_input.multigrid = algoParams.multigrid;

if ( (algoParams.estimate_R2) && (length(imDataParams.TE(:))>3) ),
    DixonQPBO_input.estimate_R2 = true;
else
    DixonQPBO_input.estimate_R2 = false;
end

DixonQPBO_input.num_B0_labels = algoParams.num_B0_labels;
DixonQPBO_input.ICM_iterations = algoParams.ICM_iterations;
DixonQPBO_input.mu = algoParams.mu;
DixonQPBO_input.R2_stepsize = algoParams.R2_stepsize;
DixonQPBO_input.max_R2 = algoParams.max_R2;
DixonQPBO_input.max_label_change = algoParams.max_label_change;
DixonQPBO_input.fine_R2_stepsize = algoParams.fine_R2_stepsize;
DixonQPBO_input.coarse_R2_stepsize = algoParams.coarse_R2_stepsize;
DixonQPBO_input.water_R2 = algoParams.water_R2;
DixonQPBO_input.InplaneOverThroughplaneVoxelsize = imDataParams.voxelSize(1) / imDataParams.voxelSize(3);
DixonQPBO_input.te_used = imDataParams.TE(1);
DixonQPBO_input.dte_used = imDataParams.TE(2) - imDataParams.TE(1);
DixonQPBO_input.use_num_echos = length(imDataParams.TE(:));

%% Call DixonQPBO
%DixonQPBO_output = DixonQPBO(DixonQPBO_input);

%% Call DixonQPBO via DixonApp
if exist('FattyRiot_DixonApp.m')~=2,
    %mfilename('fullpath')
    [pathstr,name,ext] = fileparts(mfilename('fullpath'));
    addpath( sprintf('%s/DixonApp/', pathstr) );
end
DixonQPBO_output = FattyRiot_DixonApp(DixonQPBO_input);

%% Build outParams structure
outParams.species(1).name = algoParams.species(1).name;
outParams.species(2).name = algoParams.species(2).name;
if strcmpi(outParams.species(2).name,'water'),
    outParams.species(2).amps = DixonQPBO_output.water;
    outParams.species(1).amps = DixonQPBO_output.fat;
else
    outParams.species(1).amps = DixonQPBO_output.water;
    outParams.species(2).amps = DixonQPBO_output.fat; 
end

outParams.fieldmap = DixonQPBO_output.fieldmap;
outParams.r2starmap = DixonQPBO_output.r2starmap;
outParams.status = DixonQPBO_output.status;

outParams.execution_time_seconds = toc(time_start);

end % end function fw_i3cm1i_3pluspoint_berglund_QPBO

function [validParams, imDataParams, algoParams] = checkParamsAndSetDefaults_QPBO( imDataParams, algoParams ),

    % sort echo times (just to be sure)
    imDataParams.TE = imDataParams.TE(:);
    [dummy,sorted_idx] = sort(imDataParams.TE);
    imDataParams.TE = imDataParams.TE(sorted_idx);
    imDataParams.images = imDataParams.images(:,:,:,:,sorted_idx);
    
    % check echo times to find largest set of evenly spaced echoes
    dTE_usec = round(1e6*diff(imDataParams.TE));
    dTE_candidates = unique([dTE_usec cumsum(dTE_usec) ]);
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
    
    if length(best_echoes)<3,
        error('Could not find at least 3 evenly spaced echoes');
    else
        imDataParams.TE = imDataParams.TE(best_echoes);
        imDataParams.images = imDataParams.images(:,:,:,:,best_echoes);
    end
    
    % check precession direction
    if ~isfield(imDataParams,'PrecessionIsClockwise'),
        imDataParams.PrecessionIsClockwise = 1;
    end
    
    if imDataParams.PrecessionIsClockwise ~= 1, 
        imDataParams.images = conj(imDataParams.images);
        imDataParams.PrecessionIsClockwise = 1;
    end   
    
    % check voxelsize
	if ~isfield(imDataParams,'voxelSize'),
        imDataParams.voxelSize = [1 1 1];
    end
    
    [nx,ny,nz,ncoils,nTE] = size(imDataParams.images);
    
    % if more than one channel, coil combine
    if ncoils > 1
        imDataParams.images = coilCombine3D(imDataParams.images);
    else
        imDataParams.images = reshape( imDataParams.images(:,:,:,1,:), [nx ny nz nTE]);
    end
    
    if algoParams.decoupled_estimation == 1,
        algoParams.decoupled_estimation = true;
    else
        algoParams.decoupled_estimation = false;
    end

    if algoParams.Fibonacci_search == 1,
        algoParams.Fibonacci_search = true;
    else
        algoParams.Fibonacci_search = false;
    end    

    if algoParams.B0_smooth_in_stack_direction == 1,
        algoParams.B0_smooth_in_stack_direction = true;
    else
        algoParams.B0_smooth_in_stack_direction = false;
    end    

    if algoParams.multigrid == 1,
        algoParams.multigrid = true;
    else
        algoParams.multigrid = false;
    end   

    if algoParams.estimate_R2 == 1,
        algoParams.estimate_R2 = true;
    else
        algoParams.estimate_R2 = false;
    end     
    
    if algoParams.verbose == 1,
        algoParams.verbose = true;
    else
        algoParams.verbose = false;
    end      
    
    % other default values
    if ~isfield(algoParams,'water_R2'),
        algoParams.water_R2  = 0.0;
    end
    
    if ~isfield(algoParams,'fat_R2s'),
        if strcmpi(algoParams.species(2).name,'water'),
            algoParams.fat_R2s  = zeros(1,length(algoParams.species(1).frequency));
        else
            algoParams.fat_R2s  = zeros(1,length(algoParams.species(2).frequency));
        end
    end
    
    validParams = 1;
    
end % end function checkParamsAndSetDefaults_QPBO

% Function: coilCombine
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

function im2 = coilCombine( im1 )

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

end % end coilCombine

% Function: coilCombine3D
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

function im2 = coilCombine3D( im1 )


[sx,sy,sz,C,N] = size(im1);

% Maintain the data type (e.g., single, double) of the input data
ims = zeros([sx,sy,sz,1,N],class(im1));
for kz=1:sz
  im2(:,:,kz,1,:) = coilCombine(im1(:,:,kz,:,:));
end

end % end coilCombine3D