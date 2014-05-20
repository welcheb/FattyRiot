%% Function name: fw_i3cm0i_3point_berglund
%%
%% Description: Fat-water separation from three complex echoes with uniform
%%              echo time spacing, using a multi-seeded region growing
%%              scheme, as described in:
%%
%% Berglund J, Johansson L, Ahlström H, Kullberg J. Three-point Dixon method enables whole-body 
%% water and fat imaging of obese subjects. Magn Reson Med. 2010, Jun;63(6):1659-1668.
%% 
%% Some properties:
%%   - Image-space
%%   - 2 species (water-fat)
%%   - Complex-fitting
%%   - Multi-peak fat (pre-calibrated)
%%   - No R2*
%%   - Independent water/fat phase
%%   - Requires 3 uniformly spaced echoes
%%
%% Input: structures imDataParams and algoParams
%%   - imDataParams.images: acquired images, array of size[nx,ny,nz,ncoils,3]
%%   - imDataParams.TE: echo times (in seconds)
%%   - imDataParams.FieldStrength: (in Tesla)
%%   - imDataParams.voxelSize: (mm x mm x mmm)
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
%%   - algoParams.c1 = 0.75; % Threshold on magnitude weight for seed points
%%   - algoParams.c2 = 0.25; % Threshold on |log(W/F)| for seed points
%%
%% Output: structure outParams
%%   - outParams.water: estimated water image, size [nx,ny,nz]
%%   - outParams.fat: estimated fat image, size [nx,ny,nz]
%%   - outParams.fieldmap: field map (in Hz, size [nx,ny,nz])
%%
%%
%% Author: Johan Berglund
%% Date created: November 16, 2011
%% Date last modified: December 20, 2011

function outParams = fw_i3cm0i_3point_berglund( imDataParams, algoParams )

%% Check validity of params, and set default algorithm parameters if not provided
[validParams,imDataParams,algoParams] = checkParamsAndSetDefaults_MultiSeedRegionGrowing( imDataParams,algoParams );
if validParams==0
  disp(['Exiting -- data not processed']);
  outParams = [];
  return;
end

%% create linear model matrix
%gyromagnetic ratio for hydrogen [MHz/T]
gyro = 42.6;
%resonance frequency vector [in radians]:
if (imDataParams.PrecessionIsClockwise>0) %if not clockwise, (most) fat frequencies will be negative
    omega_p = -2*pi*gyro*imDataParams.FieldStrength*(algoParams.species(1).frequency-algoParams.species(2).frequency);
else %if clockwise, (most) fat frequencies will be positive
    omega_p = +2*pi*gyro*imDataParams.FieldStrength*(algoParams.species(1).frequency-algoParams.species(2).frequency);
end

a=exp(complex(0,imDataParams.TE(1:3)'*omega_p))*algoParams.species(2).relAmps';
A=[ones(3,1) a];

%% DH*: Combine multiple coils (replaced dummyCoilCombine)
if size(imDataParams.images,4) > 1
  S = (permute(coilCombine3D( imDataParams.images(:,:,:,:,1:3) ),[1 2 3 5 4]));
else
  S = (permute( imDataParams.images(:,:,:,:,1:3),[1 2 3 5 4]));
end

%% get two field map phasor candidates in each voxel
[bA bB] = getPhasorCandidates(S,A);

%% choose one candidate in each voxel through multi-seeded region growing:
b = regionGrow(A,bA,bB,algoParams.c1,algoParams.c2,getMagnitudeWeight(S),imDataParams.voxelSize);

%% remove field map "phase error":
S(:,:,:,2)=S(:,:,:,2)./b;
S(:,:,:,3)=S(:,:,:,3)./b.^2;

%% find least squares estimates for water and fat
[W,F]=SolveLS(S,A);

%% Put results in outParams structure
outParams.species(1).name = 'water';
outParams.species(2).name = 'fat';
outParams.species(1).amps = W;
outParams.species(2).amps = F;
if (imDataParams.PrecessionIsClockwise<=0)
  outParams.fieldmap = +angle(b)/(2*pi*diff(imDataParams.TE(1:2)));
else
  outParams.fieldmap = -angle(b)/(2*pi*diff(imDataParams.TE(1:2)));
end

  
  