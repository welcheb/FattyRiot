function [W,F] = FattyRiot_Berglund2010(S,DataPar,MethodPar)
%ASR 
%Analytical water/fat separation with a Safest-first Region-growing
%scheme. Three-point Dixon method described in detail in:
%Berglund J, Johansson L, Ahlström H, Kullberg J. 'Three-point Dixon method 
%enables whole-body water and fat imaging of obese subjects'.
%Magn Reson Med 2010; 63:1659-68.
%
%This code is written by Johan Berglund, Uppsala University, Sweden.

% gyro = 42.6;                                                             % gyromagnetic ratio for hydrogen [MHz/T]
gyro = 42.58;                                                              % gyromagnetic ratio for hydrogen [MHz/T]
omega = 2*pi*gyro*DataPar.B0*(MethodPar.CS_F-MethodPar.CS_W);              % resonance frequency vector [in radians]
alfa=MethodPar.CS_F_weights/sum(MethodPar.CS_F_weights);                   % normalize weights

%t=DataPar.TE1+DataPar.dTE*(0:2);                                          % echo time vector
t=[DataPar.TE1 DataPar.TE2 DataPar.TE3];                                   % echo time vector

a=exp(complex(0,t'*omega))*alfa';                                          % fat offset vector
mw=MagnitudeWeight(S);                                                     % get magnitude weight image
[b1 b2] = getEstimates(S,a);                                               % get two error phasor candidates in each voxel

mw = single(mw);
b1 = single(b1);
b2 = single(b2);

b = RegionGrowing(a,b1,b2,MethodPar.c1,MethodPar.c2,mw,DataPar.voxelsize); % choose one candidate in each voxel through multi-seeded region growing

S(:,:,:,2)=S(:,:,:,2)./b;                                                  % linearize signal equation
S(:,:,:,3)=S(:,:,:,3)./b.^2;

[W,F]=SolveLS(S,a);                                                        % find least squares estimates for water and fat

clear mw b1 b2 b;

end

function [b1 b2] = getEstimates(S, a)
%GETESTIMATES Calculates two alternative estimates of b equivalent to eq 4

b0 = S(:,:,:,2)*(a(1)-a(3))./(2*S(:,:,:,1)*(a(2)-a(3)));
db = sqrt(b0.*b0-S(:,:,:,3)*(a(1)-a(2))./(S(:,:,:,1)*(a(2)-a(3))));
b1=b0+db;
b2=b0-db;
b1(isnan(b1))=1;    %make sure all voxels are defined
b2(isnan(b2))=-1;   %make sure all voxels are defined

clear b0 db;

end

function mw = MagnitudeWeight(S)
%MAGNITUDEWEIGHT Calculates magnitude weight images as defined in eq 6

m=sum(abs(S),4);                                                           % Magnitude according to eq 5
m=m/max(max(max(m)));                                                      % Normalize magnitude prior to Otsu's method
mt = graythresh(m);                                                        % Find threshold using Otsu's method
mf = mean(mean(m(m>mt)));                                                  % Mean value of foreground voxels
mw=1.0./(1+exp(3*(m-mt)/(mt-mf)));                                         % Calculate magnitude weight according to eq 6

clear mf mt m;

end

function b = RegionGrowing(a,b1,b2,c1,c2,mw,voxelsize)
%REGIONGROWING Multi-seeded safest-first region growing

res=true(size(b1));                                                        % true=b1, false=b2
determined=false(size(b1));                                                % keeps track on determined voxels

%Find seed points
fg=mw>c1;                                                                  % foreground voxels have magnitude weight >c1
%|log(R1)| , R1 calculated equivalently to eq. 7
r1=abs(log((a(1)*(a(2)-a(3))*b2-a(3)*(a(1)-a(2))*b1)./((a(1)-a(2))*b1-(a(2)-a(3))*b2)));
%|log(R2)| , R2 calculated equivalently to eq. 8
r2=abs(log((a(1)*(a(2)-a(3))*b1-a(3)*(a(1)-a(2))*b2)./((a(1)-a(2))*b2-(a(2)-a(3))*b1)));
determined(fg)=or(r1(fg)<c2,r2(fg)<c2);                                    % voxels with mw>c1 and |logR|<c2 are seeds
res(determined)=r1(determined)<r2(determined);                             % pick solution with smallest |logR| in seeds

%Region growing
status = zeros(size(b1), 'uint8');                                         % create status image; 0=undetermined, 1=b1, 2=b2
status(and(determined,res))=1;
status(and(determined,~res))=2;

b_angle = FattyRiot_RG(status, angle(b1), angle(b2), mw, voxelsize);       % Call subroutine RG written in c++
b = exp(complex(0,b_angle));

clear res determined fg r1 r2 status b_angle;

end

function [W,F] = SolveLS(S,a)
%SOLVELS Finds least squares solution of water and fat

A=[ones(3,1,'single') a];                                                  % Linear model matrix
Ainv = pinv(A);                                                            % Pseudoinverse of A
nx=size(S,1); ny=size(S,2); nz=size(S,3);

W=reshape(reshape(S,[],3)*Ainv(1,:).',[nx ny nz]);                         % W according to eq 12
F=reshape(reshape(S,[],3)*Ainv(2,:).',[nx ny nz]);                         % F according to eq 12

end


