% DixonApp  Binary MRF energy minimization on non-submodular graphs
%
%   [] = FattyRiot_DixonApp()
%

% v1.3 2013/03/16 welcheb

function output = FattyRiot_DixonApp(input)

[DixonApp_parent_folder, name, ext] = fileparts(mfilename('fullpath'));

DixonApp_vtk_folder = sprintf('%s/vtk', DixonApp_parent_folder);
DixonApp_vtk_prefix_input  = 'DixonApp_Input';
DixonApp_vtk_prefix_output = 'DixonApp_Output';

switch(computer)
    case {'PCWIN','PCWIN64'},
        %error( sprintf('DixonApp is not available for computer type: %s', computer) );
        path_to_DixonApp_executable = sprintf('%s/PC/DixonApp_PC.exe', DixonApp_parent_folder);
        system_command_string_prefix = sprintf('"%s"', path_to_DixonApp_executable);
	case {'GLNX86','GLNXA64'},
        %error( sprintf('DixonApp is not available for computer type: %s', computer) );
        path_to_DixonApp_executable = sprintf('%s/LINUX/DixonApp_LINUX.exe', DixonApp_parent_folder);
        path_to_DixonApp_libraries  = sprintf('%s/LINUX/lib', DixonApp_parent_folder);
        system_command_string_prefix = sprintf('export LD_LIBRARY_PATH="%s"; "%s"', path_to_DixonApp_libraries, path_to_DixonApp_executable);
    case {'MACI','MACI64'},
        path_to_DixonApp_executable = sprintf('%s/MAC/DixonApp_MAC.exe', DixonApp_parent_folder);
        path_to_DixonApp_libraries  = sprintf('%s/MAC/lib', DixonApp_parent_folder);
        system_command_string_prefix = sprintf('export DYLD_LIBRARY_PATH="%s"; "%s"', path_to_DixonApp_libraries, path_to_DixonApp_executable);
    case {'SOL64'},
        error( sprintf('DixonApp is not available for computer type: %s', computer) );
    otherwise
        error( sprintf('DixonApp is not available for computer type: %s', computer) );
end

%% Clean vtk folder
delete(sprintf('%s/*.vtk', DixonApp_vtk_folder));

%% Write vtk files
%vtk_max = 2^15-1;
vtk_max = 2^14-1;
[nx ny nz ne] = size(input.images);
images_max = max([ max(abs(real(input.images(:)))) max(abs(imag(input.images(:)))) ]);
for e=1:ne,
    % VTK supports single and double, DixonApp.app appears to want int16
    image_re = int16( real(input.images(:,:,:,e)) * (vtk_max/images_max) );
    image_im = int16( imag(input.images(:,:,:,e)) * (vtk_max/images_max) );
    vtk_filename_re = sprintf('%s/%s_echo%04d_re.vtk', DixonApp_vtk_folder, DixonApp_vtk_prefix_input, e);
    vtk_filename_im = sprintf('%s/%s_echo%04d_im.vtk', DixonApp_vtk_folder, DixonApp_vtk_prefix_input, e);
    writeVTK(image_re,vtk_filename_re);
    writeVTK(image_im,vtk_filename_im);
end

%% Create system_command_string
system_command_string = system_command_string_prefix;

% strings (3)
system_command_string = sprintf('%s "%s"', system_command_string, DixonApp_vtk_folder);
system_command_string = sprintf('%s "%s"', system_command_string, DixonApp_vtk_prefix_input);
system_command_string = sprintf('%s "%s"', system_command_string, DixonApp_vtk_prefix_output);

% flags (6)
system_command_string = sprintf('%s "%d"', system_command_string, input.verbose);
system_command_string = sprintf('%s "%d"', system_command_string, input.estimate_R2);
system_command_string = sprintf('%s "%d"', system_command_string, input.decoupled_estimation);
system_command_string = sprintf('%s "%d"', system_command_string, input.Fibonacci_search);
system_command_string = sprintf('%s "%d"', system_command_string, input.multigrid);
system_command_string = sprintf('%s "%d"', system_command_string, input.B0_smooth_in_stack_direction);

% floats (11)
system_command_string = sprintf('%s "%.5e"', system_command_string, input.field_strength);
system_command_string = sprintf('%s "%.5e"', system_command_string, input.te_used);
system_command_string = sprintf('%s "%.5e"', system_command_string, input.dte_used);
system_command_string = sprintf('%s "%.5e"', system_command_string, input.water_chemical_shift);
system_command_string = sprintf('%s "%.5e"', system_command_string, input.mu);
system_command_string = sprintf('%s "%.5e"', system_command_string, input.R2_stepsize);
system_command_string = sprintf('%s "%.5e"', system_command_string, input.max_R2);
system_command_string = sprintf('%s "%.5e"', system_command_string, input.max_label_change);
system_command_string = sprintf('%s "%.5e"', system_command_string, input.InplaneOverThroughplaneVoxelsize);
system_command_string = sprintf('%s "%.5e"', system_command_string, input.fine_R2_stepsize);
system_command_string = sprintf('%s "%.5e"', system_command_string, input.coarse_R2_stepsize);
system_command_string = sprintf('%s "%.5e"', system_command_string, input.water_R2);

% integers (4)
system_command_string = sprintf('%s "%d"', system_command_string, input.num_B0_labels);
system_command_string = sprintf('%s "%d"', system_command_string, input.ICM_iterations);
system_command_string = sprintf('%s "%d"', system_command_string, input.use_num_echos);
num_fat_peaks = length(input.fat_chemical_shifts);
system_command_string = sprintf('%s "%d"', system_command_string, num_fat_peaks );

% fat_chemical_shifts (num_fat_peaks floats)
for idx=1:num_fat_peaks,
    system_command_string = sprintf('%s "%.5e"', system_command_string, input.water_chemical_shift + input.fat_chemical_shifts(idx) ); % ppm for DixonApp is difference between water and peak
end

% fat_amps (num_fat_peaks floats)
for idx=1:num_fat_peaks,
    system_command_string = sprintf('%s "%.5e"', system_command_string, input.fat_amps(idx) );
end

% fat_R2s (num_fat_peaks floats)
for idx=1:num_fat_peaks,
    system_command_string = sprintf('%s "%.5e"', system_command_string, input.fat_R2s(idx) );
end

%disp(system_command_string);

%% Call DixonApp executable
[status,result] = system(system_command_string);

if input.verbose,
    disp( sprintf('status = %d', status) );
    disp(result);
end

%% Check status returned by exectuable and throw error if necessary
if (status~=0),
    error_str = sprintf('FattyRiot: ERROR: FattyRiot_DixonApp.m - returned status is %d - check DixonApp executable, dynamic FLTK libraries, and X11 display for LINUX', status);
    disp(error_str);
    error(error_str);
end

%% Read vtk files
wat_re = double( readVTK( sprintf('%s/%s22_wat_re.vtk', DixonApp_vtk_folder, DixonApp_vtk_prefix_output) ) ) * (images_max/vtk_max);
wat_im = double( readVTK( sprintf('%s/%s24_wat_im.vtk', DixonApp_vtk_folder, DixonApp_vtk_prefix_output) ) ) * (images_max/vtk_max);
output.water = complex(wat_re, wat_im);

fat_re = double( readVTK( sprintf('%s/%s23_fat_re.vtk', DixonApp_vtk_folder, DixonApp_vtk_prefix_output) ) ) * (images_max/vtk_max);
fat_im = double( readVTK( sprintf('%s/%s25_fat_im.vtk', DixonApp_vtk_folder, DixonApp_vtk_prefix_output) ) ) * (images_max/vtk_max);
output.fat = complex(fat_re, fat_im);

fieldmap_degrees = double( readVTK( sprintf('%s/%s6_fieldmap.vtk', DixonApp_vtk_folder, DixonApp_vtk_prefix_output) ) );
output.fieldmap = fieldmap_degrees;

if input.estimate_R2,
    output.r2starmap = double( readVTK( sprintf('%s/%s7_R2star.vtk', DixonApp_vtk_folder, DixonApp_vtk_prefix_output) ) );
else
    output.r2starmap = zeros(nx,ny,nz);
end

output.status = double( readVTK( sprintf('%s/%s21_status.vtk', DixonApp_vtk_folder, DixonApp_vtk_prefix_output) ) );

%% Clean vtk folder
delete(sprintf('%s/*.vtk', DixonApp_vtk_folder));

end

function [] = writeVTK(vol,vtkfile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: writeVTK(vol,vtkfile)
%
%   vol:     The 3D matrix to be saved to file
%   vtkfile: The output filename (string)
%   notes:   Only writes binary STRUCTURED_POINTS
%  
% Erik Vidholm 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dimensions
volinfo = whos('vol');

sz = volinfo.size;
X = sz(1); Y = sz(2); Z = 1;
if( length(sz) == 3 )
  Z = sz(3);
end

% open file (OBS! big endian format)
fid = fopen(vtkfile,'w','b');

% write header
fprintf(fid, '%s\n', '# vtk DataFile Version 3.0');
fprintf(fid, '%s\n', 'created by writeVTK (Matlab implementation by Erik Vidholm)');
fprintf(fid, '%s\n', 'BINARY');  
fprintf(fid, '%s\n', 'DATASET STRUCTURED_POINTS');  
fprintf(fid, '%s%d%c%d%c%d\n', 'DIMENSIONS ', X, ' ', Y, ' ', Z);
fprintf(fid, '%s%f%c%f%c%f\n', 'ORIGIN ', 0.0, ' ', 0.0, ' ', 0.0); 
fprintf(fid, '%s%f%c%f%c%f\n', 'SPACING ', 1.0, ' ', 1.0, ' ', 1.0); 
fprintf(fid, '%s%d\n', 'POINT_DATA ', X*Y*Z);

tp = volinfo.class;
if( strcmp(tp, 'uint8') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data unsigned_char');
elseif( strcmp(tp, 'int8') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data char');
elseif( strcmp(tp, 'uint16') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data unsigned_short');
elseif( strcmp(tp, 'int16') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data short');
elseif( strcmp(tp, 'uint32') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data unsigned_int');
elseif( strcmp(tp, 'int32') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data int');
elseif( strcmp(tp, 'single') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data float');
elseif( strcmp(tp, 'double') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data double');
end

fprintf(fid, '%s\n', 'LOOKUP_TABLE default');

% write data as binary
fwrite(fid,vol,tp);

% close file
fclose(fid);

end

function V = readVTK(vtkfile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: V = readVTK(vtkfile)
%
%   V:       The matrix to be stored
%   vtkfile: The filename
%   notes:   Only reads binary STRUCTURED_POINTS
%
% Erik Vidholm 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = 0;

% open file (OBS! big endian format)
fid = fopen(vtkfile,'r','b');

if( fid == -1 )
  return
end

fgetl(fid); % # vtk DataFile Version x.x
fgetl(fid); % comments
fgetl(fid); % BINARY
fgetl(fid); % DATASET STRUCTURED_POINTS

s = fgetl(fid); % DIMENSIONS NX NY NZ
sz = sscanf(s, '%*s%d%d%d').';

fgetl(fid); % ORIGIN OX OY OZ
fgetl(fid); % SPACING SX SY SZ
fgetl(fid); % POINT_DATA NXNYNZ

s = fgetl(fid); % SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
svstr = sscanf(s, '%s', 1);
dtstr = sscanf(s, '%*s%*s%s');

if( strcmp(svstr,'SCALARS') > 0 )
  fgetl(fid); % the lookup table
  if( strcmp(dtstr,'unsigned_char') > 0 ) 
    % read data
    V = fread(fid,prod(sz),'*uint8');
    V = reshape(V,sz);
  elseif( strcmp(dtstr,'char') > 0 )
    % read data
    V = fread(fid,prod(sz),'*int8');
    V = reshape(V,sz);
  elseif( strcmp(dtstr,'unsigned_short') > 0 )
    % read data
    V = fread(fid,prod(sz),'*uint16');
    V = reshape(V,sz);
  elseif( strcmp(dtstr,'short') > 0 )
    % read data
    V = fread(fid,prod(sz),'*int16');
    V = reshape(V,sz);
  elseif( strcmp(dtstr,'unsigned_int') > 0 )
    % read data
    V = fread(fid,prod(sz),'*uint32');
    V = reshape(V,sz);
  elseif( strcmp(dtstr,'int') > 0 )
    % read data
    V = fread(fid,prod(sz),'*int32');
    V = reshape(V,sz);
  elseif( strcmp(dtstr,'float') > 0 )
    % read data
    V = fread(fid,prod(sz),'*single');
    V = reshape(V,sz);
  elseif( strcmp(dtstr,'double') > 0 )
    % read data
    V = fread(fid,prod(sz),'*double');
    V = reshape(V,sz);
  end
  
elseif( strcmp(svstr,'VECTORS') > 0 )
  if( strcmp(dtstr,'float') > 0 ) 
    % read data
    V = fread(fid,3*prod(sz),'*single');
    V = reshape(V,[3 sz]);
    V = permute(V,[2 3 4 1]);
  end
end

fclose(fid);

end
