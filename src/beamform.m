clc
clear all
close all
%% Execution
% Beamforming and plotting
[img1, x1, z1] = beamform(1, 'rect', false, 60); % Fixed Focus

[img2, x2, z2] = beamform(1, 'rect', true, 60); % Dynamic Focus - Rectangular

[img3, x3, z3] = beamform(1, 'hann', true, 60); % Dynamic Focus - Hann

[img4, x4, z4] = beamform(1, 'hamming', true, 60); % Dynamic Focus - Hamming

[img5, x5, z5] = beamform(1, 'blackmann', true, 60); % Dynamic Focus - Blackmann

% Calculation of FWHM
[depth1, fwh1] = fwhm_calc(img1, x1, z1); % Fixed Focus
[depth2, fwh2] = fwhm_calc(img2, x2, z2); % Dynamic Focus - Rectangular
[depth3, fwh3] = fwhm_calc(img3, x3, z3); % Dynamic Focus - Hann
[depth4, fwh4] = fwhm_calc(img4, x4, z4); % Dynamic Focus - Hamming
[depth5, fwh5] = fwhm_calc(img5, x5, z5); % Dynamic Focus - Blackmann

% Plotting FWHM
figure;
hold on
plot(depth1,fwh1,"Color","k","LineWidth",3)
plot(depth2,fwh2,"Color",[0.25, 0.25, 0.25],"LineStyle","--","LineWidth",3)
plot(depth3,fwh3,"Color",[0.50, 0.50, 0.50],"LineWidth",3)
plot(depth4,fwh4,"Color",[0.37, 0.37, 0.37],"LineStyle",":","LineWidth",3)
plot(depth5,fwh5,"Color",[0.65, 0.65, 0.65],"LineStyle","-.","LineWidth",3)
x = xlabel('Axial distance [mm]');
y = ylabel("FWHM [mm]");
lgd = legend("Fixed-Rectangular","Dynamic-Rectangular","Dynamic-Hann","Dynamic-Hamming","Dynamic-Blackman","Location","northwest");
lgd.FontSize = 16;
x.FontSize = 16;
y.FontSize = 16;
grid

% Beamforming of other data sets
beamform(2, 'hann', true, 60); % Phantom measured

beamform(3, 'hann', true, 60); % In-vivo

%% Functions
function [img, x, z] = beamform(id_data, id_apod, id_focus, db_level)
%beamform    Create the image resulted by beamforming ultrasound data of a convex array
%
%USAGE
%   [img, x, z] = beamform(id_data)
%   [img, x, z] = beamform(id_data, id_apod, id_focus, db_level)
%
%INPUT ARGUMENTS
% id_data  : integer specifying selection of data to beamform
%            1  = simulated data
%            2  = phantom data
%            3  = in-vivo data
% id_apod  : string specifying window type for apodization (default, type = 'hamming')
%            'rect'      = rectangular window (no apodization)
%            'hann'      = hann window
%            'hamming'   = hamming window 
%            'blackmann' = blackmann window
% id_focus : boolean specifying approach to focusing (default, id_focus = true)
%            true  = dynamic receive focusing
%            false = fixed receive focusing
% db_level : integer specifying db level for compression (default, db_level = 60)
%
%OUTPUT ARGUMENTS
%   img : resulting image [1024 x 1024]
%     x : vector containing lateral coordinates [1024 x 1]
%     z : vector containing axial coordinates [1024 x 1]
%
%   beamform(...) plots the image in a new figure.
%% CHECK INPUT ARGUMENTS
% Check for proper input arguments
if nargin < 1 || nargin > 4
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 4 || isempty(db_level);  db_level  = 60; end
if nargin < 3 || isempty(id_focus);  id_focus  = true; end
if nargin < 2 || isempty(id_apod);   id_apod   = 'hamming';  end

%% Initializing requirements for data

switch id_data % Determining folder path for obtaining data
    case 1 % Simulated data for convex array
        file_path = "phantom_matrix_sim_files\duplex\B_mode\seq_0001\";
    case 2 % Measured data for convex array and matrix phantom
        file_path = "phantom_matrix_files\duplex\B_mode\seq_0001\";
    case 3 % Measured in-vivo data for convex array
        file_path = "convex_array_files\B_mode\seq_0001\";
    otherwise
        error('Id "%d" is not supported',id_data);
end

% Loading first file to initialize some general parameters with its size
tag = file_path + "elem_data_em0" + num2str(1,"%.3d") + ".mat"; % General format for files
load(tag)

%% Setting general parameters

n = size(samples,1); % Number of data points recorded by each transducer per emission

m = 129; % Total amount of emissions
c = 1491; % Speed of sound (m/s)
fs = 17.5E6; % Sampling frequency (Hz)

n_elements = 192; % Number of elements (transducers)
width = 0.30E-3; % Width of one element (m)
kerf = 0.03E-3; % Kerf (gap) between elements (m)
radius = 60.3E-3; % Convex radius (m)

%% Defining transducer geometry
% Calculating the element pitch
pitch = width + kerf;
% Obtain the total length of the convex element arc
total_arc_length = n_elements * pitch;
% Calculate the total angle which corresponds to the ratio of the length of
% the element arc and the convex radius
total_angle = total_arc_length / radius;

% Initializing variables to calculate the relative position of each element
theta_i = zeros(n_elements,1); % Variable that stores the angles of each element
r_i = zeros(n_elements,2); % Variable that stores the positions (x, z) of each element

% Cycle to calculate positions of each element
for j = 1:1:n_elements
    theta_i(j) = -total_angle/2 + (j-0.5) * (total_angle/192); % Calculate the angle corresponding to the i'th element
    r_i(j,1) = radius * sin(theta_i(j)); % Obtain x coordinate for the i'th element
    r_i(j,2) = radius * cos(theta_i(j)); % Obtain z coordinate for the i'th element
end

%% Focus preparation
% In this section it is decided if the approach is going to be the one
% which has delay oriented calculation or the one which finds each point
% individually

if id_focus % Dynamic receive focusing
    r_rf = (1:1:n) / fs * c + radius; % Creation of vector for different radius of focal points
else % Fixed receive focusing
    r_rf = zeros(1,n) + radius + 40E3; % Creation of vector for fixed focal point (40 mm)
end

%% Apodization preparation
% Obtaining the apodization window for the fixed active probe elements
n_w = 64; % Active probe elements
w = win_apod(n_w, id_apod); % Calculation of window

%% Beamforming

% Initializing variables
s_t = zeros(n,m); % Matrix with beamformed data
r_f = [0,0]; % Vector with position of focal point
r_c = [0,0]; % Vector with position of beam reference point

for j = 1:1:m % Iterating through number of emissions
    % Loading j'th file
    tag = file_path + "elem_data_em0" + num2str(j,"%.3d") + ".mat"; % General format for files
    load(tag)
    % Standarizing format, particullarly for cases when data is in uint16
    samples = double(samples);
    % Center data around cero, baseline removal
    baseline = mean(samples);
    samples = samples - baseline;

    % Obtaining position of beam reference point according to the emmision
    theta_rf = mean(theta_i(j:j+63)); % Obtain mean angle of the 64 emitting elements
    r_c(1) = radius * sin(theta_rf); % Obtain x coordinate for the j'th beam reference
    r_c(2) = radius * cos(theta_rf); % Obtain z coordinate for the j'th beam reference
    
    for k = 1:1:n % Iterating through all the focal points previously defined
        % Obtaining position of focal point
        r_f(1) = r_rf(k) * sin(theta_rf); % Obtain x coordinate for the k'th focal point
        r_f(2) = r_rf(k) * cos(theta_rf); % Obtain z coordinate for the k'th focal point

        for i = j:j+n_w-1 % Iterating through active elements, window of 64 elements per emission
            
            % Calculating time delay for current element
            d_fi = norm(r_f - r_i(i,:)); % Calculating the distance from the focal point to the i'th element
            % The next is only done when fixed focusing, due to making the
            % delay relative to another element
            d_fc = norm(r_f - r_c) * ~id_focus; % Calculating the distance from the focal point to the beam reference
            tau_i = (d_fi - d_fc) / c; % Converting distance to time
            sample_delay = tau_i * fs + k * ~id_focus; % Converting time to sample index
            % The term "k" is added in each iteration to effectively cycle
            % through the different samples when fixed receiving focal
            % point is selected
                
            % Due to the limitation in resolution by to the sampling
            % frequency a linear interpolation approach was made to
            % compensate and obtain a better precision
            
            % Start of the interpolation
            idx_floor = floor(sample_delay); % Truncate to lower sample index
            idx_ceil = ceil(sample_delay); % Truncate to upper sample index
            frac = sample_delay - idx_floor; % Obtain decimal value of the sample index

            if idx_ceil <= size(samples,1) % Sanity check to always obtain real values
                if idx_floor == idx_ceil % In this case there is no fraction, so it landed on an integer index
                    temp = samples(idx_floor, i); % Obtain value of the signal corresponding to that sample
                else % If there exists a fraction, it means the actual point is located in between samples
                    % The equation of a straight line is followed to
                    % make an approximation for the actual value
                    % corresponding to that particular time
                    temp = (1-frac)*samples(idx_floor,i) + frac*samples(idx_ceil,i); % Obtain value of the signal corresponding to that sample
                end
                s_t(k, j) = s_t(k, j) + temp * w(i-j+1); % Multiplying data by apodization window and adding the line of the signal to the whole sum
            end
        end
    end
end

%% Envelope and compression of signal
env_data = zeros(size(s_t)); % Initializing variable for envelope calculation

for j = 1:1:m % Calculating for each line of signal
    env_data(:,j) = abs(hilbert(s_t(:,j))); % Obtaining the envelope through hilbert transform
end

% Realizing the compression of the signal
norm_RFdata = env_data / max(env_data(:)); % Normalizing it by its max value
floor_linear = 10^(-db_level / 20); % Setting the floor for compression
compressed_signal = max(norm_RFdata, floor_linear); % Compressing the signal
% Converting it to logarithmic scale
compressed_dB = 20 * log10(compressed_signal);

%% Interpolation
% Obtaining the distance vector in radial shape according to the speed of
% sound, the sampling frequency and the amount of samples
r = radius + (1:size(compressed_dB,1)) / fs * c / 2;

% Obtaining the vector with the angles of each emission beam, which
% corresponds to the angles of the calculated signal lines
for j = 1:1:129 % Calculating for each line of signal
    theta_rf(j) = mean(theta_i(j:j+63)); % Obtain mean angle of the 64 emitting elements
end

[Theta, R] = meshgrid(theta_rf, r); % Obtaining the grid for the interpolation in polar coordinates

% Defining the cartesian coordinates
x_l = r(end) * sin(theta_rf(1)); % Max lateral distance
z_i = radius * cos(theta_rf(1)); % Min axial distance (Since first element of transducer)
z_o = r(end); % Max axial distance

x = linspace(-x_l, x_l, 1024); % Vector containing lateral coordinates
z = linspace(z_i, z_o, 1024); % Vector containing axial coordinates

[X, Z] = meshgrid(x,z); % Obtaining the grid for the interpolation in cartesian coordinates

[Tho, Ro] = cart2pol(X,Z); % Transforming cartesian coordinates to polar

img = interp2(Theta, R, compressed_dB, -Tho + pi/2, Ro); % Performing interpolation

img(isnan(img)) = -60; % Setting pixels without value to the floor level

x = x * 1000; % Adjusting lateral vector to milimeter scale
z = (z - z_i) * 1000; % Adjusting axial vector to milimeter scale and adjusting the values so it starts at the height of the first element

%% Plotting Image

figure; % Plotting the image
imagesc(x, z, img);
axis equal; % Adjusting the scale
axis image;
xlabel('Lateral distance [mm]');
ylabel('Axial distance [mm]');
colormap(gray(128))
colorbar; % Displaying the colorbar
end

function w = win_apod(N, type)
%win_apod    Create apodization window
%
%USAGE
%       w = win_apod(N, type)
%
%INPUT ARGUMENTS
%    N : integer specifying selection size of window
% type : string specifying window type for apodization
%            'rect'      = rectangular window (no apodization)
%            'hann'      = hann window
%            'hamming'   = hamming window 
%            'blackmann' = blackmann window
%
%OUTPUT ARGUMENTS
%    w : resulting apodization window [N x 1]
%
%% Modify parameter according to window type
switch lower(type)
    case 'rect' % No apodization applied
        a = 1; 
        b = 0; 
        c = 0;
    case 'hann'
        a = 0.5; 
        b = 0.5; 
        c = 0;
    case 'hamming'
        a = 0.54; 
        b = 0.46; 
        c = 0;
    case 'blackmann'
        a = 0.42; 
        b = 0.5; 
        c = 0.08;
    otherwise
        error('Window type "%s" is not supported',lower(type));
end

% Sample indices
n = (0:N-1)';

% Generalized form of cosine windows
w = a - b * cos(2*pi*n/(N-1)) + c * cos(4*pi*n/(N-1));
end

function [depths, fwhm] = fwhm_calc(img, x, z)
%fwhm_calc    Calculate Full Width at Half Maximum of wire targets
%
%USAGE
%       [depths, fwhm] = fwhm_calc(img, x, z)
%
%INPUT ARGUMENTS
%   img : beamformed image [N x N]
%     x : lateral coordinate vector [1 x N]  
%     z : axial coordinate vector [1 x N]
%
%OUTPUT ARGUMENTS
% depths : axial positions of detected wires [M x 1]
%   fwhm : lateral resolution measurements at each depth [M x 1]
%
%% INITIALIZE PARAMETERS
N = size(img,1); % Get number of axial samples

% Calculate -6 dB threshold for half maximum in logarithmic scale
half_dB = 20 * log10(0.5);

%% DETECT WIRE PEAKS
% Find local maxima along center lateral line to identify wire positions
[~, idx] = findpeaks(img(:,N/2),"MinPeakDistance",20);

M = length(idx); % Number of detected wires

% Initialize output arrays
fwhm = zeros(M,1);
depths = zeros(M,1);

%% MEASURE FWHM FOR EACH WIRE
for i = 1:M
    % Extract lateral slice at current wire depth
    slice = img(idx(i),:);

    % Find peak magnitude and position in lateral direction
    [mag_center, idx_center] = max(slice);
    
    % Calculate -6 dB threshold relative to peak
    threshold = mag_center + half_dB;
    
    %% FIND THRESHOLD CROSSINGS
    % Scan right from peak to find where signal drops below threshold
    idx_right = idx_center;
    while idx_right < length(slice) && slice(idx_right) > threshold
        idx_right = idx_right + 1;
    end
    
    % Scan left from peak to find where signal drops below threshold 
    idx_left = idx_center;
    while idx_left > 1 && slice(idx_left) > threshold
        idx_left = idx_left - 1;
    end
    
    %% INTERPOLATE CROSSING POSITIONS
    % Linear interpolation for right crossing between adjacent samples
    t = (threshold - slice(idx_right-1)) / (slice(idx_right) - slice(idx_right-1));
    x_right = x(idx_right-1) + t * (x(idx_right) - x(idx_right-1));

    % Linear interpolation for left crossing between adjacent samples
    t = (threshold - slice(idx_left+1)) / (slice(idx_left) - slice(idx_left+1));
    x_left = x(idx_left+1) + t * (x(idx_left) - x(idx_left+1));

    %% CALCULATE FWHM
    % Compute lateral distance between -6 dB points
    fwhm(i) = abs(x_right - x_left);

    % Store corresponding depth of current wire
    depths(i) = z(idx(i));
end
end