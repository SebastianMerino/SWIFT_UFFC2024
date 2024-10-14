
% This example demonstrates the use of k-Wave for the reconstruction of a
% two-dimensional photoacoustic wave-field recorded  over a linear array of
% sensor elements  The sensor data is simulated using kspaceFirstOrder2D
% and reconstructed using kspaceLineRecon. It builds on the Homogeneous
% Propagation Medium and Heterogeneous Propagation Medium examples.

function [] = simulationFocusedRTx(BaseDir)

addpath(genpath(pwd))
simuNames = {'inclusion1','inclusion2','inclusion3'};
% mkdir(BaseDir)

% medium parameters
c0              = 1540;     % sound speed [m/s]
rho0            = 1000;     % density [kg/m^3]

% source parameters
source_f0       = 6.66e6;   % source frequency [Hz]
source_amp      = 1e6;      % source pressure [Pa]
source_cycles   = 3.5;      % number of toneburst cycles
source_focus    = 20e-3;    % focal length [m]
element_pitch   = 0.3e-3;   % pitch [m]
element_width   = 0.25e-3;  % width [m]
focal_number    = 2;
nLines          = 96;

% grid parameters
grid_size_x     = 50e-3;    % [m]
grid_size_y     = 40e-3;    % [m]

% transducer position
translation     = [-20e-3, 0];
rotation        = 0;

% computational parameters
DATA_CAST       = 'gpuArray-single'; % set to 'single' or 'gpuArray-single'
ppw             = 6;        % number of points per wavelength, 4 to 8
depth           = 40e-3;    % imaging depth [m]
cfl             = 0.3;      % CFL number, could be 0.3 or 0.5

%% For looping simulations

for iSim = 1:length(simuNames)
    %% GRID

    % calculate the grid spacing based on the PPW and F0
    dx = c0 / (ppw * source_f0);   % [m]
    
    % compute the size of the grid
    Nx = roundEven(grid_size_x / dx);
    Ny = roundEven(grid_size_y / dx);
    
    % create the computational grid
    kgrid = kWaveGrid(Nx, dx, Ny, dx);
    
    % create the time array
    t_end           = depth*2/c0;     % [s];    % total compute time [s]
    kgrid.makeTime(c0, cfl, t_end);

    %% MEDIUM
    rz = kgrid.x - translation(1);
    rx = kgrid.y;

    switch iSim
        case 1
            background_std = 0.004;
            background_alpha = 0.5;       % [dB/(MHz^y cm)]
            layer_std = 0.004;
            layer_alpha = 1.0;
        case 2
            background_std = 0.002;
            background_alpha = 0.5;       % [dB/(MHz^y cm)]
            layer_std = 0.008;
            layer_alpha = 1.0;
        case 3
            background_std = 0.008;
            background_alpha = 0.5;       % [dB/(MHz^y cm)]
            layer_std = 0.002;
            layer_alpha = 1.0;
    end

    medium = addRegionSimu([],c0,rho0,background_std,...
        background_alpha,ones(Nx,Ny));

    cz = 20e-3; cx = 0; r = 7e-3;
    maskNodule = (rz-cz).^2 + (rx-cx).^2 < r^2;
    % maskLayer = rz>2e-2;
    medium = addRegionSimu(medium,c0,rho0,layer_std,...
        layer_alpha,maskNodule);

    medium.alpha_power = 1;
    medium.alpha_mode = 'no_dispersion';
    medium.sound_speed_ref = 1540;

    %% SOURCE
    aperture = source_focus/focal_number;
    element_num = floor(aperture/element_pitch);
    
    % set indices for each element
    ids = (0:element_num-1) - (element_num-1)/2;
    
    % set time delays for each element to focus at source_focus
    time_delays = -(sqrt((ids .* element_pitch).^2 + source_focus.^2) - source_focus) ./ c0;
    time_delays = time_delays - min(time_delays);
    
    % create time varying source signals (one for each physical element)
    source_sig = source_amp .* toneBurst(1/kgrid.dt, source_f0, ...
        source_cycles, 'SignalOffset', round(time_delays / kgrid.dt));
    
    % create empty kWaveArray
    karray = kWaveArray('BLITolerance', 0.05, 'UpsamplingRate', 10);
    
    % add rectangular elements
    for ind = 1:element_num
        
        % set element y position
        y_pos = 0 - (element_num * element_pitch/2 - element_pitch/2) + (ind-1) * element_pitch;
        
        % add element (set rotation angle to match the global rotation angle)
        karray.addRectElement([0, y_pos], element_width/4, element_width, rotation);
    end

    %% For Looping
    yCords = ( (0:nLines-1) - (nLines-1)/2 )* element_pitch; % Lateral cord of each element
    bf_data_final = zeros(kgrid.Nt ,nLines);

    for iLine = 1:nLines
        %% SOURCE  
        % move the array
        translation(2) = yCords(iLine);
        karray.setArrayPosition(translation, rotation)
    
        % assign binary mask
        source.p_mask = karray.getArrayBinaryMask(kgrid);
            
        % assign source signals
        source.p = karray.getDistributedSourceSignal(kgrid, source_sig);
        clc
        disp(['Line: ',num2str(iLine),' of ',num2str(nLines)]);
    
        %% SENSOR
        
        % set sensor mask to record central plane
        sensor.mask = karray.getArrayBinaryMask(kgrid);
        sensor.directivity_size = 10*kgrid.dx;
        sensor.directivity_angle = zeros(size(sensor.mask));
        
        %% SIMULATION
        
        % set input options
        input_args = {...
            'PMLInside', false, ...
            'PMLSize', 'auto', ... 
            'DataCast', DATA_CAST, ...
            'DataRecast', true, ...
            'PlotSim', false};
        
        % MATLAB CPU
        sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
            input_args{:});

                % combine sensor data
        combined_sensor_data = karray.combineSensorData(kgrid, sensor_data);

        % beamforming
        fs = 1/kgrid.dt;
        bf_data = BFangle(combined_sensor_data',element_num,fs,c0,element_pitch,...
            'rect',focal_number,0);
        index = floor(element_num/2)+1;
        bf_data_final(:,iLine) = bf_data(:,index);

    end

    %% VISUALISATION
    % offset = 60;
    offset = 400;
    
    axAxis = 0:size(bf_data_final,1)-1; axAxis = axAxis*kgrid.dt*c0/2;
    latAxis = 0:nLines-1; latAxis = latAxis-mean(latAxis); latAxis = latAxis *element_pitch;
    
    rf = bf_data_final(offset:end,:);
    z = axAxis(offset:end);
    x = latAxis;
    
    % plot the pressure field 
%     Bmode = db(hilbert(rf));
%     Bmode = Bmode - max(Bmode(:));
%     figure;
%     imagesc(1e3 * x, 1e3 * z, Bmode, [-50,-5]);
%     xlabel('z-position [mm]');
%     ylabel('x-position [mm]');
%     axis image;
%     title('Intensity');
%     cb = colorbar;
%     title(cb, '[dB]');
%     colormap gray

    save(fullfile(BaseDir,['rf_',simuNames{iSim},'.mat']),...
        'rf','x','z','fs','medium');

end

end