% ====================================================================== %
% Script to fin the acoustic enhancement posterior to a colloid nodule.
% For Journal.
% ====================================================================== %
setup,
dataDir = ['C:\Users\sebas\Documents\Data\Attenuation' ...
    '\ID316V2\06-08-2023-Generic'];
refDir = ['C:\Users\sebas\Documents\Data\Attenuation' ...
    '\ID544V2\06-08-2023-Generic'];
resultsDir = 'C:\Users\sebas\Documents\Data\Attenuation\JournalResults\24-09-18';

% targetDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\phantoms\ID316V2\06-08-2023-Generic'];
% refDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\phantoms\ID544V2\06-08-2023-Generic'];
% resultsDir = 'C:\Users\smerino.C084288\Pictures\JOURNAL\24-02-20\BS_8_12';

rawFiles = dir([dataDir,'\*.rf']);
targetFiles = dir([dataDir,'\*.mat']);
targetFiles = targetFiles(end-2:end);
if ~exist("resultsDir","dir"); mkdir(resultsDir); end


%% Loading case
iAcq = 1;
load(fullfile(dataDir,targetFiles(iAcq).name));
dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]
sam1 = RF(:,:,1);


BmodeFull = db(hilbert(sam1));
BmodeFull = BmodeFull - max(BmodeFull(:));

% Manual cropping
dynRange = [-50,-10];
figure('Units','centimeters', 'Position',[5 5 15 15]),
imagesc(x,z,BmodeFull,dynRange); axis image; 
colormap gray; clim(dynRange);
hb2=colorbar; ylabel(hb2,'dB')
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
ylim([0.1 3.5])

% confirmation = '';
% while ~strcmp(confirmation,'Yes')
%     rect = getrect;
%     confirmation = questdlg('Sure?');
%     if strcmp(confirmation,'Cancel')
%         disp(rect)
%         break
%     end
% end
% close,


%% Cropping and finding sample sizes
% Region for attenuation imaging
% x_inf = rect(1); x_sup = rect(1)+rect(3);
% z_inf = rect(2); z_sup = rect(2)+rect(4);
x_inf = 0; x_sup = 5;
z_inf = 0.1; z_sup = 3.5;

% Limits for ACS estimation
ind_x = x_inf <= x & x <= x_sup;
ind_z = z_inf <= z & z <= z_sup;
roi = ind_x.*ind_z';
x = x(ind_x);
z = z(ind_z);
sam1 = sam1(ind_z,ind_x);

% Plot region of interest B-mode image
Bmode = db(hilbert(sam1));
Bmode = Bmode - max(Bmode(:));

dynRange = [-50,0];
attRange = [0.3,1.7];
%attRange = [0,1]; % Just for 13 acq
bsRange = [-2 2];



%% ACOUSTIC ENHANCEMENT

[pxx,fpxx] = pwelch(sam1,300,250,512,fs);

figure('Units','centimeters', 'Position',[5 5 18 6]),
tiledlayout(1,2)
nexttile, imagesc(x,z,Bmode,dynRange)
axis image
colormap(gray)
colorbar('westoutside')
title('Bmode')

nexttile,plot(fpxx/1e6,mean(pxx,2))
title('Spectrum')
xlabel('f [MHz]')
grid on


%%
fc = 7.5E6; freqTol = 0.5e6;
[bFilt,aFilt] = butter(3,[fc-freqTol fc+freqTol]/fs*2, "bandpass");
samFilt = filtfilt(bFilt,aFilt,sam1);
[pxx,fpxx] = pwelch(samFilt,300,250,512,fs);

BmodeFilt = db(hilbert(samFilt));
BmodeFilt = BmodeFilt - max(BmodeFilt(:));

figure('Units','centimeters', 'Position',[5 5 18 6]),
tiledlayout(1,2)
nexttile, imagesc(x,z,BmodeFilt,[-50 0])
axis image
colormap(gray)
colorbar('westoutside')
title('Bmode')

z0 = 3; zf = 3.4;
% switch iAcq
%     case 1
%         z0 = 1.8; zf = 1.9;
%         % z0 = 2; zf = 2.2;
%     case 2
%         z0 = 2.7; zf = 3;
%     case 3
%         z0 = 1.5; zf = 2.5;
%     otherwise
% 
%         % z0 = 2.2; zf = 2.5;
% end
yline(z0, 'b--','LineWidth',1.5)
yline(zf, 'b--','LineWidth',1.5)

nexttile,plot(fpxx/1e6,mean(pxx,2))
title('Spectrum')
xlabel('f [MHz]')
grid on

%%
[~,Z] = meshgrid(x,z);
mask = Z>z0 & Z<zf;

x0Inc = 1.6; xfInc = 2.2;
x0Out = 0.25; xfOut = 0.85;
x0Out2 = 2.95; xfOut2 = 3.55;

% switch iAcq
%     case 1
%         x0Inc = 1.9; xfInc = 2.1;
%         x0Out = 0.1; xfOut = 0.9;
%     case 2
%         x0Inc = 2.3; xfInc = 3;
%         x0Out = 0.6; xfOut = 1.6;
%     case 3
%         x0Inc = 1.4; xfInc = 2.1;
%         % x0Out = 2.5; xfOut = 3.5;
%         x0Out = 2.3; xfOut = 3.2;
%     otherwise
% 
% end
BmodeFilt(~mask) = NaN;
latProfile = median(BmodeFilt,"omitmissing");
figure('Units','centimeters', 'Position',[5 5 9 6]),
plot(x,latProfile)
grid on
axis tight
xlabel('Lateral [cm]')
title('Acoustic enhancement')
ylim([-45 -15])

xline(x0Inc, 'g--', 'LineWidth',2)
xline(xfInc, 'g--', 'LineWidth',2)
xline(x0Out, 'r--', 'LineWidth',2)
xline(xfOut, 'r--', 'LineWidth',2)
xline(x0Out2, 'r--', 'LineWidth',2)
xline(xfOut2, 'r--', 'LineWidth',2)

%%
incDiameter = 1.9;
% switch iAcq
%     case 1
%         incDiameter = 1;
%     case 2
%         incDiameter = 1.3; % before 1.2
%     case 3
%         incDiameter = 0.85;
%     otherwise
%         % 
%         % 
% end
underInclusion = mean(latProfile(x>x0Inc & x <xfInc))
ousideInclusion = mean(latProfile((x>x0Out & x <xfOut)|(x>x0Out2 & x <xfOut2)))
acEnhancement = underInclusion - ousideInclusion

attDiff = acEnhancement/2/incDiameter/fc*1E6
