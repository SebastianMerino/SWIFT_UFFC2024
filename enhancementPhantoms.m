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
for iAcq = 1:3
load(fullfile(dataDir,targetFiles(iAcq).name));
dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]
sam1 = RF(:,:,1);


BmodeFull = db(hilbert(sam1));
BmodeFull = BmodeFull - max(BmodeFull(:));


%% Cropping and finding sample sizes
% Region for attenuation imaging
x_inf = 0; x_sup = 5;
z_inf = 0.1; z_sup = 4;

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

% figure('Units','centimeters', 'Position',[5 5 18 6]),
% tiledlayout(1,2)
% nexttile, imagesc(x,z,Bmode,dynRange)
% axis image
% colormap(gray)
% colorbar('westoutside')
% title('Bmode')
% 
% nexttile,plot(fpxx/1e6,mean(pxx,2))
% title('Spectrum')
% xlabel('f [MHz]')
% grid on


%%
fc = 7.5E6; freqTol = 0.5e6;
z0 = 2.9; zf = 3.5;
x0Inc = 1.5; xfInc = 2.05;
x0Out = 0.25; xfOut = 0.70;
x0Out2 = 2.95; xfOut2 = 3.40;
% x0Out = 0.25; xfOut = 0.80;
% x0Out2 = 2.85; xfOut2 = 3.45;
% x0Inc = 1.25; xfInc = 1.85;


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

yline(z0, 'b--','LineWidth',1.5)
yline(zf, 'b--','LineWidth',1.5)
xline(x0Inc, 'g--', 'LineWidth',2)
xline(xfInc, 'g--', 'LineWidth',2)
xline(x0Out, 'r--', 'LineWidth',2)
xline(xfOut, 'r--', 'LineWidth',2)
xline(x0Out2, 'r--', 'LineWidth',2)
xline(xfOut2, 'r--', 'LineWidth',2)

nexttile,plot(fpxx/1e6,mean(pxx,2))
title('Spectrum')
xlabel('f [MHz]')
grid on

%%
[~,Z] = meshgrid(x,z);
mask = Z>z0 & Z<zf;

BmodeFilt(~mask) = NaN;
latProfile = mean(BmodeFilt,"omitmissing");
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
underInclusion = mean(latProfile(x>x0Inc & x <xfInc));
ousideInclusion = mean(latProfile((x>x0Out & x <xfOut)|(x>x0Out2 & x <xfOut2)));
acEnhancement = underInclusion - ousideInclusion
attDiff = acEnhancement/2/incDiameter/fc*1E6

end