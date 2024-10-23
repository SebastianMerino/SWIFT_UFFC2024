% ====================================================================== %
% Script to fin the acoustic enhancement posterior to a colloid nodule.
% For Journal.
% ====================================================================== %
clear,
baseDir = 'C:\Users\sebas\Documents\Data\Attenuation\Thyroid_Data_PUCP_UTD';
refsDir = 'C:\Users\sebas\Documents\Data\Attenuation\REFERENCES';
resultsDir = 'C:\Users\sebas\Documents\Data\Attenuation\JournalResults\24-09-18';

% resultsDir = 'C:\Users\sebas\Pictures\Journal2024\24-02-29\';
if (~exist(resultsDir,"dir")), mkdir(resultsDir); end

T = readtable('params.xlsx');


blocksize = 8;     % Block size in wavelengths
fixedBW = true;
ratio = db2mag(-30);
freq_L = 3.5e6; freq_H = 8e6;
overlap_pc      = 0.8;
ratio_zx        = 12/8;

%% Loading case
for iAcq = 1:3;
patient = num2str(T.patient(iAcq));
samPath = fullfile(baseDir,patient,[patient,'-',T.sample{iAcq},'.rf']);
refDir = fullfile(refsDir,T.reference{iAcq});

out =lectura_OK(samPath);
sam1 = out.RF(:,:,1);
fs = out.fs;
fc = out.fc;
x = out.x; z = out.z;
fprintf("\n Selecting acq. no. %i, patient %s\n",iAcq,patient);


% Manual cropping
dx = x(2)-x(1);
dz = z(2)-z(1);
xFull = x*1e2; % [cm]
zFull = z*1e2; % [cm]

BmodeFull = db(hilbert(sam1));
BmodeFull = BmodeFull - max(BmodeFull(:));

% Manual cropping
dynRange = [-50,-10];
figure('Units','centimeters', 'Position',[5 5 15 15]),
imagesc(xFull,zFull,BmodeFull,dynRange); axis image; 
colormap gray; clim(dynRange);
hb2=colorbar; ylabel(hb2,'dB')
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
ylim([0.1 3.5])
title(patient)

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
x_inf = 0; x_sup = 4;
z_inf = 0.1; z_sup = 3.5;

% Limits for ACS estimation
ind_x = x_inf <= xFull & xFull <= x_sup;
ind_z = z_inf <= zFull & zFull <= z_sup;
roi = ind_x.*ind_z';
x = xFull(ind_x);
z = zFull(ind_z);
sam1 = sam1(ind_z,ind_x);

% Plot region of interest B-mode image
Bmode = db(hilbert(sam1));
Bmode = Bmode - max(Bmode(:));

dynRange = [-50,0];
attRange = [0.3,1.7];
%attRange = [0,1]; % Just for 13 acq
bsRange = [-2 2];



%% ACOUSTIC ENHANCEMENT

fs = 40000000;
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
fc = 8E6;
freqTol = 0.5e6;
[bFilt,aFilt] = butter(2,[fc-freqTol fc+freqTol]/fs*2, "bandpass");
samFilt = filtfilt(bFilt,aFilt,sam1);
[pxx,fpxx] = pwelch(samFilt,300,250,512,fs);

BmodeFilt = db(hilbert(samFilt));
BmodeFilt = BmodeFilt - max(BmodeFilt(:));

figure('Units','centimeters', 'Position',[5 5 18 6]),
tiledlayout(1,2)
nexttile, imagesc(x,z,BmodeFilt,[-70 0])
axis image
colormap(gray)
colorbar('westoutside')
title('Bmode')

switch iAcq
    case 1
        z0 = 1.8; zf = 1.9;
    case 2
        z0 = 2.7; zf = 3;
    case 3
        z0 = 1.6; zf = 2.5;
end
switch iAcq
    case 1
        x0Inc = 1.5; xfInc = 3;
        x0Out = 0.2; xfOut = 0.9;
    case 2
        x0Inc = 2.3; xfInc = 3;
        x0Out = 0.6; xfOut = 1.6;
    case 3
        x0Inc = 1.4; xfInc = 2.1;
        x0Out = 2.3; xfOut = 3.2;
    otherwise
        
end
yline(z0, 'b--','LineWidth',1.5)
yline(zf, 'b--','LineWidth',1.5)
xline(x0Inc, 'g--', 'LineWidth',2)
xline(xfInc, 'g--', 'LineWidth',2)
xline(x0Out, 'r--', 'LineWidth',2)
xline(xfOut, 'r--', 'LineWidth',2)

nexttile,plot(fpxx/1e6,mean(pxx,2))
title('Spectrum')
xlabel('f [MHz]')
grid on

%%
[~,Z] = meshgrid(x,z);
mask = Z>z0 & Z<zf;


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

%%
switch iAcq
    case 1
        incDiameter = 1;
    case 2
        incDiameter = 1.3; % before 1.2
    case 3
        incDiameter = 0.85;
end
disp("Thyroid #"+iAcq)
underInclusion = mean(latProfile(x>x0Inc & x <xfInc))
ousideInclusion = mean(latProfile(x>x0Out & x <xfOut))
acEnhancement = underInclusion - ousideInclusion
attDiff = acEnhancement/2/incDiameter/fc*1E6
end