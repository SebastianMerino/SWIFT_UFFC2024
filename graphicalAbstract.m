% ====================================================================== %
% Script for plotting figures and results for simulation on circular
% inclusions, one resembling thyroid tissue
% Created on Jan 31, 2024
% ====================================================================== %

setup,
% baseDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\simulations_processed\inc_journal'];
% resultsDir = 'C:\Users\smerino.C084288\Pictures\JOURNAL\24-02-20\BS_8_12';

targetDir = ['C:\Users\sebas\Documents\Data\' ...
    'Attenuation\Simulation\24_04_04_inc'];
refDir = ['C:\Users\sebas\Documents\Data\' ...
    'Attenuation\Simulation\24_04_25_ref'];
resultsDir = ['C:\Users\sebas\Documents\Data\Attenuation\' ...
    'UFFC2024Results\graphical_abstract'];
if (~exist(resultsDir,"dir")), mkdir(resultsDir); end

targetFiles = dir([targetDir,'\rf*.mat']);
refFiles = dir([refDir,'\rf*.mat']);
%% Generating cropped data
% SETTING PARAMETERS
blocksize = 8;     % Block size in wavelengths
ratio_zx        = 12/8;

freq_L = 3e6; freq_H = 8e6; % GOOD
% freq_L = 3e6; freq_H = 9e6; % Also GOOD

overlap_pc      = 0.8;
referenceAtt    = 0.6;

% Weights SWTV
aSNR = 1; bSNR = 0.1;
desvMin = 15;

% Weight parameters
muB = 10^3; muC = 10^0;
ratioCutOff = 10;
order = 5;
reject = 0.1;
extension = 3;

% groundTruthThyroid = [0.6,1.5];
% groundTruthNodule = [1.2,0.8];
groundTruthThyroid = [0.8,1.5];
groundTruthNodule = [1.5,0.8];

%attRange = [0.6 1.7];
attRange = [0.4,1.1];

% Plotting
dynRange = [-60,0];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

%% For looping
iAcq = 2;

load(fullfile(targetDir,targetFiles(iAcq).name));

muBswift = 10^3.5; muCswift = 10^1;

fprintf("Acquisition no. %i, patient %s\n",iAcq,targetFiles(iAcq).name);
dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]

sam1 = rf(:,:,1);

%% Cropping and finding sample sizes
% Region for attenuation imaging
x_inf = -1.5; x_sup = 1.5;
z_inf = 0.4; z_sup = 3.7;

% Limits for ACS estimation
ind_x = x_inf <= x & x <= x_sup;
ind_z = z_inf <= z & z <= z_sup;
x = x(ind_x);
z = z(ind_z);
sam1 = sam1(ind_z,ind_x);

% Wavelength size
c0 = 1540;
wl = c0/mean([freq_L freq_H]);   % Wavelength (m)

% Lateral samples
wx = round(blocksize*wl*(1-overlap_pc)/dx);  % Between windows
nx = round(blocksize*wl/dx);                 % Window size
x0 = 1:wx:length(x)-nx;
x_ACS = x(1,x0+round(nx/2));
n  = length(x0);

% Axial samples
wz = round(blocksize*wl*(1-overlap_pc)/dz * ratio_zx); % Between windows
% nz = 2*round(blocksize*wl/dz /2); % Window size
nz = 2*round(blocksize*wl/dz /2 * ratio_zx); % Window size
L = (nz/2)*dz*100;   % (cm)
z0p = 1:wz:length(z)-nz;
z0d = z0p + nz/2;
z_ACS = z(z0p+ nz/2);
m  = length(z0p);

%% BW and spectrogram
% ratio = db2mag(-50);
% [pxx,fpxx] = pwelch(sam1-mean(sam1),500,400,500,fs);
% meanSpectrum = mean(pxx,2);
% figure,plot(fpxx/1e6,meanSpectrum)
% [freq_L,freq_H] = findFreqBand(fpxx, meanSpectrum, ratio);
% xline([freq_L,freq_H]/1e6)
% xlim([0 15])
% xlabel('Frequency [MHz]')
% ylabel('Magnitude')
% grid on

% Frequency samples
NFFT = 2^(nextpow2(nz/2)+2);
band = (0:NFFT-1)'/NFFT * fs;   % [Hz] Band of frequencies
rang = band > freq_L & band < freq_H ;   % useful frequency range
f  = band(rang)*1e-6; % [MHz]
p = length(f);

% Plot region of interest B-mode image
% dynRange = [-40 -10];
Bmode = db(hilbert(sam1));
Bmode = Bmode - max(Bmode(:));
% figure, imagesc(x,z,Bmode);
% axis image; colormap gray; clim(dynRange);
% hb2=colorbar; ylabel(hb2,'dB')
% xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');


fprintf('\nFrequency range: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
fprintf('Blocksize in wavelengths: %i\n',blocksize)
fprintf('Blocksize x: %.2f mm, z: %.2f mm\n',nx*dx*1e3,nz*dz*1e3)
fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
fprintf('Region of interest rows: %i, col: %i\n\n',m,n);

%% Generating Diffraction compensation
% Windows for spectrum
windowing = tukeywin(nz/2,0.25);
windowing = windowing*ones(1,nx);

% For looping
Nref = length(refFiles);

% Memory allocation
Sp_ref = zeros(m,n,p);
Sd_ref = zeros(m,n,p);
compensation = zeros(m,n,p,Nref);

for iRef = 1:Nref %Nref
    load(fullfile(refDir,refFiles(iRef).name),"rf","medium");
    acs_mean = medium.alpha_coeff(1,1);
    att_ref = acs_mean*(f.^medium.alpha_power)/NptodB;
    att_ref_map = repmat(reshape(att_ref,[1 1 p]),m,n,1);

    samRef = rf;
    samRef = samRef(ind_z,ind_x); % Cropping
    for jj=1:n
        for ii=1:m
            xw = x0(jj) ;   % x window
            zp = z0p(ii);
            zd = z0d(ii);

            sub_block_p = samRef(zp:zp+nz/2-1,xw:xw+nx-1);
            sub_block_d = samRef(zd:zd+nz/2-1,xw:xw+nx-1);
            [tempSp,~] = spectra(sub_block_p,windowing,0,nz/2,NFFT);
            [tempSd,~] = spectra(sub_block_d,windowing,0,nz/2,NFFT);

            Sp_ref(ii,jj,:) = (tempSp(rang));
            Sd_ref(ii,jj,:) = (tempSd(rang));
        end
    end
    compensation(:,:,:,iRef) = log(Sp_ref) - log(Sd_ref) - 4*L*att_ref_map;
end

compensation = mean(compensation,4);


%% Spectrum
Sp = zeros(m,n,p);
Sd = zeros(m,n,p);
%cuec = zeros([NFFT,1]);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = sam1(zp:zp+nz/2-1,xw:xw+nx-1);
        sub_block_d = sam1(zd:zd+nz/2-1,xw:xw+nx-1);

        [tempSp,~] = spectra(sub_block_p,windowing,0,nz/2,NFFT);
        [tempSd,~] = spectra(sub_block_d,windowing,0,nz/2,NFFT);
%         if (ii==floor(m/2))
%             cuec = cuec + tempSp;
%         end
        Sp(ii,jj,:) = (tempSp(rang));
        Sd(ii,jj,:) = (tempSd(rang));
    end
end

%% Setting Up

% System of equations
b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];
tol = 1e-3;
clear mask
mask = ones(m,n,p);


[X,Z] = meshgrid(x_ACS,z_ACS);
cz = 2; cx = 0;
rInc = 0.8;
inc = (Z-cz).^2 + (X-cx).^2 > rInc^2;
% thyroidMask = ~noduleMask;

c1x = 0; c1z = 2;
roiL = 0.7; roiD = 0.5;
roiLz = 1.2;
x0mask = c1x - roiL/2; 
z0mask = c1z - roiLz/2;

sldLine = reshape(b.*inc,m*n,p);
sldLine = sldLine(inc(:),:);
sldLine = mean(sldLine)'*NptodB;
% sldLine = sldLine(30,:)'*NptodB;

%fit1 = f\sldLine;
fit2 = [f ones(length(f),1)]\sldLine;
figure('Units','points', 'Position',[100 100 120 280]),

tl = tiledlayout(5,1, "Padding","tight", 'TileSpacing','tight');

t1 = nexttile([3,1]);
imagesc(x,z,Bmode,dynRange)
axis image
axis off
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
%xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(t1,gray)
% colorbar
title('B-mode')
%subtitle(' ')
%c = colorbar('Location', 'westoutside');
% c.Label.String = 'dB';
% fontsize(gcf,8,'points')

nexttile([2,1]);
plot(f,sldLine, 'LineWidth',1),
hold on,
%plot(f,fit1*f, '--')
plot(f,fit2(1)*f + fit2(2), 'k--')
hold off,
grid on,
xlim([freq_L,freq_H]/1e6),
ylim([1 3.5]),
% xlabel('Freq. [MHz]')
% ylabel('Att. [dB/cm]')
title('SLD')
%legend({'SLD','Fit 1', 'Fit 2'}, 'Location','northwest')

%plot(f,sldLine)
%grid on
fontsize(gcf,18,'points')

%% SWIFT
% First iteration
[Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muBswift,muCswift,m,n,tol,mask(:));
BTVL1 = reshape(Bn*NptodB,m,n);
CTVL1 = reshape(Cn*NptodB,m,n);

% Weight map
w = (1-reject)*(abs(CTVL1)<ratioCutOff)+reject;
wExt = movmin(w,extension);

% Weight matrices and new system
W = repmat(wExt,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;

% Second iteration
tic
[Bn,~,ite] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBswift,muCswift,m,n,tol,mask(:),w);
exTime = toc;
BSWIFT = reshape(Bn*NptodB,m,n);

%% Plotting
figure('Units','points', 'Position',[100 100 100 120]);
imagesc(x_ACS,z_ACS,BTVL1, attRange)
colormap(turbo)
axis image
title('ACS')
fontsize(gcf,18,'points')
axis off
figure('Units','points', 'Position',[100 100 100 120]);
imagesc(x_ACS,z_ACS,CTVL1, [-20 20])
colormap(parula)
axis image
axis off
title('\DeltaBCS')
fontsize(gcf,18,'points')

figure('Units','points', 'Position',[100 100 100 120]);
imagesc(x_ACS,z_ACS,wExt, [0 1])
colormap(parula)
axis image
title('Weights')
axis off
% c = colorbar;
fontsize(gcf,18,'points')

figure('Units','points', 'Position',[100 100 100 120]);
imagesc(x_ACS,z_ACS,BSWIFT, attRange)
colormap(turbo)
axis image
axis off
title('ACS - SWIFT')
% c = colorbar;
fontsize(gcf,18,'points')

%%

save_all_figures_to_directory(resultsDir,'abstract','svg');
close all

