setup,

% dataDir = 'C:\Users\sebas\Documents\Data\Attenuation\Simulation\24_04_04_layered';
% refDir = 'C:\Users\sebas\Documents\Data\Attenuation\Simulation\24_04_25_ref';
% resultsDir = 'C:\Users\sebas\Documents\Data\Attenuation\UFFC2024results\sim_layered';
dataDir = 'P:\smerino\simulation_acs\rf_data\24_04_04_layered';
refDir = 'P:\smerino\simulation_acs\rf_data\24_04_25_ref';
resultsDir = 'P:\smerino\UFFC2024results\simulation';

[~,~] = mkdir(resultsDir);

targetFiles = dir([dataDir,'\rf*.mat']);
refFiles = dir([refDir,'\rf*.mat']);

tableName = 'simuLayered.xlsx';

%% Generating cropped data
% SETTING PARAMETERS
blocksize = 8;     % Block size in wavelengths
freq_L = 3.5e6; freq_H = 8.5e6;

overlap_pc      = 0.8;
ratio_zx        = 12/8;

% Weight parameters
ratioCutOff = 10;
reject = 0.1;
extension = 3;

% SWTV
aSNR = 5; bSNR = 0.09;
desvMin = 15;

% Plotting
dynRange = [-50,0];
attRange = [0.4,1.1];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

% GT
groundTruthTop = [0.5,0.5,0.5,0.75,0.75];
groundTruthBottom = [1,1,1,0.75,0.75];

% Region for attenuation imaging
x_inf = -1.5; x_sup = 1.5;
z_inf = 0.4; z_sup = 3.7;

iAcq = 1;
%% For looping

for iAcq = 1:4
load(fullfile(dataDir,targetFiles(iAcq).name));

fprintf("Acquisition no. %i, patient %s\n",iAcq,targetFiles(iAcq).name);
dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]
z = z-3.5*medium.sound_speed_ref/6.66e6*100/4;

sam1 = rf(:,:,1);

%%
switch iAcq
    case 1
        % Regularization
        muBtv = 10^4; muCtv = 10^4;
        muBswtv = 10^3; muCswtv = 10^2;
        muBswift = 10^3.5; muCswift = 10^2;

    case 2
        % Regularization
        muBtv = 10^4; muCtv = 10^4;
        muBswtv = 10^3; muCswtv = 10^0;
        muBswift = 10^4; muCswift = 10^1.5;

    case 3
        % Regularization
        muBtv = 10^4; muCtv = 10^4;
        muBswtv = 10^3; muCswtv = 10^0.5;
        muBswift = 10^4; muCswift = 10^1.5;
    case 4
        % Regularization
        muBtv = 10^4; muCtv = 10^4;
        muBswtv = 10^3.5; muCswtv = 10^0;
        muBswift = 10^4; muCswift = 10^0;
end
%% Cropping and finding sample sizes

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
nz = 2*round(blocksize*wl/dz /2 * ratio_zx); % Window size
L = (nz/2)*dz*100;   % (cm)
z0p = 1:wz:length(z)-nz;
z0d = z0p + nz/2;
z_ACS = z(z0p+ nz/2);
m  = length(z0p);

%% BW and spectrogram
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
fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);

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

for iRef = 1:Nref
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
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = sam1(zp:zp+nz/2-1,xw:xw+nx-1);
        sub_block_d = sam1(zd:zd+nz/2-1,xw:xw+nx-1);

        [tempSp,~] = spectra(sub_block_p,windowing,0,nz/2,NFFT);
        [tempSd,~] = spectra(sub_block_d,windowing,0,nz/2,NFFT);
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

% Creating reference
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);

attIdeal = ones(size(Z));
attIdeal(Z<=2) = groundTruthTop(iAcq);
attIdeal(Z>2) = groundTruthBottom(iAcq);

top = Zq < 1.9; % 0.1 cm interface
bottom = Zq > 2.1;
%% TV
tic
[Bn,Cn,ite] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
exTime = toc;
BRTV = reshape(Bn*NptodB,m,n);
CRTV = reshape(Cn*NptodB,m,n);

axialTV = mean(BRTV,2);

AttInterp = interp2(X,Z,BRTV,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
r.rmseTop = sqrt(mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,"omitnan"));
r.rmseBottom = sqrt(mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
r.method = 'TV';
r.sample = iAcq;
r.ite = ite;
r.exTime = exTime;
MetricsTV(iAcq) = r;
%% SWTV
% Calculating SNR
envelope = abs(hilbert(sam1));
SNR = zeros(m,n);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = envelope(zp:zp+nz/2-1,xw:xw+nx-1);
        sub_block_d = envelope(zd:zd+nz/2-1,xw:xw+nx-1);
        
        temp = [sub_block_p(:);sub_block_d(:)];
        SNR(ii,jj) = mean(temp)/std(temp);
    end
end

% Calculating weights
SNRopt = sqrt(1/(4/pi - 1));
desvSNR = abs(SNR-SNRopt)/SNRopt*100;
wSNR = aSNR./(1 + exp(bSNR.*(desvSNR - desvMin)));

% Method
tic
[Bn,Cn,ite] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muBswtv,muCswtv, ...
    m,n,tol,mask(:),wSNR);
exTime = toc;
BRSWTV = reshape(Bn*NptodB,m,n);
CRSWTV = reshape(Cn*NptodB,m,n);

axialSWTV = mean(BRSWTV,2);
AttInterp = interp2(X,Z,BRSWTV,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
r.rmseTop = sqrt(mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,"omitnan"));
r.rmseBottom = sqrt(mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
r.method = 'SWTV';
r.sample = iAcq;
r.ite = ite;
r.exTime = exTime;
MetricsSWTV(iAcq) = r;

%% SWIFT
% First iteration
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muBswift,muCswift,m,n,tol,mask(:));
bscMap = reshape(Cn*NptodB,m,n);

% Weight map
w = (1-reject)*(abs(bscMap)<ratioCutOff)+reject;
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

axialWFR = mean(BSWIFT,2);
AttInterp = interp2(X,Z,BSWIFT,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
r.rmseTop = sqrt(mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,"omitnan"));
r.rmseBottom = sqrt(mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
r.method = 'SWIFT';
r.sample = iAcq;
r.ite = ite;
r.exTime = exTime;
MetricsWFR(iAcq) = r;


%% Plotting
figure('Units','centimeters', 'Position',[5 5 22 4]);
tiledlayout(1,5, "Padding","tight", 'TileSpacing','compact');

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
c = colorbar(t1, 'westoutside');
c.Label.String = '[db]';
title('B-mode')
ylabel('Axial [cm]')
xlabel('Lateral [cm]')
% hold on 
% xline(2.05, 'w--', 'LineWidth',1.5)
% hold off

t2 = nexttile;
imagesc(x_ACS,z_ACS,attIdeal,attRange)
xlabel('Lateral [cm]'), % ylabel('Axial [cm]')
colormap(t2,turbo)
axis image
title('Ideal')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRTV, attRange)
colormap(t1,turbo)
axis image
title('RSLD')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')
% hold on 
% yline(2, 'w--', 'LineWidth',1.5)
% hold off

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRSWTV, attRange)
colormap(t1,turbo)
axis image
title('SWTV-ACE')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')
% hold on 
% yline(2, 'w--', 'LineWidth',1.5)
% hold off

t4 = nexttile; 
imagesc(x_ACS,z_ACS,BSWIFT, attRange)
colormap(t4,turbo)
axis image
title('SWIFT')
c = colorbar;
c.Label.String = 'ACS [db/cm/MHz]';
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')
% hold on 
% yline(2, 'w--', 'LineWidth',1.5)
% hold off

fontsize(gcf,8,'points')
%%
figure('Units','centimeters', 'Position',[5 5 10 4])
tl = tiledlayout(1,2, "Padding","tight");

t2 = nexttile; 
imagesc(x_ACS,z_ACS,bscMap, bsRange)
colormap(t2,parula)
axis image
title('TVL1')
c = colorbar;
c.Label.String = '\DeltaBSC [dB/cm]';
xlabel('Lateral [cm]')
ylabel('Axial [cm]')

t3 = nexttile; 
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t3,parula)
axis image
title('Weights')
c = colorbar;
xlabel('Lateral [cm]')
ylabel('Axial [cm]')
% c.Label.String = '[a.u.]';

fontsize(gcf,8,'points')
%% Lateral and axial profiles
lineColors = [0.635 0.078 0.184; 0.466 0.674 0.188; 0.301 0.745 0.933];
lw = 2;

figure('Units','centimeters', 'Position',[5 5 14 4])
tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')
nexttile,
plot(z_ACS, axialTV, ':', 'LineWidth',lw, 'Color',lineColors(1,:) ),
hold on
plot(z_ACS, axialSWTV, '-.', 'LineWidth',lw, 'Color',lineColors(2,:) ),
% plot(z_ACS, axialTVL1, 'b:', 'LineWidth',1.5, 'Color',lineColors(2,:) ),
plot(z_ACS, axialWFR, '-', 'LineWidth',lw, 'Color',lineColors(3,:) ),
plot(z_ACS,mean(attIdeal,2), '--', 'Color', [0.2 0.2 0.2])
hold off
grid on
ylim([0 1.5])
xlim([z_ACS(1) z_ACS(end)])
%title('Axial profile')
xlabel('Axial [cm]')
ylabel('ACS [dB/cm/MHz]')
if iAcq==1
    legend({'RSLD','SWTV-ACE','SWIFT'}, 'Location','northwest') ;
end
%%
save_all_figures_to_directory(resultsDir,['simLayered',num2str(iAcq),'Fig'],'svg');
close all

end

%%
results1 = struct2table(MetricsTV);
results2 = struct2table(MetricsSWTV);
results4 = struct2table(MetricsWFR);

T = [results1;results2;results4];
writetable(T,fullfile(resultsDir,tableName),...
     'WriteRowNames',true);

%%
