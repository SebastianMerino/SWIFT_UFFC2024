clear,clc

dataDir = 'C:\Users\sebas\Documents\Data\Attenuation\Simulation\24_04_04_inc';
refDir = 'C:\Users\sebas\Documents\Data\Attenuation\Simulation\24_04_25_ref';
resultsDir = 'C:\Users\sebas\Documents\Data\Attenuation\JournalResults\ablation';

[~,~] = mkdir(resultsDir);
targetFiles = dir([dataDir,'\rf*.mat']);
refFiles = dir([refDir,'\rf*.mat']);
tableName = 'simuInc.xlsx';

%%
blocksize = 8;     % Block size in wavelengths
freq_L = 3.5e6; freq_H = 8.5e6; % original 3.3-8.7s
overlap_pc      = 0.8;
ratio_zx        = 12/8;

% New simu
referenceAtt    = 0.6;
groundTruthBack = [0.5,0.5,0.5];
groundTruthInc = [1,1,1];

% Weight parameters
muB = 10^3; muC = 10^0;
ratioCutOff = 10;
order = 5;
reject = 0.1;
extension = 3;

% SWTV
aSNR = 5; bSNR = 0.09;
desvMin = 15;

% Plotting
dynRange = [-40,0];
attRange = [0.4,1.1];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

% Region for attenuation imaging
x_inf = -1.5; x_sup = 1.5;
z_inf = 0.4; z_sup = 3.7;
%% Setting up

for iAcq = 2:3
% Regularization
switch iAcq
    case 1
        muBtv = 10^3.5; muCtv = 10^3.5;
        muBswtv = 10^3; muCswtv = 10^3;
        muBtvl1 = 10^3.5; muCtvl1 = 10^2;
        muBswift = 10^3.5; muCswift = 10^2;
    case 2
        muBtv = 10^2.5; muCtv = 10^2.5;
        muBswtv = 10^2.5; muCswtv = 10^0;
        muBtvl1 = 10^3.5; muCtvl1 = 10^1;
        muBswift = 10^3.5; muCswift = 10^1;
    case 3
        muBtv = 10^3.5; muCtv = 10^3.5;
        muBswtv = 10^2.5; muCswtv = 10^0;
        muBtvl1 = 10^3.5; muCtvl1 = 10^1;
        muBswift = 10^3.5; muCswift = 10^1;
end

load(fullfile(dataDir,targetFiles(iAcq).name));

fprintf("Acquisition no. %i, patient %s\n",iAcq,targetFiles(iAcq).name);
dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]

sam1 = rf(:,:,1);
Bmode = db(hilbert(sam1));
dynRange = [-50,0];

%% Cropping and finding sample sizes
% Limits for ACS estimation
ind_x = x_inf <= x & x <= x_sup;
ind_z = z_inf <= z & z <= z_sup;
x = x(ind_x);
z = z(ind_z);
sam1 = sam1(ind_z,ind_x);
Bmode = Bmode(ind_z,ind_x);
Bmode = Bmode - max(Bmode(:));

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

%% Spectrum
% BW from spectrogram
% [pxx,fpxx] = pwelch(sam1-mean(sam1),500,400,500,fs);
% meanSpectrum = mean(pxx,2);
% [freq_L,freq_H] = findFreqBand(fpxx, meanSpectrum, 0.1);
% meanSpectrum = db(meanSpectrum./max(meanSpectrum));
% figure,plot(fpxx/1e6,meanSpectrum)
% xline([freq_L,freq_H]/1e6)
% xlabel('Frequency [MHz]')
% ylabel('Magnitude')
% xlim([0 15])
% grid on

% Frequency samples
NFFT = 2^(nextpow2(nz/2)+2);
band = (0:NFFT-1)'/NFFT * fs;   % [Hz] Band of frequencies
rang = band > freq_L & band < freq_H ;   % useful frequency range
f  = band(rang)*1e-6; % [MHz]
p = length(f);

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
% compensation = repmat(mean(compensation,3),1,1,p);

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

%% Setting up
% System of equations
b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];
tol = 1e-3;
clear mask
mask = ones(m,n,p);

% Creating masks and ideal map
rInc = 0.7; c1x = 0; c1z = 2.05;
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
inclusion = (Xq.^2 + (Zq-c1z).^2)<= (rInc-0.1)^2;
back = (Xq.^2 + (Zq-c1z).^2) >= (rInc+0.1)^2;
attIdeal = ones(size(Xq))*groundTruthBack(iAcq);
attIdeal((Xq.^2 + (Zq-c1z).^2)<= rInc^2) = groundTruthInc(iAcq);

inclusionACS = (X.^2 + (Z-c1z).^2)<= rInc^2;
attIdealACS = ones(size(X))*groundTruthBack(iAcq);
attIdealACS(inclusionACS) = groundTruthInc(iAcq); %incl = inclusion


%% TVL1
[Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muBswift,muCswift,m,n,tol,mask(:));
BTVL1 = reshape(Bn*NptodB,m,n);
CTVL1 = reshape(Cn*NptodB,m,n);

AttInterp = interp2(X,Z,BTVL1,Xq,Zq);
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.meanInc = mean(AttInterp(inclusion),"omitnan");
r.stdInc = std(AttInterp(inclusion),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
r.biasInc = mean( AttInterp(inclusion) - groundTruthInc(iAcq),"omitnan");
r.rmseBack = sqrt(mean( (AttInterp(back) - groundTruthBack(iAcq)).^2,"omitnan"));
r.rmseInc = sqrt(mean( (AttInterp(inclusion) - groundTruthInc(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdInc^2);
r.method = 'TVL1';
MetricsTVL1(iAcq) = r;


%% Weight map
% Weight map
w = (1-reject)*(abs(CTVL1)<ratioCutOff)+reject;
wExt = movmin(w,extension);

% Weight matrices and new system
W = repmat(wExt,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;

%% Weighted fidelity term
[Bn,~] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBswift,muCswift,m,n,tol, ...
    mask(:),ones(m,n));
BWFid = reshape(Bn*NptodB,m,n);


AttInterp = interp2(X,Z,BWFid,Xq,Zq);
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.meanInc = mean(AttInterp(inclusion),"omitnan");
r.stdInc = std(AttInterp(inclusion),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
r.biasInc = mean( AttInterp(inclusion) - groundTruthInc(iAcq),"omitnan");
r.rmseBack = sqrt(mean( (AttInterp(back) - groundTruthBack(iAcq)).^2,"omitnan"));
r.rmseInc = sqrt(mean( (AttInterp(inclusion) - groundTruthInc(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdInc^2);
r.method = 'WFid';
MetricsWFid(iAcq) = r;

%% Weighted regularization term
[Bn,~] = optimAdmmWeightedTvTikhonov(A1,A2,b(:),muBswift,muCswift,m,n,tol, ...
    mask(:),w);
BWReg = reshape(Bn*NptodB,m,n);


AttInterp = interp2(X,Z,BWReg,Xq,Zq);
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.meanInc = mean(AttInterp(inclusion),"omitnan");
r.stdInc = std(AttInterp(inclusion),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
r.biasInc = mean( AttInterp(inclusion) - groundTruthInc(iAcq),"omitnan");
r.rmseBack = sqrt(mean( (AttInterp(back) - groundTruthBack(iAcq)).^2,"omitnan"));
r.rmseInc = sqrt(mean( (AttInterp(inclusion) - groundTruthInc(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdInc^2);
r.method = 'WReg';
MetricsWReg(iAcq) = r;


%% SWIFT
[Bn,~] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBswift,muCswift,m,n,tol,mask(:),w);
BSWIFT = reshape(Bn*NptodB,m,n);


AttInterp = interp2(X,Z,BSWIFT,Xq,Zq);
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.meanInc = mean(AttInterp(inclusion),"omitnan");
r.stdInc = std(AttInterp(inclusion),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
r.biasInc = mean( AttInterp(inclusion) - groundTruthInc(iAcq),"omitnan");
r.rmseBack = sqrt(mean( (AttInterp(back) - groundTruthBack(iAcq)).^2,"omitnan"));
r.rmseInc = sqrt(mean( (AttInterp(inclusion) - groundTruthInc(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdInc^2);
r.method = 'SWIFT';
MetricsSWIFT(iAcq) = r;


%% Plotting
figure('Units','centimeters', 'Position',[5 5 24 4]);
tiledlayout(1,6, "Padding","tight", 'TileSpacing','compact');

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
c = colorbar(t1, 'westoutside');
c.Label.String = 'dB';
title('B-mode')
ylabel('Axial [cm]')
xlabel('Lateral [cm]')
% hold on 
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1)
% hold off

t2 = nexttile;
imagesc(x,z,attIdeal,attRange)
xlabel('Lateral [cm]'), % ylabel('Axial [cm]')
colormap(t2,turbo)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
title('Ideal')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BTVL1, attRange)
colormap(t1,turbo)
axis image
title('L1 norm')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')
% hold on 
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1, 'EdgeColor',[1 1 1])
% hold off

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BWFid, attRange)
colormap(t1,turbo)
axis image
title('W Fidelity')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')
% hold on 
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1, 'EdgeColor',[1 1 1])
% hold off
% 
t1 = nexttile; 
imagesc(x_ACS,z_ACS,BWReg, attRange)
colormap(t1,turbo)
axis image
title('W Regularization')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')
% hold on 
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1, 'EdgeColor',[1 1 1])
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
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1, 'EdgeColor',[1 1 1])
% hold off

fontsize(gcf,8,'points')

%% SWIFT2
% Weight matrices and new system
W = repmat(wExt,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;

% Second iteration
[Bn,~] = optimAdmmWeightedTv(A1w,A2w,bw,muBswift,muCswift,m,n,tol,mask(:),w);
BSWIFT2 = reshape(Bn*NptodB,m,n);


AttInterp = interp2(X,Z,BSWIFT2,Xq,Zq);
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.meanInc = mean(AttInterp(inclusion),"omitnan");
r.stdInc = std(AttInterp(inclusion),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
r.biasInc = mean( AttInterp(inclusion) - groundTruthInc(iAcq),"omitnan");
r.rmseBack = sqrt(mean( (AttInterp(back) - groundTruthBack(iAcq)).^2,"omitnan"));
r.rmseInc = sqrt(mean( (AttInterp(inclusion) - groundTruthInc(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdInc^2);
r.method = 'SWIFT2';
MetricsSWIFT2(iAcq) = r;


figure('Units','centimeters', 'Position',[5 5 16 4]);
tiledlayout(1,4, "Padding","tight", 'TileSpacing','compact');

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
c = colorbar(t1, 'westoutside');
c.Label.String = 'dB';
title('B-mode')
ylabel('Axial [cm]')
xlabel('Lateral [cm]')
% hold on 
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1)
% hold off

t2 = nexttile;
imagesc(x,z,attIdeal,attRange)
xlabel('Lateral [cm]'), % ylabel('Axial [cm]')
colormap(t2,turbo)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
title('Ideal')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BSWIFT, attRange)
colormap(t1,turbo)
axis image
title('SWIFT')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')
% hold on 
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1, 'EdgeColor',[1 1 1])
% hold off

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BSWIFT2, attRange)
colormap(t1,turbo)
axis image
title('SWIFT TV')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')
% hold on 
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1, 'EdgeColor',[1 1 1])
% hold off
% 


fontsize(gcf,8,'points')
%%
cbPlacement = 'eastoutside';
figure('Units','centimeters', 'Position',[5 5 11 8])
tl = tiledlayout(2,2, "Padding","compact", 'TileSpacing','compact');
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
c = colorbar(t1, cbPlacement);
c.Label.String = 'dB';
title('Bmode')
xlabel('Lateral [cm]')
ylabel('Axial [cm]')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,CTVL1, bsRange)
colormap(t2,parula)
axis image
title('\DeltaBSC')
c = colorbar(t2, cbPlacement);
c.Label.String = 'dB/cm';
xlabel('Lateral [cm]')
ylabel('Axial [cm]')

t3 = nexttile([1,2]); 
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t3,parula)
axis image
title('Weights')
c = colorbar(t3, cbPlacement);
xlabel('Lateral [cm]')
ylabel('Axial [cm]')
fontsize(gcf,9,'points')


% %%
% axialTV = mean(BRTV(:,41:49),2);
% axialSWTV = mean(BSWTV(:,41:49),2);
% axialTVL1 = mean(BTVL1(:,41:49),2);
% axialSWIFT = mean(BSWIFT(:,41:49),2);
% 
% lateralTV = mean(BRTV(24:26,:),1);
% lateralSWTV = mean(BSWTV(24:26,:),1);
% lateralTVL1 = mean(BTVL1(24:26,:),1);
% lateralSWIFT = mean(BSWIFT(24:26,:),1);
% 
% %% Lateral and axial profiles
% lineColors = [0.635 0.078 0.184; 0.466 0.674 0.188; 0.301 0.745 0.933];
% 
% figure('Units','centimeters', 'Position',[5 5 14 4])
% tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')
% % figure('Units','centimeters', 'Position',[5 5 12 12])
% nexttile,
% plot(z_ACS, axialTV, ':', 'LineWidth',1.5, 'Color',lineColors(1,:) ),
% hold on
% plot(z_ACS, axialSWTV, '-.', 'LineWidth',1.5, 'Color',lineColors(2,:) ),
% % plot(z_ACS, axialTVL1, 'b:', 'LineWidth',1.5, 'Color',lineColors(2,:) ),
% plot(z_ACS, axialSWIFT, '-', 'LineWidth',1.5, 'Color',lineColors(3,:) ),
% plot(z_ACS,mean(attIdealACS(:,41:49),2), '--', 'Color', [0.2 0.2 0.2])
% hold off
% grid on
% ylim([0 1.5])
% xlim([z_ACS(1) z_ACS(end)])
% %title('Axial profile')
% ylabel('ACS [dB/cm/MHz]')
% xlabel('Axial [cm]')
% 
% nexttile,
% plot(x_ACS, lateralTV, ':', 'LineWidth',1.5, 'Color',lineColors(1,:) ),
% hold on
% plot(x_ACS, lateralSWTV, '-.', 'LineWidth',1.5, 'Color',lineColors(2,:) ),
% % plot(z_ACS, lateralTVL1, 'b:', 'LineWidth',1.5, 'Color',lineColors(2,:) ),
% plot(x_ACS, lateralSWIFT, '-', 'LineWidth',1.5, 'Color',lineColors(3,:) ),
% plot(x_ACS,mean(attIdealACS(24:26,:),1), 'k--')
% hold off
% grid on
% ylim([0 1.5])
% xlim([x_ACS(1) x_ACS(end)])
% %title('Lateral profile')
% xlabel('Lateral [cm]')

%%
save_all_figures_to_directory(resultsDir,['ablationSim',num2str(iAcq),'Fig'],'svg');
close all

end


%%
results1 = struct2table(MetricsTVL1);
results2 = struct2table(MetricsWFid);
results3 = struct2table(MetricsWReg);
results4 = struct2table(MetricsSWIFT);
results5 = struct2table(MetricsSWIFT2);



T = [results1;results2;results3;results4;results5];
writetable(T,fullfile(resultsDir,tableName),...
     'WriteRowNames',true);
