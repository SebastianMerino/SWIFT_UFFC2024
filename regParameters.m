clear,clc

% dataDir = 'C:\Users\sebas\Documents\Data\Attenuation\Simulation\24_04_04_inc';
% refDir = 'C:\Users\sebas\Documents\Data\Attenuation\Simulation\24_04_25_ref';
% resultsDir = 'C:\Users\sebas\Documents\Data\Attenuation\UFFC2024results\reg';
dataDir = 'P:\smerino\simulation_acs\rf_data\24_04_04_inc';
refDir = 'P:\smerino\simulation_acs\rf_data\24_04_25_ref';
resultsDir = 'P:\smerino\UFFC2024results\reg';

[~,~] = mkdir(resultsDir);
targetFiles = dir([dataDir,'\rf*.mat']);
refFiles = dir([refDir,'\rf*.mat']);
tableName = 'simuInc.xlsx';

%%
blocksize = 8;     % Block size in wavelengths
freq_L = 3.5e6; freq_H = 8.5e6; % original 3.3-8.7s
overlap_pc      = 0.8;
ratio_zx        = 12/8;
tol = 1e-3;

% New simu
referenceAtt    = 0.6;
groundTruthBack = [0.5,0.5,0.5];
groundTruthInc = [1,1,1];

% Weight parameters
ratioCutOff = 10;
reject = 0.1;
extension = 3;

% SWTV
aSNR = 5; bSNR = 0.09;
desvMin = 15;

% Plotting
attRange = [0.4,1.1];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

% Region for attenuation imaging
x_inf = -1.5; x_sup = 1.5;
z_inf = 0.4; z_sup = 3.7;

iAcq = 2;

%% Setting up


load(fullfile(dataDir,targetFiles(iAcq).name));

fprintf("Acquisition no. %i, patient %s\n",iAcq,targetFiles(iAcq).name);
dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]
z = z-3.5*medium.sound_speed_ref/6.66e6*100*0.5;

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
clear mask
mask = ones(m,n,p);


c1x = 0; c1z = 2;
roiL = 0.7; roiD = 0.5;
roiLz = 1.2;
rInc = 0.7;
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);

attIdeal = getIdealAcsInc([c1x,c1z],rInc,groundTruthBack(iAcq), ...
    groundTruthInc(iAcq),x_ACS,z_ACS,[blocksize blocksize*ratio_zx]*wl*100);
inc = (Xq.^2 + (Zq-2).^2)<= (rInc-0.1)^2;
back = (Xq.^2 + (Zq-2).^2) >= (rInc+0.1)^2;

%% TV
disp('RSLD')
muRange = 10.^(0:0.25:6);
rmseBack = zeros(size(muRange));
rmseInc = zeros(size(muRange));
cnr = zeros(size(muRange));
for iMu = 1:length(muRange)
    mutv = muRange(iMu);
    tic
    [Bn,~] = AlterOpti_ADMM(A1,A2,b(:),mutv,mutv,m,n,tol,mask(:));
    toc
    BRTV = reshape(Bn*NptodB,m,n);

    % imagesc(x_ACS,z_ACS,BRTV, attRange)
    % colormap(turbo)
    % axis image
    % title('RSLD')
    % % ylabel('Axial [cm]')
    % xlabel('Lateral [cm]')
    pause(0.01)

    AttInterp = interp2(X,Z,BRTV,Xq,Zq);
    r.meanBack = mean(AttInterp(back),"omitnan");
    r.stdBack = std(AttInterp(back),"omitnan");
    r.meanInc = mean(AttInterp(inc),"omitnan");
    r.stdInc = std(AttInterp(inc),"omitnan");
    r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
    r.biasInc = mean( AttInterp(inc) - groundTruthInc(iAcq),"omitnan");
    r.rmseBack = sqrt(mean( (AttInterp(back) - groundTruthBack(iAcq)).^2,"omitnan"));
    r.rmseInc = sqrt(mean( (AttInterp(inc) - groundTruthInc(iAcq)).^2,"omitnan"));
    rmseBack(iMu) = r.rmseBack ;
    rmseInc(iMu) = r.rmseInc ;
    cnr(iMu) = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdInc^2);
end
save(fullfile(resultsDir,'rsld.mat'),"rmseBack","rmseInc","cnr","muRange")

figure,
semilogx(muRange,(rmseBack/0.5 + rmseInc)/2)
grid on
xlabel('\mu')
ylabel('RMSE')

figure,
semilogx(muRange,cnr)
grid on
xlabel('\mu')
ylabel('CNR')


%% SWTV
disp('SWTV')

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

% Finding optimal reg parameters
muB = 10.^(1.5:0.25:4.5);
muC = 10.^(-1.5:0.25:2.5);
rmseBack = zeros(length(muC),length(muB));
rmseInc = zeros(length(muC),length(muB));
cnr = zeros(length(muC),length(muB));

for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        tic
        [Bn,~] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muB(mmB),muC(mmC),...
        m,n,tol*10,mask(:),wSNR);
        toc
        pause(0.01)
        BSWTV = reshape(Bn*NptodB,m,n);

        AttInterp = interp2(X,Z,BSWTV,Xq,Zq);
        r.meanBack = mean(AttInterp(back),"omitnan");
        r.stdBack = std(AttInterp(back),"omitnan");
        r.meanInc = mean(AttInterp(inc),"omitnan");
        r.stdInc = std(AttInterp(inc),"omitnan");
        r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
        r.biasInc = mean( AttInterp(inc) - groundTruthInc(iAcq),"omitnan");
        r.rmseBack = sqrt(mean( (AttInterp(back) - groundTruthBack(iAcq)).^2,"omitnan"));
        r.rmseInc = sqrt(mean( (AttInterp(inc) - groundTruthInc(iAcq)).^2,"omitnan"));
        % rmse(mmC,mmB) = (r.rmseBack + r.rmseInc)/2;
        rmseBack(mmC,mmB) = r.rmseBack ;
        rmseInc(mmC,mmB) = r.rmseInc ;
        cnr(mmC,mmB) = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdInc^2);
    end
end
save(fullfile(resultsDir,'swtv.mat'),"rmseBack","rmseInc","cnr","muB","muC")

figure,
imagesc(log10(muB),log10(muC),(rmseBack/0.5 + rmseInc)/2)
colorbar,
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')

figure,
imagesc(log10(muB),log10(muC),cnr)
colorbar,
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')


%% SWIFT
disp('SWIFT')

muB = 10.^(1.5:0.25:4.5);
muC = 10.^(-1.5:0.25:2.5);
rmseBack = zeros(length(muC),length(muB));
rmseInc = zeros(length(muC),length(muB));
cnr = zeros(length(muC),length(muB));

for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        muBswift = muB(mmB);
        muCswift = muC(mmC);
        
        tic
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
        [Bn,~,~] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBswift,muCswift,m,n,tol,mask(:),w);
        BSWIFT = reshape(Bn*NptodB,m,n);
        toc
        pause(0.01)

        AttInterp = interp2(X,Z,BSWIFT,Xq,Zq);
        r.meanBack = mean(AttInterp(back),"omitnan");
        r.stdBack = std(AttInterp(back),"omitnan");
        r.meanInc = mean(AttInterp(inc),"omitnan");
        r.stdInc = std(AttInterp(inc),"omitnan");
        r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
        r.biasInc = mean( AttInterp(inc) - groundTruthInc(iAcq),"omitnan");
        r.rmseBack = sqrt(mean( (AttInterp(back) - groundTruthBack(iAcq)).^2,"omitnan"));
        r.rmseInc = sqrt(mean( (AttInterp(inc) - groundTruthInc(iAcq)).^2,"omitnan"));
        % rmse(mmC,mmB) = (r.rmseBack + r.rmseInc)/2;
        rmseBack(mmC,mmB) = r.rmseBack ;
        rmseInc(mmC,mmB) = r.rmseInc ;
        cnr(mmC,mmB) = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdInc^2);
    end
end
%%
save(fullfile(resultsDir,'swift.mat'),"rmseBack","rmseInc","cnr","muB","muC")
% load(fullfile(resultsDir,'swift.mat'),"rmse","muB","muC")

figure,
imagesc(log10(muB),log10(muC),(rmseBack/0.5 + rmseInc)/2)
colorbar,
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image

figure,
imagesc(log10(muB),log10(muC),cnr)
colorbar,
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
%%
save_all_figures_to_directory(resultsDir,['sim',num2str(iAcq),'Fig']);
close all,
%%

rsld = load(fullfile(resultsDir,'rsld.mat'));
swtv = load(fullfile(resultsDir,'swtv.mat'));
swift = load(fullfile(resultsDir,'swift.mat'));

figure("Units","centimeters", "Position",[5 5 18 5]) 
tiledlayout(1,3, "TileSpacing","compact", "Padding","compact")
% nexttile([1,2])
% nexttile, axis off
rsld.nrmse = (rsld.rmseBack/0.5 + rsld.rmseInc)/2;
nexttile,
plot(log10(rsld.muRange),rsld.nrmse)
xlabel('log_{10}(\mu)')
ylabel('NRMSE')
grid on
xlim([0.5, 6])
ylim([0.05 1])
title('RSLD')

swtv.nrmse = (swtv.rmseBack/0.5 + swtv.rmseInc)/2;
nexttile,
imagesc(log10(swtv.muB),log10(swtv.muC),swtv.nrmse, [0.1 1])
% colorbar,
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWTV-ACE')
ylim([-1.4 2.6])

swift.nrmse = (swift.rmseBack/0.5 + swift.rmseInc)/2;
nexttile,
imagesc(log10(swift.muB),log10(swift.muC),swift.nrmse, [0.1 1])
c = colorbar;
c.Label.String = 'NRMSE';
xlabel('log_{10}(\mu_B)')
% ylabel('log_{10}(\mu_C)')
axis image
title('SWIFT')
ylim([-1.4 2.6])
hold on 
contour(log10(swift.muB),log10(swift.muC),swift.nrmse<min(swtv.nrmse(:)),...
    1,'w--', 'LineWidth',1.5)
% contour(log10(swift.muB),log10(swift.muC),swift.nrmse<(swtv.nrmse),...
%     1,'w--', 'LineWidth',1.5)
hold off


%%
figure("Units","centimeters", "Position",[5 5 18 5]) 
tiledlayout(1,3, "TileSpacing","compact", "Padding","compact")
% nexttile([1,2])
% nexttile, axis off
nexttile,
plot(log10(rsld.muRange),rsld.cnr)
xlabel('log_{10}(\mu)')
ylabel('CNR')
grid on
xlim([0.5, 6])
ylim([0.05 1])
title('RSLD')

nexttile,
imagesc(log10(swtv.muB),log10(swtv.muC),swtv.cnr, [0,7])
% colorbar,
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWTV-ACE')
ylim([-1.4 2.6])

nexttile,
imagesc(log10(swift.muB),log10(swift.muC),swift.cnr, [0,7])
c = colorbar;
c.Label.String = 'CNR';
xlabel('log_{10}(\mu_B)')
% ylabel('log_{10}(\mu_C)')
axis image
title('SWIFT')
ylim([-1.4 2.6])
hold on 
contour(log10(swift.muB),log10(swift.muC),swift.cnr>max(swtv.cnr(:)),...
    1,'w--', 'LineWidth',1.5)
% contour(log10(swift.muB),log10(swift.muC),swift.cnr>(swtv.cnr),...
%     1,'w--', 'LineWidth',1.5)
hold off

% save_all_figures_to_directory(resultsDir,'regFinal','svg');
% close all