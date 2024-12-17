% ======================================================================
% ======================================================================
%% PHANTOMSSS
clear, clc
% dataDir = ['C:\Users\sebas\Documents\Data\Attenuation' ...
%     '\ID316V2\06-08-2023-Generic'];
% refDir = ['C:\Users\sebas\Documents\Data\Attenuation' ...
%     '\ID544V2\06-08-2023-Generic'];
% resultsDir = 'C:\Users\sebas\Documents\Data\Attenuation\UFFC2024results\phantoms';

dataDir = ['C:\Users\smerino.C084288\Documents\Datasets\' ...
    'Attenuation\phantoms\ID316V2\06-08-2023-Generic'];
refDir = ['C:\Users\smerino.C084288\Documents\Datasets\' ...
    'Attenuation\phantoms\ID544V2\06-08-2023-Generic'];
resultsDir0 = 'P:\smerino\UFFC2024results\reg';
[~,~] = mkdir(resultsDir0);

rawFiles = dir([dataDir,'\*.rf']);
targetFiles = dir([dataDir,'\*.mat']);
targetFiles = targetFiles(end-2:end);
if ~exist("resultsDir","dir"); mkdir(resultsDir); end
tableName = 'phantoms.xlsx';

%% Constants
blocksize = 8;     % Block size in wavelengths
freq_L = 2.5e6; freq_H = 7.5e6;
freq_C = (freq_L + freq_H)/2;

overlap_pc      = 0.8;
ratio_zx        = 12/8;
x_inf = 0.1; x_sup = 3.8;
z_inf = 0.2; z_sup = 3.5;
NptodB = log10(exp(1))*20;

% Weight parameters
ratioCutOff = 10;
reject = 0.1;
extension = 3;

% SWTV
aSNR = 5; bSNR = 0.09;
desvMin = 15;

groundTruthTargets = [0.97,0.95,0.95,0.55];

% Plotting constants
dynRange = [-50,0];
attRange = [0.2,1.2];

tol = 1e-3;

c1x = 1.95; c1z = 1.93;
roiL = 1; roiD = 0.6;
roiLz = 1.5;
%% For looping each phantom

for iAcq = 1:3
resultsDir = fullfile(resultsDir0,"P"+iAcq);
mkdir(resultsDir)

switch iAcq
    % Optimal reg for BS 8x12, circular ROIs
    case 1
        muBtv = 10^3.5; muCtv = 10^3.5;
        muBswtv = 10^3; muCswtv = 10^2.5;
        muBtvl1 = 10^3.5; muCtvl1 = 10^2;
        muBswift = 10^3.5; muCswift = 10^2;
    case 2
        muBtv = 10^3; muCtv = 10^3;
        muBswtv = 10^3; muCswtv = 10^0;
        muBtvl1 = 10^3.5; muCtvl1 = 10^1;
        muBswift = 10^3.5; muCswift = 10^1;
    case 3
        muBtv = 10^3.5; muCtv = 10^3.5;
        muBswtv = 10^3; muCswtv = 10^0;
        muBtvl1 = 10^3.5; muCtvl1 = 10^1;
        muBswift = 10^3.5; muCswift = 10^1;
end

switch iAcq
    case 1
        c1x = 1.8; c1z = 1.9;
    case 2
        c1x = 1.95; c1z = 1.95;
    case 3
        c1x = 1.85; c1z = 1.9;

end

%%
fprintf("Phantom no. %i, %s\n",iAcq,targetFiles(iAcq).name);
load(fullfile(dataDir,targetFiles(iAcq).name));

dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]

sam1 = RF(:,:,1);

%% Cropping and finding sample sizes

% Limits for ACS estimation
ind_x = x_inf <= x & x <= x_sup;
ind_z = z_inf <= z & z <= z_sup;
roi = ind_x.*ind_z';
x = x(ind_x);
z = z(ind_z);
sam1 = sam1(ind_z,ind_x);

% Wavelength size
c0 = 1540;
wl = c0/mean(freq_C);   % Wavelength (m)

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

% Frequency samples
NFFT = 2^(nextpow2(nz/2)+2);
band = (0:NFFT-1)'/NFFT * fs;   % [Hz] Band of frequencies
rang = band > freq_L & band < freq_H ;   % useful frequency range
f  = band(rang)*1e-6; % [MHz]
p = length(f);

% Plot region of interest B-mode image
Bmode = db(hilbert(sam1));
Bmode = Bmode - max(Bmode(:));

fprintf('\nFrequency range: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
fprintf('Blocksize in wavelengths: %i\n',blocksize)
fprintf('Blocksize x: %.2f mm, z: %.2f mm\n',nx*dx*1e3,nz*dz*1e3)
fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);



%% Generating Diffraction compensation
if true %iAcq == 1
    % Generating references
    att_ref = 0.53*f/NptodB; % From phantom especifications
    att_ref_map = zeros(m,n,p);
    for jj=1:n
        for ii=1:m
            att_ref_map(ii,jj,:) = att_ref;
        end
    end
    
    % Windows for spectrum
    windowing = hamming(nz/2);
    windowing = windowing*ones(1,nx);
    
    % For looping
    refFiles = dir([refDir,'\*.mat']);
    Nref = length(refFiles);
    
    % Memory allocation
    Sp_ref = zeros(m,n,p,Nref);
    Sd_ref = zeros(m,n,p,Nref);
    for iRef = 1:Nref
        out = load([refDir,'\',refFiles(iRef).name]);
        samRef = out.RF;
        samRef = samRef(ind_z,ind_x); % Cropping
        % figure,imagesc(db(hilbert(samRef)))
        for jj=1:n
            for ii=1:m
                xw = x0(jj) ;   % x window
                zp = z0p(ii);
                zd = z0d(ii);
    
                sub_block_p = samRef(zp:zp+nz/2-1,xw:xw+nx-1);
                sub_block_d = samRef(zd:zd+nz/2-1,xw:xw+nx-1);
                [tempSp,~] = spectra(sub_block_p,windowing,0,nz/2,NFFT);
                [tempSd,~] = spectra(sub_block_d,windowing,0,nz/2,NFFT);
    
                Sp_ref(ii,jj,:,iRef) = (tempSp(rang));
                Sd_ref(ii,jj,:,iRef) = (tempSd(rang));
            end
        end
    end
    
    Sp = mean(Sp_ref,4); Sd = mean(Sd_ref,4);
    compensation = ( log(Sp) - log(Sd) ) - 4*L*att_ref_map;
    
    % Liberating memory to avoid killing my RAM
    clear Sp_ref Sd_ref
end

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

%% ROI selection
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
rInc = 0.95;
inc = ((Xq-c1x).^2 + (Zq-c1z).^2)<= (rInc-0.1)^2;
back = ((Xq-c1x).^2 + (Zq-c1z).^2) >= (rInc+0.1)^2;

x0mask = c1x - roiL/2; 
z0mask = c1z - roiLz/2;
% [back,inc] = getRegionMasks(x,z,c1x,c1z,roiL,roiD,roiLz);

% Setting up
b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
mask = ones(m,n,p);

%% RSLD-TV
disp('RSLD')
muRange = 10.^(0:0.25:7);
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
    r.meanInc = mean(AttInterp(inc),"omitnan");
    r.stdInc = std(AttInterp(inc),"omitnan");
    r.meanBack = mean(AttInterp(back),"omitnan");
    r.stdBack = std(AttInterp(back),"omitnan");
    r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
    r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
    r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
        "omitnan") );
    r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
        "omitnan") );
    r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
    rmseBack(iMu) = r.rmseBack ;
    rmseInc(iMu) = r.rmseInc ;
    cnr(iMu) = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdInc^2);
end
save(fullfile(resultsDir,'rsld.mat'),"rmseTop","rmseBottom","cnr","muRange")


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


%% British Columbia Approach

% Weights
envelope = abs(hilbert(sam1));
SNR = zeros(m,n);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = envelope(zp:zp+nz/2-1,xw:xw+nx-1);
        sub_block_d = envelope(zd:zd+nz/2-1,xw:xw+nx-1);

        temp = [sub_block_p(:) sub_block_d(:)];
        SNR(ii,jj) = mean(temp)/std(temp);
    end
end
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
        r.meanInc = mean(AttInterp(inc),"omitnan");
        r.stdInc = std(AttInterp(inc),"omitnan");
        r.meanBack = mean(AttInterp(back),"omitnan");
        r.stdBack = std(AttInterp(back),"omitnan");
        r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
        r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
        r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
            "omitnan") );
        r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
            "omitnan") );
        r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
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
        r.meanInc = mean(AttInterp(inc),"omitnan");
        r.stdInc = std(AttInterp(inc),"omitnan");
        r.meanBack = mean(AttInterp(back),"omitnan");
        r.stdBack = std(AttInterp(back),"omitnan");
        r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
        r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
        r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
            "omitnan") );
        r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
            "omitnan") );
        r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
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
end
