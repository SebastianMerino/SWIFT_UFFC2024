%% Clinical colloid case
setup,
baseDir = 'C:\Users\sebas\Documents\Data\Attenuation\Thyroid_Data_PUCP_UTD';
refsDir = 'C:\Users\sebas\Documents\Data\Attenuation\REFERENCES';
resultsDir = 'C:\Users\sebas\Documents\Data\Attenuation\JournalResults\clinical';

tableName = 'clinical.xlsx';
T = readtable('params.xlsx');
if (~exist(resultsDir,"dir")), mkdir(resultsDir); end

%% Params
blocksize = 8;     % Block size in wavelengths
fixedBW = true;
ratio = db2mag(-30);
freq_L = 3.5e6; freq_H = 8e6;
overlap_pc      = 0.8;
ratio_zx        = 12/8;

% Weight parameters
ratioCutOff = 10;
reject = 0.1;
extension = 3;

% SWTV
aSNR = 5; bSNR = 0.09;
desvMin = 15;

% Plotting constants
dynRange = [-50,0];
attRange = [0.2,2];
bsRange = [-15 15];

iAcq = 3;

%% Loading case
for iRoi = [2,1]
patient = num2str(T.patient(iAcq));
class = T.class(iAcq);
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

if iRoi == 1
    % rect = [1.03; 0.49; 1.6; 1.69]; % Previous rectangle
    rect = [1.03; 0.49; 1.6; 1.65]; % Previous rectangle
    muBtv = 10^2.5; muCtv = 10^2.5;
    muBswtv = 10^2.5; muCswtv = 10^-0.5;
    muBswift = 10^3; muCswift = 10^0.5;
else
    % rect = [2.63; 0.49; 1.6; 1.69]; % Previous rectangle
    rect = [2.63; 0.49; 1.6; 1.65]; % Previous rectangle
    muBtv = 10^3; muCtv = 10^3;
    muBswtv = 10^3; muCswtv = 10^0;
    muBswift = 10^3.5; muCswift = 10^1;
end
% hold on
% rectangle('Position',rect)
% hold off
%% Cropping and finding sample sizes
% Region for attenuation imaging
x_inf = rect(1); x_sup = rect(1)+rect(3);
z_inf = rect(2); z_sup = rect(2)+rect(4);

% Limits for ACS estimation
ind_x = x_inf <= xFull & xFull <= x_sup;
ind_z = z_inf <= zFull & zFull <= z_sup;
roi = ind_x.*ind_z';
x = xFull(ind_x);
z = zFull(ind_z);
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

% Generating references
att_ref = attenuation_phantoms_Np(f, 4, []); % CAMBIAESTO
att_ref_map = zeros(m,n,p);
for jj=1:n
    for ii=1:m
        att_ref_map(ii,jj,:) = att_ref;
    end
end

% Windows for spectrum
% windowing = tukeywin(nz/2,0.25);
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

%% Setup

b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
tol = 1e-3;
mask = ones(m,n,p);
NptodB = log10(exp(1))*20;

[X,Z] = meshgrid(xFull,zFull);
roiACS{iRoi} = X >= x_ACS(1) & X <= x_ACS(end) & ...
    Z >= z_ACS(1) & Z <= z_ACS(end);


%% RSLD-TV
[Bn,~] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
BRTV = reshape(Bn*NptodB,m,n);

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
[Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muBswtv,muCswtv,...
m,n,tol,mask(:),wSNR);
BRSWTV = reshape(Bn*NptodB,m,n);
CRSWTV = reshape(Cn*NptodB,m,n);

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
[Bn,~] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBswift,muCswift,m,n,tol,mask(:),w);
BSWIFT = reshape(Bn*NptodB,m,n);

%% Metrics
if iRoi == 1
    load(fullfile('newMasks',[patient,'.mat']));
else
    maskThyroid = true(size(BmodeFull));
    maskNodule = maskThyroid;
end


[X,Z] = meshgrid(xFull,zFull);
[Xq,Zq] = meshgrid(x_ACS,z_ACS);
maskNoduleACS = interp2(X,Z,maskNodule,Xq,Zq, 'nearest');
maskThyroidACS = interp2(X,Z,maskThyroid,Xq,Zq, 'nearest');

se = strel('diamond',1);
maskThyroidACS = imerode(maskThyroidACS,se);
maskNoduleACS = imerode(maskNoduleACS,se);
%figure, imagesc(x_ACS,z_ACS,maskThyroidACS|maskNoduleACS)

patCol(iRoi) = {patient}; 
classCol(iRoi) = {class};
dataCols(iRoi,:) = [mean(BRTV(maskNoduleACS)), std(BRTV(maskNoduleACS)),...
    mean(BRSWTV(maskNoduleACS)), std(BRSWTV(maskNoduleACS)),...
    mean(BSWIFT(maskNoduleACS)), std(BSWIFT(maskNoduleACS)), ...
    mean(BRTV(maskThyroidACS)), std(BRTV(maskThyroidACS)),...
    mean(BRSWTV(maskThyroidACS)), std(BRSWTV(maskThyroidACS)),...
    mean(BSWIFT(maskThyroidACS)), std(BSWIFT(maskThyroidACS))]; 

fprintf("D ACS, TV: %.2f\n",...
    mean(BRTV(maskThyroidACS)) - mean(BRTV(maskNoduleACS)));
fprintf("D ACS, SWIFT: %.2f\n",...
    mean(BSWIFT(maskThyroidACS)) - mean(BSWIFT(maskNoduleACS)));


%% Saving
imageData.x = x_ACS;
imageData.z = z_ACS;
imageData.roi = roi;
imageData.TV = BRTV;
imageData.SWTV = BRSWTV;
imageData.SWIFT = BSWIFT;
dataRoi{iRoi} = imageData;

end

%% PLOTTING IMAGES
alpha_img = 0.7;
zlim1 = 0.21; zlim2 = 2.46;
xlim1 = 0.95; xlim2 = 3.8;

load(fullfile('newMasks',[patient,'.mat']));
[X,Z] = meshgrid(xFull,zFull);
[Xq,Zq] = meshgrid(x_ACS,z_ACS);
maskNoduleACS = interp2(X,Z,maskNodule,Xq,Zq, 'nearest');
maskThyroidACS = interp2(X,Z,maskThyroid,Xq,Zq, 'nearest');
se = strel('diamond',1);
maskThyroidACS = imerode(maskThyroidACS,se);
maskNoduleACS = imerode(maskNoduleACS,se);



figure('Units','centimeters', 'Position',[5 5 24 4.6])
tiledlayout(1,4, 'TileSpacing','compact', 'Padding','compact')
t1 = nexttile;
imagesc(xFull,zFull,BmodeFull,dynRange); axis image; 
title('B-mode')
ylim([0.1, 3])
hold on
contour(xFull,zFull,roiACS{1},1,'w--')
contour(xFull,zFull,roiACS{2},1,'w--')
hold off
xlabel('Lateral [cm]')
ylabel('Axial [cm]')
hBm = colorbar;
hBm.Label.String = 'dB';
hBm.Location = 'westoutside';

[X,Z] = meshgrid(x,z);
roiSub = X >= x_ACS(1) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);

t2 = nexttile;
iRoi = 1;
[~,hB,hColor] = imOverlayInterp(BmodeFull,dataRoi{iRoi}.TV,dynRange,attRange,alpha_img,...
    dataRoi{iRoi}.x,dataRoi{iRoi}.z,dataRoi{iRoi}.roi,xFull,zFull);
% Interpolation
iRoi = 2;
[X,Z] = meshgrid(dataRoi{iRoi}.x,dataRoi{iRoi}.z);
[Xq,Zq] = meshgrid(xFull,zFull);
imgInterp = interp2(X,Z,dataRoi{iRoi}.TV,Xq,Zq);
emptyRegion = isnan(imgInterp);
newRoi = ~emptyRegion & dataRoi{iRoi}.roi;
% Overlap
hold on;
hF = imagesc(dataRoi{iRoi}.x,dataRoi{iRoi}.z,imgInterp,attRange);
set(hF,'XData',get(hB,'XData'),'YData',get(hB,'YData'))
alphadata = alpha_img.*(newRoi);
set(hF,'AlphaData',alphadata);

contour(xFull,zFull,roiACS{1},1,'w--')
contour(xFull,zFull,roiACS{2},1,'w--')
contour(xFull,zFull,maskThyroid & roiACS{1},1,'w--')
hold off

ylim([0.1, 3])
xlabel('Lateral [cm]'), % ylabel('Axial [cm]')
title('RSLD')
colorbar off

t2 = nexttile;
iRoi = 1;
[~,hB,hColor] = imOverlayInterp(BmodeFull,dataRoi{iRoi}.SWTV,dynRange,attRange,alpha_img,...
    dataRoi{iRoi}.x,dataRoi{iRoi}.z,dataRoi{iRoi}.roi,xFull,zFull);
% Interpolation
iRoi = 2;
[X,Z] = meshgrid(dataRoi{iRoi}.x,dataRoi{iRoi}.z);
[Xq,Zq] = meshgrid(xFull,zFull);
imgInterp = interp2(X,Z,dataRoi{iRoi}.SWTV,Xq,Zq);
emptyRegion = isnan(imgInterp);
newRoi = ~emptyRegion & dataRoi{iRoi}.roi;
% Overlap
hold on;
hF = imagesc(dataRoi{iRoi}.x,dataRoi{iRoi}.z,imgInterp,attRange);
set(hF,'XData',get(hB,'XData'),'YData',get(hB,'YData'))
alphadata = alpha_img.*(newRoi);
set(hF,'AlphaData',alphadata);

contour(xFull,zFull,roiACS{1},1,'w--')
contour(xFull,zFull,roiACS{2},1,'w--')
contour(xFull,zFull,maskThyroid & roiACS{1},1,'w--')
hold off

ylim([0.1, 3])
xlabel('Lateral [cm]'), % ylabel('Axial [cm]')
title('SWTV-ACE')
colorbar off

nexttile,
iRoi = 1;
[~,hB,hColor] = imOverlayInterp(BmodeFull,dataRoi{iRoi}.SWIFT,dynRange,attRange,alpha_img,...
    dataRoi{iRoi}.x,dataRoi{iRoi}.z,dataRoi{iRoi}.roi,xFull,zFull);
% Interpolation
iRoi = 2;
[X,Z] = meshgrid(dataRoi{iRoi}.x,dataRoi{iRoi}.z);
[Xq,Zq] = meshgrid(xFull,zFull);
imgInterp = interp2(X,Z,dataRoi{iRoi}.SWIFT,Xq,Zq);
emptyRegion = isnan(imgInterp);
newRoi = ~emptyRegion & dataRoi{iRoi}.roi;
% Overlap
hold on;
hF = imagesc(dataRoi{iRoi}.x,dataRoi{iRoi}.z,imgInterp,attRange);
set(hF,'XData',get(hB,'XData'),'YData',get(hB,'YData'))
alphadata = alpha_img.*(newRoi);
set(hF,'AlphaData',alphadata);

contour(xFull,zFull,roiACS{1},1,'w--')
contour(xFull,zFull,roiACS{2},1,'w--')
contour(xFull,zFull,maskThyroid & roiACS{1},1,'w--')
hold off
ylim([0.1, 3])
xlabel('Lateral [cm]'), % ylabel('Axial [cm]')
title('SWIFT')

hColor.Layout.Tile = 'east';
hColor.Label.String = 'ACS [dB/cm/MHz]';
% colorbar off
colormap(t1,'gray')
fontsize(gcf,9,'points')


%% Table
infoTable = table(patCol',classCol',...
          'VariableNames',{'patient','type'});
dataTable = array2table(dataCols,...
    'VariableNames',{ ...
    'nod-TV-mean','nod-TV-std', ...
    'nod-SWTV-mean','nod-SWTV-std', ...
    'nod-SWIFT-mean','nod-SWIFT-std', ...
    'thy-TV-mean','thy-TV-std', ...
    'thy-SWTV-mean','thy-SWTV-std', ...
    'thy-SWIFT-mean','thy-SWIFT-std', ...
    });

writetable([infoTable,dataTable],fullfile(resultsDir,'clinical2Rois.xlsx'),...
     'WriteRowNames',true);

save_all_figures_to_directory(resultsDir,'specialFig','svg');
close all