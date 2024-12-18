setup,

% dataDir = 'C:\Users\sebas\Documents\Data\Attenuation\Simulation\24_11_05_multiInc';
% refDir = 'C:\Users\sebas\Documents\Data\Attenuation\Simulation\24_04_25_ref';
% resultsDir = 'C:\Users\sebas\Documents\Data\Attenuation\UFFC2024results\sim_multiInc';
dataDir = 'P:\smerino\simulation_acs\rf_data\24_11_05_multiInc';
refDir = 'P:\smerino\simulation_acs\rf_data\24_04_25_ref';
resultsDir = 'P:\smerino\UFFC2024results\simulation';

[~,~] = mkdir(resultsDir);
targetFiles = dir([dataDir,'\rf*.mat']);
refFiles = dir([refDir,'\rf*.mat']);
tableName = 'simuInc.xlsx';

%%
blocksize = 8;     % Block size in wavelengths
freq_L = 3.5e6; freq_H = 8.5e6; % original 3.3-8.7s
overlap_pc      = 0.8;
ratio_zx        = 3/2;
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
dynRange = [-40,0];
attRange = [0.4,1.1];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

% Region for attenuation imaging
x_inf = -1.5; x_sup = 1.5;
z_inf = 0.4; z_sup = 3.7;

iAcq = 3;

%% Setting up

% Regularization

muBtv = 10^4; muCtv = 10^4;
muBswtv = 10^3; muCswtv = 10^0;
muBswift = 10^3.5; muCswift = 10^1;

out = load(fullfile(dataDir,targetFiles(iAcq).name));
x = out.x; z = out.z;  fs = out.fs;

fprintf("Acquisition no. %i, patient %s\n",iAcq,targetFiles(iAcq).name);
dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]
z = z-0.05;

sam1 = out.rf(:,:,1);
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

if true; % iAcq == 1
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

%% Setting up
% System of equations
b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];
clear mask
mask = ones(m,n,p);

% Ideal maps
% [Xq,Zq] = meshgrid(x/1e2,z/1e2);
xMedium = out.rx(1,:)*100;
zMedium = out.rz(:,1)*100;
attIdeal = out.medium.alpha_coeff;
%% Creating masks
% Creating masks and ideal map
Xq = out.rx*100; Zq = out.rz*100;
rInc = 0.3; c1x = -0.6; c1z = 1.1;
inc1 = ((Xq-c1x).^2 + (Zq-c1z).^2)<= (rInc-0.1)^2;
rInc = 0.5; c1x = 0.6; c1z = 1.8;
inc2 = ((Xq-c1x).^2 + (Zq-c1z).^2)<= (rInc-0.1)^2;
rInc = 0.4; c1x = -0.4; c1z = 2.8;
inc3 = ((Xq-c1x).^2 + (Zq-c1z).^2)<= (rInc-0.1)^2;
back = ones(size(Xq)) & ~(inc1|inc2|inc3);
se = strel('disk',50,8);
back = imerode(back,se);

% figure, 
% imagesc(xMedium,zMedium,attIdeal, attRange)
% hold on
% contour(xMedium,zMedium,back,1,'w--')
% contour(xMedium,zMedium,inc1,1,'w--')
% contour(xMedium,zMedium,inc2,1,'w--')
% contour(xMedium,zMedium,inc3,1,'w--')
% hold off
% axis image
% colormap turbo

%% TV
tic
[Bn,Cn,ite] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
exTime = toc;
BRTV = reshape(Bn*NptodB,m,n);
CRTV = reshape(Cn*NptodB,m,n);

fprintf('\nExecution time: %.4f\n',exTime)
fprintf('Number of iterations: %d\n',ite)
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
BSWTV = reshape(Bn*NptodB,m,n);
CSWTV = reshape(Cn*NptodB,m,n);
fprintf('\nExecution time: %.4f\n',exTime)
fprintf('Number of iterations: %d\n',ite)

%% SWIFT
% muBswift = 10^3.3; muCswift = 10^0.8;
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
fprintf('\nExecution time: %.4f\n',exTime)
fprintf('Number of iterations: %d\n',ite)

%% Plotting
figure('Units','centimeters', 'Position',[5 5 13 8]);
tiledlayout(2,3, "Padding","tight", 'TileSpacing','compact');

t1 = nexttile([1,2]);
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
hold on
contour(xMedium,zMedium,back,1,'w--')
contour(xMedium,zMedium,inc1,1,'w--')
contour(xMedium,zMedium,inc2,1,'w--')
contour(xMedium,zMedium,inc3,1,'w--')
hold off
t2 = nexttile;
imagesc(xMedium,zMedium,attIdeal, attRange)
xlabel('Lateral [cm]'), % ylabel('Axial [cm]')
colormap(t2,turbo)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
title('Ideal')
c = colorbar;
c.Label.String = 'ACS [db/cm/MHz]';

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRTV, attRange)
colormap(t1,turbo)
axis image
title('RSLD')
ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BSWTV, attRange)
colormap(t1,turbo)
axis image
title('SWTV-ACE')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')

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
% contour(xMedium,zMedium,back,1,'w--')
% contour(xMedium,zMedium,inc1,1,'w--')
% contour(xMedium,zMedium,inc2,1,'w--')
% contour(xMedium,zMedium,inc3,1,'w--')
% hold off
fontsize(gcf,8,'points')

%% Metrics
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(xMedium,zMedium);
rsldInt = interp2(X,Z,BRTV,Xq,Zq);
swtvInt = interp2(X,Z,BSWTV,Xq,Zq);
swiftInt = interp2(X,Z,BSWIFT,Xq,Zq);

T = struct2table([
    getMetrics(rsldInt,back,1,'TV','back');
    getMetrics(rsldInt,inc1,0.5,'TV','inc1');
    getMetrics(rsldInt,inc2,0.5,'TV','inc2');
    getMetrics(rsldInt,inc3,1,'TV','inc3');
    getMetrics(swtvInt,back,1,'SWTV','back');
    getMetrics(swtvInt,inc1,0.5,'SWTV','inc1');
    getMetrics(swtvInt,inc2,0.5,'SWTV','inc2');
    getMetrics(swtvInt,inc3,1,'SWTV','inc3');
    getMetrics(swiftInt,back,1,'SWIFT','back');
    getMetrics(swiftInt,inc1,0.5,'SWIFT','inc1');
    getMetrics(swiftInt,inc2,0.5,'SWIFT','inc2');
    getMetrics(swiftInt,inc3,1,'SWIFT','inc3') ]);

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
% hold on
% contour(xMedium,zMedium,back,1,'w--')
% contour(xMedium,zMedium,inc1,1,'w--')
% contour(xMedium,zMedium,inc2,1,'w--')
% contour(xMedium,zMedium,inc3,1,'w--')
% hold off

t2 = nexttile; 
imagesc(x_ACS,z_ACS,bscMap, bsRange)
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

%%
save_all_figures_to_directory(resultsDir,['simMultiInc',num2str(iAcq),'Fig'],'svg');
% close all

writetable(T,fullfile(resultsDir,'multiInc.xlsx'))

%%
function r = getMetrics(image,mask,gt,method,region)
r.mean = mean(image(mask),'all','omitmissing');
r.std = std(image(mask),[],'all','omitmissing');
r.nrmse = sqrt(mean((image(mask) - gt).^2,'all','omitmissing'))/gt;
r.nbias = mean(image(mask) - gt,'all','omitmissing')/gt;
r.method = method;
r.region = region;
end