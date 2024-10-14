% ====================================================================== %
% Script to create arbitrary masks for homogeneous ROIs in clinical data. 
% Created on March 25, 2024
% ====================================================================== %
setup,
baseDir = 'C:\Users\sebas\Documents\Data\Attenuation\Thyroid_Data_PUCP_UTD';
refsDir = 'C:\Users\sebas\Documents\Data\Attenuation\REFERENCES';
resultsDir = 'newMasks';

T = readtable('params.xlsx');
if (~exist(resultsDir,"dir")), mkdir(resultsDir); end


dynRange = [-50 0];

%% Loading case FULL VERSION
for iAcq = 1:height(T)

patient = num2str(T.patient(iAcq));
samPath = fullfile(baseDir,patient,[patient,'-',T.sample{iAcq},'.rf']);
refDir = fullfile(refsDir,T.reference{iAcq});

%%
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

figure('Units','centimeters', 'Position',[3 5 35 15]),
tiledlayout(1,2)
nexttile,
orig = imread(fullfile(baseDir,patient,[patient,'-',T.sample{iAcq},'.png']));
image(orig)
axis image
nexttile,
imagesc(xFull,zFull,BmodeFull,dynRange); axis image; 
colormap gray; clim(dynRange);
hb2=colorbar; ylabel(hb2,'dB')
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
ylim([0.1 3.5])

confirmation = '';
while ~strcmp(confirmation,'Yes')
    h = drawfreehand('Multiclick',true);
    confirmation = questdlg('Sure?');
    if strcmp(confirmation,'Cancel')
        break
    elseif strcmp(confirmation,'No')
        delete(h)
    end
end
regionMask = createMask(h);

[X,Z] = meshgrid(xFull,zFull);
x_inf = min(X(regionMask));
x_sup = max(X(regionMask));
z_inf = min(Z(regionMask));
z_sup = max(Z(regionMask));


rectangle('Position', [x_inf, z_inf, x_sup-x_inf, z_sup-z_inf], ...
    'LineStyle','--', 'LineWidth',2, 'EdgeColor','w')


%%
maskNodule = regionMask;
maskThyroid = ~regionMask; % for now

save(fullfile(resultsDir,[patient,'.mat']),...
    "maskNodule",'maskThyroid')


end

