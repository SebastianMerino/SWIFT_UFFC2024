%% S1 to S3

clear rsldNrmse swtvNrmse swiftNrmse


for ii = 1:3;
    resultsDir = "P:\smerino\UFFC2024results\reg\S"+ii;
    
    rsld = load(fullfile(resultsDir,'rsld.mat'));
    swtv = load(fullfile(resultsDir,'swtv.mat'));
    swift = load(fullfile(resultsDir,'swift.mat'));
    
    rsldNrmse(:,:,ii) = (rsld.rmseTop/0.5 + rsld.rmseBottom)/2;
    swtvNrmse(:,:,ii) = (swtv.rmseTop/0.5 + swtv.rmseBottom)/2;
    swiftNrmse(:,:,ii) = (swift.rmseTop/0.5 + swift.rmseBottom)/2;
    
    rsldCnr(:,:,ii) = rsld.cnr;
    swtvCnr(:,:,ii) = swtv.cnr;
    swiftCnr(:,:,ii) = swift.cnr;
end

%% S4 
ii = 4;
resultsDir = "P:\smerino\UFFC2024results\reg\S"+ii;
    
rsld = load(fullfile(resultsDir,'rsld.mat'));
swtv = load(fullfile(resultsDir,'swtv.mat'));
swift = load(fullfile(resultsDir,'swift.mat'));

rsldNrmse(:,:,ii) = (rsld.rmseTop/0.5 + rsld.rmseBottom)/2;
swtvNrmse(:,:,ii) = (swtv.rmseTop/0.5 + swtv.rmseBottom)/2;
swiftNrmse(:,:,ii) = (swift.rmseTop/0.5 + swift.rmseBottom)/2;

rsldCnr(:,:,ii) = rsld.cnr;
swtvCnr(:,:,ii) = swtv.cnr;
swiftCnr(:,:,ii) = swift.cnr;

%% S5 to S7
for ii = 5:7
    resultsDir = "P:\smerino\UFFC2024results\reg\S"+(ii);
    
    rsld = load(fullfile(resultsDir,'rsld.mat'));
    swtv = load(fullfile(resultsDir,'swtv.mat'));
    swift = load(fullfile(resultsDir,'swift.mat'));
    
    rsldNrmse(:,:,ii) = (rsld.rmseBack/0.5 + rsld.rmseInc)/2;
    swtvNrmse(:,:,ii) = (swtv.rmseBack/0.5 + swtv.rmseInc)/2;
    swiftNrmse(:,:,ii) = (swift.rmseBack/0.5 + swift.rmseInc)/2;
    
    rsldCnr(:,:,ii) = rsld.cnr;
    swtvCnr(:,:,ii) = swtv.cnr;
    swiftCnr(:,:,ii) = swift.cnr;
end

%%
selectedCases = [2,3,6,7];
rsldNrmse = rsldNrmse(:,:,selectedCases);
swtvNrmse = swtvNrmse(:,:,selectedCases);
swiftNrmse = swiftNrmse(:,:,selectedCases);
rsldCnr = rsldCnr(:,:,selectedCases);
swtvCnr = swtvCnr(:,:,selectedCases);
swiftCnr = swiftCnr(:,:,selectedCases);

%%
figure("Units","centimeters", "Position",[5 5 18 5]) 
tiledlayout(1,3, "TileSpacing","compact", "Padding","compact")
nexttile,
plot(log10(rsld.muRange),mean(rsldNrmse,3))
xlabel('log_{10}(\mu)')
ylabel('NRMSE')
grid on
xlim([0.5, 6])
ylim([0.05 1])
title('RSLD')

nexttile,
imagesc(log10(swtv.muB),log10(swtv.muC),mean(swtvNrmse,3), [0.1 1])
% c = colorbar;
% c.Label.String = 'NRMSE';
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWTV')
ylim([-1.4 2.6])
hold on
minNrmse = min(mean(swtvNrmse,3),[],'all');

nexttile,
imagesc(log10(swift.muB),log10(swift.muC),mean(swiftNrmse,3), [0.1 1])
c = colorbar;
c.Label.String = 'NRMSE';
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWIFT')
hold on
for ii = 1:4
    valid = swiftNrmse(:,:,ii)<minNrmse;
    % valid = swiftNrmse(:,:,ii)<min(swtvNrmse(:,:,ii),[],'all');
    contourf(log10(swift.muB),log10(swift.muC),valid,[1 2], ...
        'w-.', "FaceAlpha",0.1, "FaceColor", [1 1 1])
end
hold off
xlim([1.5 4.5])
ylim([-1.4 2.5])


%%
figure("Units","centimeters", "Position",[5 5 18 5]) 
tiledlayout(1,3, "TileSpacing","compact", "Padding","compact")
nexttile,
plot(log10(rsld.muRange),mean(rsldCnr(:,:,[3,4]),3))
xlabel('log_{10}(\mu)')
ylabel('CNR')
grid on
xlim([0.5, 6])
ylim([0.05 1])
title('RSLD')

nexttile,
imagesc(log10(swtv.muB),log10(swtv.muC),mean(swtvCnr,3), [0 8])
%c = colorbar;
% c.Label.String = 'CNR';
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWTV')
ylim([-1.4 2.6])
hold on
maxCnr = max(mean(swtvCnr,3),[],'all');

nexttile,
imagesc(log10(swift.muB),log10(swift.muC),mean(swiftCnr,3), [0 8])
c = colorbar;
c.Label.String = 'CNR';
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWIFT')
hold on
for ii = 1:4
    valid = swiftCnr(:,:,ii)>maxCnr;
    % valid = swiftCnr(:,:,ii)<min(swtvCnr(:,:,ii),[],'all');
    contourf(log10(swift.muB),log10(swift.muC),valid,[1 2], ...
        'w-.', "FaceAlpha",0.1, "FaceColor", [1 1 1])
end
hold off
xlim([1.5 4.5])
ylim([-1.4 2.5])

save_all_figures_to_directory('P:\smerino\UFFC2024results\reg','regFinal','svg');
close all