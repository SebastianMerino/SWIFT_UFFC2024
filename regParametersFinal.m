%% S1 to S3
setup,
% clear rsldNrmse swtvNrmse swiftNrmse


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

rsldNrmse(:,:,ii) = (rsld.rmseTop/0.75 + rsld.rmseBottom/0.75)/2;
swtvNrmse(:,:,ii) = (swtv.rmseTop/0.75 + swtv.rmseBottom/0.75)/2;
swiftNrmse(:,:,ii) = (swift.rmseTop/0.75 + swift.rmseBottom/0.75)/2;

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

%% P1 to P3
gtAcs = [0.97,0.95,0.95,0.55];

for ii =1:3
    resultsDir = "P:\smerino\UFFC2024results\reg\P"+ii;
    
    rsld = load(fullfile(resultsDir,'rsld.mat'));
    swtv = load(fullfile(resultsDir,'swtv.mat'));
    swift = load(fullfile(resultsDir,'swift.mat'));
    
    rsldNrmse(:,:,ii+7) = (rsld.rmseBack/gtAcs(4) + rsld.rmseInc/gtAcs(ii))/2;
    swtvNrmse(:,:,ii+7) = (swtv.rmseBack/gtAcs(4) + swtv.rmseInc)/gtAcs(ii)/2;
    swiftNrmse(:,:,ii+7) = (swift.rmseBack/gtAcs(4) + swift.rmseInc/gtAcs(ii))/2;
    
    rsldCnr(:,:,ii+7) = rsld.cnr;
    swtvCnr(:,:,ii+7) = swtv.cnr;
    swiftCnr(:,:,ii+7) = swift.cnr;
end

%% ALL CASES
muRangeLog = log10(rsld.muRange);
muBRangeLog = log10(swtv.muB);
muCRangeLog = log10(swtv.muC);

for ii = 1:10
    figure("Units","centimeters", "Position",[5 5 18 5]) 
    tiledlayout(1,3, "TileSpacing","compact", "Padding","compact")
    nexttile,
    plot(muRangeLog,rsldNrmse(:,:,ii))
    xlabel('log_{10}(\mu)')
    ylabel('NRMSE')
    grid on
    xlim([0, 7])
    ylim([0.05 1.1])
    title('RSLD')
    
    nexttile,
    imagesc(muBRangeLog,muCRangeLog,swtvNrmse(:,:,ii), [0.05 1.1])
    xlabel('log_{10}(\mu_B)')
    ylabel('log_{10}(\mu_C)')
    axis image
    title('SWTV-ACE')
    % ylim([-1.4 2.6])
    
    nexttile,
    imagesc(muBRangeLog,muCRangeLog,swiftNrmse(:,:,ii), [0.05 1.1])
    c = colorbar;
    c.Label.String = 'NRMSE';
    xlabel('log_{10}(\mu_B)')
    ylabel('log_{10}(\mu_C)')
    axis image
    title('SWIFT')
    % ylim([-1.4 2.6])
    hold on 
    valid = swiftNrmse(:,:,ii)<min(swtvNrmse(:,:,ii),[],'all');
    contour(muBRangeLog,muCRangeLog,valid,...
        1,'w--', 'LineWidth',1.5)
    hold off
    
    
    
    figure("Units","centimeters", "Position",[5 5 18 5]) 
    tiledlayout(1,3, "TileSpacing","compact", "Padding","compact")
    nexttile,
    plot(muRangeLog,rsldCnr(:,:,ii))
    xlabel('log_{10}(\mu)')
    ylabel('CNR')
    grid on
    xlim([0, 7])
    ylim([0.05 5])
    title('RSLD')
    
    nexttile,
    imagesc(muBRangeLog,muCRangeLog,swtvCnr(:,:,ii), [0,8])
    % colorbar,
    xlabel('log_{10}(\mu_B)')
    ylabel('log_{10}(\mu_C)')
    axis image
    title('SWTV-ACE')
    % ylim([-1.4 2.6])
    
    nexttile,
    imagesc(muBRangeLog,muCRangeLog,swiftCnr(:,:,ii), [0,8])
    c = colorbar;
    c.Label.String = 'CNR';
    xlabel('log_{10}(\mu_B)')
    % ylabel('log_{10}(\mu_C)')
    axis image
    title('SWIFT')
    % ylim([-1.4 2.6])
    hold on 
    valid = swiftCnr(:,:,ii)>max(swtvCnr(:,:,ii),[],'all');
    contour(muBRangeLog,muCRangeLog,valid,...
        1,'w--', 'LineWidth',1.5)
    hold off

    pause(0.5)
    save_all_figures_to_directory('P:\smerino\UFFC2024results\reg', ...
        "sam"+ii+"fig");
    close all
end


%% Region of improved performance, cases with changes in BSC and ACS
selectedCases = [2,3,6,7,9,10];

figure("Units","centimeters", "Position",[5 5 18 5]) 
tiledlayout(1,3, "TileSpacing","compact", "Padding","compact")
nexttile,
plot(log10(rsld.muRange),mean(rsldNrmse,3,"omitmissing"))
xlabel('log_{10}(\mu)')
ylabel('NRMSE')
grid on
xlim([0.5, 7])
ylim([0.05 1.1])
title('RSLD')

nexttile,
imagesc(muBRangeLog,muCRangeLog,mean(swtvNrmse,3), [0.1 1])
% c = colorbar;
% c.Label.String = 'NRMSE';
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWTV')
xlim([1.5 5.5])
ylim([-1.5 2.5])
hold on

nexttile,
imagesc(muBRangeLog,muCRangeLog,mean(swiftNrmse,3), [0.1 1])
c = colorbar;
c.Label.String = 'NRMSE';
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWIFT')
hold on
for ii = selectedCases
    valid = swiftNrmse(:,:,ii)<min(swtvNrmse(:,:,ii),[],'all');
    contourf(muBRangeLog,muCRangeLog,valid,[1 2], ...
        'w-.', "FaceAlpha",0.1, "FaceColor", [1 1 1])
end
hold off
xlim([1.5 5.5])
ylim([-1.5 2.5])


figure("Units","centimeters", "Position",[5 5 18 5]) 
tiledlayout(1,3, "TileSpacing","compact", "Padding","compact")
nexttile,
plot(log10(rsld.muRange),mean(rsldCnr,3,"omitmissing"))
xlabel('log_{10}(\mu)')
ylabel('CNR')
grid on
xlim([0.5, 7])
ylim([0.05, 2.5])
title('RSLD')

nexttile,
imagesc(muBRangeLog,muCRangeLog,mean(swtvCnr(:,:,[1:3,5:10]),3), [0 8])
%c = colorbar;
% c.Label.String = 'CNR';
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWTV')
xlim([1.5 5.5])
ylim([-1.5 2.5])
hold on

nexttile,
imagesc(muBRangeLog,muCRangeLog,mean(swiftCnr(:,:,[1:3,5:10]),3), [0 8])
c = colorbar;
c.Label.String = 'CNR';
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWIFT')
hold on
for ii = selectedCases
    valid = swiftCnr(:,:,ii)>max(swtvCnr(:,:,ii),[],'all');
    contourf(muBRangeLog,muCRangeLog,valid,[1 2], ...
        'w-.', "FaceAlpha",0.1, "FaceColor", [1 1 1])
end
hold off
xlim([1.5 5.5])
ylim([-1.5 2.5])

%%  Level curves, Half of cases
minNrmse = 0.25;
maxCnr = 1.8;
selectedCases = 1:2:10;

figure("Units","centimeters", "Position",[5 5 18 5]) 
tiledlayout(1,3, "TileSpacing","compact", "Padding","compact")
nexttile,
plot(log10(rsld.muRange),mean(rsldNrmse,3))
xlabel('log_{10}(\mu)')
ylabel('NRMSE')
grid on
xlim([0.5, 7])
ylim([0.05 1.1])
title('RSLD')

nexttile,
imagesc(muBRangeLog,muCRangeLog,mean(swtvNrmse,3), [0.1 1])
% c = colorbar;
% c.Label.String = 'NRMSE';
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWTV')
xlim([1.5 5.5])
ylim([-1.5 2.5])
hold on

nexttile,
imagesc(muBRangeLog,muCRangeLog,mean(swiftNrmse,3), [0.1 1])
c = colorbar;
c.Label.String = 'NRMSE';
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWIFT')
hold on
for ii = selectedCases
    valid = swiftNrmse(:,:,ii)<minNrmse;
    contourf(muBRangeLog,muCRangeLog,valid,[1 2], ...
        'w-.', "FaceAlpha",0.1, "FaceColor", [1 1 1])
end
hold off
xlim([1.5 5.5])
ylim([-1.5 2.5])


figure("Units","centimeters", "Position",[5 5 18 5]) 
tiledlayout(1,3, "TileSpacing","compact", "Padding","compact")
nexttile,
plot(log10(rsld.muRange),mean(rsldCnr,3,"omitmissing"))
xlabel('log_{10}(\mu)')
ylabel('CNR')
grid on
xlim([0.5, 7])
ylim([0.05 2.5])
title('RSLD')

nexttile,
imagesc(muBRangeLog,muCRangeLog,mean(swtvCnr(:,:,[1:3,5:10]),3), [0 10])
%c = colorbar;
% c.Label.String = 'CNR';
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWTV')
xlim([1.5 5.5])
ylim([-1.5 2.5])
hold on

nexttile,
imagesc(muBRangeLog,muCRangeLog,mean(swiftCnr(:,:,[1:3,5:10]),3), [0 10])
c = colorbar;
c.Label.String = 'CNR';
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWIFT')
hold on
for ii = selectedCases
    valid = swiftCnr(:,:,ii)>maxCnr;
    contourf(muBRangeLog,muCRangeLog,valid,[1 2], ...
        'w-.', "FaceAlpha",0.1, "FaceColor", [1 1 1])
end
hold off
xlim([1.5 5.5])
ylim([-1.5 2.5])


%% 
save_all_figures_to_directory('P:\smerino\UFFC2024results\reg','regFinal','svg');
close all