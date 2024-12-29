setup,
% clear rsldNrmse swtvNrmse swiftNrmse

%% S6 and S7
for ii = 1:2
    resultsDir = "P:\smerino\UFFC2024results\reg\S"+(ii+5);
    
    rsld = load(fullfile(resultsDir,'rsld.mat'));
    swtv = load(fullfile(resultsDir,'swtv.mat'));
    swift = load(fullfile(resultsDir,'swift.mat'));
    
    % rsldNrmse(:,:,ii) = (rsld.rmseBack/0.5 + rsld.rmseInc)/2;
    % swtvNrmse(:,:,ii) = (swtv.rmseBack/0.5 + swtv.rmseInc)/2;
    % swiftNrmse(:,:,ii) = (swift.rmseBack/0.5 + swift.rmseInc)/2;
    rsldNrmse(:,:,ii) = rsld.rmseInc;
    swtvNrmse(:,:,ii) = swtv.rmseInc;
    swiftNrmse(:,:,ii) = swift.rmseInc;

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
    
    % rsldNrmse(:,:,ii+1) = (rsld.rmseBack/gtAcs(4) + rsld.rmseInc/gtAcs(ii))/2;
    % swtvNrmse(:,:,ii+1) = (swtv.rmseBack/gtAcs(4) + swtv.rmseInc)/gtAcs(ii)/2;
    % swiftNrmse(:,:,ii+1) = (swift.rmseBack/gtAcs(4) + swift.rmseInc/gtAcs(ii))/2;
    rsldNrmse(:,:,ii+1) = rsld.rmseInc/gtAcs(ii);
    swtvNrmse(:,:,ii+1) = swtv.rmseInc/gtAcs(ii);
    swiftNrmse(:,:,ii+1) = swift.rmseInc/gtAcs(ii);

    rsldCnr(:,:,ii+1) = rsld.cnr;
    swtvCnr(:,:,ii+1) = swtv.cnr;
    swiftCnr(:,:,ii+1) = swift.cnr;
end

%%  Level curves, Half of cases
selectedCases = 1:4;
color1 = [0 0.4470 0.7410];

meanRsldNrmse = mean(rsldNrmse(:,:,selectedCases),3, "omitmissing");
meanSwtvNrmse = mean(swtvNrmse(:,:,selectedCases),3, "omitmissing");
meanSwiftNrmse = mean(swiftNrmse(:,:,selectedCases),3, "omitmissing");
stdRsldNrmse = std(rsldNrmse(:,:,selectedCases),[],3, "omitmissing");
stdSwtvNrmse = std(swtvNrmse(:,:,selectedCases),[],3, "omitmissing");
stdSwiftNrmse = std(swiftNrmse(:,:,selectedCases),[],3, "omitmissing");

figure("Units","centimeters", "Position",[5 5 18 5]) 
tiledlayout(1,3, "TileSpacing","compact", "Padding","compact")
nexttile,
errorbar(log10(rsld.muRange),meanRsldNrmse,stdRsldNrmse,...
    'CapSize',1.5, 'Color',color1)
hold on
plot(log10(rsld.muRange),meanRsldNrmse, 'Color',color1)
hold off
xlabel('log_{10}(\mu)')
ylabel('NRMSE')
grid on
xlim([0.5, 7])
ylim([0.05,0.9])
title('RSLD')

nexttile,
imagesc(log10(swtv.muB),log10(swtv.muC),meanSwtvNrmse, [0.1 0.95])
% c = colorbar;
% c.Label.String = 'NRMSE';
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWTV')
xlim([1.5 5.5])
ylim([-1.5 2.5])
hold on
contour(log10(swtv.muB),log10(swtv.muC),meanSwtvNrmse,[0.19, 10], ...
    '--', 'LineWidth',1.5, "EdgeColor",[0 0 0]);
contour(log10(swtv.muB),log10(swtv.muC),meanSwtvNrmse,[0.3, 10], ...
    '--', 'LineWidth',1.5, "EdgeColor",[0.3 0.3 0.3]);
% [~,ii] = min(meanSwtvNrmse(:));
% [i,j] = ind2sub(size(meanSwtvNrmse),ii);
% plot(log10(swtv.muB(j)),log10(swtv.muC(i)),'kx')
hold off
xlim([1.5 5.5]-0.5)
ylim([-1.5 2.5]-0.5)


nexttile,
imagesc(log10(swift.muB),log10(swift.muC),meanSwiftNrmse, [0.1 0.95])
c = colorbar;
c.Label.String = 'NRMSE';
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWIFT')
hold on
contour(log10(swift.muB),log10(swift.muC),meanSwiftNrmse,[0.19, 10], ...
    '--', 'LineWidth',1.5, "EdgeColor",[0 0 0]);
contour(log10(swift.muB),log10(swift.muC),meanSwiftNrmse,[0.3, 10], ...
    '--', 'LineWidth',1.5, "EdgeColor",[0.3 0.3 0.3]);
hold off
xlim([1.5 5.5])
ylim([-1.5 2.5])


%%
minNrmse = 0.3;
colors = lines(4);
lineSpec = {":","--","-","-."};
figure("Units","centimeters", "Position",[5 5 18 5.1]) 
tiledlayout(1,2, "TileSpacing","compact", "Padding","compact")

nexttile,
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWTV')
xlim([1.5 5.5])
ylim([-1.5 2.5])
hold on
for ii = selectedCases
    fprintf('SWTV NRMSE %i\n',ii)
    valid = swtvNrmse(:,:,ii)<minNrmse;
    contourf(log10(swtv.muB),log10(swtv.muC),-swtvNrmse(:,:,ii),-[minNrmse 0], ...
        lineSpec{ii}, "FaceAlpha",0.1, "FaceColor",colors(ii,:), ...
        "EdgeColor",colors(ii,:), 'LineWidth',1.5)
end
xlim([1.5 5.5]-0.5)
ylim([-1.5 2.5]-0.5)
clim([0.1 0.95])
hold off
ax = gca;
ax.YDir = "reverse";
grid on

nexttile,
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWIFT')
hold on
for ii = selectedCases
    fprintf('SWIFT NRMSE %i\n',ii)
    plot(0,0,lineSpec{ii}, 'LineWidth',1.5)
    contourf(log10(swift.muB),log10(swift.muC),-swiftNrmse(:,:,ii),-[minNrmse 0], ...
        lineSpec{ii}, "FaceAlpha",0.1, "FaceColor",colors(ii,:), ...
        "EdgeColor",colors(ii,:), 'LineWidth',1.5)
end
xlim([1.5 5.5])
ylim([-1.5 2.5])
clim([0.1 0.95])
hold off
ax = gca;
ax.YDir = "reverse";
grid on
legend('S6','','P1','','P2','','P3', 'Location','eastoutside')

%% CNR
meanRsldCnr = mean(rsldCnr(:,:,selectedCases),3, "omitmissing");
meanSwtvCnr = mean(swtvCnr(:,:,selectedCases),3, "omitmissing");
meanSwiftCnr = mean(swiftCnr(:,:,selectedCases),3, "omitmissing");
stdRsldCnr = std(rsldCnr(:,:,selectedCases),[],3, "omitmissing");
stdSwtvCnr = std(swtvCnr(:,:,selectedCases),[],3, "omitmissing");
stdSwiftCnr = std(swiftCnr(:,:,selectedCases),[],3, "omitmissing");

figure("Units","centimeters", "Position",[5 5 18 5]) 
tiledlayout(1,3, "TileSpacing","compact", "Padding","compact")
nexttile,
errorbar(log10(rsld.muRange),meanRsldCnr,stdRsldCnr,...
    'CapSize',1.5, 'Color',color1)
hold on
plot(log10(rsld.muRange),meanRsldCnr, 'Color',color1)
hold off
xlabel('log_{10}(\mu)')
ylabel('CNR')
grid on
xlim([0 7])
ylim([0,5])
title('RSLD')

nexttile,
imagesc(log10(swtv.muB),log10(swtv.muC),meanSwtvCnr, [0 5])
% c = colorbar;
% c.Label.String = 'NRMSE';
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWTV')
xlim([1.5 5.5])
ylim([-1.5 2.5])
hold on
contour(log10(swtv.muB),log10(swtv.muC),meanSwtvCnr,[1, 10], ...
    '--', 'LineWidth',1.5, "EdgeColor",[0 0 0]);
contour(log10(swtv.muB),log10(swtv.muC),meanSwtvCnr,[2, 10], ...
    '--', 'LineWidth',1.5, "EdgeColor",[0.3 0.3 0.3]);
hold off
xlim([1.5 5.5]-0.5)
ylim([-1.5 2.5]-0.5)


nexttile,
imagesc(log10(swift.muB),log10(swift.muC),meanSwiftCnr, [0 5])
c = colorbar;
c.Label.String = 'CNR';
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWIFT')
hold on
contour(log10(swift.muB),log10(swift.muC),meanSwiftCnr,[1, 10], ...
    '--', 'LineWidth',1.5, "EdgeColor",[0 0 0]);
contour(log10(swift.muB),log10(swift.muC),meanSwiftCnr,[2, 10], ...
    '--', 'LineWidth',1.5, "EdgeColor",[0.3 0.3 0.3]);
hold off
xlim([1.5 5.5])
ylim([-1.5 2.5])


%%
maxCnr = 1;
colors = lines(4);
figure("Units","centimeters", "Position",[5 5 18 5.1]) 
tiledlayout(1,2, "TileSpacing","compact", "Padding","compact")

nexttile,
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWTV')
xlim([1.5 5.5])
ylim([-1.5 2.5])
hold on
for ii = selectedCases
    fprintf('SWTV CNR %i\n',ii)
    contourf(log10(swtv.muB),log10(swtv.muC),swtvCnr(:,:,ii),[maxCnr 100], ...
        lineSpec{ii}, "FaceAlpha",0.1, "FaceColor",colors(ii,:), ...
        "EdgeColor",colors(ii,:), 'LineWidth',1.5)
end
xlim([1.5 5.5]-0.5)
ylim([-1.5 2.5]-0.5)
clim([0.1 0.95])
hold off
ax = gca;
ax.YDir = "reverse";
grid on

nexttile,
xlabel('log_{10}(\mu_B)')
ylabel('log_{10}(\mu_C)')
axis image
title('SWIFT')
hold on
for ii = selectedCases
    fprintf('SWIFT CNR %i\n',ii)
    plot(0,0,lineSpec{ii}, 'LineWidth',1.5)
    contourf(log10(swift.muB),log10(swift.muC),swiftCnr(:,:,ii),[maxCnr 100], ...
        lineSpec{ii}, "FaceAlpha",0.1, "FaceColor",colors(ii,:), ...
        "EdgeColor",colors(ii,:), 'LineWidth',1.5)
end
xlim([1.5 5.5])
ylim([-1.5 2.5])
clim([0.1 0.95])
hold off
ax = gca;
ax.YDir = "reverse";
grid on
legend('S6','','P1','','P2','','P3', 'Location','eastoutside')

%% 
save_all_figures_to_directory('P:\smerino\UFFC2024results\reg','regFinalEstaSi','svg');
close all