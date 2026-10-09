clc;
clear;
close all;

%% 1. Settings & Global Formatting
dataFolder = 'D:\DATA\IMD DATA\RAINFALL ( RESL. 0.25-0.25)';
alpha = 0.05;

% --- APPLY GLOBAL PLOT FORMATTING ---
% Forces all plots to have a white background, Times New Roman font, size 10, in black.
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesColor', 'W');
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 10);
set(groot, 'defaultAxesXColor', 'k');
set(groot, 'defaultAxesYColor', 'k');
set(groot, 'defaultAxesZColor', 'k');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultTextFontSize', 10);
set(groot, 'defaultTextColor', 'k');
set(groot, 'defaultLegendFontName', 'Times New Roman');
set(groot, 'defaultLegendFontSize', 10);
set(groot, 'defaultLegendTextColor', 'W');
set(groot, 'defaultColorbarFontName', 'Times New Roman');
set(groot, 'defaultColorbarFontSize', 10);
set(groot, 'defaultColorbarColor', 'k');

%% Output folder
outputFolder = 'D:\matlad\results\50 ANNUAL RF TA';
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end
fprintf('Output folder: %s\n', outputFolder);

%% 2. Define six target points
pointNames = {'POINT-1'; 'POINT-2'; 'POINT-3'; 'POINT-4'; 'POINT-5'; 'POINT-6'};
pointLon = [70.50; 70.75; 70.25; 70.50; 70.75; 70.25];
pointLat = [21.50; 21.50; 21.50; 21.75; 21.75; 21.75];
numPoints = length(pointNames);

%% 3. Locate NetCDF files
if ~isfolder(dataFolder)
    error('Folder not found: %s', dataFolder);
end
allFiles = dir(fullfile(dataFolder, '**', '*.nc'));
if isempty(allFiles)
    error('No NetCDF files found.');
end
fileNames = {allFiles.name};
keep = ~cellfun('isempty', regexp(fileNames, '^RF25_ind\d{4}_rfp25\.nc$', 'once'));
allFiles = allFiles(keep);

%% 4. Extract years & Keep one file per year
fileYears = NaN(length(allFiles), 1);
for k = 1:length(allFiles)
    token = regexp(allFiles(k).name, 'RF25_ind(\d{4})_rfp25\.nc', 'tokens', 'once');
    fileYears(k) = str2double(token{1});
end
[fileYears, uniqueIndex] = unique(fileYears, 'stable');
allFiles = allFiles(uniqueIndex);
[fileYears, sortIndex] = sort(fileYears);
allFiles = allFiles(sortIndex);
nYears = length(fileYears);

%% 5. POP-UP: Ask user for number of years to analyze
prompt = {sprintf('Total available years: %d\n\nEnter the number of years to use for trend analysis (most recent years will be selected):', nYears)};
dlgtitle = 'Data Subset Selection';
dims = [1 60];
definput = {num2str(nYears)}; % Default is all available years
answer = inputdlg(prompt, dlgtitle, dims, definput);

% Handle user cancellation
if isempty(answer)
    error('Execution cancelled by the user.');
end

userYears = round(str2double(answer{1}));

% Validate input
if isnan(userYears) || userYears <= 0
    error('Invalid input. Please enter a positive number.');
end

% Subset the data if the user requested fewer years than available
if userYears < nYears
    fprintf('\nSubsetting data to the most recent %d years...\n', userYears);
    allFiles = allFiles(end-userYears+1:end);
    fileYears = fileYears(end-userYears+1:end);
    nYears = length(fileYears);
else
    fprintf('\nUsing all available %d years (or user entered a number larger than available).\n', nYears);
end

%% Check continuous yearly coverage
expectedYears = (fileYears(1):fileYears(end))';
missingYears = setdiff(expectedYears, fileYears);
if isempty(missingYears)
    fprintf('Year continuity check: PASS (%d to %d)\n', ...
        fileYears(1), fileYears(end));
else
    warning('Missing years: %s', mat2str(missingYears'));
end
if nYears ~= length(expectedYears)
   warning(['Number of files differs from expected continuous years. ', ...
        'Check duplicate or missing annual files.']);
end

% Data Set Summary
startDateStr = sprintf('0101%04d', fileYears(1));
endDateStr = sprintf('3112%04d', fileYears(end));
fprintf('\n--- Data Set Summary ---\n');
fprintf('Files used: %d\n', nYears);
fprintf('Start date: %s\n', startDateStr);
fprintf('End date: %s\n', endDateStr);
fprintf('------------------------\n');

%% 6. Read spatial grid
firstFile = fullfile(allFiles(1).folder, allFiles(1).name);
info = ncinfo(firstFile);
variableNames = {info.Variables.Name};
if ismember('lon', variableNames)
    lon = double(ncread(firstFile, 'lon'));
    lat = double(ncread(firstFile, 'lat'));
else
    lon = double(ncread(firstFile, 'LONGITUDE'));
    lat = double(ncread(firstFile, 'LATITUDE'));
end

%% 7. Find nearest grid cell
lonIndex = zeros(numPoints, 1);
latIndex = zeros(numPoints, 1);
selectedLon = zeros(numPoints, 1);
selectedLat = zeros(numPoints, 1);
for p = 1:numPoints
    [~, latIndex(p)] = min(abs(lat - pointLat(p)));
    [~, lonIndex(p)] = min(abs(lon - pointLon(p)));
    selectedLat(p) = lat(latIndex(p));
    selectedLon(p) = lon(lonIndex(p));
end

%% 8. Extract Annual Rainfall ONLY for the 6 Points
pointAnnualRain = NaN(nYears, numPoints);
fprintf('\nCalculating annual rainfall...\n');
for k = 1:nYears
    currentFile = fullfile(allFiles(k).folder, allFiles(k).name);
    info = ncinfo(currentFile);
    variableNames = {info.Variables.Name};
    if ismember('rf', variableNames)
        rain = double(ncread(currentFile, 'rf'));
    else
        rain = double(ncread(currentFile, 'RAINFALL'));
    end
    pointAnnualRain(k, :) = NaN; % Initialize for the current year
    if size(rain,1) ~= length(lon) || size(rain,2) ~= length(lat)
        error('Rainfall grid size differs in file: %s', allFiles(k).name);
    end
    for p = 1:numPoints
        dailyRain = squeeze(rain(lonIndex(p), latIndex(p), :));
        dailyRain(dailyRain == -999) = NaN; 
        if sum(~isnan(dailyRain)) > 0
            pointAnnualRain(k, p) = sum(dailyRain, 'omitnan');
        end
    end
end

%% 9. Apply Trend and Homogeneity Tests
% Preallocate metric arrays
pointMKZ = NaN(numPoints, 1); pointMKP = NaN(numPoints, 1);
pointMKH = NaN(numPoints, 1); pointSenSlope = NaN(numPoints, 1);
pointSRC_P = NaN(numPoints, 1); pointSRC_Rho = NaN(numPoints, 1);
pointLR_Slope = NaN(numPoints, 1); pointLR_P = NaN(numPoints, 1);
pointKS_P = NaN(numPoints, 1); pointWilcoxon_P = NaN(numPoints, 1);
pointPettitt_CP = NaN(numPoints, 1); pointPettitt_P = NaN(numPoints, 1);
pointTTest_P = NaN(numPoints, 1); pointTFPW_MK_P = NaN(numPoints, 1);

% Preallocate time-series arrays for plotting
all_sqmk_u = NaN(nYears, numPoints);
all_sqmk_uprime = NaN(nYears, numPoints);
all_cusum = NaN(nYears, numPoints);

% Prepare half-splits for ITA
halfLength = floor(nYears / 2);
ita_x = NaN(halfLength, numPoints);
ita_y = NaN(halfLength, numPoints);

for p = 1:numPoints
    y = pointAnnualRain(:, p);
    t = (1:nYears)';
    valid = ~isnan(y);
    
    if sum(valid) >= 3
        y_valid = y(valid);
        t_valid = t(valid);
        
        % Split valid data in half for tests requiring it
        half_idx = floor(length(y_valid) / 2);
        y1 = y_valid(1:half_idx);
        y2 = y_valid(end-half_idx+1:end); 
        
        % Mann-Kendall & Sen's Slope
        [h, p_val, S, Z] = mk_test_corrected(y_valid, alpha);
        pointMKH(p) = h; pointMKP(p) = p_val; pointMKZ(p) = Z;
        pointSenSlope(p) = sens_slope_corrected(y_valid);
        
        % Spearman Rank Correlation
        [rho, src_p] = corr(t_valid, y_valid, 'Type', 'Spearman');
        pointSRC_Rho(p) = rho; pointSRC_P(p) = src_p;
        
        % Linear Regression
        mdl = fitlm(t_valid, y_valid);
        pointLR_Slope(p) = mdl.Coefficients.Estimate(2);
        pointLR_P(p) = mdl.Coefficients.pValue(2);
        
        % Homogeneity Tests
        [~, pointKS_P(p)] = kstest2(y1, y2);
        pointWilcoxon_P(p) = ranksum(y1, y2);
        [~, pointTTest_P(p)] = ttest2(y1, y2);
        
        % Pettitt's Change-Point Test
        [cp_idx, p_val_pettitt] = pettitt_test(y_valid);
        valid_years = fileYears(valid); 
        pointPettitt_CP(p) = valid_years(cp_idx); 
        pointPettitt_P(p) = p_val_pettitt;
        
        % Trend-Free Pre-Whitening (TFPW) followed by MK Test
        y_tfpw = tfpw_filter(y_valid);
        [~, pointTFPW_MK_P(p), ~, ~] = mk_test_corrected(y_tfpw, alpha);
        
        % Sequential MK (store full series)
        [u_sqmk, uprime_sqmk] = sqmk_test(y_valid);
        all_sqmk_u(valid, p) = u_sqmk;
        all_sqmk_uprime(valid, p) = uprime_sqmk;
        
        % CUSUM (store full series)
        all_cusum(valid, p) = cusum_test(y_valid);
        
        % Innovative Trend Analysis using valid annual rainfall values
        itaHalfLength = floor(length(y_valid) / 2);
        ita_y1 = y_valid(1:itaHalfLength);
        ita_y2 = y_valid(end-itaHalfLength+1:end);
        ita_x(1:itaHalfLength, p) = sort(ita_y1);
        ita_y(1:itaHalfLength, p) = sort(ita_y2);
    end
end

%% 10. Display Console Output
fprintf('\n====================================================\n');
fprintf('TREND RESULTS FOR SIX POINTS\n');
fprintf('====================================================\n');
for p = 1:numPoints
    if pointMKH(p) == 1
        mkDecision = 'Significant trend at 5% level';
    elseif pointMKH(p) == 0
        mkDecision = 'Not significant at 5% level';
    else
        mkDecision = 'Not calculated';
    end
    if pointMKZ(p) > 0
        trendDirection = 'Increasing';
    elseif pointMKZ(p) < 0
        trendDirection = 'Decreasing';
    else
        trendDirection = 'No direction';
    end
    fprintf('\n%s\n', pointNames{p});
    fprintf('Selected grid: Latitude = %.4f, Longitude = %.4f\n', ...
        selectedLat(p), selectedLon(p));
    fprintf('MK Z-statistic: %.4f\n', pointMKZ(p));
    fprintf('MK p-value: %.4f\n', pointMKP(p));
    fprintf('MKH (0 = no significant trend, 1 = significant): %d\n', ...
        pointMKH(p));
    fprintf('Trend direction: %s\n', trendDirection);
    fprintf('MK conclusion: %s\n', mkDecision);
    fprintf('Sen''s slope: %.4f mm/year\n', pointSenSlope(p));
    fprintf('Spearman rho: %.4f | p-value: %.4f\n', ...
        pointSRC_Rho(p), pointSRC_P(p));
    fprintf('Linear-regression slope: %.4f mm/year | p-value: %.4f\n', ...
        pointLR_Slope(p), pointLR_P(p));
    fprintf('Pettitt change year: %.0f | p-value: %.4f\n', ...
        pointPettitt_CP(p), pointPettitt_P(p));
    fprintf('KS homogeneity p-value: %.4f\n', pointKS_P(p));
    fprintf('Wilcoxon rank-sum p-value: %.4f\n', pointWilcoxon_P(p));
    fprintf('Two-sample t-test p-value: %.4f\n', pointTTest_P(p));
    fprintf('TFPW-MK p-value: %.4f\n', pointTFPW_MK_P(p));
end
fprintf('\n====================================================\n');

%% 11. Plot Maps & Charts
% MAP PLOTS
figure('Name', 'Spatial Trends');
lonLimits = [min(selectedLon) - 0.25, max(selectedLon) + 0.25];
latLimits = [min(selectedLat) - 0.25, max(selectedLat) + 0.25];
subplot(2,1,1);
scatter(selectedLon, selectedLat, 1200, pointSenSlope, 'square', 'filled', 'MarkerEdgeColor', 'k');
colorbar; colormap(subplot(2,1,1), 'parula');
title('Sen''s Slope of Annual Rainfall'); xlim(lonLimits); ylim(latLimits); grid on;
subplot(2,1,2);
scatter(selectedLon, selectedLat, 1200, pointMKH, 'square', 'filled', 'MarkerEdgeColor', 'k');
colorbar; colormap(subplot(2,1,2), [0 0.5 1; 1 1 0]); clim([0 1]);
title('Significant MK Trends (5% Level)'); xlim(lonLimits); ylim(latLimits); grid on;

% --- SQMK PLOTS FOR ALL 6 POINTS ---
figure('Name', 'Sequential Mann-Kendall Test', 'Position', [100, 100, 1200, 700]);
for p = 1:numPoints
    subplot(2, 3, p);
    valid_p = ~isnan(pointAnnualRain(:, p));
    plot(fileYears(valid_p), all_sqmk_u(valid_p, p), 'b-', 'LineWidth', 1.5); hold on;
    plot(fileYears(valid_p), all_sqmk_uprime(valid_p, p), 'r--', 'LineWidth', 1.5);
    yline(1.96, 'k:'); yline(-1.96, 'k:'); % 95% confidence bounds
    title(sprintf('SQMK: %s', pointNames{p}));
    xlabel('Year'); ylabel('Test Statistic (u)');
    grid on;
    if p == 1
        legend('Progressive u(t)', 'Retrograde u''(t)', '95% Confidence Level', 'Location', 'best');
    end
end

% --- CUSUM PLOTS FOR ALL 6 POINTS ---
figure('Name', 'Cumulative Sum (CUSUM)', 'Position', [150, 150, 1200, 700]);
for p = 1:numPoints
    subplot(2, 3, p);
    valid_p = ~isnan(pointAnnualRain(:, p));
    plot(fileYears(valid_p), all_cusum(valid_p, p), 'k-', 'LineWidth', 1.5);
    title(sprintf('CUSUM: %s', pointNames{p}));
    xlabel('Year'); ylabel('CUSUM');
    grid on;
end

% ITA PLOTS
figure('Name', 'Innovative Trend Analysis (ITA)', 'Position', [200, 200, 1200, 700]);
for p = 1:numPoints
    subplot(2, 3, p);
    scatter(ita_x(:, p), ita_y(:, p), 30, 'b', 'filled'); hold on;
    maxVal = max([ita_x(:, p); ita_y(:, p)]) * 1.1;
    if isnan(maxVal); maxVal = 1; end % Fallback for empty data
    plot([0, maxVal], [0, maxVal], 'k-', 'LineWidth', 1.5);
    plot([0, maxVal], [0, maxVal*1.1], 'r--', 'LineWidth', 1);
    plot([0, maxVal], [0, maxVal*0.9], 'r--', 'LineWidth', 1);
    title(sprintf('ITA: %s', pointNames{p}));
    xlabel('First Half Sorted (mm)'); ylabel('Second Half Sorted (mm)');
    xlim([0 maxVal]); ylim([0 maxVal]); grid on;
end

%% 12. Save Results
% One combined table: all points, all years
allAnnualTable = table(fileYears, ...
    pointAnnualRain(:,1), ...
    pointAnnualRain(:,2), ...
    pointAnnualRain(:,3), ...
    pointAnnualRain(:,4), ...
    pointAnnualRain(:,5), ...
    pointAnnualRain(:,6), ...
    'VariableNames', { ...
        'Year', ...
        'POINT_1_AnnualRainfall_mm', ...
        'POINT_2_AnnualRainfall_mm', ...
        'POINT_3_AnnualRainfall_mm', ...
        'POINT_4_AnnualRainfall_mm', ...
        'POINT_5_AnnualRainfall_mm', ...
        'POINT_6_AnnualRainfall_mm'});
writetable(allAnnualTable, ...
    fullfile(outputFolder, 'All_Points_Annual_Rainfall.csv'));

% One combined table: trend and homogeneity results for all points
trendSummaryTable = table( ...
    pointNames, ...
    pointLat, ...
    pointLon, ...
    selectedLat, ...
    selectedLon, ...
    pointMKZ, ...
    pointMKP, ...
    pointMKH, ...
    pointSenSlope, ...
    pointSRC_Rho, ...
    pointSRC_P, ...
    pointLR_Slope, ...
    pointLR_P, ...
    pointKS_P, ...
    pointWilcoxon_P, ...
    pointPettitt_CP, ...
    pointPettitt_P, ...
    pointTTest_P, ...
    pointTFPW_MK_P, ...
    'VariableNames', { ...
        'Point', ...
        'RequestedLatitude', ...
        'RequestedLongitude', ...
        'SelectedLatitude', ...
        'SelectedLongitude', ...
        'MK_Z', ...
        'MK_P', ...
        'MKH_Significant_At_5Percent', ...
        'Sen_Slope_mm_per_year', ...
        'Spearman_Rho', ...
        'Spearman_P', ...
        'Linear_Slope_mm_per_year', ...
        'Linear_P', ...
        'KS_Homogeneity_P', ...
        'Wilcoxon_P', ...
        'Pettitt_Change_Year', ...
        'Pettitt_P', ...
        'TTest_P', ...
        'TFPW_MK_P'});
writetable(trendSummaryTable, ...
    fullfile(outputFolder, 'All_Points_Trend_Summary.csv'));

% Separate outputs for each point
for p = 1:numPoints
    annualTable = table( ...
        fileYears, ...
        pointAnnualRain(:, p), ...
        'VariableNames', { ...
            'Year', ...
            'AnnualRainfall_mm'});
    annualFile = fullfile(outputFolder, ...
        sprintf('%s_Annual_Rainfall.csv', pointNames{p}));
    writetable(annualTable, annualFile);
    
    trendTable = table( ...
        pointNames(p), ...
        pointLat(p), ...
        pointLon(p), ...
        selectedLat(p), ...
        selectedLon(p), ...
        pointMKZ(p), ...
        pointMKP(p), ...
        pointMKH(p), ...
        pointSenSlope(p), ...
        pointSRC_Rho(p), ...
        pointSRC_P(p), ...
        pointLR_Slope(p), ...
        pointLR_P(p), ...
        pointKS_P(p), ...
        pointWilcoxon_P(p), ...
        pointPettitt_CP(p), ...
        pointPettitt_P(p), ...
        pointTTest_P(p), ...
        pointTFPW_MK_P(p), ...
        'VariableNames', { ...
            'Point', ...
            'RequestedLatitude', ...
            'RequestedLongitude', ...
            'SelectedLatitude', ...
            'SelectedLongitude', ...
            'MK_Z', ...
            'MK_P', ...
            'MKH_Significant_At_5Percent', ...
            'Sen_Slope_mm_per_year', ...
            'Spearman_Rho', ...
            'Spearman_P', ...
            'Linear_Slope_mm_per_year', ...
            'Linear_P', ...
            'KS_Homogeneity_P', ...
            'Wilcoxon_P', ...
            'Pettitt_Change_Year', ...
            'Pettitt_P', ...
            'TTest_P', ...
            'TFPW_MK_P'});
    trendFile = fullfile(outputFolder, ...
        sprintf('%s_Trend_Summary.csv', pointNames{p}));
    writetable(trendTable, trendFile);
end

% Save MATLAB variables for later use
save(fullfile(outputFolder, ...
    'Rainfall_Trend_AllPoints_Results.mat'), ...
    'fileYears', 'pointNames', 'pointLon', 'pointLat', ...
    'selectedLon', 'selectedLat', 'pointAnnualRain', ...
    'pointMKZ', 'pointMKP', 'pointMKH', 'pointSenSlope', ...
    'pointSRC_Rho', 'pointSRC_P', ...
    'pointLR_Slope', 'pointLR_P', ...
    'pointKS_P', 'pointWilcoxon_P', ...
    'pointPettitt_CP', 'pointPettitt_P', ...
    'pointTTest_P', 'pointTFPW_MK_P');
fprintf('\nResults successfully saved in:\n%s\n', outputFolder);
fprintf('Files created:\n');
fprintf('1. All_Points_Annual_Rainfall.csv\n');
fprintf('2. All_Points_Trend_Summary.csv\n');
fprintf('3. POINT-1 to POINT-6 annual rainfall CSV files\n');
fprintf('4. POINT-1 to POINT-6 trend summary CSV files\n');
fprintf('5. Rainfall_Trend_AllPoints_Results.mat\n');

%% Helper Functions
function [h, p, S, Z] = mk_test_corrected(x, alpha)
    x = x(:); x = x(~isnan(x)); n = length(x);
    if n < 3; h = NaN; p = NaN; S = NaN; Z = NaN; return; end
    S = 0;
    for i = 1:n-1; for j = i+1:n; S = S + sign(x(j) - x(i)); end; end
    uniqueValues = unique(x); tieCorrection = 0;
    for k = 1:length(uniqueValues)
        tieCount = sum(x == uniqueValues(k));
        if tieCount > 1; tieCorrection = tieCorrection + tieCount * (tieCount - 1) * (2*tieCount + 5); end
    end
    varianceS = (n*(n-1)*(2*n+5) - tieCorrection) / 18;
    if varianceS == 0; Z = 0; p = 1; h = 0; return; end
    if S > 0; Z = (S - 1) / sqrt(varianceS);
    elseif S < 0; Z = (S + 1) / sqrt(varianceS);
    else; Z = 0; end
    p = erfc(abs(Z) / sqrt(2)); h = p < alpha;
end

function slope = sens_slope_corrected(x)
    x = x(:); x = x(~isnan(x)); n = length(x);
    if n < 2; slope = NaN; return; end
    slopes = NaN(n*(n-1)/2, 1); counter = 0;
    for i = 1:n-1; for j = i+1:n; counter = counter + 1; slopes(counter) = (x(j) - x(i)) / (j - i); end; end
    slope = median(slopes, 'omitnan');
end

function [u, uprime] = sqmk_test(x)
    x = x(:); 
    n = length(x); 
    u = zeros(n, 1); 
    uprime = zeros(n, 1);
    
    % Progressive series u(t)
    t_sum = 0; 
    for i = 2:n
        P = 0;
        for j = 1:i-1
            if x(i) > x(j)
                P = P + 1; 
            end
        end
        t_sum = t_sum + P; 
        
        E_t = i * (i - 1) / 4;
        Var_t = i * (i - 1) * (2 * i + 5) / 72;
        
        if Var_t > 0
            u(i) = (t_sum - E_t) / sqrt(Var_t); 
        else
            u(i) = 0; 
        end
    end
    
    % Retrograde series u'(t)
    x_rev = flipud(x);
    t_sum_rev = 0; 
    for i = 2:n
        P = 0;
        for j = 1:i-1
            if x_rev(i) > x_rev(j)
                P = P + 1; 
            end
        end
        t_sum_rev = t_sum_rev + P; 
        
        E_t = i * (i - 1) / 4;
        Var_t = i * (i - 1) * (2 * i + 5) / 72;
        
        if Var_t > 0
            uprime(n - i + 1) = -(t_sum_rev - E_t) / sqrt(Var_t); 
        else
            uprime(n - i + 1) = 0; 
        end
    end
    uprime(1) = 0; 
end

function [cp_year_idx, p_value] = pettitt_test(x)
    x = x(:); n = length(x); U = zeros(n-1, 1);
    for t = 1:n-1
        sum_sign = 0;
        for i = 1:t; for j = t+1:n; sum_sign = sum_sign + sign(x(i) - x(j)); end; end
        U(t) = sum_sign;
    end
    [K, cp_year_idx] = max(abs(U));
    p_value = 2 * exp((-6 * K^2) / (n^3 + n^2));
    p_value = min(p_value, 1);
end

function x_pw = tfpw_filter(x)
    x = x(:); n = length(x); t = (1:n)';
    beta = sens_slope_corrected(x);
    Y = x - beta .* t;
    r1 = corr(Y(1:end-1), Y(2:end), 'Type', 'Pearson');
    if abs(r1) > 0.1 
        Y_prime = zeros(n, 1); Y_prime(1) = Y(1);
        for i = 2:n; Y_prime(i) = Y(i) - r1 * Y(i-1); end
        x_pw = Y_prime + beta .* t;
    else
        x_pw = x;
    end
end

function cusum_stat = cusum_test(x)
    x = x(:); x_mean = mean(x, 'omitnan'); cusum_stat = cumsum(x - x_mean);
end