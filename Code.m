%% ========================================================================
%  ECOSYSTEM RESPIRATION (ER) PARAMETER ESTIMATION AND GLOBAL EXTRAPOLATION
%  This script processes FLUXNET site data, performs non-linear fitting, 
%  trains/uses ML models, and predicts global ER using CMIP6 climate data.
%  ========================================================================

clc; clear all; 
%% 1. Configuration and Data Initialization
site_n = 369;
monthDays = [31 28 31 30 31 30 31 31 30 31 30 31];
filelist = dir('D:\BaiduSyncdisk\FLUX_DATA\*.csv');
plot_p = 0; % Plot flag for fitting functions

% Pre-allocate storage matrices for site parameters and statistics
para = zeros(10, 3, site_n); 
statistic_data = zeros(10, 5, site_n);
environment_factors = zeros(10, 49, site_n);

% Pre-allocate storage for residual data (sites/years failing the QC)
para_rest = zeros(10, 3, site_n);
environment_factors_rest = zeros(10, 49, site_n);

% Pre-allocate flux data storage
ER_flux_net = zeros(366, 10, site_n);
Ta_flux_net = zeros(366, 10, site_n);
GPP_flux_net = zeros(366, 10, site_n);
ER_flux_net_mean = zeros(366, 10, site_n);
Ta_flux_net_mean = zeros(366, 10, site_n);
GPP_flux_net_mean = zeros(366, 10, site_n);

h = waitbar(0, 'Processing Site Data...');

%% 2. Loop Through FLUXNET Sites
for i = 1:site_n
    waitbar(i/site_n, h, sprintf('Site %i of %i processing...', i, site_n));
    
    % Construct filename and read data
    filename = fullfile(filelist(i).folder, filelist(i).name);
    Daily_data = csvread(filename, 1, 0);
    [num, txt, raw] = xlsread(filename);
    
    % Identify column indices for required variables
    Ernumber = find(ismember(txt, {'RECO_DT_VUT_REF'}) == 1);
    Ernumber_gpp = find(ismember(txt, {'GPP_DT_VUT_REF'}) == 1);
    Ernumber2 = find(ismember(txt, {'SW_IN_F'}) == 1);
    Ernumber3 = find(ismember(txt, {'LW_IN_F'}) == 1);
    Ernumber4 = find(ismember(txt, {'VPD_F'}) == 1);
    Ernumber5 = find(ismember(txt, {'WS_F'}) == 1);
    Ernumber6 = find(ismember(txt, {'P_F'}) == 1);
    Ernumber7 = find(ismember(txt, {'LE_CORR'}) == 1);
    Ernumber8 = find(ismember(txt, {'PA_F'}) == 1);
    Ernumber9 = find(ismember(txt, {'NETRAD'}) == 1);
    Ernumber10 = find(ismember(txt, {'PPFD_IN'}) == 1);
    Ernumber11 = find(ismember(txt, {'CO2_F_MDS'}) == 1);
    Ernumber12 = find(ismember(txt, {'TS_F_MDS_1'}) == 1);
    Ernumber13 = find(ismember(txt, {'SWC_F_MDS_1'}) == 1);
    Ernumber14 = find(ismember(txt, {'G_F_MDS'}) == 1);
    Ernumber15 = find(ismember(txt, {'H_CORR'}) == 1);

    % Assign variables and handle missing data with -9999 placeholder
    Time = Daily_data(:,1);
    Ta = Daily_data(:,2);
    Er = Daily_data(:, Ernumber);
    GPP = Daily_data(:, Ernumber_gpp);
    SW = Daily_data(:, Ernumber2);
    LW = Daily_data(:, Ernumber3);
    VPD = Daily_data(:, Ernumber4);
    WS = Daily_data(:, Ernumber5);
    LE = Daily_data(:, Ernumber7);
    P_F = Daily_data(:, Ernumber6);
    PA = Daily_data(:, Ernumber8);
    
    % Conditional checks for optional variables
    if isempty(Ernumber9), NETRAD = zeros(length(Ta),1)-9999; else, NETRAD = Daily_data(:, Ernumber9); end
    if isempty(Ernumber10), PPFD = zeros(length(Ta),1)-9999; else, PPFD = Daily_data(:, Ernumber10); end
    if isempty(Ernumber11), CO2_F_MDS = zeros(length(Ta),1)-9999; else, CO2_F_MDS = Daily_data(:, Ernumber11); end
    if isempty(Ernumber12), TS_F_MDS_1 = zeros(length(Ta),1)-9999; else, TS_F_MDS_1 = Daily_data(:, Ernumber12); end
    if isempty(Ernumber13), SWC_F_MDS_1 = zeros(length(Ta),1)-9999; else, SWC_F_MDS_1 = Daily_data(:, Ernumber13); end
    if isempty(Ernumber14), G_F_MDS = zeros(length(Ta),1)-9999; else, G_F_MDS = Daily_data(:, Ernumber14); end
    if isempty(Ernumber15), H_CORR = zeros(length(Ta),1)-9999; else, H_CORR = Daily_data(:, Ernumber15); end

    Environment_f = [Ta SW LW VPD WS LE P_F PA NETRAD PPFD CO2_F_MDS TS_F_MDS_1 SWC_F_MDS_1 G_F_MDS H_CORR];
    Time_Ta_Er_Gpp = [Time Ta Er GPP];
    
    % Calculate Year details
    Year_amount = fix(length(Time)/365);
    First_year = fix(Time(1)/10000);
    last_day = 1;
    j = 1;

    %% 3. Annual Processing Loop
    for year = First_year:(First_year + Year_amount - 1)
        % Determine leap year and year length
        is_leap = (mod(year, 4) == 0 && mod(year, 100) ~= 0) || (mod(year, 400) == 0);
        days_in_yr = 365 + is_leap;
        
        % Ensure index doesn't exceed data length
        if (last_day + days_in_yr - 1) > size(Time_Ta_Er_Gpp, 1), break; end
        
        tpye_year = Time_Ta_Er_Gpp(last_day:(last_day + days_in_yr - 1), :);
        tpye_year_env = Environment_f(last_day:(last_day + days_in_yr - 1), :);
        
        % Clean environmental data
        tpye_year_env(tpye_year_env < -9000) = NaN;
        env_mean = nanmean(tpye_year_env, 1);
        env_max = max(tpye_year_env, [], 1);
        env_min = min(tpye_year_env, [], 1);
        ann_precip = sum(P_F(last_day:(last_day + days_in_yr - 1)));
        
        % Monthly GPP Analysis
        curr_GPP = tpye_year(:, 4);
        curr_GPP(curr_GPP == -9999) = NaN;
        Gppp_month = zeros(1,12);
        d_ptr = 1;
        for m = 1:12
            Gppp_month(m) = nanmean(curr_GPP(d_ptr : min(d_ptr + monthDays(m)-1, length(curr_GPP))));
            d_ptr = d_ptr + monthDays(m);
        end
        
        envir_all = [env_mean, mean(Gppp_month), env_max, max(Gppp_month), env_min, ann_precip, min(Gppp_month)];
        
        % Filter invalid ER data for binning
        Taa = tpye_year(:, 2); Err = tpye_year(:, 3);
        invalid_mask = (Err == -9999);
        Taa(invalid_mask) = []; Err(invalid_mask) = [];
        
        % Create Bins for Temperature-Response curve
        if ~isempty(Taa)
            minT = floor(min(Taa)); maxT = ceil(max(Taa));
            bins = maxT - minT - 2;
            movingT = []; movingER = []; movingGPP = [];
            
            for u = 1:bins
                T_target = minT + (u-1);
                pts = find(Taa <= (T_target+3) & Taa >= T_target);
                if ~isempty(pts)
                    movingT(u) = mean(Taa(pts));
                    movingER(u) = mean(Err(pts));
                    movingGPP(u) = mean(curr_GPP(pts));
                end
            end
            
            % Curve Fitting
            addpath 'D:\BaiduSyncdisk\F2015_JIEYASUO'
            try
                [reg1, reg2] = createFit_gpp(movingT, movingER, plot_p);
                
                % Quality Control (QC): R2 > 0.4 and realistic Peak Temperature
                if reg2.rsquare >= 0.4 && reg1.c < 1.2 * max(movingT) && reg1.c > -5 && reg2.dfe > 3
                    para(j, :, i) = [reg1.b, reg1.c, reg1.d];
                    statistic_data(j, :, i) = [reg2.sse, reg2.rsquare, reg2.dfe, reg2.adjrsquare, reg2.rmse];
                    environment_factors(j, :, i) = envir_all;
                else
                    para_rest(j, :, i) = [reg1.b, reg1.c, reg1.d];
                    environment_factors_rest(j, :, i) = envir_all;
                end
            catch
                % Fit failed
            end
        end
        
        last_day = last_day + days_in_yr;
        j = j + 1;
    end
    clear num txt raw Daily_data Time_Ta_Er_Gpp Environment_f
end
close(h);

%% 4. Data Aggregation and Smoothing
% Consolidate parameters into 2D arrays for ML training
% [Code for flattening site-year-parameter matrices omitted for brevity but logic follows reshaping para(j,:,i)]

%% 5. Global Extrapolation using CMIP6 Data
% Goal: Use site-trained parameters to estimate global ER from 2015-2100 (SSP585)
filelist_cmip = dir('F:\CMIP6 data\CLM\585\*.nc');
lati = 288; longti = 192;
para_b1_86 = zeros(lati, longti, 9, 10);

% Loop through CMIP6 decade files
for year_n = 1:9
    for year10 = 1:10
        if year_n == 9 && year10 > 6, continue; end % SSP585 limit
        
        % Load variables (Humidity, Precip, Radiation, Wind, Temperature)
        % Note: ncread logic remains specific to your local file indexing
        % Calculate VPD globally
        % VPD = 0.611 * exp(17.27 * (Ta - 273.15) ./ (Ta - 273.15 + 237.3)) .* (1 - Hurs/100) * 10;
        
        % Prepare grid-cell features for ML Model prediction
        % Features: means, max, min of T, SW, LW, VPD, WS, Precip, GPP
        
        parfor r = 1:lati
            for c = 1:longti
                % Extract local time series and calculate stats
                % Predict Topt (Optimal Temperature) using trained model
                % map_b1(r, c) = trainedModel_1318_b_e.predictFcn(temp_envir);
            end
        end
        % Save spatial maps for statistical analysis
    end
end

%% 6. Spatial Mapping and Plotting
% Use m_map toolbox for visualization
% [m_proj, m_pcolor, m_coast etc. logic as per your script]

%% 7. Statistical Analysis by Climate Zones (Koppen-Geiger)
% Load Koppen-Geiger classification and calculate zonal mean/std
% Compare Topt and ERmax differences between periods (e.g., late vs early century)

fprintf('Processing complete.\n');