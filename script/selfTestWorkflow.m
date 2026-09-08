function selfTestWorkflow
%SELFTESTWORKFLOW Run a small end-to-end test of the site-year OOF workflow.
rng(3); numberOfSites=6; numberOfDays=365; n=numberOfSites*numberOfDays;
dates=datetime(2019,1,1)+caldays(0:numberOfDays-1);
daily=table;
daily.site_id=repelem("S"+string((1:numberOfSites)'),numberOfDays);
allDates=repmat(dates',numberOfSites,1);
daily.TIMESTAMP=str2double(string(allDates,'yyyyMMdd'));
season=repmat(12+13*sin(2*pi*(1:numberOfDays)'/365),numberOfSites,1);
daily.TA_F=season+randn(n,1);
daily.VPD_F=max(0,8+0.4*daily.TA_F+randn(n,1));
daily.P_F=max(0,randn(n,1)+1);
daily.GPP_DT_VUT_REF=max(0.1,5+0.15*daily.TA_F+randn(n,1));
x=daily.TA_F-12;
daily.RECO_DT_VUT_REF=exp(1+0.08*x-0.004*x.^2).*(1+0.05*randn(n,1));
daily.location_lat=repelem(linspace(30,55,numberOfSites)',numberOfDays);
daily.location_long=repelem(linspace(-120,-70,numberOfSites)',numberOfDays);
daily.igbp=repmat("GRA",n,1); daily.data_hub=repmat("TEST",n,1);
temporaryDirectory=string(tempname); mkdir(temporaryDirectory);
cleanup=onCleanup(@() rmdir(temporaryDirectory,'s'));
dailyFile=fullfile(temporaryDirectory,'daily.csv'); writetable(daily,dailyFile);
results=refitPredictBetaSimulateER(dailyFile,fullfile(temporaryDirectory,'output'),NumTrees=5,Seed=9);
assert(height(results.fits)==numberOfSites,'Unexpected number of fitted site-years.');
assert(all(isfinite(results.oofPredictions.predicted_beta0)),'OOF beta prediction failed.');
assert(height(results.annualER)==numberOfSites,'OOF annual ER aggregation failed.');
modelDirectory=fullfile(temporaryDirectory,'model');
bundle=trainGlobalBetaModel(dailyFile, ...
    fullfile(temporaryDirectory,'output','all_site_year_beta_fits.csv'), ...
    modelDirectory,NumTrees=5,Seed=13);
latitude=[35;45]; longitude=[-110;-90;-70];
spatialOffset=reshape([-2 0 2 1 3 5],1,2,3);
taDaily=reshape(12+13*sin(2*pi*(1:365)'/365),365,1,1)+spatialOffset;
vpdDaily=max(0,8+0.4*taDaily);
precipDaily=ones(365,2,3);
gppMonthly=ones(12,2,3)*5;
forcing=struct('taDaily',taDaily,'vpdDaily',vpdDaily, ...
    'precipDaily',precipDaily,'gppMonthly',gppMonthly);
projection=projectDiermYear(forcing,bundle,latitude,longitude);
assert(isequal(size(projection.hot_days),[2 3]),'Projection grid dimensions are incorrect.');
assert(all(isfinite(projection.hot_days),'all'),'Synthetic projection produced missing HOT days.');
fprintf('DIERM end-to-end synthetic workflow passed.\n');
end
