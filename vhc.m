function [lagData, timeLags, counts, edges] = ...
    vhc(data, fs, dT, preview, edgeMin, edgeMax, binCount)

%Check inputs
narginchk(2,7);

%Data check
if (~isvector(data) || isscalar(data))
    error("VHC input data is not a vector.");
end

%Sampling frequency check
if(~isscalar(fs))
    error("Sampling frequency fs must be a scalar.");    
end
frameT = 1/fs;

%Time lags, if dT exists and is a vector tranaslate from seconds to samples
%with rounding, otherwise use the default
if isempty(dT)
    sampleLags = 1:1:floor(length(data)/2);
else
    if ~isvector(dT)
        warning("dT is not a vector, using default time lags");
        sampleLags = 1:1:floor(length(data)/2);
    else
        %convert time lags to number of samples
        sampleLags = round(dT*fs);
        %check and remove any time lags <= 0
        len = length(sampleLags);
        sampleLags = sampleLags(sampleLags>0);
        if(len > length(sampleLags))
            warning("All time lags <= 0 are not being considered.");
        end
        %check and remove any time lags greater than the size of the data
        %set
        len = length(sampleLags);
        sampleLags = sampleLags(sampleLags<(length(data)-1));
        if(len > length(sampleLags))
            warning("All time lags greater than the length of the data set are not being considered.");
        end
    end
end
timeLags = sampleLags/fs;

%Find time lag differences
numLags = length(sampleLags);
minLag = sampleLags(1);
lagData = NaN([numLags (length(data)-minLag)]);

%Progress string set up
fprintf("Calculating time lag data: 00%%");
lastprint = 0;

%Calculation of time lag data
for vIndex = 1:numLags
    lag = sampleLags(vIndex);
    for hIndex = 1:(length(data)-lag)
        lagData(vIndex, hIndex) = data(hIndex) - data(hIndex+lag);
    end
    
    %Progress string checks
    if(floor(vIndex/numLags*100)) > lastprint
        lastprint = floor(vIndex/numLags*100);
        fprintf("\b\b\b%02d%%", lastprint);
    end
end
%Progress string finish
fprintf("\b\b\b\bDone\n");

%Find mean and std for default bin range
sz = size(lagData);
m = 0;
s = 0;
numEl = 0;
for vIndex = 1:sz(1)
    for hIndex = 1:sz(2)
        temp = lagData(vIndex, hIndex);
        if(~isnan(temp))
            m = m + temp;
            s = s + temp^2;
            numEl = numEl + 1;
        end
    end
end
m = m/numEl;
s = s/numEl;
s = s - m^2;
s = sqrt(abs(s));

%check args for bin edges, else use default of 
%50 bins ranging from mean-4*std to mean+4*std
if nargin > 3 && ~isempty(preview) && islogical(preview)
    printPreview = lower(preview);
else
    printPreview = false;
end

if nargin <= 4 || isempty(edgeMin)
    edgeMin = m-4*s;
end

if nargin <= 5 || isempty(edgeMax)
    edgeMax = m+4*s;
end

if nargin > 6 && ~isempty(binCount) && isscalar(binCount)
    bins = floor(binCount);
else
    bins = 50;
end

edges = edgeMin:(edgeMax-edgeMin)/bins:edgeMax;
counts = zeros(numLags, length(edges)-1);

%progress string setup
lastprint = 0;
fprintf("Calculating distributions: 00%%");
for vIndex = 1:numLags
    %The probability normalization of histcounts is roughly 
    %valI = countI/numel
    %which deflates later timelags as they are padded with nans
    %we fix this by passing only the non-nan values into histcounts
    [counts(vIndex,:) ~] = histcounts(lagData(vIndex,1:(end - sampleLags(vIndex)+minLag)), ...
                                edges, 'Normalization', 'Probability');
    
    %Check if we need to update the progress string
    if(floor(vIndex/numLags*100)) > lastprint
        lastprint = floor(vIndex/numLags*100);
        fprintf("\b\b\b%02d%%", lastprint);
    end
end
%finish the progress string
fprintf("\b\b\b\bDone\n");

if printPreview
    figure
    surf(edges(1:end-1), timeLags, counts);
    xlabel("d(t)-d(t+lag)");
    ylabel("Lag (s)");
    zlabel("Probability");
    shading interp;
end
return
end


