function [timeLags, counts, edges] = ...
    microvhc(data, fs, dT, preview, edgeMin, edgeMax, binCount)

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

%Time lags, if dT exists and is a vector tranaslate from seconds to samples
%with rounding, otherwise use the default
if isempty(dT)
    sampleLags = 1:1:floor(length(data)/2);
elseif ~isvector(dT)
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
timeLags = sampleLags/fs;

%Find time lag differences
numLags = length(sampleLags);
minLag = sampleLags(1);

if nargin <= 4 || isempty(edgeMin)
    needEdgeMin = true;
else
    needEdgeMin = false;
end

if nargin <= 5 || isempty(edgeMax)
    needEdgeMax = true;
else
    needEdgeMax = false;
end


if needEdgeMin || needEdgeMax
    m = 0;
    s = 0;
    numEl = 0;
    %Progress string set up
    fprintf("Finding edge(s): 00%%");
    lastprint = 0;
    %Calculation of time lag data
    for vIndex = 1:numLags
        lag = sampleLags(vIndex);
        for hIndex = 1:(length(data)-lag)
            temp = data(hIndex) - data(hIndex+lag);
            m = m + temp;
            s = s + temp^2;
            numEl = numEl + 1;
        end
        
        %Progress string checks
        if(floor(vIndex/numLags*100)) > lastprint
            lastprint = floor(vIndex/numLags*100);
            fprintf("\b\b\b%02d%%", lastprint);
        end
    end
    m = m/numEl;
    s = s/numEl;
    s = s - m^2;
    s = sqrt(abs(s));
    
    %Progress string finish
    fprintf("\b\b\b\bDone\n");
end

if needEdgeMin
    edgeMin = m - 4*s;
end
if needEdgeMax
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
    lag = sampleLags(vIndex);
    lagData = NaN([1 size(data)]);
    for hIndex = 1:(length(data)-lag)
        lagData(1,hIndex) = data(hIndex) - data(hIndex+lag);
    end
    %The probability normalization of histcounts is roughly
    %valI = countI/numel
    %which deflates later timelags as they are padded with nans
    %we fix this by passing only the non-nan values into histcounts
    [counts(vIndex,:) ~] = histcounts(lagData(1:(end - sampleLags(vIndex))), ...
        edges, 'Normalization', 'Probability');
    
    %Check if we need to update the progress string
    if(floor(vIndex/numLags*100)) > lastprint
        lastprint = floor(vIndex/numLags*100);
        fprintf("\b\b\b%02d%%", lastprint);
    end
end
%finish the progress string
fprintf("\b\b\b\bDone\n");


if nargin > 3 && ~isempty(preview) && islogical(preview)
    printPreview = preview;
else
    printPreview = false;
end

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