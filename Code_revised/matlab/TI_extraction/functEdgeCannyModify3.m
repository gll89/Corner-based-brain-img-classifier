function [TI, TM, thresh] = functEdgeCannyModify3(varargin);
%      The Canny method finds edges by looking for local maxima of the
%      gradient of I. The gradient is calculated using the derivative of a
%      Gaussian filter. The method uses two thresholds, to detect strong
%      and weak edges, and includes the weak edges in the output only if
%      they are connected to strong edges. This method is therefore less
%      likely than the others to be "fooled" by noise, and more likely to
%      detect true weak edges.
%
%     Modification:
%     Change the while-operation in the function of selectThresholds (197-210) into a
%     matrix operation without loop (212-224)


    [a,method,thresh,sigma,thinning,H,kx,ky] = parse_inputs(varargin{:});
    if ~isa(a,'double') && ~isa(a,'single')
        a = im2single(a);
    end
    [m,n] = size(a);

    % Magic numbers
    PercentOfPixelsNotEdges = .8; % Used for selecting thresholds
    ThresholdRatio = .3;          % Low thresh is this fraction of the high.
    
    % Calculate gradients using a derivative of Gaussian filter
    [dx, dy] = smoothGradient(a, sigma);
    
    % Calculate Magnitude of Gradient
    magGrad = hypot(dx, dy);
    
    % Normalize for threshold selection
    magmax = max(magGrad(:));
    if magmax > 0
        magGrad = magGrad / magmax;
    end
    
    % Determine Hysteresis Thresholds
    [lowThreshArr, highThreshArr] = selectThresholds(thresh, magGrad, PercentOfPixelsNotEdges, ThresholdRatio, mfilename);
  
    %% extract multilayer texture image
    BW = zeros(size(a));
    TI = zeros(size(a));
    TM = zeros(size(a));
%     BW = logical(BW);
    leng = length(lowThreshArr);
%     disp([num2str(leng) 'layers']);
    for i = 1: leng
        lowThresh = lowThreshArr(i);
        highThresh = highThreshArr(i);
        % Perform Non-Maximum Suppression Thining and Hysteresis Thresholding of Edge
        % Strength
        e = thinAndThreshold(dx, dy, magGrad, lowThresh, highThresh);
        thresh = [lowThresh highThresh];
        if i == 1
            newLayerInd = find (e == 1);
        else
            e = double(e);
%             intersectInd = find (e&BW == 1);
            newLayerInd = find(xor(e, BW) == 1);
        end
        BW = e;
        TI(newLayerInd) = 1/i;
        TM(newLayerInd) = i;
%         figure, imshow(TI); 
%         title(num2str(i));
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : parse_inputs
%
function [I,Method,Thresh,Sigma,Thinning,H,kx,ky] = parse_inputs(varargin)
% OUTPUTS:
%   I      Image Data
%   Method Edge detection method
%   Thresh Threshold value
%   Sigma  standard deviation of Gaussian
%   H      Filter for Zero-crossing detection
%   kx,ky  From Directionality vector

narginchk(1,5)

I = varargin{1};

validateattributes(I,{'numeric','logical'},{'real','nonsparse','2d'},mfilename,'I',1);

% Defaults
Method    = 'sobel';
Direction = 'both';
Thinning  = true;

methods    = {'canny','canny_old','prewitt','sobel','marr-hildreth','log','roberts','zerocross'};
directions = {'both','horizontal','vertical'};
options    = {'thinning','nothinning'};

% Now parse the nargin-1 remaining input arguments

% First get the strings - we do this because the interpretation of the
% rest of the arguments will depend on the method.
nonstr = [];   % ordered indices of non-string arguments
for i = 2:nargin
    if ischar(varargin{i})
        str = lower(varargin{i});
        j = find(strcmp(str,methods));
        k = find(strcmp(str,directions));
        l = find(strcmp(str,options));
        if ~isempty(j)
            Method = methods{j(1)};
            if strcmp(Method,'marr-hildreth')
                error(message('images:removed:syntax','EDGE(I,''marr-hildreth'',...)','EDGE(I,''log'',...)')) 
            end
        elseif ~isempty(k)
            Direction = directions{k(1)};
        elseif ~isempty(l)
            if strcmp(options{l(1)},'thinning')
                Thinning = true;
            else
                Thinning = false;
            end
        else
            error(message('images:edge:invalidInputString', varargin{ i }))
        end
    else
        nonstr = [nonstr i]; %#ok<AGROW>
    end
end

% Now get the rest of the arguments
[Thresh,Sigma,H,kx,ky] = images.internal.parseNonStringInputsEdge(varargin,Method,Direction,nonstr);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : smoothGradient
%
function [GX, GY] = smoothGradient(I, sigma)

% Create an even-length 1-D separable Derivative of Gaussian filter

% Determine filter length
filterLength = 8*ceil(sigma);
n = (filterLength - 1)/2;
x = -n:n;

% Create 1-D Gaussian Kernel
c = 1/(sqrt(2*pi)*sigma);
gaussKernel = c * exp(-(x.^2)/(2*sigma^2));

% Normalize to ensure kernel sums to one
gaussKernel = gaussKernel/sum(gaussKernel);

% Create 1-D Derivative of Gaussian Kernel
derivGaussKernel = gradient(gaussKernel);

% Normalize to ensure kernel sums to zero
negVals = derivGaussKernel < 0;
posVals = derivGaussKernel > 0;
derivGaussKernel(posVals) = derivGaussKernel(posVals)/sum(derivGaussKernel(posVals));
derivGaussKernel(negVals) = derivGaussKernel(negVals)/abs(sum(derivGaussKernel(negVals)));

% Compute smoothed numerical gradient of image I along x (horizontal)
% direction. GX corresponds to dG/dx, where G is the Gaussian Smoothed
% version of image I.
GX = imfilter(I, gaussKernel', 'conv', 'replicate');
GX = imfilter(GX, derivGaussKernel, 'conv', 'replicate');

% Compute smoothed numerical gradient of image I along y (vertical)
% direction. GY corresponds to dG/dy, where G is the Gaussian Smoothed
% version of image I.
GY = imfilter(I, gaussKernel, 'conv', 'replicate');
GY  = imfilter(GY, derivGaussKernel', 'conv', 'replicate');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : selectThresholds
%
function [lowThreshArr, highThreshArr] = selectThresholds(thresh, magGrad, PercentOfPixelsNotEdges, ThresholdRatio, ~)

[m,n] = size(magGrad);

% Select the thresholds
if isempty(thresh)
    counts=imhist(magGrad, 64);
%  figure(1), plot(1:length(counts), counts);
    countCum = cumsum(counts);  % column vector
%  figure(2), plot(1:length(countCum), countCum);
    highThresh = find(countCum> PercentOfPixelsNotEdges*m*n,...
        1,'first') / 64;
    lowThresh = ThresholdRatio*highThresh;
    %% added code
    highThreshArr = [highThresh];
    lowThreshArr = [lowThresh];
    timeN = 0.1*1;
    
%     tmp = [PercentOfPixelsNotEdges];
%     tmpPercent = PercentOfPixelsNotEdges - timeN;
%     while tmpPercent >0
%         tmp = [tmp, tmpPercent];
%         highTmp = find(countCum > tmpPercent*m*n,...
%         1,'first') / 64;
%         if highTmp ~= highThresh
%             lowTmp = ThresholdRatio*highTmp;
%             highThreshArr = [highThreshArr; highTmp];
%             lowThreshArr = [lowThreshArr; lowTmp];
%             highThresh = highTmp;
%         end
%         tmpPercent = tmpPercent-timeN;
%     end
    
    leng = floor(PercentOfPixelsNotEdges/timeN);
    tmpPercent = zeros(1, leng);
    tmpPercent = PercentOfPixelsNotEdges : -timeN : 0.01;
    tmpPercent1 = (m*n) .* tmpPercent ;
    highElement = bsxfun(@ge, countCum', tmpPercent1');   
    %highElement is a size(countCum,1)*size(tmpPercent, 2) matrix, 
    %It is a binary matrix, returing 1 if countCum(i)>=tmpPercent(j),
    %countCum is a column vector, and tmpPercent is a row vector
    [elemS, indS] = sort(highElement, 2, 'descend');
    highTmp1 = indS(:, 1);
    highTmp2 = unique(highTmp1, 'stable');   %remove the duplicate elements
    highThreshArr = highTmp2./64;
    lowThreshArr = ThresholdRatio .* highThreshArr;
    
    
elseif length(thresh)==1
    highThresh = thresh;
    if thresh>=1
        error(message('images:edge:thresholdMustBeLessThanOne'))
    end
    lowThresh = ThresholdRatio*thresh;
elseif length(thresh)==2
    lowThresh = thresh(1);
    highThresh = thresh(2);
    if (lowThresh >= highThresh) || (highThresh >= 1)
        error(message('images:edge:thresholdOutOfRange'))
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : thinAndThreshold
%
function H = thinAndThreshold(dx, dy, magGrad, lowThresh, highThresh)
% Perform Non-Maximum Suppression Thining and Hysteresis Thresholding of
% Edge Strength
    
% We will accrue indices which specify ON pixels in strong edgemap
% The array e will become the weak edge map.

E = cannyFindLocalMaxima(dx,dy,magGrad,lowThresh);

if ~isempty(E)
    [rstrong,cstrong] = find(magGrad>highThresh & E);
    
    if ~isempty(rstrong) % result is all zeros if idxStrong is empty
        H = bwselect(E, cstrong, rstrong, 8);
    else
        H = false(size(E));
    end
else
    H = false(size(E));
end
