function [CRS, COMPILED] = PIV_CORRELATION(data)
% PIV_CORRELATION is used to calculate x- and y- spatial velocity
% correlations for a timelapse movie

% Input: PIV matrix with the form [x-position y-position x-velocity
% y-velocity Frame#]

% Output: Matrix with frame data of the form [Distance Mean-Ix Mean-Iy
% StdErr-Ix StdErr-Iy Frame#] and a matrix averaged across the frames of the
% same form (without Frame# column of course)

% Bin size: please update accordingly
bSize = 25;

% Initilize CRS
CRS = [];

%% Calculate correlation curves for each frame
for k = min(data(:,5)):max(data(:,5))
    % Get current frame's data
    X = data(data(:,5)==k,1);
    Y = data(data(:,5)==k,2);
    Vx = data(data(:,5)==k,3);
    Vy = data(data(:,5)==k,4);
    
    % Generate symmetrical matrices to generate all possible correlations
    X = toeplitz(X);
    Y = toeplitz(Y);
    Vx = toeplitz(Vx);
    Vy = toeplitz(Vy);
    
    % Normalize velocities
    Vx = Vx./abs(Vx);
    Vx(isnan(Vx)) = 0;
    Vy = Vy./abs(Vy);
    Vy(isnan(Vy)) = 0;
    
    % Calculate distance matrix
    D = sqrt(bsxfun(@minus,X(:,2:end),X(:,1)).^2+...
        bsxfun(@minus,Y(:,2:end),Y(:,1)).^2);
    clear X Y
        
    % Calculate correlations
    Ix = bsxfun(@times,Vx(:,2:end),Vx(:,1));
    Iy = bsxfun(@times,Vy(:,2:end),Vy(:,1));
    clear Vx Vy
    
    % Rearrange matrixes into arrays
    D = D(:);
    Ix = Ix(:);
    Iy = Iy(:);
    
    % Reduce the amount of data via binning
    % Define bin edges
    bEdges = 0:bSize:max(D)+bSize;
    [~,Dbins] = histc(D,bEdges);
    clear D
    % Move edge values to center of bins
    bEdges = bEdges+bSize/2;
    
    % Initialize matrix for compiling data
    A = zeros(length(bEdges),6);
    
    % Calculate mean and stderr for each bin
    for c = 1:length(bEdges)
        A(c,2:5) = [mean(Ix(Dbins==c)) mean(Iy(Dbins==c))...
            std(Ix(Dbins==c))/sqrt(length(Ix(Dbins==c)))...
            std(Iy(Dbins==c))/sqrt(length(Iy(Dbins==c)))];
    end
    clear Ix Iy
    
    % Add distance info
    A(:,1) = bEdges;
    
    % Remove any values with SE = 0
    A = A(A(:,4)~=0 & A(:,5)~=0 &...
        ~isnan(A(:,4)) & ~isnan(A(:,5)),:);
    
    % Add frame identification
    A(:,6) = k;
    
    % Compile master
    CRS = [CRS; A];
    clear A
    
    clc
    disp([int2str((k)/max(data(:,5))*100) ' %']);
end

% Compile correlations across all frames
% Define bin edges
bEdges = 0:bSize:max(CRS(:,1))+bSize;
[~,Dbins] = histc(CRS(:,1),bEdges);
clear D
% Move edge values to center of bins
bEdges = bEdges+bSize/2;

% Initialize matrix for compiling data
COMPILED = zeros(length(bEdges),5);

% Calculate mean and stderr for each bin
for c = 1:length(bEdges)
    COMPILED(c,2:5) = [mean(CRS(Dbins==c,2)) mean(CRS(Dbins==c,3))...
        std(CRS(Dbins==c,2))/sqrt(length(CRS(Dbins==c,2)))...
        std(CRS(Dbins==c,3))/sqrt(length(CRS(Dbins==c,3)))];
end

% Add distance info
COMPILED(:,1) = bEdges;

% Remove any values with SE = 0
COMPILED = COMPILED(COMPILED(:,4)~=0 & COMPILED(:,5)~=0 &...
    ~isnan(COMPILED(:,4)) & ~isnan(COMPILED(:,5)),:);

end

