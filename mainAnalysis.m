function [] = mainAnalysis()
%% leadingEdgeAnalysis Overview
%{

Requirements:
-Folder containing single timelapse experiment
-Original DIC images, in .TIF format (original required for timestamps)

Function workflow:
-User input and data import
++Get directory & files
++Rotate, segment, and crop movies
-Edge finder

%}

close all

%% Establish necessary paths
addpath(genpath('../'));

%% User input and data import
% Import files from specified directory
fprintf('Please locate folder\n');
tifFolder = uigetdir;
fileType = '.TIF';
% Data = stack of 2D images, datetime = vector of time stamps (days, see
% datenum)
[data, datetime] = readOriginalTif(tifFolder, fileType);
if data == 0
    return
end
clc

% Rotate image to satisfy edge identification code
% Cells should migrate left to right across the screen for this code
imshow3D(data);
data = rotateImage(data);
clc
close
% Segment movie to meet appropriate criteria
imshow3D(data);
[data, datetime] = segmentMovie(data, datetime);
% Convert datetime to elapsed time (rather than absolute), in hours
datetime = (datetime-datetime(1))*24;
close
clc
% Crop movie segments to meet appropriate criteria and remove portions that
% cannot be analyzed
data = cropMovie(data);
clc
% Get base filename
tog = input('Please enter filename for output: ','s');

%% Edge Finder (complements of MEB)
% Use covariance filter to create mask of cells
[y,x,z] = size(data);
mask = zeros(y,x,z); % Binary mask
points = cell(1,z); % Edge mask only

% Edge ID
mask = edgeID(data,y,x,z);
% Isolate points along leading edge
points = edgeBoundary(mask,x,y,z);

%% PIV Analysis (adopted from PIVlab)
[PIV_vectors, PIV_theta] = PIVAngles(data,mask);
% Write PIV to folder
PIVfolder = fullfile(tifFolder,'PIV_images');
mkdir_no_err(PIVfolder);
for k = 1:z
    figure('Visible','off')
    imshow(imadjust(data(:,:,k)))
    hold on
    quiver(PIV_vectors(:,3),PIV_vectors(:,4),...
        PIV_vectors(:,1),PIV_vectors(:,2),'LineWidth',2);
    saveas(gcf,fullfile(PIVfolder, ['PIV_' int2str(k) '.png'])]);
end

%% Curve Fitting to Leading Edge
coeffs = edgeFit(points,z);

%% Calculate Parameters
% Call function to perform analysis and calculate parameters
[dA_orig, dA_line, fitLength, contourLength, lineTheta, r2] =...
    params(mask, points, coeffs, datetime,x,y,z);

%% Align PIV_theta with lineTheta
% Normalize each frame with mean PIV_theta
for k = 1:length(lineTheta)-1
    PIV_theta(PIV_theta(:,2) == k,1) = ...
        PIV_theta(PIV_theta(:,2) == k,1)-...
        mean(PIV_theta(PIV_theta(:,2) == k,1));
end

%% Write params to file
% Plot parameters for viewing
% dA
figure('Visible','off');
subplot(4,1,1)
title('Plots of dA/(l dt) params');
plot(repmat(datetime(2:end)',1,2),...
    [(dA_orig./fitLength(1:end-1))',...
    (dA_orig./contourLength(1:end-1))']);
ylabel('(um/hr)');  axis([0 inf -25 50]);
legend('Original (Fit Normal)',...
    'Original (Contour Nrm)','Location','NorthEast');
set(gca, 'XMinorTick','on','YMinorTick','on');
grid on
grid minor

% Fit error
subplot(4,1,2)
ax = plotyy(datetime,fitLength./contourLength,datetime,r2);
ylabel(ax(1),'Cntr Ratio');
ylabel(ax(2),'Avg Err: R^2');
axis(ax(1),[0 inf 0 1]);
% Figure out ideal range of error
axis(ax(2),[0 inf mean(r2)-2*std(r2) mean(r2)+2*std(r2)]);
set(gca, 'XMinorTick','on','YMinorTick','on');

% Line orientation
subplot(4,1,3)
plot(datetime,lineTheta);
ylabel('Line Ornt (deg)');
xlabel('Time elapsed (hr)');
axis([0 inf -60 60]);
set(gca, 'XMinorTick','on','YMinorTick','on');
grid on
grid minor

% PIV vector angles
subplot(4,1,4)
histogram(PIV_theta(:,1),36);
xlabel('Degrees')
axis([-180 180 0 inf])
saveas(gcf, fullfile(tifFolder,[tog '-analysis.jpg']));
saveas(gcf, fullfile(tifFolder,[tog '-analysis.svg']));

% Write variables to a file for future use
save(fullfile(tifFolder, [tog '-analysis.mat']),...
    'datetime','mask','points','coeffs',...
    'dA_orig', 'dA_line', 'fitLength', 'contourLength',...
    'lineTheta', 'r2','PIV_vectors','PIV_theta');

%% Create movie visualization
% Create visualization
fprintf('Creating movie visualization...');
tic
visual = visualize(data,mask,coeffs, PIV_vectors);
% Write to file
v = VideoWriter(fullfile(tifFolder,[tog '.avi']));
v.FrameRate = 15;
open(v)
writeVideo(v,im2uint8(visual));
close(v)
fprintf('DONE.\n');
toc

end

%% edgeBounbary
function [points] = edgeBoundary(bw_edge,x,y,z)
% This function will isolate points along the leading edge boundary
fprintf('Edge boundary isolation started...');
tic
points = cell(1,z);
% Loop through each image
for k = 1:z
    temp = bwboundaries(bw_edge(:,:,k));
    % Remove anything along the image border
    temp = temp{1}(temp{1}(:,1) ~= 1,:);
    temp = temp(temp(:,2) ~= 1,:);
    temp = temp(temp(:,1) ~= y,:);
    temp = temp(temp(:,2) ~= x,:);
    points{k} = temp;
end
fprintf('DONE.\n');
toc
end

%% edgeFit
function [coeffs] = edgeFit(points,z)
% Fits a linear curve to the leading edge of each frame
fprintf('Leading edge curve fitting started...');
tic
coeffs = zeros([z,2]);
% Loop through each image
for b = 1:z
    coeffs(b,:) = polyfit(points{b}(:,1),points{b}(:,2),1);
end
fprintf('DONE.\n');
toc
end

%% visualize
function [visual] = visualize(data, hilite, coeffs, PIV_vectors)
visual = zeros([size(data,1) size(data,2) 3 size(data,3)]);
for k = 1:size(data,3)
    y = 1:size(data,1);
    % Make highlight
    visual(:,:,:,k) = ...
        create_highlighted_image(normalize_image(data(:,:,k)),...
        thicken_perimeter(bwperim(hilite(:,:,k)),3),...
        'mix_percent',0.5);
%     f = figure('Visible', 'off');
%     axes('Visible','off');
%     imshow(visual(:,:,:,k));
%     hold on;
%     temp = PIV_vectors(PIV_vectors(:,5)==k,1:4);
%     quiver(temp(:,3),temp(:,4),temp(:,1),temp(:,2),'LineWidth',2);
%     temp = getframe(f);
%     [visual(:,:,:,k),~] = frame2im(temp);
    % Make line
    x = polyval(coeffs(k,:),y);
    x = round(x);
    % Get rid of points outside the image matrix
    y = y(x>1 & x<size(data,2));
    x = x(x>1 & x<size(data,2));
    % Overlay line on RGB
    % Green
    visual(sub2ind(size(visual),[y y y],[x-1 x x+1],...
        ones(1, 3*length(x)),repmat(k,[1 3*length(x)])))...
        = intmax('uint8');
end
end

%% params
function [dA_orig, dA_line, fitLength, contourLength, lineTheta, r2] =...
    params(mask, points, coeffs, datetime,x,y,z)
% dA_orig will return a vector of the change in area between frames divided
% by both the time elapsed and the average line fit length
% dA_orig will be in units of um/s

fprintf('Calculating analysis paramters...');
tic

% For VivaView 20x, conversion is 1px = 0.32um
conv = 0.32;
dA_orig = zeros(1,z-1);
dA_line = zeros(1,z-1);
lineTheta = zeros(1,z);
r2 = lineTheta;

% Initialize starting frame for variables
y1 = 1:y;
% Make line
x1 = polyval(coeffs(1,:),y1);
x1 = round(x1);
% Get rid of points outside the image matrix
y1 = y1(x1>=1 & x1<=x);
x1 = x1(x1>=1 & x1<=x);
% Make mask from line
im1 = zeros(y,x);
for k = 1:size(y1,2)
    im1(y1(k),1:x1(k)) = 1;
end
% Adjust for cutoff
if y1(end) < y && x1(end) == x
    im1(y1(end)+1:y,:) = 1;
end
if y1(1) > 1 && x1(end) == x
    im1(1:y1(1),:) = 1;
end

% Start fitLength & contourLength
fitLength = zeros(1,z);
contourLength = fitLength;
fitLength(1) = sqrt((x1(end)-x1(1))^2+(y1(end)-y1(1))^2);
contourLength(1) = sum(sqrt(...
    (points{1}(2:end,1)-points{1}(1:end-1,1)).^2+...
    (points{1}(2:end,2)-points{1}(1:end-1,2)).^2));

% Start lineTheta
lineTheta(1) = atan(coeffs(1,1))*180/pi;

% Start r^2
r2(1) = mean((x1-polyval(coeffs(1,:),y1)).^2);

% Begin looping through all frames
for k = 1:z-1
    % Make new line
    y1 = 1:y;
    x1 = polyval(coeffs(k+1,:),y1);
    x1 = round(x1);
    % Get rid of points outside the image matrix
    y1 = y1(x1>=1 & x1<=x);
    x1 = x1(x1>=1 & x1<=x);
    
    % Update fitLength and contourLength
    fitLength(k+1) = sqrt((x1(end)-x1(1))^2+(y1(end)-y1(1))^2);
    contourLength(k+1) = sum(sqrt(...
        (points{k+1}(2:end,1)-points{k+1}(1:end-1,1)).^2+...
        (points{k+1}(2:end,2)-points{k+1}(1:end-1,2)).^2));
    
	% Update lineTheta
    lineTheta(k+1) = atan(coeffs(k+1,1))*180/pi;   
    
    % Update r^2
    r2(k+1) = mean((x1-polyval(coeffs(k+1,:),y1)).^2);
    
    % Calculate dA param from edge finder
    dA_orig(k) = sum(sum(mask(:,:,k+1)-mask(:,:,k)))*conv/...
        ((datetime(k+1)-datetime(k)));
    
    % Calculate dA param from line fits
    % First, make new image
    im2 = zeros(y,x);
    for ii = 1:size(y1,2)
        im2(y1(ii),1:x1(ii)) = 1;
    end
    
    % Adjust for cutoff
    if y1(end) < y && x1(end) == x
        im2(y1(end)+1:y,:) = 1;
    end
    if y1(1) > 1 && x1(end) == x
        im2(1:y1(1),:) = 1;
    end
    
    % Calculate dA from images (can use same l_avg)
    dA_line(k) = sum(sum(im2-im1))*conv/...
        ((datetime(k+1)-datetime(k)));
    
    im1 = im2;
end

fprintf('DONE.\n');
toc

end