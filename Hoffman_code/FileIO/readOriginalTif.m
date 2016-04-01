function [tiff_stack, datetime] = readOriginalTif(originalTifFolder, fileType)
% This function identifies .TIF files and forms a stack in the natural
% ascending order

%originalTifFolder = folder containing desired image files
%fileType = .tif, .TIF, .tiff, etc.

fprintf('File import initiated...');
tic;

files = dir(fullfile(originalTifFolder, ['*' fileType]));

if isempty(files)
    fprintf('Files could not be located.\n');
    tiff_stack = 0;
    datetime = 0;
    return;
end

% Sort files in correct numerical order
files = sort_nat({files.name});

tiff_stack = imread(fullfile(originalTifFolder,files{1}));
% Read in timestamp information
temp_date = imfinfo(fullfile(originalTifFolder,files{1}));
datetime = datenum(temp_date.DateTime,'yyyymmdd HH:MM:SS.FFF');

for ii = 2 : length(files)
    temp_tiff = imread(fullfile(originalTifFolder,files{ii}));
    temp_date = imfinfo(fullfile(originalTifFolder,files{ii}));
    datetime(ii) = datenum(temp_date.DateTime,'yyyymmdd HH:MM:SS.FFF');
    tiff_stack = cat(3 , tiff_stack, temp_tiff);
end

fprintf('DONE.\n');
toc;

end