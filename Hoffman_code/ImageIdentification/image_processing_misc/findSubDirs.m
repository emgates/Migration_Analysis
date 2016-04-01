function [folders] = findSubDirs(parentDir, target)
% Seek out any target folders downstream of the parentDir and return the
% full paths of all targets as a structured cell array

% Target should obey REGEXP rules

% Determine home directory
if ispc
    home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
else
    home = getenv('HOME');
end

% Make sure REGEXP only finds where target is the terminal folder i9n the
% path
target = [target '[;:]'];

% Generate long string of all possible sub-directory paths
subPaths = genpath(parentDir);

% Find hits for the target
target_start = regexp(subPaths, target);
% Find all home directory listings
home_start = strfind(subPaths, home);

folders = cell(length(target_start),1);

% Break the filenames into cells
for k = 1:length(target_start)
    [start, index] = max(home_start(home_start<target_start(k)));
    if index ~= length(home_start)
        folders{k}.name = subPaths(start:(home_start(index+1)-2));
    else
        folders{k}.name = subPaths(start:(end-1));
    end
end

end

