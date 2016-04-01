function [] = padFiles(myFolder, fileType, padLength)
% Rename files to be padded to padLength

myFiles = dir(fullfile(myFolder, ['*' fileType]));

for id = 1:length(myFiles)
    % Get the file name (minus the extension)
    [~, f] = fileparts(myFiles(id).name);
    
    for k = length(f):-1:1
        if isnan(str2double(f(k)))
            break
        end
    end
    
    movefile(fullfile(myFolder, myFiles(id).name),...
        fullfile(myFolder,...
        sprintf([f(1:k) '%0' num2str(padLength) 'd' fileType], str2double(f(k+1:end)))));
    
end

end