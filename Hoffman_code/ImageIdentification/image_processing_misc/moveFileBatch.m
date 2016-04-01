function [] = moveFileBatch(myFolder, fileID, destination)

allFiles = dir(fullfile(myFolder,fileID));
mkdir_no_err(fullfile(myFolder,destination));
for m = 1:length(allFiles)
    movefile(fullfile(myFolder,allFiles(m).name),...
        fullfile(myFolder,destination));
end

end