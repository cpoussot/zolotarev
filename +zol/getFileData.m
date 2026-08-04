function out = getFileData(fileData)

allFiles    = dir(fileData);
allNames    = {allFiles.name};
%fileName    = strrep(fileName,'_',' ');
k = 0;
for i = 1:numel(allNames)
    if length(allNames{i})>6
        if strcmp('data_r',allNames{i}(1:6))
            k       = k+1;
            out{k}  = allNames{i};
        end
    end
end
