function [] = leadingEdgeCSV()
% Let's compile data into a CSV!

data_lead = cell(0);
data_PIV = cell(0);
datetime = [];
contourLength = [];
fitLength = [];
lineTheta = [];
r2 = [];
headers_lead = {'Date','Sample','Type',...
    'Contour Length','Fit Length','Fit Orientation','Fit Error (r-sq)',...
    'dA.dt Line','dA.dt Edge',...
    'Time Difference (hr)'};
headers_PIV = {'Date','Sample','Type','X-Velocity','Y-Velocity',...
    'X-Position','Y-Position','Frame'};

% User input
clc
tog(1) = {input('Date: ','s')};
tog(2) = {input('Sample: ','s')};
tog(3) = {input('Type: ','s')};

while ~isequal(tog(1),{''})
    
    % Get directory & file .mat
    myDirectory = uigetdir();
    myFile = dir(fullfile(myDirectory,'*.mat'));
    myFile = myFile(1).name;
    
    % Load files
    load(fullfile(myDirectory,myFile),'contourLength','fitLength',...
        'dA_line','dA_orig','datetime','lineTheta','r2',...
        'PIV_vectors');
    
    t = datetime(2:end)-datetime(1:end-1);
    % Convert to um/frame
    PIV_vectors(:,1:2) = PIV_vectors(:,1:2)*0.53;
    % Convert vectors to micron/hr
    for k = 1:size(PIV_vectors,1)
        PIV_vectors(k,1:2) = PIV_vectors(k,1:2)/t(PIV_vectors(k,5));
    end
    
    n_lead = length(dA_line);
    n_PIV = size(PIV_vectors,1);
    % Organize into columns
    data_lead = [data_lead ...
        [repmat(tog(1),1,n_lead); repmat(tog(2),1,n_lead);repmat(tog(3),1,n_lead);...
        num2cell(contourLength(2:end)); num2cell(fitLength(2:end));...
        num2cell(lineTheta(2:end)); num2cell(r2(2:end));...
        num2cell(dA_line); num2cell(dA_orig); num2cell(t)]];
    data_PIV = [data_PIV ...
        [repmat(tog(1),1,n_PIV); repmat(tog(2),1,n_PIV);repmat(tog(3),1,n_PIV);...
        num2cell(PIV_vectors')]];
    
    % User input
    clc
    tog(1) = {input('Date: ','s')};
    tog(2) = {input('Sample: ','s')};
    tog(3) = {input('Type: ','s')};
end

data_lead = data_lead';
data_PIV = data_PIV';

% Write leading edge
fid = fopen('leading_edge.csv','wt');
fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',headers_lead{:});
if fid>0
    for k=1:size(data_lead,1)
        fprintf(fid,'%s,%s,%s,%f,%f,%f,%f,%f,%f,%f\n',data_lead{k,:});
    end
fclose(fid);
end

% Write PIV
fid = fopen('PIV_data.csv','wt');
fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s\n',headers_PIV{:});
if fid>0
    for k=1:size(data_PIV,1)
        fprintf(fid,'%s,%s,%s,%f,%f,%d,%d,%d\n',data_PIV{k,:});
    end
fclose(fid);
end

end

