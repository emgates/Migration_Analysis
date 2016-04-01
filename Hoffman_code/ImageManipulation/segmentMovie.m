function [data_new, datatime_new] = segmentMovie(data, datetime)
%SEGMENTMOVIE prompts the user to segment a long image stack

err_tog = 'n';

while strcmp(err_tog,'n')
    err_tog = input('Perform segmentation? (y/n): ','s');

    if strcmp(err_tog, 'y')
        intervals(1) = input('First frame in segment?:');
        intervals(2) = input('Last frame in segment?:');
        
        if intervals(1) < intervals(2) && intervals(1) >= 1 &&...
               intervals(2) <= size(data,3) && intervals(1) < intervals(2)
            data_new = data(:,:,intervals(1):intervals(2));
            datatime_new = datetime(intervals(1):intervals(2));
            return;
        else
            fprintf('Problem with your choice. Try again...');
            err_tog = 'n';
        end
    elseif strcmp(err_tog, 'n')
        intervals = [1 length(datetime)];
        data_new = data(:,:,intervals(1):intervals(2));
        datatime_new = datetime(intervals(1):intervals(2));
        return;
    else
        fprintf('Problem with your choice. Try again...');
        err_tog = 'n';
    end
    
end
end