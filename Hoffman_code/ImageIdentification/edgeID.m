function [hilite] = edgeID(data,y,x,z)
%EDGEID Takes in a stack of 2D images and performs a standard deviation
%filter based on COV_thresh to identify leading edge of cells imaged with
%DIC

hilite = zeros(y,x,z);
ii = 1;

while ii == 1
    %Determine best COV_thresh value
    COV_thresh = input('Input covariance threshold (0.014 is a good start): ');
    
    base = data(:,:,ii);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate filtered images to conduct segmentation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    std_filt = stdfilt(base,getnhood(strel('disk',15))); %15 for vivaview
    
    I_filt = fspecial('disk',15); %15 for vivaview
    mean_filt = imfilter(base,I_filt,'replicate');
    
    cov_image = double(std_filt) ./ double(mean_filt);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Segment and clean the image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cell_mask = cov_image > COV_thresh;
    
    %Pick out the largest object, we will assume it is the cell region.
    %Then fill in any holes missed by the COV filter.
    cell_mask = filter_to_largest_object(cell_mask);
    
    %The fairly large search region used for in stdfilt has a tendency to
    %expand the boundary around the cells. This will remove that 15 pixel
    %boundary and then re-filter to pick out the largest object
    cell_mask(1:end,1) = 1;
    cell_mask(1,1:end) = 1;
    cell_mask(end,1:end)=1;
    cell_mask = imfill(cell_mask,'holes');
    cell_mask = imclose(cell_mask,strel('disk',15));
    cell_mask = imfill(cell_mask,'holes');
    cell_mask = imerode(cell_mask,strel('disk',15));
    cell_mask = filter_to_largest_object(cell_mask);

    hilite(:,:,ii) = cell_mask;
    
    imshow(create_highlighted_image(normalize_image(data(:,:,ii)),...
        thicken_perimeter(bwperim(hilite(:,:,ii)),3),...
        'mix_percent',0.5),'InitialMagnification',50);
    tog = input('Accept covariance threshold? (y/n): ','s');
    close;
    if strcmp(tog,'y')
        % Reset COV_thresh
        ii = 2;
    end    
end

fprintf('Edge identification started...');
tic

parfor ii = 2:z
    
    base = data(:,:,ii);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate filtered images to conduct segmentation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    std_filt = stdfilt(base,getnhood(strel('disk',15))); %15 for vivaview
    
    I_filt = fspecial('disk',15); %15 for vivaview
    mean_filt = imfilter(base,I_filt,'replicate');
    
    cov_image = double(std_filt) ./ double(mean_filt);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Segment and clean the image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cell_mask = cov_image > COV_thresh;
    
    %Pick out the largest object, we will assume it is the cell region.
    %Then fill in any holes missed by the COV filter.
    cell_mask = filter_to_largest_object(cell_mask);
    
    %The fairly large search region used for in stdfilt has a tendency to
    %expand the boundary around the cells. This will remove that 15 pixel
    %boundary and then re-filter to pick out the largest object
    cell_mask(1:end,1) = 1;
    cell_mask(1,1:end) = 1;
    cell_mask(end,1:end)=1;
    cell_mask = imfill(cell_mask,'holes');
    cell_mask = imclose(cell_mask,strel('disk',15));
    cell_mask = imfill(cell_mask,'holes');
    cell_mask = imerode(cell_mask,strel('disk',15));
    cell_mask = filter_to_largest_object(cell_mask);

    hilite(:,:,ii) = cell_mask;
end

fprintf('DONE.\n');
toc

end