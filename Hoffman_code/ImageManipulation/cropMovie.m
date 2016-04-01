function [data] = cropMovie(data)
%CROPMOVIE Displays movie segment with imshow3D and has user draw rectangle
%around region of interest, cropping the movie to that space.

imshow3D(data);
tog = input('Would you like to crop this movie? (y/n): ', 's');
if strcmp(tog,'y')
    rect = getrect;
    rect= round(rect);
    close
    data = data(rect(2):(rect(2)+rect(4)),rect(1):(rect(1)+rect(3)),:);
end
close

end

