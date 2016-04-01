function [IM_out] = rotateImage(IM_in)
%ROTATEIMAGE User will visualize stack with imshow3D and determine amount
%to rotate stack counterclockwise

fprintf('Cells should migrate left to right across the screen.\n');
tog = input('Does this image need to be rotated? (y/n):','s');

if strcmp(tog,'n')
    IM_out = IM_in;
    return;
elseif strcmp(tog,'y')
    clear tog;
    tog = input('Turns to rotate counterclockwise (1-3):');
    
    if mod(tog,1) == 0 && tog > 0 && tog < 4
        IM_out = rot90(IM_in,tog);
    else
        fprintf('Error: incorrect input\n');
        IM_out = rotateImage(IM_in);
        return;
    end
    
else
    fprintf('User input not recognized (try y/n)\n');
    IM_out = rotateImage(IM_in);
    return;
end

end