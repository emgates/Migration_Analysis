function label_filtered = labelpropopen(label_image,prop_name,min_val)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.addRequired('label_image',@(x)isnumeric(x) || islogical(x));
i_p.addRequired('prop_name',@ischar);
i_p.addRequired('min_val',@isnumeric);

i_p.parse(label_image,prop_name,min_val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

props = regionprops(label_image,prop_name);

filtered_image = ismember(label_image,find([props.(prop_name)] >= min_val));

label_filtered = label_image.*filtered_image;

end