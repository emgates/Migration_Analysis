function filtered_image = bwpropopen(bw_image,prop_name,min_val,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.addRequired('bw_image',@(x)isnumeric(x) || islogical(x));
i_p.addRequired('prop_name',@ischar);
i_p.addRequired('min_val',@isnumeric);

i_p.addParamValue('connectivity',8,@(x)x == 4 || x == 8);

i_p.parse(bw_image,prop_name,min_val,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bw_label = bwlabel(bw_image,i_p.Results.connectivity);
props = regionprops(bw_label,prop_name);

filtered_image = ismember(bw_label,find([props.(prop_name)] >= min_val));

end