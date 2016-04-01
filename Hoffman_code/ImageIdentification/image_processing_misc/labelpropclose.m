function [label_filtered,varargout] = labelpropclose(label_image,prop_name,max_val,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.addRequired('label_image',@(x)isnumeric(x) || islogical(x));
i_p.addRequired('prop_name',@ischar);
i_p.addRequired('max_val',@isnumeric);

i_p.addOptional('return_filtered_ids',0);

i_p.parse(label_image,prop_name,max_val,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
label_image = double(label_image);

props = regionprops(label_image,prop_name);

if (i_p.Results.return_filtered_ids)
    id_nums = nonzeros(unique(label_image(:)));
    passed_ids = find([props.(prop_name)] <= max_val);
    not_passed_ids = setdiff(id_nums,passed_ids);
    varargout{1} = not_passed_ids;
end

filtered_image = ismember(label_image,find([props.(prop_name)] <= max_val));

label_filtered = label_image.*filtered_image;

end