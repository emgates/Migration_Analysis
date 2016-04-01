function largestObj = filter_to_largest_object(labelImage,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'filter_to_largest_object';

i_p.addRequired('labelImage',@(x)isnumeric(x) || islogical(x));
i_p.addParameter('connectivity',8,@(x)isnumeric);

i_p.parse(labelImage);

i_p.parse(labelImage,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (max(labelImage(:)) == 0)
    error('Are you sure you gave me a labeled image? I''m only seeing zeros');
elseif (max(labelImage(:)) == 1)
    labelImage = bwlabel(labelImage,i_p.Results.connectivity);
end

props = regionprops(labelImage,'Area');

largestObj = ismember(labelImage, find([props.Area] == max([props.Area])));

end
