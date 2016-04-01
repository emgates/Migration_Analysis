function image = normalize_image(image,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;

i_p.addRequired('image',@isnumeric);

i_p.addParameter('quantiles',@isnumeric);
i_p.addParameter('limits',@isnumeric);

i_p.addParameter('debug',0,@(x)(isnumeric(x) && (x == 0 || x == 1) || islogical(x)));

i_p.parse(image,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image = double(image);

if (not(any(strcmp(i_p.UsingDefaults,'quantiles'))))
    limits = quantile(image(:),i_p.Results.quantiles);
elseif (not(any(strcmp(i_p.UsingDefaults,'limits'))))
    limits = i_p.Results.limits;
else
    limits = [min(image(:)),max(image(:))];
end

image = (image - limits(1))/range(limits);

image(image > 1) = 1;
image(image < 0) = 0;
