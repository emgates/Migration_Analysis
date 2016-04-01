function [vectors, theta] = PIVAngles(data, hilite)
% PIVangles has been adapted from PIVlab (PIVlab_commandline.m)

fprintf('PIV analysis started...');
amount = size(data,3);
tic

%% Standard PIV Settings
s = cell(10,2); % To make it more readable, let's create a "settings table"
%Parameter                       %Setting           %Options
s{1,1}= 'Int. area 1';           s{1,2}=256;         % window size of first pass
s{2,1}= 'Step size 1';           s{2,2}=128;         % step of first pass
s{3,1}= 'Subpix. finder';        s{3,2}=1;          % 1 = 3point Gauss, 2 = 2D Gauss
s{4,1}= 'Mask';                  s{4,2}=[];         % EVAN: We will replace this with hilite during call to PIV_fftmulti
s{5,1}= 'ROI';                   s{5,2}=[];         % Region of interest: [x,y,width,height] in pixels, may be left empty
s{6,1}= 'Nr. of passes';         s{6,2}=2;          % 1-4 nr. of passes
s{7,1}= 'Int. area 2';           s{7,2}=128;         % second pass window size
s{8,1}= 'Int. area 3';           s{8,2}=64;         % third pass window size
s{9,1}= 'Int. area 4';           s{9,2}=32;         % fourth pass window size
s{10,1}='Window deformation';    s{10,2}='*linear'; % '*spline' is more accurate, but slower

%% PIV analysis loop
x=cell(amount-1,1);
y=x;
u=x;
v=x;
typevector=x; %typevector will be 1 for regular vectors, 0 for masked areas

%% PIV analysis loop:
for i=1:amount-1
    image1=data(:,:,i); % read images
    image2=data(:,:,i+1);
    [x{i}, y{i}, u{i}, v{i}, typevector{i}] = piv_FFTmulti (image1,image2,s{1,2},s{2,2},s{3,2},s{4,2},s{5,2},s{6,2},s{7,2},s{8,2},s{9,2},s{10,2});
    clc
    disp([int2str((i+1)/amount*100) ' %']);
    
    % Graphical output (disable to improve speed)
    
%     imagesc(double(image1)+double(image2));colormap('gray');
%     hold on
%     quiver(x{i}, y{i}, u{i}, v{i},'g','AutoScaleFactor', 1.5);
%     hold off;
%     axis image;
%     title(['Frame ' int2str(i)],'interpreter','none')
%     set(gca,'xtick',[],'ytick',[])
%     drawnow;

end

%% PIV postprocessing loop
% Settings
umin = -25; % minimum allowed u velocity
umax = 25; % maximum allowed u velocity
vmin = -25; % minimum allowed v velocity
vmax = 25; % maximum allowed v velocity
stdthresh=6; % threshold for standard deviation check
epsilon=0.15; % epsilon for normalized median test
thresh=3; % threshold for normalized median test

u_filt=cell(amount-1,1);
v_filt=u_filt;
typevector_filt=u_filt;
for PIVresult=1:size(x,1)
    u_filtered=u{PIVresult,1};
    v_filtered=v{PIVresult,1};
    typevector_filtered=typevector{PIVresult,1};
    %vellimit check
    u_filtered(u_filtered<umin)=NaN;
    u_filtered(u_filtered>umax)=NaN;
    v_filtered(v_filtered<vmin)=NaN;
    v_filtered(v_filtered>vmax)=NaN;
    % stddev check
    meanu=nanmean(nanmean(u_filtered));
    meanv=nanmean(nanmean(v_filtered));
    std2u=nanstd(reshape(u_filtered,size(u_filtered,1)*size(u_filtered,2),1));
    std2v=nanstd(reshape(v_filtered,size(v_filtered,1)*size(v_filtered,2),1));
    minvalu=meanu-stdthresh*std2u;
    maxvalu=meanu+stdthresh*std2u;
    minvalv=meanv-stdthresh*std2v;
    maxvalv=meanv+stdthresh*std2v;
    u_filtered(u_filtered<minvalu)=NaN;
    u_filtered(u_filtered>maxvalu)=NaN;
    v_filtered(v_filtered<minvalv)=NaN;
    v_filtered(v_filtered>maxvalv)=NaN;
    % normalized median check
    %Westerweel & Scarano (2005): Universal Outlier detection for PIV data
    [J,I]=size(u_filtered);
    medianres=zeros(J,I);
    normfluct=zeros(J,I,2);
    b=1;
    for c=1:2
        if c==1; velcomp=u_filtered;else;velcomp=v_filtered;end %#ok<*NOSEM>
        for i=1+b:I-b
            for j=1+b:J-b
                neigh=velcomp(j-b:j+b,i-b:i+b);
                neighcol=neigh(:);
                neighcol2=[neighcol(1:(2*b+1)*b+b);neighcol((2*b+1)*b+b+2:end)];
                med=median(neighcol2);
                fluct=velcomp(j,i)-med;
                res=neighcol2-med;
                medianres=median(abs(res));
                normfluct(j,i,c)=abs(fluct/(medianres+epsilon));
            end
        end
    end
    info1=(sqrt(normfluct(:,:,1).^2+normfluct(:,:,2).^2)>thresh);
    u_filtered(info1==1)=NaN;
    v_filtered(info1==1)=NaN;

    typevector_filtered(isnan(u_filtered))=2;
    typevector_filtered(isnan(v_filtered))=2;
    typevector_filtered(typevector{PIVresult,1}==0)=0; %restores typevector for mask
    
    %Interpolate missing data
    u_filtered=inpaint_nans(u_filtered,4);
    v_filtered=inpaint_nans(v_filtered,4);
    
    u_filt{PIVresult,1}=u_filtered;
    v_filt{PIVresult,1}=v_filtered;
    typevector_filt{PIVresult,1}=typevector_filtered;
end

% Scrub data, using hilite as mask
vectors = scrub(u,v,x,y,hilite);

% First column are angles, second column indicates frame
theta = [atan2(vectors(:,2),vectors(:,1))*180/pi vectors(:,5)];
fprintf('DONE.\n')
toc
end

function [vectors] = scrub(u, v, x, y, hilite)
% Remove unwanted regions from dataset, keep frame and x- y- position data
count = 1;
for k = 1:size(u,1)
    for a = 1:size(u{k},1)
        for b = 1:size(u{k},2)
            if hilite(y{k}(a,b),x{k}(a,b),k) == 1
                vectors(count,:) = [u{k}(a,b) v{k}(a,b) x{k}(a,b) y{k}(a,b) k];
                count = count+1;
            end
        end
    end
end
end