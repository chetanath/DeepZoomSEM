clearvars;
close all;

fov=350;

dowriteEDSmaps=false; % write EDS maps

increaseborderfactor=1.5; % increase area by this fraction

overlap_fraction=0.1; % 10% = 0.1

fovx=fov;
fovy=fov;

% fovx=256; % For High-Resolution Map, in microns
% 
% fovy=fovx*870/1114; %=fovy; normally || =fovx*870/1114; for EDAX scans (Eltanin residue mounts)


if dowriteEDSmaps
    
    overlap_fraction_x=overlap_fraction+((512-506)/512);
    overlap_fraction_y=overlap_fraction+((512-423)/512);
else

    overlap_fraction_x=overlap_fraction;
    overlap_fraction_y=overlap_fraction;
    
end


dowritepythonfiles=true; % write files for python scripting

docuts=true;% Switch to cut by WD (before outlier removal)
WD_range=0.5e-3;

dooutliers=true; % Switch to remove outliers

doremovepoints=false; % Switch to particular points (after the two prior cuts)


% Must save the ImageSnapper XML file for "Focus Map", enter its filename here:
imagesnapperfocusmapfile='C:\Users\Tescan\Documents\LST12313_FocusMap.xml';
fprintf('Analyzing %s\n',imagesnapperfocusmapfile);

A=xml2struct(imagesnapperfocusmapfile); % https://www.mathworks.com/matlabcentral/fileexchange/28518-xml2struct

samplestruct=A.ImageSnapperProject.Samples.PointSample;
%samplestruct=A.ImageSnapperProject.Samples.RectangleSample; % For Rectangle
WD_initial=str2double(samplestruct.Attributes.WD);

% WD_min=WD_initial-WD_range;
% WD_max=WD_initial+WD_range;


%samplestruct=A.ImageSnapperProject.Samples.RectangleSample; % For Rectangle
%samplestruct=A.ImageSnapperProject.Samples.CircleSample; % For Circle



nsamples=numel(samplestruct);

if nsamples>1
    
    ZStage_=cell(nsamples,1);
    
    for ii=1:nsamples
        ZStage_{ii}=samplestruct(ii).Attributes.Z;
    end
    
    if numel(unique(ZStage_))~=1
        disp('Different Z Stage values in XML!!!');
        return
    else
        ZStage=ZStage_{1};
        samplename=samplestruct(1).Attributes.Name;
    end
    
else % only one sample
    
    ZStage=samplestruct.Attributes.Z;
    samplename=samplestruct.Attributes.Name;
    
end

surfaceshape='poly22'; % poly11=plane, poly22=surface quadratic

focusimagedir=[A.ImageSnapperProject.Settings.Attributes.Path '\' samplename '\'];
parentpath=fileparts(A.ImageSnapperProject.Settings.Attributes.Path);

imagepath=[parentpath '\HighRes\']; % Output directory for high-res images

if (~exist(imagepath, 'dir'))
    mkdir(imagepath);
end

imageformat=A.ImageSnapperProject.Settings.Attributes.ImageFormat;


imagedir=dir([focusimagedir '*hdr']);
nfiles=numel(imagedir);

WD=zeros(nfiles,1);
YStage=zeros(nfiles,1);
XStage=zeros(nfiles,1);


for ii=1:nfiles
    
    fid = fopen([imagedir(ii).folder '/' imagedir(ii).name]);
    dd = textscan(fid,'%s');
    dd=dd{1};
    fclose(fid);
    
    wdstr='WD=';
    ss=strfind(dd,wdstr);
    index = false(1, numel(ss));
    for k = 1:numel(ss)
        if numel(ss{k} == 1)==0
            index(k) = 0;
        else
            index(k) = 1;
        end
    end
    ll=dd{index};
    WD(ii)=str2double(ll((numel(wdstr)+1):end));
    
    xstr='StageX=';
    ss=strfind(dd,xstr);
    index = false(1, numel(ss));
    for k = 1:numel(ss)
        if numel(ss{k} == 1)==0
            index(k) = 0;
        else
            index(k) = 1;
        end
    end
    ll=dd{index};
    XStage(ii)=str2double(ll((numel(xstr)+1):end));
    
    ystr='StageY=';
    ss=strfind(dd,ystr);
    index = false(1, numel(ss));
    for k = 1:numel(ss)
        if numel(ss{k} == 1)==0
            index(k) = 0;
        else
            index(k) = 1;
        end
    end
    ll=dd{index};
    YStage(ii)=str2double(ll((numel(ystr)+1):end));
    
end

%rr=randi(numel(XStage),15,1);

%convexhullindices=convhull(XStage(rr),YStage(rr));
shrinkfactor=0.5; % default = 0.5
boundaryindices=boundary(XStage,YStage,shrinkfactor);
%%
% Add boundary indices:
%boundaryindices=[boundaryindices;[58,59,60,61,62,63,64,65,68,66,67,70,26,25,24,23,22,21,19,18,17,59,16,15,14,13]'];

% xv=(min(XStage)-range(XStage)*0.1):1e-6*fovx*(1-overlap_fraction):(max(XStage)+range(XStage)*0.1);
% yv=(min(YStage)-range(YStage)*0.1):1e-6*fovy*(1-overlap_fraction):(max(YStage)+range(YStage)*0.1);

gridspacex=(1e-6*fovx*(1-overlap_fraction_x));
gridspacey=(1e-6*fovy*(1-overlap_fraction_y));
%inPoints = polygrid(XStage(boundaryindices),YStage(boundaryindices),mean([gridspacex,gridspacey])^-2);
%xv = unique(inPoints(:,1));
xv = [min(XStage)-(10:-1:1)*gridspacex, min(XStage):gridspacex:max(XStage), max(XStage)+(1:10)*gridspacex]+gridspacex/2;
if dowriteEDSmaps
    xv=xv+1e-6*fovx*0.19;
end
%yv = unique(inPoints(:,2));
yv = [min(YStage)-(10:-1:1)*gridspacey, min(YStage):gridspacey:max(YStage), max(YStage)+(1:10)*gridspacey]+gridspacey/2;
thresh_dist=increaseborderfactor*sqrt(2)*0.5*mean([gridspacex,gridspacey]);
%gridvec=[];
R=zeros(numel(yv),numel(xv));
for ii=1:numel(xv)
    for jj=1:numel(yv)
        d_min = p_poly_dist(xv(ii), yv(jj), XStage(boundaryindices), YStage(boundaryindices));
        if d_min<=thresh_dist || inpolygon(xv(ii),yv(jj),XStage(boundaryindices),YStage(boundaryindices))
            R(jj,ii)=1;
            %gridvec=[gridvec [Cx(ii);Cy(jj)]];
        end
    end
end

[rows, columns] = find(R);
topRow = min(rows);
bottomRow = max(rows);
leftColumn = min(columns);
rightColumn = max(columns);
R = R(topRow:bottomRow, leftColumn:rightColumn);
yv=yv(topRow:bottomRow);
xv=xv(leftColumn:rightColumn);


nxv=numel(xv);
nyv=numel(yv);

xx1=interp1(xv,1:nxv,XStage(boundaryindices)');
yy1=interp1(yv,1:nyv,YStage(boundaryindices)');


% xx1=interp1(xv,1:nxv,XStage(rr(convexhullindices))');
% yy1=interp1(yv,1:nyv,YStage(rr(convexhullindices))');
% xx1=interp1(xv,1:nxv,XStage(boundaryindices)');
% yy1=interp1(yv,1:nyv,YStage(boundaryindices)');
% 
% R_=poly2mask(xx1,yy1,nyv,nxv);
% Rr = imdilate(R_, strel('disk', 3));
% %Rr=imresize(R_,increaseborderfactor);
% center_=size(R_)/2;
% centerr=size(Rr)/2;
% Rr=imtranslate(Rr,fliplr(center_-centerr)+[0,1]);
% if increaseborderfactor~=1
%     %R=Rr(cornercrop(1):cornercrop(1)-1+size(R_,1),cornercrop(2):cornercrop(2)-1+size(R_,2));
%     R=Rr(1:size(R_,1),1:size(R_,2));
% else
%     R=R_;
% end
figure;
imagesc(R); axis image; colormap(gray)
axis on
hold on
plot(xx1,yy1,'-r','LineWidth',2)
hold off

figure;
plot(XStage,YStage,'.',XStage(boundaryindices),YStage(boundaryindices),'*'); axis image
%%
if docuts
    
    WD_min=mean(WD)-WD_range;
    WD_max=mean(WD)+WD_range;

    
    cut_index=WD>WD_min & WD<WD_max;
    
    WD=WD(cut_index);
    XStage=XStage(cut_index);
    YStage=YStage(cut_index);
    
end



if doremovepoints
    
    badpointindices=[1,170];
    
    WD(badpointindices)=[];
    XStage(badpointindices)=[];
    YStage(badpointindices)=[];
    
end

fprintf('%d files, %d outliers\n',nfiles,nfiles-numel(WD));


sf=fit([XStage,YStage],-WD,surfaceshape);
residual=-WD-sf([XStage,YStage]);

if dooutliers
    
    %[~,hampel_index] = hampel(WD); % outler detection and removal (out of focus image removal)
    outlier_index = isoutlier(residual);
    
    WD=WD(~outlier_index);
    XStage=XStage(~outlier_index);
    YStage=YStage(~outlier_index);
    
end

fid=fopen('Poly.txt','w');
fprintf(fid,'%.15f,%.15f,%.15f,%.15f,%.15f,%.15f',sf.p00,sf.p10,sf.p01,sf.p20,sf.p11,sf.p02);
fclose(fid);

sfi=scatteredInterpolant(XStage,YStage,-WD,'linear','nearest');

[xq,yq] = ndgrid(XStage,YStage);

figure
plot(sf,[XStage,YStage],-WD)
shading interp
title('Curve Fit');
figure
plot(sf,[XStage,YStage],-WD)
hold on
surf(xq,yq,sfi(xq,yq))
colormap bone
shading interp
hold off
title('Interp Fit');


n_vec=[sf.p10;sf.p01;-1];
n_vec=-n_vec./norm(n_vec);
qfac=5*std(-WD);
if strcmp(surfaceshape,'poly11') % plot surface normal only for plane
    hold on
    quiver3(mean(XStage),mean(YStage),sf.p00+sf.p10*mean(XStage) + sf.p01*mean(YStage),qfac*n_vec(1),qfac*n_vec(2),qfac*n_vec(3),'Color','red','LineWidth',4)
    hold off
end
daspect([1 1 0.2])
%axis equal

ptsxyz=zeros(nxv*nyv,3);
coordn=zeros(nxv*nyv,2);
mask=zeros(nxv*nyv,1);


for jj=1:nyv
    % This will flip the direction every other line to minimize stage
    % movement:
    if mod(jj,2)==1
        iiloopvec=1:nxv;
    else
        iiloopvec=nxv:-1:1;
    end
    %
    for ii=1:nxv
        %
        % For surface fit:
         ptsxyz(ii+(jj-1)*nxv,:)=[xv(iiloopvec(ii)),yv(jj),-sf(xv(iiloopvec(ii)),yv(jj))];
        % For interp fit:
        % ptsxyz(ii+(jj-1)*nxv,:)=[xv(iiloopvec(ii)),yv(jj),-sfi(xv(iiloopvec(ii)),yv(jj))];
        %
        coordn(ii+(jj-1)*nxv,:)=[iiloopvec(ii),jj];
        mask(ii+(jj-1)*nxv)=R(jj,iiloopvec(ii));
        
    end
end


%if ~dowriteEDSmaps
    
fprintf('%d images\n',sum(mask));
    
res=writeImageSnapper(ptsxyz,ZStage,fovx,coordn,imagepath,imageformat,dowritepythonfiles,mask);
    
%end

if dowriteEDSmaps
    

    
    
    [~,CurrentPOS]=system('python C:\Users\Tescan\Documents\TescanGetWD.py');
    CurrentPOS=textscan(CurrentPOS, '%s', 'Delimiter','\n' );
    CurrentWD=CurrentPOS{1}{1};
    CurrentX=CurrentPOS{1}{2};
    CurrentY=CurrentPOS{1}{3};
    CurrentZ=CurrentPOS{1}{4};
    CurrentRotation=CurrentPOS{1}{5};
    CurrentTilt=CurrentPOS{1}{6};
    
    outputfile='MultiMaps.tmt';
    magnification=190836/fovx;
    stagetiltdegrees=CurrentTilt;
    stagerotationdegrees=CurrentRotation;
    amptime=0.96;
    ev_per_channel=5;
    mapmatrix=4;
    nPoints=512;
    dwelltime_microseconds=50;
    nLines=nPoints*fovy/fovx;
    total_map_time_minutes=80;
    
    
    total_map_time_microseconds=total_map_time_minutes*60*1e6;
    nFrames=total_map_time_microseconds./(nPoints.*nLines*dwelltime_microseconds);
    
    nFramesR=2.^ceil(log2(nFrames));
    
    fprintf('%d images\n',sum(mask));
    fprintf('Set FRAMES to %d\n',nFramesR)
    fprintf('Set to Current WD to %s\n',CurrentWD)
    fprintf('Set View Field (Tescan) to %d microns\n',fov)
    fprintf('Set AMP TIME and IMAGE SIZE\n')
    
    fprintf('%.0f hours = %.1f days\n',(nFramesR/nFrames).*sum(mask)*total_map_time_minutes/60,(nFramesR/nFrames).*sum(mask)*total_map_time_minutes/(24*60));
    
    MapPresetDwell=dwelltime_microseconds*1e-6;
        
    fmask=find(mask);
    
    % --------------------
    % Skip these maps because no BSE image:
    %load('IGD.mat','IGD'); fprintf('\nREMOVING %d X-RAY MAPS WITH IGD\n\n',sum(~IGD)); 
    %fmask=fmask(logical(IGD));
    % --------------------
    
    %  ONLY SELECT SOME OF THE MAPS:
%     middleind=round(size(ptsxyz,1)/2);
%     ind00=middleind-49:middleind+50;
%     ind00_complement=setdiff(fmask,ind00);
    
    nskip=0; % Skip this many maps
    [status,Ztotal]=writeEDAXMultifieldMaps(outputfile,magnification,ptsxyz(fmask,3),ptsxyz(fmask,1),ptsxyz(fmask,2),...
        ZStage,stagetiltdegrees,stagerotationdegrees,amptime,dwelltime_microseconds,ev_per_channel,mapmatrix,MapPresetDwell,nLines,nPoints,nskip,CurrentWD,CurrentZ);
    figure
    plot(Ztotal)
    xlabel('Map #')
    ylabel('Z (mm)')
    
end

fclose('all');

