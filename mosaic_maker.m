close all
clearvars;
delete(gcp('nocreate'));
%cols_array=10:98;
%; %32:41;
tic

pixelcountresize=25e6;
pixelcountresizewidth=10e3; % width in pixels

transformtype='affine';
transformtypev='affine';
mincontrastvalue=0.05; % default = 0.2
minqualityvalue=0.1	; % default = 0.1
maxdistancesep=0.1; % default = 1.5
% https://www.mathworks.com/help/vision/examples/feature-based-panoramic-image-stitching.html
outliertformerror_thresh=1;


dooutlierremoval=1;
dosavefigurematches=0;
dosavefigurematchesv=0;
dopause=0;



imgfolder='~/Documents/Acfer-TEMP2/';
samplename=fileparts(imgfolder);


D=dir([imgfolder '*.png']);
allnames={D.name};
nfilestotal=numel(allnames);
row_t=zeros(nfilestotal,1);
col_t=zeros(nfilestotal,1);
for ii=1:nfilestotal
    [~,fname_,~]=fileparts(allnames{ii});
    st_=str2double(strsplit(fname_,'_'));
    row_t(ii)=st_(1);
    col_t(ii)=st_(2);
end

colt=unique(sort(col_t));

% OLD:
%cols_array=colt';

%NEW: 
% colt_split=floor(numel(colt)/2);
% cols_array=[colt(colt_split+1:end)',fliplr(colt(1:colt_split)')];

% DEBUG:
%colt=[40:50]';
 colt_split_index=floor(numel(colt)/2);
 colt_split=colt(colt_split_index);
 cols_array=[colt(colt_split_index+1:end)',fliplr(colt(1:colt_split_index)')];
%colt_split=-inf;
%cols_array=colt';

ncols=numel(colt);

nrows=zeros(ncols,1);

rowsforcol=cell(ncols,1);

for ii=1:ncols
ind=col_t==colt(ii);  
nrows(ii)=sum(ind);
rowsforcol{ii}=row_t(ind);
end


%converttext=sprintf('/opt/local/bin/convert '); % sprintf([imgfolder '%03d_' col '.' imgext],rowsforcol{qq}(1))]);
%converttext0=sprintf('/opt/local/bin/convert ');
converttext0=sprintf('convert ');


converttexts=[];
converttext1=sprintf(' -background transparent -layers merge +repage strip%03d-%03d.png',min(cols_array),max(cols_array));

%tformtprev=affine2d(eye(3));

maxrows=[];
for qq=cols_array
maxrows=max([maxrows,rowsforcol{colt==qq}(end)]);
end


brightnessratio=nan(max(col_t),maxrows);
brightnessratiov=nan(max(col_t),maxrows);

errortform=zeros(max(col_t),maxrows);

rowmax=0;
for qq=cols_array
    rowmax=max([rowmax;rowsforcol{colt==qq}]);
end

colmax=max(cols_array);

matchedPointsPrev_v=cell(colmax,rowmax);
matchedPoints_v=cell(colmax,rowmax);
matchedPointsPrev=cell(colmax,rowmax);
matchedPoints=cell(colmax,rowmax);

matchedvcount=zeros(colmax,rowmax);
matchedcount=zeros(colmax,rowmax);
maxc_index_t=zeros(colmax,1);
maxvc_index_t=zeros(colmax,1);
imgshiftx=nan(colmax,rowmax);
imgshifty=nan(colmax,rowmax);


%poolobj = parpool;

for qq=cols_array
%for qq=1:ncols
    tic
    
    col=sprintf('%03d',qq);
 %   disp(col);
    nf=nrows(colt==qq);
    
 
     % DEBUGGING:
    %rowsforcol{colt==qq}=[16:50]';%[24:26]';

    rows_array=rowsforcol{colt==qq}'; %1:nf;

    imsz=2048;

xlim=zeros(nf,2);
ylim=zeros(nf,2);
It=zeros(nf,imsz,imsz,'uint16');
mincountpoints=2;
imgext='png';

       
% Tt=nan(nf,3,3);
% Tt(1,:,:)=eye(3);

I1=imread([imgfolder sprintf(['%03d_' col '.' imgext],rows_array(1))]);

width=size(I1,2);
height=size(I1,1);
overlap=0.20;
searchfac=1;

searchwidth=round(width*overlap*searchfac);
searchheight=round(height*overlap*searchfac);

It(1,:,:)=I1;


inliercount=0;

fprintf('Col %s Rows %03d – %03d\n',col,rows_array(1),rows_array(end))


%parfor kk=rows_array %2:nf
for kk=rows_array %2:nf


% imshow(I1); hold on; plot(points); hold off;

imname2=sprintf(['%03d_%03d.' imgext],kk,qq);

if qq>colt_split
imname2v=sprintf(['%03d_%03d.' imgext],kk,qq-1);
searchROI_1=[1 height-searchheight width searchheight];
searchROI_2=[1 1 width searchheight];
else
imname2v=sprintf(['%03d_%03d.' imgext],kk,qq+1);
searchROI_2=[1 height-searchheight width searchheight];
searchROI_1=[1 1 width searchheight];
end

I2=imread([imgfolder imname2]);
imname1=[];

if isfile([imgfolder imname1]) && kk>rows_array(1)
imname1=sprintf(['%03d_' col '.' imgext],kk-1);

I1=imread([imgfolder imname1]);
pointsPrevious  = detectFASTFeatures(I1,'ROI',[1 1 searchwidth height],'MinQuality',minqualityvalue,'MinContrast',mincontrastvalue); %,'MinContrast',0.01); %,'MinContrast',0.01); %detectSURFFeatures(I1,'ROI',[1 1 searchwidth height],'MetricThreshold',1,'NumScaleLevels',6,'NumOctaves',1);
[featuresPrevious, pointsPrevious]  = extractFeatures(I1,  pointsPrevious, 'Method','Block');

points = detectFASTFeatures(I2,'ROI',[width-searchwidth 1 searchwidth height],'MinQuality',minqualityvalue,'MinContrast',mincontrastvalue); %,'MinContrast',0.01); %,'MetricThreshold',100,'NumScaleLevels',6);
[features, points]  = extractFeatures(I2,  points, 'Method','Block');
indexPairs = matchFeatures(features, featuresPrevious, 'Unique', true);

disp(['Features = ' num2str(size(indexPairs,1))])

if size(indexPairs,1)>0
matchedPoints{qq,kk}  = points(indexPairs(:,1));
matchedPointsPrev{qq,kk} = pointsPrevious(indexPairs(:,2));

matchedcount(qq,kk)=matchedPoints{qq,kk}.Count;
else
    matchedcount(qq,kk)=0;
end

if matchedcount(qq,kk)>3
    [~, inlierDistorted, inlierOriginal, ~] = estimateGeometricTransform(matchedPoints{qq,kk}, matchedPointsPrev{qq,kk}, transformtype); %, 'MaxDistance', round(height/100));
    bbo=convhull(double(inlierOriginal.Location(:,1)),double(inlierOriginal.Location(:,2)));
    bbd=convhull(double(inlierDistorted.Location(:,1)),double(inlierDistorted.Location(:,2)));
    BW1=poly2mask(double(inlierOriginal.Location(bbo,1)),double(inlierOriginal.Location(bbo,2)),width,height);
    BW2=poly2mask(double(inlierDistorted.Location(bbd,1)),double(inlierDistorted.Location(bbd,2)),width,height);
    brightnessratio(kk,qq)=mean(double(I2(BW2)))./mean(double(I1(BW1)));
end

    if dosavefigurematches && matchedcount(qq,kk)>3
       
    f=figure('Visible','off');
    showMatchedFeatures(I1,I2,inlierOriginal,inlierDistorted);
    title([imname1 ' <--> ' imname2])
    saveas(f,[imname1(1:end-4) '-' imname2(1:end-4) '_match.png'])
    close(f);

        
    end
    
else
    
        matchedcount(qq,kk)=0;

end


if isfile([imgfolder imname2v]) && qq~=cols_array(1)

I1_=imread([imgfolder imname2v]);

pointsv_ = detectFASTFeatures(I1_,'ROI',searchROI_1,'MinQuality',minqualityvalue,'MinContrast',mincontrastvalue); %,'MinContrast',0.01); %,'MetricThreshold',100,'NumScaleLevels',6);    
[featuresv_, pointsv_]  = extractFeatures(I1_,  pointsv_, 'Method','Block');

pointsv = detectFASTFeatures(I2,'ROI',searchROI_2,'MinQuality',minqualityvalue,'MinContrast',mincontrastvalue); %,'MinContrast',0.01); %,'MetricThreshold',100,'NumScaleLevels',6);    
[featuresv, pointsv]  = extractFeatures(I2,  pointsv, 'Method','Block');

%imshow(I1_); hold on; plot(pointsv_); hold off;
%imshow(I1); hold on; plot(pointsv); hold off;


indexPairsv = matchFeatures(featuresv, featuresv_, 'Unique', true);

if numel(indexPairsv)>0
        
        matchedPoints_v{qq,kk}=pointsv(indexPairsv(:,1));
        matchedPointsPrev_v{qq,kk}=pointsv_(indexPairsv(:,2));
        
        matchedvcount(qq,kk)=matchedPointsPrev_v{qq,kk}.Count;
                
else

    matchedvcount(qq,kk)=0;
    
end


if matchedvcount(qq,kk)>3
    
    [~, inlierDistortedv, inlierOriginalv, ~] = estimateGeometricTransform(matchedPoints_v{qq,kk}, matchedPointsPrev_v{qq,kk}, transformtype); %, 'MaxDistance', round(height/100));
    bbo=convhull(double(inlierOriginalv.Location(:,1)),double(inlierOriginalv.Location(:,2)));
    bbd=convhull(double(inlierDistortedv.Location(:,1)),double(inlierDistortedv.Location(:,2)));
    BW1=poly2mask(double(inlierOriginalv.Location(bbo,1)),double(inlierOriginalv.Location(bbo,2)),width,height);
    BW2=poly2mask(double(inlierDistortedv.Location(bbd,1)),double(inlierDistortedv.Location(bbd,2)),width,height);
    brightnessratiov(kk,qq)=mean(double(I2(BW2)))./mean(double(I1_(BW1)));
    
end
    

    if dosavefigurematchesv && matchedvcount(qq,kk)>3
       
    %tft_(qq-cols_array(1))=tf_;
    f=figure('Visible','off');
    showMatchedFeatures(I1_,I2,inlierOriginalv,inlierDistortedv);
    title('Matched inlier points');
    %num2str(inlierOriginal.Count);
    title([imname1 ' <--> ' imname2])
    saveas(f,[imname2v(1:end-4) '-' imname2(1:end-4) '_vmatch.png'])
    close(f);

        
    end
    
else
    matchedvcount(qq,kk)=0;
end

% matchedPointsPrev_v{53,52}

end


end

%delete(poolobj);





% figure;
% showMatchedFeatures(I1,I2,matchedOriginal,matchedDistorted);
% title('Putatively matched points (including outliers)');


 %  if exist('tformtprev')
     %   num2str(tformtprev__.T)
     
%      if qq==3
%          disp('here');
%      end
%%

% nn=round(0.5*numel(rowsforcol{colt==cols_array(1)}));
% startrow=rowsforcol{colt==cols_array(1)}(nn);
% maxc_index_t(cols_array(1))=startrow+1;
% startrow=rowsforcol{colt==cols_array(1)}(1);
% maxc_index_t(cols_array(1))=startrow+1;
% Tt(cols_array(1),startrow,:,:)=eye(3);
isoutlierall=zeros(max(col_t),maxrows);
nsigmaoutlier=10;

close all

if strcmp(transformtype,'affine') || strcmp(transformtype,'similarity')
tformt(max(col_t),maxrows)=affine2d(eye(3));
%tfti(max(col_t),maxrows)=affine2d(eye(3));
end

if strcmp(transformtype,'projective')
tformt(max(col_t),maxrows)=projective2d(eye(3));
%tfti(max(col_t),maxrows)=projective2d(eye(3));
end

Tt=nan(max(col_t),maxrows,3,3);



for qq=cols_array
    
    if qq<=colt_split
        
        qqp=qq+1;

    else
        
        qqp=qq-1;
        
    end
    
    [~,maxc_index_t(qq)]=max(matchedcount(qq,:));
        
    if qq==cols_array(1)   
               
    Tt(qq,maxc_index_t(qq)-1,:,:)=eye(3);
    tformt(qq,maxc_index_t(qq)-1).T=eye(3);
    loopvector=[maxc_index_t(qq):rowsforcol{colt==qq}(end),(maxc_index_t(qq)-2):-1:rowsforcol{colt==qq}(1)];
    
    else
        
        [~,maxvc_index_t(qq)]=max(matchedvcount(qq,:).*(~isoutlierall(qqp,:)));
        loopvector=[maxvc_index_t(qq):rowsforcol{colt==qq}(end),(maxvc_index_t(qq)-1):-1:rowsforcol{colt==qq}(1)];
    
    end
    

    
    % sum(sum(isnan(squeeze(Tt(qq,kk-1,:,:)))))
    
    
    
    %loopvector=rowsforcol{colt==qq}'
    
  for kk=loopvector  %[maxc_index_t(qq):rowsforcol{colt==qq}(end)] %rowsforcol{colt==qq}' %startrow+1:rows_array(end)
  
  ll=[];
  mp_ll=[];    
  
if (kk+1)<=maxrows && ((qq==cols_array(1) && kk<maxc_index_t(qq)) ||  (qq~=cols_array(1) && kk<maxvc_index_t(qq))) && ~isnan(sum(sum(Tt(qq,kk+1,:,:),3),4))
if (matchedcount(qq,kk+1)>=mincountpoints)
    ll=transformPointsForward(tformt(qq,kk+1),matchedPoints{qq,kk+1}.Location);    
    mp_ll=matchedPointsPrev{qq,kk+1}.Location;
end
end

if kk>1 && ((matchedcount(qq,kk)>=mincountpoints) && ((qq==cols_array(1) && kk>=maxc_index_t(qq)) ||  (qq~=cols_array(1) && kk>=maxvc_index_t(qq))) && ~isnan(sum(sum(Tt(qq,kk-1,:,:),3),4)))
    ll=transformPointsForward(tformt(qq,kk-1),matchedPointsPrev{qq,kk}.Location);    
    mp_ll=matchedPoints{qq,kk}.Location;
end

  
 llp=[];
 mp_llp=[];
 
 if matchedvcount(qq,kk)>=mincountpoints && ~isoutlierall(qqp,kk)
    llp=transformPointsForward(tformt(qqp,kk),matchedPointsPrev_v{qq,kk}.Location);
    mp_llp=matchedPoints_v{qq,kk}.Location;
 end
   
 fprintf('(%d,%d): %d,%d / (%d,%d), Outlier: %d|%d\n',qq,kk,size(ll,1),size(llp,1),matchedcount(qq,kk),matchedvcount(qq,kk),isoutlierall(qq,kk),isoutlierall(qqp,kk))

 
 if (numel([mp_ll;mp_llp])+numel([ll;llp])) > 0
 
                if matchedvcount(qq,kk)>=mincountpoints
                     transformtype=transformtypev;
                 end
                 
     
[tformt(qq,kk), inlierDistorted, inlierOriginal, status] = estimateGeometricTransform(...
    [mp_ll;mp_llp], [ll;llp], transformtype, 'Confidence', 99.999, 'MaxNumTrials', 1e5, 'MaxDistance', maxdistancesep); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
%[tfti(qq,kk), ~, ~, ~] = estimateGeometricTransform(...
%     [ll;llp], [mp_ll;mp_llp], transformtype, 'Confidence', 99.999, 'MaxNumTrials', 1e5, 'MaxDistance', maxdistancesep); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));



if numel(inlierDistorted)>=mincountpoints
U=transformPointsForward(tformt(qq,kk),inlierDistorted);
errortform(qq,kk)=sum(sqrt((U(:,1)-inlierOriginal(:,1)).^2+(U(:,2)-inlierOriginal(:,2)).^2))./size(U,1);
else
errortform(qq,kk)=inf;
brightnessratio(kk,qq)=nan;
end

Tt(qq,kk,:,:)=tformt(qq,kk).T;

% sprintf('%03d,%03d\n',kk,qq) 
% num2str(tformt(qq,kk).T)

     
 end

% num2str(tformtprev.T)
    

kkprev=kk;

end
    


if dooutlierremoval

vv=rowsforcol{colt==qq};    


[~,outlierindex11] = hampel(Tt(qq,vv,1,1),nf,nsigmaoutlier);
[~,outlierindex12] = hampel(Tt(qq,vv,1,2),nf,nsigmaoutlier);
[~,outlierindex21] = hampel(Tt(qq,vv,2,1),nf,nsigmaoutlier);
[~,outlierindex22] = hampel(Tt(qq,vv,2,2),nf,nsigmaoutlier);
[~,outlierindex31] = hampel(Tt(qq,vv,3,1),nf,nsigmaoutlier);
[~,outlierindex32] = hampel(Tt(qq,vv,3,2),nf,nsigmaoutlier);

isnanindex = isnan(sum(sum(Tt(qq,vv,:,:),4),3));

outliertformerror = errortform(qq,vv) > outliertformerror_thresh; %3e-3;

vvo=vv((outlierindex11 | outlierindex12 | outlierindex21 | outlierindex22 | outlierindex31 | outlierindex32) & ~isnanindex);
if numel(vvo)>0
fprintf('Outliers: ');    
for bb=1:numel(vvo)
fprintf(['%03d_%03d.' imgext ' '],vvo(bb),qq)
end
fprintf('\n');
end

outlierindex= isnanindex | outlierindex11 | outlierindex12 | outlierindex21 | outlierindex22 | outlierindex31 | outlierindex32 | outliertformerror;

isoutlierall(qq,vv(outlierindex))=ones(sum(outlierindex),1);
Tt(qq,vv(outlierindex),:,:)=nan(1,sum(outlierindex),3,3);

fprintf('Col %d, Outliers/Total: %d/%d\n',qq, sum(outlierindex),numel(outlierindex))

end

imgshiftx_=nan(max(rowsforcol{colt==qq}),1);
imgshifty_=nan(max(rowsforcol{colt==qq}),1);

tf_=affine2d;

for kk=rowsforcol{colt==qq}' 

if sum(sum(isfinite(Tt(qq,kk,:,:))))==9
    tf_.T=squeeze(Tt(qq,kk,:,:));
    [imgshiftx_(kk),imgshifty_(kk)]=transformPointsForward(tf_,0,0);
else
    imgshiftx_(kk)=nan;
    imgshifty_(kk)=nan;
end

end

%figure
plot(imgshiftx_,imgshifty_,'s'); title(['Column ' num2str(qq)])

end

T2=Tt;
for ii=1:3
for jj=1:2
    T2(:,:,ii,jj)=inpaint_nans(Tt(:,:,ii,jj),1);
end
end
T2(:,:,1,3)=zeros(size(T2,1),size(T2,2));
T2(:,:,2,3)=zeros(size(T2,1),size(T2,2));
T2(:,:,3,3)=ones(size(T2,1),size(T2,2));


for qq=cols_array
for kk=rowsforcol{colt==qq}'

 tformt(qq,kk).T=squeeze(T2(qq,kk,:,:));
 
end
end


xlim=zeros(max(col_t),maxrows,2);
ylim=zeros(max(col_t),maxrows,2);
for qq=cols_array
for kk=rowsforcol{colt==qq}'
    [xlim(qq,kk,:), ylim(qq,kk,:)] = outputLimits(tformt(qq,kk),  [1 width], [1 height]);
end
end

% Find the minimum and maximum output limits
xMin = min([1; xlim(:)]);
xMax = max([width; xlim(:)]);

yMin = min([1; ylim(:)]);
yMax = max([height; ylim(:)]);

% Width and height of panorama.
widthp  = round(xMax - xMin);
heightp = round(yMax - yMin);


fid=fopen(sprintf('imageinformation_%03d_%03d.txt',min(cols_array),max(cols_array)),'w+');

fprintf(fid,'%d,%d,%d,%d\n',xMin,xMax,yMin,yMax);

converttexts=[];
for qq=cols_array
for kk=rowsforcol{colt==qq}' %2:nf

fname=[imgfolder sprintf(['%03d_%03d.' imgext],kk,qq)];

[imgshiftx(qq,kk),imgshifty(qq,kk)]=transformPointsForward(tformt(qq,kk),0,0);

converttexts=[converttexts sprintf(' \\( %s -alpha set -virtual-pixel transparent +distort AffineProjection ''%7.6f,%7.6f,%7.6f,%7.6f,%7.6f,%7.6f'' \\)',fname,tformt(qq,kk).T(1,1),tformt(qq,kk).T(1,2),tformt(qq,kk).T(2,1),tformt(qq,kk).T(2,2),tformt(qq,kk).T(3,1),tformt(qq,kk).T(3,2))]; %,sxv,rxv,ryv,syv,txv,tyv)];

fprintf(fid,'%s,%7.6f,%7.6f,%7.6f,%7.6f,%7.6f,%7.6f\n',fname,tformt(qq,kk).T(1,1),tformt(qq,kk).T(1,2),tformt(qq,kk).T(2,1),tformt(qq,kk).T(2,2),tformt(qq,kk).T(3,1),tformt(qq,kk).T(3,2));
end
end

fclose(fid);


cc=distinguishable_colors(numel(cols_array));
cc2=distinguishable_colors(maxrows);
figure('visible','on');
set(gcf,'Position',[0,0,2000,2000])
hold on
for jj=1:maxrows
  plot(imgshiftx(:,jj),imgshifty(:,jj),'-','Color',cc2(jj,:)); 
end
for jj=cols_array
  plot(imgshiftx(jj,:),imgshifty(jj,:),'-','Color',cc(jj-min(cols_array)+1,:)); 
end
for ii=cols_array
indout=logical(isoutlierall(ii,rowsforcol{colt==ii}));
plot(imgshiftx(ii,rowsforcol{colt==ii}(~indout)),imgshifty(ii,rowsforcol{colt==ii}(~indout)),'s','Color',cc(ii-min(cols_array)+1,:),'MarkerFaceColor',cc(ii-min(cols_array)+1,:)); 
plot(imgshiftx(ii,rowsforcol{colt==ii}(indout)),imgshifty(ii,rowsforcol{colt==ii}(indout)),'s','Color',cc(ii-min(cols_array)+1,:),'MarkerFaceColor','w'); 
if maxvc_index_t(ii)>0
plot(imgshiftx(ii,maxvc_index_t(ii)),imgshifty(ii,maxvc_index_t(ii)),'O','Color',cc(ii-min(cols_array)+1,:),'MarkerFaceColor',cc(ii-min(cols_array)+1,:),'MarkerSize',7);
end
text(max(imgshiftx(ii,rowsforcol{colt==ii}))+width,mean(imgshifty(ii,rowsforcol{colt==ii})),sprintf('%03d',ii))
end
plot(imgshiftx(cols_array(1),maxc_index_t(cols_array(1))),imgshifty(cols_array(1),maxc_index_t(cols_array(1))),'p','Color','k','MarkerFaceColor','k','MarkerSize',11);
hold off
xlabel('Image X (pixels)')
ylabel('Image Y (pixels)')
box on
mosaiclocationsfname=sprintf('mosaic_locations_%03d_%03d.png',min(cols_array),max(cols_array));
axis equal
set(gca,'TickLength',[0 0])
print('-depsc',[mosaiclocationsfname(1:end-3) 'eps']);
[~,~]=system(['epstopdf ' mosaiclocationsfname(1:end-3) 'eps ' mosaiclocationsfname(1:end-3) 'pdf']);

figure('visible','on');
hold on
for ii=cols_array
plot(rowsforcol{colt==ii},errortform(ii,rowsforcol{colt==ii}),'-O','Color',cc(ii-min(cols_array)+1,:),'MarkerFaceColor',cc(ii-min(cols_array)+1,:)); 
end
ax=axis;
hold off
xlabel('Image Number')
ylabel('Transform Error')
box on
set(gca, 'YScale', 'log')
tformfname=sprintf('tform_errors_%03d_%03d.png',min(cols_array),max(cols_array));

print('-dpng','-r300',tformfname);

imwrite(uint8(255*matchedcount./max(matchedcount(:))),sprintf('matchedcount_%03d_%03d.png',min(cols_array),max(cols_array)));
imwrite(uint8(255*matchedvcount./max(matchedvcount(:))),sprintf('matchedvcount_%03d_%03d.png',min(cols_array),max(cols_array)));
imwrite(uint8(255*(brightnessratio./max(brightnessratio(:))).^1),sprintf('brightnessratio_%03d_%03d.png',min(cols_array),max(cols_array)));
imwrite(uint8(255*(brightnessratiov./max(brightnessratiov(:))).^1),sprintf('brightnessratiov_%03d_%03d.png',min(cols_array),max(cols_array)));


%%
if dopause
pause
end

%imcall=['echo ' converttexts converttext1 '| xargs ' converttext0]
imcall=[converttext0 converttexts converttext1]; % num2str(heightp)

save(sprintf('mosaic_maker_%03d-%03d.mat',min(cols_array),max(cols_array)))

%system(imcall);
imresizecall=sprintf('/usr/local/bin/vipsthumbnail strip%03d-%03d.png --size %d',min(cols_array),max(cols_array),pixelcountresizewidth);
%system(imresizecall);

fid = fopen('imcall.sh','wt');
fprintf(fid, '%s\n\n%s','#!/bin/bash',imcall);
fclose(fid);

toc

% system(['python ~/Documents/affinetransform.py' samplename]);

