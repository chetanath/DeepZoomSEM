close all
clearvars;
delete(gcp('nocreate'));
%cols_array=10:98;
%; %32:41;

pixelcountresize=25e6;
pixelcountresizewidth=10e3; % width in pixels

transformtype='affine';
transformtypev='affine';
nsigmaoutlier=inf;
mincontrastvalue=0.05; % default = 0.2
minqualityvalue=0.1	; % default = 0.1
maxdistancesep=0.1; % default = 1.5
% https://www.mathworks.com/help/vision/examples/feature-based-panoramic-image-stitching.html
outliertformerror_thresh=1;


dooutlierremoval=1;
dosavefigurematches=0;
dosavefigurematchesv=0;
dosavefigurematchesstrip=0;
dowritestrips=0;
dopause=0;
dopointtransform=1;



imgfolder='~/Documents/Acfer-TEMP/';
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

if strcmp(transformtype,'affine') || strcmp(transformtype,'similarity')
tformtprev_t(max(col_t),maxrows)=affine2d(eye(3));
end

if strcmp(transformtype,'projective')
tformtprev_t(max(col_t),maxrows)=projective2d(eye(3));
end

Tt=nan(max(col_t),maxrows,3,3);

isoutlierall=zeros(max(col_t),maxrows,3,3);

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
maxvc_index_t=zeros(colmax,1);
imgshiftx=zeros(colmax,rowmax);
imgshifty=zeros(colmax,rowmax);


poolobj = parpool;

for qq=cols_array
%for qq=1:ncols
    tic
    
    col=sprintf('%03d',qq);
    disp(col);
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

sprintf('Rows : %03d -- %03d',rows_array(1),rows_array(end))


parfor kk=rows_array %2:nf
%for kk=rows_array %2:nf


% imshow(I1); hold on; plot(points); hold off;

imname2=sprintf(['%03d_%03d.' imgext],kk,qq);

%imname_total(qq,kk)=imname2;

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

if kk>rows_array(1)
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

    if dosavefigurematches && matchedcount(qq,kk)>0
       
    f=figure('Visible','off');
    showMatchedFeatures(I1,I2,inlierOriginal,inlierDistorted);
    title([imname1 ' <--> ' imname2])
    saveas(f,[imname1(1:end-4) '-' imname2(1:end-4) '_match.png'])
    close(f);

        
    end
    
% [~, inlierDistorted, inlierOriginal, ~] = estimateGeometricTransform(...
%     matchedPoints, matchedPointsPrev, transformtype, 'Confidence', 99.999, 'MaxNumTrials', 1e5, 'MaxDistance', round(width/100));

% showMatchedFeatures(I1,I2,inlierOriginal,inlierDistorted);


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
        
%         if qq==53 && kk==52
%             disp('stop')
%         end
                
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
    

    if dosavefigurematchesv && matchedvcount(qq,kk)>0
       
    [~, inlierDistortedv, inlierOriginalv, ~] = estimateGeometricTransform(matchedPoints_v{qq,kk}, matchedPointsPrev_v{qq,kk}, transformtype, 'Confidence', 99.99, 'MaxNumTrials', 1e5, 'MaxDistance', maxdistancesep); %, 'MaxDistance', round(height/100));
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

delete(poolobj);





% figure;
% showMatchedFeatures(I1,I2,matchedOriginal,matchedDistorted);
% title('Putatively matched points (including outliers)');


 %  if exist('tformtprev')
     %   num2str(tformtprev__.T)
     
%      if qq==3
%          disp('here');
%      end
%%
cols_array_middle=floor(mean(cols_array));
rows_array_middle=floor(mean(rows_array));

% nn=round(0.5*numel(rowsforcol{colt==cols_array(1)}));
% startrow=rowsforcol{colt==cols_array(1)}(nn);
% maxvc_index_t(cols_array(1))=startrow+1;
% startrow=rowsforcol{colt==cols_array(1)}(1);
% maxvc_index_t(cols_array(1))=startrow+1;
% Tt(cols_array(1),startrow,:,:)=eye(3);

for qq=cols_array
        
    if qq==cols_array(1)   
            
    [~,maxvc_index_t(qq)]=max(matchedcount(qq,:));
    Tt(qq,maxvc_index_t(qq)-1,:,:)=eye(3);
    mnoffset=2;
    
    else
        
    [~,maxvc_index_t(qq)]=max(matchedvcount(qq,:));
    mnoffset=1;

        
    end
    
    if qq<=colt_split
        
        qqp=qq+1;

    else
        
        qqp=qq-1;
        
    end
    
    % sum(sum(isnan(squeeze(Tt(qq,kk-1,:,:)))))
    
    
    loopvector=[maxvc_index_t(qq):rowsforcol{colt==qq}(end),(maxvc_index_t(qq)-mnoffset):-1:rowsforcol{colt==qq}(1)];
    
    kkprev=loopvector(1)-1;
    
  for kk=loopvector %[maxvc_index_t(qq):rowsforcol{colt==qq}(end)] %rowsforcol{colt==qq}' %startrow+1:rows_array(end)
    
  if matchedcount(qq,kk)>=mincountpoints && ~isnan(sum(sum(Tt(qq,kkprev,:,:),3),4)) && (qq==cols_array(1) || (qq~=cols_array(1) && kk>maxvc_index_t(qq)) ) %&& errortform(qq,kk-1)<outliertformerror_thresh %&& (qq==cols_array(1) || (qq>cols_array(1) && kk>rows_array(startrowindex)+1)) % sum(sum(isnan(squeeze(Tt(qq,kk-1,:,:)))))==0 %&& (qq==cols_array(1) || (qq>cols_array(1) && kk>rows_array(startrowindex)+1))
    ll=transformPointsForward(tformtprev_t(qq,kkprev),matchedPointsPrev{qq,kk}.Location);
    mp_ll=matchedPoints{qq,kk}.Location;
  else
     ll=[];
     mp_ll=[];
  end
    
 if matchedvcount(qq,kk)>=mincountpoints && ~isoutlierall(qqp,kk)
    llp=transformPointsForward(tformtprev_t(qqp,kk),matchedPointsPrev_v{qq,kk}.Location);
    mp_llp=matchedPoints_v{qq,kk}.Location;
    else
        llp=[];
        mp_llp=[];
 end
 
    
 
 if (numel([mp_ll;mp_llp])+numel([ll;llp])) > 0
 
                if matchedvcount(qq,kk)>=mincountpoints
                     transformtype=transformtypev;
                 end
                 
     
[tformtprev_t(qq,kk), inlierDistorted, inlierOriginal, status] = estimateGeometricTransform(...
    [mp_ll;mp_llp], [ll;llp], transformtype, 'Confidence', 99.999, 'MaxNumTrials', 1e5, 'MaxDistance', maxdistancesep); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));

if numel(inlierDistorted)>=mincountpoints
U=transformPointsForward(tformtprev_t(qq,kk),inlierDistorted);
errortform(qq,kk)=sum(sqrt((U(:,1)-inlierOriginal(:,1)).^2+(U(:,2)-inlierOriginal(:,2)).^2))./size(U,1);
else
errortform(qq,kk)=inf;
brightnessratio(kk,qq)=nan;
end

Tt(qq,kk,:,:)=tformtprev_t(qq,kk).T;

% sprintf('%03d,%03d\n',kk,qq) 
% num2str(tformtprev_t(qq,kk).T)

     
 end

% num2str(tformtprev.T)
    
% Tt(kk,:,:)=Tt(kk-1,:,:);
%

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

outlierindex= isnanindex | outlierindex11 | outlierindex12 | outlierindex21 | outlierindex22 | outlierindex31 | outlierindex32 | outliertformerror;

isoutlierall(qq,vv(outlierindex))=ones(sum(outlierindex),1);
Tt(qq,vv(outlierindex),:,:)=nan(1,sum(outlierindex),3,3);

fprintf('Outliers/Total: %d/%d\n',sum(outlierindex),numel(outlierindex))

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

 tformtprev_t(qq,kk).T=squeeze(T2(qq,kk,:,:));
 
end
end

% BEGIN OLD OUTLIER REMOVAL
% nanindex=isnan(Tt(qq,vv,3,1));
% 
% outlierfrac=sum(outlierindex)./numel(outlierindex);
% 
% nanfrac=sum(nanindex)./numel(outlierindex);
% 
% nonanoutlierfrac=sum(outlierindex & ~nanindex)./numel(outlierindex);
% 
% disp(['Outliers (All) (NaNs) = ' num2str(outlierfrac) ' , ' num2str(nanfrac)]);
% 
% if sum(~nanindex)>1 && sum(~outlierindex)>1
%     
% for kk=1:3
%     for jj=1:2
% %    pp=polyfit(vv(~outlierindex),squeeze(Tt(ii,vv(~outlierindex),3,1))',1);    polyval(pp,vv(outlierindex))
%     Tt(qq,vv(outlierindex),kk,jj)=interp1(vv(~outlierindex),Tt(qq,vv(~outlierindex),kk,jj),vv(outlierindex),'linear','extrap');
%     end
% end
% Tt(qq,vv,:,3)=repmat([0,0,1],numel(vv),1);
% 
% if dopointtransform    
% for jj=1:numel(vv)
%     tformtprev_t(qq,vv(jj)).T=squeeze(Tt(qq,vv(jj),:,:));
% end 
% else
%     tformtprev_t(qq,vv(1)).T=eye(3);
%     for jj=2:numel(vv)
%     tformtprev_t(qq,vv(jj)).T=squeeze(Tt(qq,vv(jj),:,:))*tformtprev_t(qq,vv(jj-1)).T;
%     end 
% end
% 
% elseif sum(~nanindex)<=1
% 
%    % Tt(ii,1,:,:)=eye(3);
%     for jj=find(nanindex)
%         Tt(qq,vv(jj),:,:)=eye(3);
%         Tt(qq,vv(jj),3,1)=-vv(jj)*imsz*(1-overlap);
%     end
% 
% if dopointtransform    
% for jj=1:numel(vv)
%     tformtprev_t(qq,vv(jj)).T=squeeze(Tt(qq,vv(jj),:,:));
% end 
% else
%     tformtprev_t(qq,vv(1)).T=eye(3);
%     for jj=2:numel(vv)
%     tformtprev_t(qq,vv(jj)).T=squeeze(Tt(qq,vv(jj),:,:))*tformtprev_t(qq,vv(jj-1)).T;
%     end 
% end
%     
% end
%     
% 
% else
%     
%     disp('Outlier Removal Off');
%     
%     end
% 
%     end
% END OLD OUTLIER REMOVAL       
      
        
    
% 
% else
%     
% Tt(kk,:,:)=Tt(kk-1,:,:);
% 
 
%     
%     tformt(kk).T = squeeze(Tt(kk,:,:)) * tformt(kk-1).T;

%num2str(squeeze(Tt(kk,:,:)))

%tform(kk).T(3,3)
%inlierDistorted.Count

    [sz1,sz2]=size(tformtprev_t);
visptext=[];

xlim=zeros(max(col_t),maxrows,2);
ylim=zeros(max(col_t),maxrows,2);
for qq=cols_array
for kk=rowsforcol{colt==qq}'
    [xlim(qq,kk,:), ylim(qq,kk,:)] = outputLimits(tformtprev_t(qq,kk),  [1 width], [1 height]);
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

%fprintf(fid,'%d,%d\n',widthp,heightp);
fprintf(fid,'%d,%d,%d,%d\n',xMin,xMax,yMin,yMax);

converttexts=[];
for qq=cols_array
for kk=rowsforcol{colt==qq}' %2:nf

        %file1='stripw128.tif';
% Row:
% file1=[imagepath sprintf(['/%0' num2str(formatc0) 'd_%0' num2str(formatc0) 'd' ext],row,jj)];
%fname=sprintf([imgfolder '%03d_%03d.' imgext],kk,qq);
fname=[imgfolder sprintf(['%03d_%03d.' imgext],kk,qq)];

[imgshiftx(qq,kk),imgshifty(qq,kk)]=transformPointsForward(tformtprev_t(qq,kk),0,0);

converttexts=[converttexts sprintf(' \\( %s -alpha set -virtual-pixel transparent +distort AffineProjection ''%7.6f,%7.6f,%7.6f,%7.6f,%7.6f,%7.6f'' \\)',fname,tformtprev_t(qq,kk).T(1,1),tformtprev_t(qq,kk).T(1,2),tformtprev_t(qq,kk).T(2,1),tformtprev_t(qq,kk).T(2,2),tformtprev_t(qq,kk).T(3,1),tformtprev_t(qq,kk).T(3,2))]; %,sxv,rxv,ryv,syv,txv,tyv)];
%converttexts=[converttexts sprintf(' \\( %s -alpha set -virtual-pixel transparent +distort AffineProjection ''%7.6f,%7.6f,%7.6f,%7.6f,%7.6f,%7.6f'' \\)',fname,tformtp(kk).T(1,1),tformtp(kk).T(1,2),tformtp(kk).T(2,1),tformtp(kk).T(2,2),tformtp(kk).T(3,1),tformtp(kk).T(3,2))]; %,sxv,rxv,ryv,syv,txv,tyv)];
%T__=tformtprev_t(qq,kk).T;
fprintf(fid,'%s,%7.6f,%7.6f,%7.6f,%7.6f,%7.6f,%7.6f\n',fname,tformtprev_t(qq,kk).T(1,1),tformtprev_t(qq,kk).T(1,2),tformtprev_t(qq,kk).T(2,1),tformtprev_t(qq,kk).T(2,2),tformtprev_t(qq,kk).T(3,1),tformtprev_t(qq,kk).T(3,2));
end
end

fclose(fid);


cc=distinguishable_colors(numel(cols_array));
figure('visible','on');
hold on
for ii=cols_array
plot(imgshiftx(ii,rowsforcol{colt==ii}),imgshifty(ii,rowsforcol{colt==ii}),'-s','Color',cc(ii-min(cols_array)+1,:)); 
plot(imgshiftx(ii,maxvc_index_t(ii)),imgshifty(ii,maxvc_index_t(ii)),'-*','Color',cc(ii-min(cols_array)+1,:));
end
ax=axis;
%plot([startrow,startrow],[ax(3),ax(4)],'k-','LineWidth',2,'Color',[0.5,0.5,0.5]);
hold off
xlabel('Image X (pixels)')
ylabel('Image Y (pixels)')
box on
mosaiclocationsfname=sprintf('mosaic_locations_%03d_%03d.png',min(cols_array),max(cols_array));

print('-dpng','-r300',mosaiclocationsfname);

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


%%
if dopause
pause
end

%imcall=['echo ' converttexts converttext1 '| xargs ' converttext0]
imcall=[converttext0 converttexts converttext1] % num2str(heightp)
system(imcall);
imresizecall=sprintf('/usr/local/bin/vipsthumbnail strip%03d-%03d.png --size %d',min(cols_array),max(cols_array),pixelcountresizewidth);
system(imresizecall);

fid = fopen('imcall.sh','wt');
fprintf(fid, '%s\n\n%s','#!/bin/bash',imcall);
fclose(fid);


% system(['python ~/Documents/affinetransform.py' samplename]);

