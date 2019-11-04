tic

% rsync -avz ~/Data/Tescan/Ryan/MurchisonThickSection/HighRes/ ~/Documents/MurchisonThickSection
% ~/Documents/MoveTescanImages.sh ~/Documents/MurchisonThickSection

close all
clearvars;
%delete(gcp('nocreate'));

%cols_array=10:98;
%; %32:41;


imgfolder='~/Documents/Acfer-TEMP2/';
%imgfolder='~/Data/deepzoom/CumberlandFalls_24_C3/';
overlap=0.20;


pixelcountresize=25e6;
pixelcountresizewidth=10e3; % width in pixels
N_Histogram_Bins=64;
offset_tolerance=inf; %0.1;
tform_tolerance=inf; %0.04;
npolynomial=2;

transformtype='similarity'; %'affine';
transformtypev=transformtype; %'projective';%'affine';
mincontrastvalue=0.001; % =0.05; default = 0.2
minqualityvalue=0.1	; % default = 0.1
maxdistancesep=1.5; % default = 1.5
outliertformerror_thresh=maxdistancesep;
minlinfitpoints=5;
confidencelevel=99.999;
searchmethod='Approximate';
maxratio=0.6;


dosavefigurematches=0;
dosavefigurematchesv=0;
dopause=0;
dodisplayfeatures=0;
doimagemagick=0;
doglobalbalance=1;
dovips=0;





pointsAreColinear = @(xy) rank(xy(2:end,:) - xy(1,:)) == 1;


% Output folder:
%imgoutfolder='~/Documents/Acfer-OUT/';
% if ~exist(imgoutfolder, 'dir')
%
%     [status, msg, msgID] = mkdir(imgoutfolder);
%
% end


% Get max and min primary beam current:
D=dir([imgfolder '*.hdr']);
nfiles=numel(D);
%emission_current_amps=zeros(nfiles,1);
[~,b,~]=fileparts(D(1).name);
imgext=b(strfind(b,'-')+1:end);

for ii=1:nfiles
    
    A_=readhdr([imgfolder D(ii).name]);
    Av=A_.value;
    emission_index=find(strcmp(A_.name,'EmissionCurrent'));
    %emission_current_amps(ii)=str2double(Av{emission_index});
    
end

% emission_current_amps_max=max(emission_current_amps);
% emission_current_amps_min=min(emission_current_amps);
% emission_current_amps_mean=mean(emission_current_amps);

samplename=fileparts(imgfolder);

D=dir([imgfolder '*.' imgext]);
imageinfo_ = imfinfo([imgfolder D(1).name]);
BitDepth=imageinfo_.BitDepth;
imsz=imageinfo_.Width;
if imsz~=imageinfo_.Height
    disp('NON-SQUARE IMAGES!!!')
    keyboard
end
allnames={D.name};
nfilestotal=numel(allnames);
col_t=zeros(nfilestotal,1);
row_t=zeros(nfilestotal,1);
for ii=1:nfilestotal
    [~,fname_,~]=fileparts(allnames{ii});
    st_=str2double(strsplit(fname_,'_'));
    col_t(ii)=st_(2);
    row_t(ii)=st_(1);
end

colt=unique(sort(col_t));
rowt=unique(sort(row_t));

%HistBinEdges=linspace(0,2^BitDepth-1,N_Histogram_Bins+1);

% OLD:
%cols_array=colt';

%NEW:
% colt_split=floor(numel(colt)/2);
% cols_array=[colt(colt_split+1:end)',fliplr(colt(1:colt_split)')];

%DEBUG:
%colt=[22:32]'; %[92:95]';

colt_split_index=0; %ceil(numel(colt)/2);
colt_split=colt(1); %colt(colt_split_index);
cols_array=[colt(colt_split_index+1:end)',fliplr(colt(1:colt_split_index)')];


ncols=numel(colt);

nrows=numel(rowt);

rowsforcol=cell(ncols,1);

for ii=1:ncols
    rowsforcol{ii}=row_t(col_t==colt(ii));
end

colsforrow=cell(nrows,1);

for ii=1:nrows
    colsforrow{ii}=col_t(row_t==rowt(ii));
end


converttext0=sprintf('convert -define registry:temporary-path=/backup/tmp');


converttexts=[];
converttext1=sprintf(' -background transparent -layers merge +repage strip%03d-%03d.png',min(cols_array),max(cols_array));

rowmax=0;
for qq=cols_array
    rowmax=max([rowmax;rowsforcol{colt==qq}]);
end

colmax=max(cols_array);

errortform=nan(colmax,rowmax,2);


matchedPointsPrev_v=cell(colmax,rowmax);
matchedPoints_v=cell(colmax,rowmax);
matchedPointsPrev=cell(colmax,rowmax);
matchedPoints=cell(colmax,rowmax);

matchedvcount=zeros(colmax,rowmax);
matchedcount=zeros(colmax,rowmax);
matchedvcount_refine=zeros(colmax,rowmax);
matchedcount_refine=zeros(colmax,rowmax);
matchmetric=cell(colmax,rowmax);
matchmetricv=cell(colmax,rowmax);
maxc_index_t=zeros(colmax,1);
maxvc_index_t=zeros(colmax,1);
imgshiftx=nan(colmax,rowmax);
imgshifty=nan(colmax,rowmax);
npt=nan(colmax,rowmax);

brightnessvalues=nan(colmax,rowmax,2);
brightnessvaluesv=nan(colmax,rowmax,2);
emission_current_amps=nan(colmax,rowmax);

delete(gcp('nocreate'));
poolobj = parpool;


for qq=cols_array
    %for qq=1:ncols
    tic
    
    %   col=sprintf('%03d',qq);
    %   disp(col);
    %   nf=nrows(colt==qq);
    
    
    % DEBUGGING:
    %rows_array=[50:51]; %[24:26]';
    rows_array=rowsforcol{colt==qq}'; %1:nf;
    
    
    % xlim=zeros(nf,2);
    % ylim=zeros(nf,2);
    % It=zeros(nf,imsz,imsz,'uint16');
    mincountpoints=2;
    
    
    % Tt=nan(nf,3,3);
    % Tt(1,:,:)=eye(3);
    
    I1=imread([imgfolder sprintf(['%03d_%03d.' imgext],rows_array(1),qq)]);
    
    width=size(I1,2);
    height=size(I1,1);
    searchfac=1;
    
    searchwidth=round(width*overlap*searchfac);
    searchheight=round(height*overlap*searchfac);
    
    It(1,:,:)=I1;
    
    
    inliercount=0;
    
    fprintf('\nCol %03d Rows %03d -- %03d',qq,rows_array(1),rows_array(end))
    
    
    %parfor kk=rows_array %2:nf
    if dodisplayfeatures
        fprintf('\nFeatures = ');
    end
    
    
    nrowsarray=numel(rows_array);
    
    matchedvcountqq=zeros(nrowsarray,1);
    matchedcountqq=zeros(nrowsarray,1);
    matchedPointsqq=cell(nrowsarray,1);
    matchedPointsPrevqq=cell(nrowsarray,1);
    matchedPoints_vqq=cell(nrowsarray,1);
    matchedPointsPrev_vqq=cell(nrowsarray,1);
    brightnessvaluesqq=nan(nrowsarray,2);
    brightnessvaluesvqq=nan(nrowsarray,2);
    emission_current_ampsqq=nan(nrowsarray,1);
    matchmetricqq=cell(nrowsarray,1);
    matchmetricvqq=cell(nrowsarray,1);
    
    parfor kk__=1:nrowsarray %2:nf
        %for kk__=1:nrowsarray %2:nf
        
        [matchedcountqq(kk__),matchedvcountqq(kk__),matchedPointsqq{kk__},matchedPointsPrevqq{kk__},matchedPoints_vqq{kk__},matchedPointsPrev_vqq{kk__},brightnessvaluesqq(kk__,:),brightnessvaluesvqq(kk__,:),emission_current_ampsqq(kk__),matchmetricqq{kk__},matchmetricvqq{kk__}] = ... % matchmetricqq(kk__)
            computeMatchedPointsMosaicMaker(kk__,rows_array,cols_array,qq,imgext,imgfolder,colt_split,height,searchheight,width,searchwidth,dosavefigurematches,dosavefigurematchesv,dodisplayfeatures,minqualityvalue,mincontrastvalue,transformtype,transformtypev,mincountpoints,pointsAreColinear,searchmethod,maxratio);
        
    end
    
    % close poolobj
    
    matchedcount(qq,rows_array)=matchedcountqq;
    matchedvcount(qq,rows_array)=matchedvcountqq;
    matchedPoints(qq,rows_array)=matchedPointsqq;
    matchedPointsPrev(qq,rows_array)=matchedPointsPrevqq;
    matchedPoints_v(qq,rows_array)=matchedPoints_vqq;
    matchedPointsPrev_v(qq,rows_array)=matchedPointsPrev_vqq;
    brightnessvalues(qq,rows_array,:)=brightnessvaluesqq;
    brightnessvaluesv(qq,rows_array,:)=brightnessvaluesvqq;
    emission_current_amps(qq,rows_array)=emission_current_ampsqq;
    matchmetric(qq,rows_array)=matchmetricqq;
    matchmetricv(qq,rows_array)=matchmetricvqq;
    
end

delete(gcp('nocreate'));


Bi=1./(emission_current_amps./nanmean(emission_current_amps(:)));
%save('brightness.mat')
%%
if doglobalbalance
    
    %[SLOPE,OFFSET]=globalbalance(Bi,brightnessvalues,brightnessvaluesv,cols_array,colt_split);
    [SLOPE,OFFSET]=globalbalance2(brightnessvalues,brightnessvaluesv,cols_array,colt_split);
    
else
    
    SLOPE=ones(colmax,rowmax);
    OFFSET=zeros(colmax,rowmax);
    
end




fprintf('\n\n')
%%

% nx=colmax;
% ny=rowmax;


close all

if strcmp(transformtype,'affine') || strcmp(transformtype,'similarity')
    tformt(colmax,rowmax)=affine2d(eye(3));
    %tfti(max(col_t),maxrows)=affine2d(eye(3));
end

if strcmp(transformtype,'projective')
    tformt(colmax,rowmax)=projective2d(eye(3));
    %tfti(max(col_t),maxrows)=projective2d(eye(3));
end

Tt=nan(colmax,rowmax,4,3,3);
Ttxy=nan(colmax,rowmax,3,3);
d1x=nan(colmax,rowmax);
d2x=nan(colmax,rowmax);
d3x=nan(colmax,rowmax);
d4x=nan(colmax,rowmax);
d1y=nan(colmax,rowmax);
d2y=nan(colmax,rowmax);
d3y=nan(colmax,rowmax);
d4y=nan(colmax,rowmax);

noutliers=inf;
itercount=0;


for qq=cols_array
    
    if qq<=colt_split
        
        qqp=qq+1;
        qqn=qq-1;
        
    else
        
        qqp=qq-1;
        qqn=qq+1;
        
    end
    
    [maxv,maxc_index_t(qq)]=max(matchedcount(qq,:));
    if maxv==0
        maxc_index_t(qq)=2;
    end
    
    
    % Tt(qq,maxc_index_t(qq)-1,:,:)=eye(3);
    %tformt(qq,maxc_index_t(qq)-1).T=eye(3);
    %loopvector=[maxc_index_t(qq):rowsforcol{colt==qq}(end),(maxc_index_t(qq)-2):-1:rowsforcol{colt==qq}(1)];
    
    if qq~=cols_array(1)
        
        [~,maxvc_index_t(qq)]=max(matchedvcount(qq,:).*~isnan(sum(sum(Tt(qqp,:,:,:),3),4)));
        %loopvector=[maxvc_index_t(qq):rowsforcol{colt==qq}(end),(maxvc_index_t(qq)-1):-1:rowsforcol{colt==qq}(1)];
        
    end
    
    for kk=rowsforcol{colt==qq}' %loopvector
        
        ll=[];
        mp_ll=[];
        
        llp=[];
        mp_llp=[];
        
        if strcmp(transformtype,'affine') || strcmp(transformtype,'similarity')
            tf_=affine2d(eye(3));
        end
        
        if strcmp(transformtype,'projective')
            tf_=projective2d(eye(3));
        end
        
        
        if strcmp(transformtypev,'affine') || strcmp(transformtypev,'similarity')
            tfv_=affine2d(eye(3));
        end
        
        if strcmp(transformtypev,'projective')
            tfv_=projective2d(eye(3));
        end
                
        if matchedcount(qq,kk)>=mincountpoints
 
            ll=matchedPointsPrev{qq,kk}.Location;
            mp_ll=matchedPoints{qq,kk}.Location;
            [~, mp_ll, ll, ~] = estimateGeometricTransform(...
                mp_ll, ll, transformtype, 'Confidence', confidencelevel, 'MaxNumTrials', 1e5, 'MaxDistance', maxdistancesep); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
            if numel(mp_ll)>=mincountpoints
                mpx=mean(mp_ll(:,1));
                mpy=mean(mp_ll(:,2));
                llx=mean(ll(:,1));
                lly=mean(ll(:,2));
                d2x(qq,kk)=-mpx+llx;
                d2y(qq,kk)=-mpy+lly;
            end
            %disp('kk-1');
        end
               
        if qqp>0 && qqp<=colmax && matchedvcount(qq,kk)>=mincountpoints
            % if qqp>0 && qqp<=colmax && matchedvcount(qq,kk)>=mincountpoints && ~isnan(sum(sum(Tt(qqp,kk,:,:)))) && ~isoutlierall(qq,kk)
            %tfv_.T=squeeze(Tt(qqp,kk,:,:));
            llp=matchedPointsPrev_v{qq,kk}.Location;
            mp_llp=matchedPoints_v{qq,kk}.Location;
            [~, mp_llp, llp, ~] = estimateGeometricTransform(...
                mp_llp, llp, transformtype, 'Confidence', confidencelevel, 'MaxNumTrials', 1e5, 'MaxDistance', maxdistancesep); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
            if numel(mp_ll)>=mincountpoints
                mpx=mean(mp_llp(:,1));
                mpy=mean(mp_llp(:,2));
                llx=mean(llp(:,1));
                lly=mean(llp(:,2));
                d3x(qq,kk)=-mpx+llx;
                d3y(qq,kk)=-mpy+lly;
            end
            %disp('qqp')
        end
        
        if numel(mp_ll) > 0
            
            [~, mp_ll, ll, ~] = estimateGeometricTransform(...
                mp_ll, ll, transformtype, 'Confidence', confidencelevel, 'MaxNumTrials', 1e5, 'MaxDistance', maxdistancesep); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
            
            mpx=mean(mp_ll(:,1));
            mpy=mean(mp_ll(:,2));
            llx=mean(ll(:,1));
            lly=mean(ll(:,2));
            
            [tformtqq, inlier1, inlier2, ~] = estimateGeometricTransform(...
                mp_ll-[mpx,mpy], ll-[llx,lly], transformtype, 'MaxNumTrials', 1e5); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
            
            [tformtqq2, ~, ~, ~] = estimateGeometricTransform(...
                ll-[llx,lly], mp_ll-[mpx,mpy], transformtype, 'MaxNumTrials', 1e5); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
                  
            Tt(qq,kk,1,:,:)=[1,0,0;0,1,0;-mpx+llx,-mpy+lly,1]*tformtqq.T;
            Tt(qq,kk,3,:,:)=[1,0,0;0,1,0;mpx-llx,mpy-lly,1]*tformtqq2.T;
            
            if numel(inlier1)>=mincountpoints
                U=transformPointsForward(tformtqq,inlier1);
                errortform(qq,kk,1)=sum((U(:,1)-inlier2(:,1)).^2+(U(:,2)-inlier2(:,2)).^2);
                matchedcount_refine(qq,kk)=size(U,1);
                %errortform(qq,kk,1)=rms(matchmetric{qq,kk})./(numel(matchmetric{qq,kk})).^1;
            end
            
        end
        
        if numel(mp_llp) > 0
            
            [~, mp_llp, llp, ~] = estimateGeometricTransform(...
                mp_llp, llp, transformtype, 'Confidence', confidencelevel, 'MaxNumTrials', 1e5, 'MaxDistance', maxdistancesep); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
            
            mpx=mean(mp_llp(:,1));
            mpy=mean(mp_llp(:,2));
            llx=mean(llp(:,1));
            lly=mean(llp(:,2));
            
            [tformtqq, inlier1, inlier2, ~] = estimateGeometricTransform(...
                mp_llp-[mpx,mpy], llp-[llx,lly], transformtype, 'MaxNumTrials', 1e5); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
            
            [tformtqq2, ~, ~, ~] = estimateGeometricTransform(...
                llp-[llx,lly], mp_llp-[mpx,mpy], transformtype, 'MaxNumTrials', 1e5); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
            
            
            Tt(qq,kk,2,:,:)=[1,0,0;0,1,0;-mpx+llx,-mpy+lly,1]*tformtqq.T;
            Tt(qq,kk,4,:,:)=[1,0,0;0,1,0;mpx-llx,mpy-lly,1]*tformtqq2.T;
            
            
            
            if numel(inlier1)>=mincountpoints
                U=transformPointsForward(tformtqq,inlier1);
                errortform(qq,kk,2)=sum((U(:,1)-inlier2(:,1)).^2+(U(:,2)-inlier2(:,2)).^2);
                matchedvcount_refine(qq,kk)=size(U,1);
                %errortform(qq,kk,2)=rms(matchmetricv{qq,kk})./(numel(matchmetricv{qq,kk})).^1;
            end
            
            
        end
        
    end
end

%%

nx=rowmax;
ny=colmax;

gg=zeros(nx*ny,nx*ny);

for ii=1:ny
    for jj=1:nx

         if ii>1
            gg(ii+(jj-1)*ny,(ii-1)+(jj-1)*ny)=1;
        end
        
        if ii<ny
            gg(ii+(jj-1)*ny,(ii+1)+(jj-1)*ny)=1;
        end
        
        if jj>1
            gg(ii+(jj-1)*ny,ii+(jj-2)*ny)=1;
        end
        
        if jj<nx
            gg(ii+(jj-1)*ny,ii+jj*ny)=1;
        end       
        
        
        
    end
    
end

G=graph(gg);

EndNodes=G.Edges.EndNodes;
G.Edges.Weight=G.Edges.Weight*inf;


%qq1=70;
%kk1=113;

%N=nearest(G,sub2ind([ny,nx],qq1,rowmax),1,'Method','unweighted'); 
%[a,b]=ind2sub([ny,nx],N);


qqo=round(median(cols_array));
kko=maxc_index_t(qqo);

pctile_thresh=99;
tform_power=1;
shift_power=1;

ET1=errortform(:,:,1)./matchedcount_refine;
ET1_p=prctile(ET1(:),pctile_thresh);
ET2=errortform(:,:,2)./matchedvcount_refine;
ET2_p=prctile(ET2(:),pctile_thresh);

ST1=abs(squeeze(Tt(:,:,1,3,1))./(imsz*(1-overlap))+1);
ST1_p=prctile(ST1(:),pctile_thresh);

ST2=abs(squeeze(Tt(:,:,2,3,2))./(imsz*(1-overlap))-1);
ST2_p=prctile(ST2(:),pctile_thresh);

isoutlier=zeros(colmax,rowmax);

for qq1=cols_array
    for kk1=rowsforcol{colt==qq1}'
        
        if kk1>1
            if (errortform(qq1,kk1,1)/matchedcount_refine(qq1,kk1))<ET1_p && abs(squeeze(Tt(qq1,kk1,1,3,1))./(imsz*(1-overlap))+1)<ST1_p
                in2a=find(EndNodes(:,1)==sub2ind([ny,nx],qq1,kk1-1) & EndNodes(:,2)==sub2ind([ny,nx],qq1,kk1));
                G.Edges.Weight(in2a)=((errortform(qq1,kk1,1)./matchedcount_refine(qq1,kk1))./ET1_p).^tform_power.*(abs(squeeze(Tt(qq1,kk1,1,3,1))./(imsz*(1-overlap))+1)./ST1_p).^shift_power; % errortform(qq1,kk1,1);
            else
                isoutlier(qq1,kk1)=1;
            end
        end
        if qq1>1
            if (errortform(qq1,kk1,2)/matchedvcount_refine(qq1,kk1))<ET2_p && abs(squeeze(Tt(qq1,kk1,2,3,2))./(imsz*(1-overlap))-1)<ST2_p
                in4a=find(EndNodes(:,1)==sub2ind([ny,nx],qq1-1,kk1) & EndNodes(:,2)==sub2ind([ny,nx],qq1,kk1));
                G.Edges.Weight(in4a)=((errortform(qq1,kk1,2)./matchedvcount_refine(qq1,kk1))./ET2_p).^tform_power.*(abs(squeeze(Tt(qq1,kk1,2,3,2))./(imsz*(1-overlap))-1)./ST2_p).^shift_power; % errortform(qq1,kk1,2);                
            else
                isoutlier(qq1,kk1)=1;
            end
        end
        
    end
end

Tt2=nan(colmax,rowmax,3,3);

for qq=cols_array
    for kk=rowsforcol{colt==qq}'
        
        [path1,d]=shortestpath(G,sub2ind([ny,nx],qq,kk),sub2ind([ny,nx],qqo,kko));
        %[path1,d]=shortestpath(G,sub2ind([ny,nx],kk,qq),sub2ind([ny,nx],kko,qqo));
        
        [x_,y_]=ind2sub([ny,nx],path1);
        %plot(x_,y_,'-s',qqo,kko,'rp',qq,kk,'kp'); axis image
        
        % [d,isoutlier(qq,kk)]
        
        if isfinite(d) && isoutlier(qq,kk)==0
            
            np=numel(path1);
            npt(qq,kk)=np;
            
            a=nan(np,1);
            b=nan(np,1);
            for ii=1:np
                [a(ii),b(ii)]=ind2sub([ny,nx], path1(np-ii+1));
                %[a(ii),b(ii)]=ind2sub([ny,nx], path1(np-ii+1));
                %disp([a(ii),b(ii)])
            end
            %disp('-------');
            diffa=diff(a);
            diffb=diff(b);
            
            T_=eye(3);
            for ii=2:np
                if diffa(ii-1)==0 && diffb(ii-1)==-1
                    T_=squeeze(Tt(a(ii),b(ii)+1,3,:,:))*T_;
                end
                if diffa(ii-1)==0 && diffb(ii-1)==1
                    T_=squeeze(Tt(a(ii),b(ii),1,:,:))*T_;
                end
                if diffa(ii-1)==-1 && diffb(ii-1)==0
                    %squeeze(Tt(a(ii)+1,b(ii),2,:,:))
                    T_=squeeze(Tt(a(ii)+1,b(ii),4,:,:))*T_;
                end
                if diffa(ii-1)==1 && diffb(ii-1)==0
                    T_=squeeze(Tt(a(ii),b(ii),2,:,:))*T_;
                end
            end
            
            Tt2(qq,kk,:,:)=T_;
            
        end
        
    end
end



T2=Tt2;
for ii=1:3
    for jj=1:2
        T2(:,:,ii,jj)=inpaint_nans(Tt2(:,:,ii,jj),0);
    end
end
T2(:,:,1,3)=0;
T2(:,:,2,3)=0;
T2(:,:,3,3)=double(single(inpaint_nans(Tt2(:,:,3,3),0)));


for qq=cols_array
    for kk=rowsforcol{colt==qq}'
        
        tformt(qq,kk).T=squeeze(T2(qq,kk,:,:));
        
    end
end


XLim=zeros(colmax,rowmax,2);
YLim=zeros(colmax,rowmax,2);
for qq=cols_array
    for kk=rowsforcol{colt==qq}'
        [XLim(qq,kk,:), YLim(qq,kk,:)] = outputLimits(tformt(qq,kk),  [1 width], [1 height]);
    end
end

% Find the minimum and maximum output limits
xMin = min([1; XLim(:)]);
xMax = max([width; XLim(:)]);

yMin = min([1; YLim(:)]);
yMax = max([height; YLim(:)]);

% Width and height of panorama.
widthp  = round(xMax - xMin);
heightp = round(yMax - yMin);

imageinfofilename=sprintf('imageinformation_%03d_%03d.txt',min(cols_array),max(cols_array));

fid=fopen(imageinfofilename,'w+');

fprintf(fid,'%d,%d,%d,%d\n',xMin,xMax,yMin,yMax);

%brightnessfactor=brightnessfactor*2;

converttexts=[];
for qq=cols_array
    for kk=rowsforcol{colt==qq}' %2:nf
        
        fname=[imgfolder sprintf(['%03d_%03d.' imgext],kk,qq)];
        
        [imgshiftx(qq,kk),imgshifty(qq,kk)]=transformPointsForward(tformt(qq,kk),0,0);
        
        converttexts=[converttexts sprintf(' \\( %s -evaluate multiply %.15f -alpha set -virtual-pixel transparent +distort AffineProjection ''%.15f,%.15f,%.15f,%.15f,%.15f,%.15f'' \\)',fname,SLOPE(qq,kk),tformt(qq,kk).T(1,1),tformt(qq,kk).T(1,2),tformt(qq,kk).T(2,1),tformt(qq,kk).T(2,2),tformt(qq,kk).T(3,1),tformt(qq,kk).T(3,2))]; %,sxv,rxv,ryv,syv,txv,tyv)];
        
        fprintf(fid,'%s,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f\n',fname,SLOPE(qq,kk),OFFSET(qq,kk),tformt(qq,kk).T(1,1),tformt(qq,kk).T(1,2),tformt(qq,kk).T(2,1),tformt(qq,kk).T(2,2),tformt(qq,kk).T(3,1),tformt(qq,kk).T(3,2));
    end
end

fclose(fid);


cc=distinguishable_colors(numel(cols_array));
cc2=distinguishable_colors(rowmax);
figure('visible','on');
set(gcf,'Position',[0,0,2000,2000])
hold on
for jj=1:rowmax
    plot(imgshiftx(:,jj),imgshifty(:,jj),'-','Color',cc2(jj,:));
end
for jj=cols_array
    plot(imgshiftx(jj,:),imgshifty(jj,:),'-','Color',cc(jj-min(cols_array)+1,:));
end
for ii=cols_array
    indout=isnan(sum(sum(Tt2(ii,rowsforcol{colt==ii},:,:),3),4)); % isoutlierall(ii,rowsforcol{colt==ii}) |
    plot(imgshiftx(ii,rowsforcol{colt==ii}(~indout)),imgshifty(ii,rowsforcol{colt==ii}(~indout)),'s','Color',cc(ii-min(cols_array)+1,:),'MarkerFaceColor',cc(ii-min(cols_array)+1,:));
    plot(imgshiftx(ii,rowsforcol{colt==ii}(indout)),imgshifty(ii,rowsforcol{colt==ii}(indout)),'s','Color',cc(ii-min(cols_array)+1,:),'MarkerFaceColor','w');
    if maxvc_index_t(ii)>0
        plot(imgshiftx(ii,maxvc_index_t(ii)),imgshifty(ii,maxvc_index_t(ii)),'O','Color',cc(ii-min(cols_array)+1,:),'MarkerFaceColor',cc(ii-min(cols_array)+1,:),'MarkerSize',7);
    end
    text(max(imgshiftx(ii,rowsforcol{colt==ii}))+width,mean(imgshifty(ii,rowsforcol{colt==ii})),sprintf('%03d',ii))
end
plot(imgshiftx(qqo,kko),imgshifty(qqo,kko),'p','Color','k','MarkerFaceColor','k','MarkerSize',11);
hold off
xlabel('Image X (pixels)')
ylabel('Image Y (pixels)')
set(gca, 'Ydir', 'reverse')
axis image
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
%imwrite(uint8(255*(brightnessratio./max(brightnessratio(:))).^1),sprintf('brightnessratio_%03d_%03d.png',min(cols_array),max(cols_array)));
%imwrite(uint8(255*(brightnessratiov./max(brightnessratiov(:))).^1),sprintf('brightnessratiov_%03d_%03d.png',min(cols_array),max(cols_array)));

[~,b]=min(abs(imgshifty(cols_array(1),:)));

origin=[cols_array(1),b];

[gridcoordsx,gridcoordsy]=meshgrid(linspace(0,width*(1-overlap)*rowmax,rowmax),linspace(0,width*(1-overlap)*colmax,colmax));
gridcoordsx=gridcoordsx(origin(1),origin(2))-gridcoordsx;
gridcoordsy=gridcoordsy-gridcoordsy(origin(1),origin(2));

grid_distance=sqrt((gridcoordsx-imgshiftx).^2 + (gridcoordsy-imgshifty).^2);

imagesc(grid_distance); axis image; colorbar
gridfname=sprintf('grid_error_%03d_%03d.png',min(cols_array),max(cols_array));
set(gca,'XDir','reverse');
print('-dpng','-r300',gridfname);
close(gcf);

imagesc(errortform(:,:,1)); axis image; colorbar
fname=sprintf('tform1_%03d_%03d.png',min(cols_array),max(cols_array));
set(gca,'XDir','reverse');
set(gca,'ColorScale','log')
print('-dpng','-r300',fname);
close(gcf);

imagesc(errortform(:,:,2)); axis image; colorbar
fname=sprintf('tform2_%03d_%03d.png',min(cols_array),max(cols_array));
axis image;
set(gca,'XDir','reverse');
set(gca,'ColorScale','log')
print('-dpng','-r300',fname);
close(gcf);

imagesc(abs(squeeze(Tt(:,:,2,3,2))./(imsz*(1-overlap))-1));
axis image;
set(gca,'XDir','reverse');
colorbar;
fname=sprintf('shifterry_%03d_%03d.png',min(cols_array),max(cols_array));
set(gca,'ColorScale','log')
print('-dpng','-r300',fname);
close(gcf);

imagesc(abs(squeeze(Tt(:,:,1,3,1))./(imsz*(1-overlap))+1)); axis image; colorbar
fname=sprintf('tform2_%03d_%03d.png',min(cols_array),max(cols_array));
set(gca,'XDir','reverse');
set(gca,'ColorScale','log')
print('-dpng','-r300',fname);
fname=sprintf('shifterrx_%03d_%03d.png',min(cols_array),max(cols_array));
print('-dpng','-r300',fname);
close(gcf);

if dopause
    pause
end

%imcall=['echo ' converttexts converttext1 '| xargs ' converttext0]
imcall=[converttext0 converttexts converttext1]; % num2str(heightp)

save(sprintf('mosaic_maker_%03d-%03d.mat',min(cols_array),max(cols_array)))

dzifilename=sprintf('deepzoom_%03d_%03d.dzi',min(cols_array),max(cols_array));
dzifilenametest='deepzoom.dzi';

pyvipscall=['python3 ~/Documents/affinetransform3.py ~/Documents/MATLAB/' imageinfofilename ' ' dzifilenametest];

if dovips
    system(pyvipscall);
else
    disp(pyvipscall);
end


if doimagemagick
    system(imcall);
    imresizecall=sprintf('/usr/local/bin/vipsthumbnail strip%03d-%03d.png --size %d',min(cols_array),max(cols_array),pixelcountresizewidth);
    system(imresizecall);
end

fid = fopen('imcall.sh','wt');
fprintf(fid, '%s\n\n%s','#!/bin/bash',imcall);
fclose(fid);

toc

%system(['python ~/Documents/affinetransform.py ' samplename]);

