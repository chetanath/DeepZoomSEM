tic

% rsync -avz ~/Data/Tescan/Ryan/MurchisonThickSection/HighRes/ ~/Documents/MurchisonThickSection
% ~/Documents/MoveTescanImages.sh ~/Documents/MurchisonThickSection

close all
clearvars;
%delete(gcp('nocreate'));

%cols_array=10:98;
%; %32:41;

%imgfolder='~/Data/Tescan/Ryan/Acfer182/HighRes/';
%imgfolder='~/Data/Tescan/Ryan/DOM143055/HighRes/';
imgfolder='~/Data/Tescan/Ryan/Acfer182/HighRes/';

pixelcountresize=25e6;
pixelcountresizewidth=10e3; % width in pixels
N_Histogram_Bins=64;
offset_tolerance=inf; %0.1;
tform_tolerance=inf; %0.04;
npolynomial=2;
mincountpoints=10;

transformtype='similarity'; %'affine';
transformtypev=transformtype; %'projective';%'affine';
mincontrastvalue=0.001; % =0.05; default = 0.2.....0.001
%minqualityvalue=0.1; % default = 0.1
maxdistancesep=1.5; % default = 1.5
outliertformerror_thresh=maxdistancesep;
minlinfitpoints=5;
searchmethod='Approximate'; %'Exhaustive'; %'Approximate';
maxratio=0.8;
numoctaves=8;
rszscale=2;
confidencelevel=99;
maxnumtrials=1e5;
maxdistance=1.5;

dooverlapfromhdr=true;
doglobalbalance=true;
dorefine=true;
dosavefigurematches=false;
dosavefigurematchesv=false;
dosavematfile=false;
dopause=false;
dodisplayfeatures=false;
doimagemagick=false;
dovips=false;





pointsAreColinear = @(xy) rank(xy(2:end,:) - xy(1,:)) == 1;


% Output folder:
%imgoutfolder='~/Documents/Acfer-OUT/';
% if ~exist(imgoutfolder, 'dir')
%
%     [status, msg, msgID] = mkdir(imgoutfolder);
%
% end

if dooverlapfromhdr
    % Get max and min primary beam current:
    D=dir([imgfolder '*.hdr']);
    nfiles=numel(D);
    nfilesread=8;
    %emission_current_amps=zeros(nfiles,1);
    [~,b,~]=fileparts(D(1).name);
    imgext=b(strfind(b,'-')+1:end);
    XStage=nan(nfiles,1);
    YStage=nan(nfiles,1);
    ZStage=nan(nfiles,1);
    WD=nan(nfiles,1);
    StageRotation=nan(nfiles,1);
    PixelSizeX=nan(nfiles,1);
    PixelSizeY=nan(nfiles,1);
    EmissionCurrent=nan(nfiles,1);
    middlefilenumber=round(nfiles/2);
    
    %for ii=1:nfilesread
    
    parfor ii=1:nfiles
        
        [XStage(ii),YStage(ii),ZStage(ii),WD(ii),StageRotation(ii),PixelSizeX(ii),PixelSizeY(ii),EmissionCurrent(ii)]=tescanreadcoords([imgfolder D(ii).name]);
        
    end
    
    
    % emission_current_amps_max=max(emission_current_amps);
    emission_current_amps_min=min(EmissionCurrent);
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
    fov_microns=round(imsz*PixelSizeX(1)*1e6);
    ss=diff(XStage);
    stage_shift_micron=ss(ss>0)*1e6;
    stage_shift_micron_mode=mode(round(stage_shift_micron));
    overlap=(fov_microns-stage_shift_micron_mode)/fov_microns;
else
    overlap=nan;
end

allnames={D.name};
col_t=zeros(nfiles,1);
row_t=zeros(nfiles,1);
for ii=1:nfiles
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




%DEBUGGING:
%colt=[90:110]'; %[92:95]';


% % colt_split_index=0; %ceil(numel(colt)/2);
% % colt_split=colt(1); %colt(colt_split_index);

cols_array=colt';

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

% DEBUGGING:
%onlyrows=99:116;
%rowsforcol=repmat({onlyrows'},numel(rowsforcol),1);


%converttext0=sprintf('convert -define registry:temporary-path=/backup/tmp');
converttext0=sprintf('/opt/local/bin/convert -define registry:temporary-path=/backup/tmp');


converttexts=[];
converttext1=sprintf(' -background transparent -layers merge +repage strip%03d-%03d.png',min(cols_array),max(cols_array));


rowmax=0;
rowmin=inf;
for qq=cols_array
    rowmax=max([rowmax;rowsforcol{colt==qq}]);
    rowmin=min([rowmin;rowsforcol{colt==qq}]);
end

colmax=max(cols_array);
colmin=min(cols_array);

errortform=nan(colmax,rowmax,4);
errortformT=nan(colmax,rowmax,2);


matchedPointsPrev_v=cell(colmax,rowmax);
matchedPoints_v=cell(colmax,rowmax);
matchedPointsPrev=cell(colmax,rowmax);
matchedPoints=cell(colmax,rowmax);
llp_t=cell(colmax,rowmax);
mp_llp_t=cell(colmax,rowmax);
ll_t=cell(colmax,rowmax);
mp_ll_t=cell(colmax,rowmax);

XL_h=nan(colmax,rowmax);
XL_l=nan(colmax,rowmax);
YL=nan(colmax,rowmax,2);
XL=nan(colmax,rowmax,2);

matchedvcount=zeros(colmax,rowmax);
matchedcount=zeros(colmax,rowmax);
matchedvcount_refine=zeros(colmax,rowmax,2);
matchedcount_refine=zeros(colmax,rowmax,2);
matchmetric=cell(colmax,rowmax);
matchmetricv=cell(colmax,rowmax);
maxc_index_t=zeros(colmax,1);
maxvc_index_t=zeros(colmax,1);
imgshiftx=nan(colmax,rowmax);
imgshifty=nan(colmax,rowmax);
npt=nan(colmax,rowmax);
tformstage=nan(colmax,rowmax,3,3);

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
    
    
    rows_array=rowsforcol{colt==qq}'; %1:nf;
    
    
    % xlim=zeros(nf,2);
    % ylim=zeros(nf,2);
    % It=zeros(nf,imsz,imsz,'uint16');
    
    
    
    % Tt=nan(nf,3,3);
    % Tt(1,:,:)=eye(3);
    
    I1=imread([imgfolder sprintf(['%03d_%03d.' imgext],rows_array(1),qq)]);
    
    %     width=size(I1,2);
    %     height=size(I1,1);
    searchfac=1;
    
    searchwidth=round(imsz*overlap*searchfac);
    searchheight=round(imsz*overlap*searchfac);
    
    %    It(1,:,:)=I1;
    
    
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
    tformstageqq=zeros(nrowsarray,3,3);
    
    parfor kk__=1:nrowsarray %2:nf
        %     for kk__=1:nrowsarray %2:nf
        
        [matchedcountqq(kk__),matchedvcountqq(kk__),matchedPointsqq{kk__},matchedPointsPrevqq{kk__},matchedPoints_vqq{kk__},matchedPointsPrev_vqq{kk__},brightnessvaluesqq(kk__,:),brightnessvaluesvqq(kk__,:),emission_current_ampsqq(kk__),matchmetricqq{kk__},matchmetricvqq{kk__}] = ... % matchmetricqq(kk__)
            computeMatchedPointsMosaicMaker(kk__,rows_array,cols_array,qq,imgext,imgfolder,imsz,  searchheight,  imsz,searchwidth,dosavefigurematches,dosavefigurematchesv,dodisplayfeatures,mincontrastvalue,transformtype,transformtypev,mincountpoints,pointsAreColinear,searchmethod,maxratio,numoctaves,rszscale,confidencelevel,maxnumtrials,maxdistance);
        %                                  (kk__,rows_array,cols_array,qq,imgext,imgfolder,height,searchheight,width,searchwidth,dosavefigurematches,dosavefigurematchesv,dodisplayfeatures,mincontrastvalue,transformtype,transformtypev,mincountpoints,pointsAreColinear,searchmethod,maxratio,numoctaves,rszscale,confidencelevel,maxnumtrials,maxdistance)
        [XStage,YStage,ZStage,WD,StageRotation,PixelSizeX,PixelSizeY,EmissionCurrent]=tescanreadcoords(sprintf([imgfolder '%03d_%03d-' imgext '.hdr'],rows_array(kk__),qq));
        tformstageqq(kk__,:,:)=[1,0,0;0,1,0;-XStage/PixelSizeX,YStage./PixelSizeY,1];
        
        
    end
    
    tformstage(qq,rows_array,:,:)=tformstageqq;
    
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
    OFFSET=globalbalance2(brightnessvalues,brightnessvaluesv,cols_array);
    figure('visible','off');
    imagesc(fliplr(OFFSET(colmin:colmax,rowmin:rowmax))); colorbar; axis image;
    print('-dpng','-r300',sprintf('OFFSET_%03d_%03d.png',min(cols_array),max(cols_array)));
    
else
    
    OFFSET=zeros(colmax,rowmax);
    
end




fprintf('\n\n')
%%

% nx=colmax;
% ny=rowmax;


close all

% if strcmp(transformtype,'affine') || strcmp(transformtype,'similarity')
%     tformt(colmax,rowmax)=affine2d(eye(3));
%     %tfti(max(col_t),maxrows)=affine2d(eye(3));
% end

% if strcmp(transformtype,'projective')
tformt(colmax,rowmax)=projective2d(eye(3));
%tfti(max(col_t),maxrows)=projective2d(eye(3));
% end

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
    
    %   if qq<=colt_split
    
    qqp=qq-1;
    qqn=qq+1;
    
    %     else
    %
    %         qqp=qq-1;
    %         qqn=qq+1;
    %
    %     end
    
    [maxv,maxc_index_t(qq)]=max(matchedcount(qq,:));
    if maxv==0
        maxc_index_t(qq)=2;
    end
    
    
    % Tt(qq,maxc_index_t(qq)-1,:,:)=eye(3);
    %tformt(qq,maxc_index_t(qq)-1).T=eye(3);
    %loopvector=[maxc_index_t(qq):rowsforcol{colt==qq}(end),(maxc_index_t(qq)-2):-1:rowsforcol{colt==qq}(1)];
    
    if qq~=cols_array(1)
        
        [~,maxvc_index_t(qq)]=max(matchedvcount(qq,:,1).*~isnan(sum(sum(Tt(qqp,:,:,:),3),4)));
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
            %             [~, mp_ll, ll, ~] = estimateGeometricTransform(...
            %                 mp_ll, ll, transformtype, 'Confidence', confidencelevel, 'MaxNumTrials', 1e5, 'MaxDistance', maxdistancesep); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
            %             if numel(mp_ll)>=mincountpoints
            %                 mpx=mean(mp_ll(:,1));
            %                 mpy=mean(mp_ll(:,2));
            %                 llx=mean(ll(:,1));
            %                 lly=mean(ll(:,2));
            %                 d2x(qq,kk)=-mpx+llx;
            %                 d2y(qq,kk)=-mpy+lly;
            %             end
            %disp('kk-1');
            
            [tformtqq, inlier1, inlier2, ~] = estimateGeometricTransform(...
                mp_ll, ll, transformtype, 'MaxNumTrials', 1e5); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
            
            mp_ll_t{qq,kk}=inlier1;
            ll_t{qq,kk}=inlier2;
            
            [tformtqq2, inlier1a, inlier2a, ~] = estimateGeometricTransform(...
                ll, mp_ll, transformtype, 'MaxNumTrials', 1e5); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
            
            Tt(qq,kk,1,:,:)=tformtqq.T;
            
            %if kk>1
            Tt(qq,kk,3,:,:)=tformtqq2.T;
            %end
            
            if numel(inlier1)>=mincountpoints
                U=transformPointsForward(tformtqq,inlier1);
                errortform(qq,kk,1)=sum((U(:,1)-inlier2(:,1)).^2+(U(:,2)-inlier2(:,2)).^2);
                matchedcount_refine(qq,kk,1)=size(U,1);
                %errortform(qq,kk,1)=rms(matchmetric{qq,kk})./(numel(matchmetric{qq,kk})).^1;
            end
            if numel(inlier1a)>=mincountpoints
                U=transformPointsForward(tformtqq2,inlier1a);
                errortform(qq,kk,3)=sum((U(:,1)-inlier2a(:,1)).^2+(U(:,2)-inlier2a(:,2)).^2);
                matchedcount_refine(qq,kk,2)=size(U,1);
                %errortform(qq,kk,1)=rms(matchmetric{qq,kk})./(numel(matchmetric{qq,kk})).^1;
            end
            
            
        end
        
        if  matchedvcount(qq,kk,1)>=mincountpoints % qqp>0 && qqp<=colmax &&
            % if qqp>0 && qqp<=colmax && matchedvcount(qq,kk)>=mincountpoints && ~isnan(sum(sum(Tt(qqp,kk,:,:)))) && ~isoutlierall(qq,kk)
            %tfv_.T=squeeze(Tt(qqp,kk,:,:));
            llp=matchedPointsPrev_v{qq,kk}.Location;
            mp_llp=matchedPoints_v{qq,kk}.Location;
            %             [~, mp_llp, llp, ~] = estimateGeometricTransform(...
            %                 mp_llp, llp, transformtype, 'Confidence', confidencelevel, 'MaxNumTrials', 1e5, 'MaxDistance', maxdistancesep); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
            %             if numel(mp_ll)>=mincountpoints
            %                 mpx=mean(mp_llp(:,1));
            %                 mpy=mean(mp_llp(:,2));
            %                 llx=mean(llp(:,1));
            %                 lly=mean(llp(:,2));
            %                 d3x(qq,kk)=-mpx+llx;
            %                 d3y(qq,kk)=-mpy+lly;
            %             end
            %disp('qqp')
            
            
            %        if numel(mp_ll) > 0
            
            %             [~, mp_ll, ll, ~] = estimateGeometricTransform(...
            %                 mp_ll, ll, transformtype, 'Confidence', confidencelevel, 'MaxNumTrials', 1e5, 'MaxDistance', maxdistancesep); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
            %
            %             mpx=mean(mp_ll(:,1));
            %             mpy=mean(mp_ll(:,2));
            %             llx=mean(ll(:,1));
            %             lly=mean(ll(:,2));
            %
            %             [tformtqq, inlier1, inlier2, ~] = estimateGeometricTransform(...
            %                 mp_ll-[mpx,mpy], ll-[llx,lly], transformtype, 'MaxNumTrials', 1e5); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
            %
            %             [tformtqq2, ~, ~, ~] = estimateGeometricTransform(...
            %                 ll-[llx,lly], mp_ll-[mpx,mpy], transformtype, 'MaxNumTrials', 1e5); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
            %
            %             Tt(qq,kk,1,:,:)=[1,0,0;0,1,0;-mpx+llx,-mpy+lly,1]*tformtqq.T;
            %             Tt(qq,kk,3,:,:)=[1,0,0;0,1,0;mpx-llx,mpy-lly,1]*tformtqq2.T;
            
            [tformtqq, inlier1, inlier2, ~] = estimateGeometricTransform(...
                mp_llp, llp, transformtype, 'MaxNumTrials', 1e5); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
            
            mp_llp_t{qq,kk}=inlier1;
            llp_t{qq,kk}=inlier2;
            
            [tformtqq2, inlier1a, inlier2a, ~] = estimateGeometricTransform(...
                llp, mp_llp, transformtype, 'MaxNumTrials', 1e5); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
            
            Tt(qq,kk,2,:,:)=tformtqq.T;
            
            %if qq>1
            Tt(qq,kk,4,:,:)=tformtqq2.T;
            %end
            
            if numel(inlier1)>=mincountpoints
                U=transformPointsForward(tformtqq,inlier1);
                errortform(qq,kk,2)=sum((U(:,1)-inlier2(:,1)).^2+(U(:,2)-inlier2(:,2)).^2);
                matchedvcount_refine(qq,kk,1)=size(U,1);
                %errortform(qq,kk,2)=rms(matchmetricv{qq,kk})./(numel(matchmetricv{qq,kk})).^1;
            end
            if numel(inlier1a)>=mincountpoints
                U=transformPointsForward(tformtqq2,inlier1a);
                errortform(qq,kk,4)=sum((U(:,1)-inlier2a(:,1)).^2+(U(:,2)-inlier2a(:,2)).^2);
                matchedvcount_refine(qq,kk,2)=size(U,1);
                %errortform(qq,kk,2)=rms(matchmetricv{qq,kk})./(numel(matchmetricv{qq,kk})).^1;
            end
        end
        
        
        
        %       end
        
        %         if numel(mp_llp) > 0
        %
        %             [~, mp_llp, llp, ~] = estimateGeometricTransform(...
        %                 mp_llp, llp, transformtype, 'Confidence', confidencelevel, 'MaxNumTrials', 1e5, 'MaxDistance', maxdistancesep); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
        
        %             mpx=mean(mp_llp(:,1));
        %             mpy=mean(mp_llp(:,2));
        %             llx=mean(llp(:,1));
        %             lly=mean(llp(:,2));
        %
        %             [tformtqq, inlier1, inlier2, ~] = estimateGeometricTransform(...
        %                 mp_llp-[mpx,mpy], llp-[llx,lly], transformtype, 'MaxNumTrials', 1e5); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
        %
        %             [tformtqq2, ~, ~, ~] = estimateGeometricTransform(...
        %                 llp-[llx,lly], mp_llp-[mpx,mpy], transformtype, 'MaxNumTrials', 1e5); %, 'MaxDistance', maxdistancesep);%;); %, 'MaxDistance', round(width/100));
        %
        %
        %             Tt(qq,kk,2,:,:)=[1,0,0;0,1,0;-mpx+llx,-mpy+lly,1]*tformtqq.T;
        %             Tt(qq,kk,4,:,:)=[1,0,0;0,1,0;mpx-llx,mpy-lly,1]*tformtqq2.T;
        
        
        
        
        
        
    end
end


%%
%

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

%%
G=digraph(gg);

EndNodes=G.Edges.EndNodes;
G.Edges.Weight=G.Edges.Weight*inf;
G1=G;

%qq1=70;
%kk1=113;

%N=nearest(G,sub2ind([ny,nx],qq1,rowmax),1,'Method','unweighted');
%[a,b]=ind2sub([ny,nx],N);






kko=inf;
rowmeanv=round(mean(cell2mat(rowsforcol)));
for ii=-round(numel(cols_array)/20):round(numel(cols_array)/20)
    qqo_=round(median(cols_array)+ii);
    
    if abs(maxc_index_t(qqo_)-rowmeanv)<abs(kko-rowmeanv)
        kko=maxc_index_t(qqo_);
        qqo=qqo_;
    end
    
end


pctile_thresh_tform=100;
pctile_thresh_shift=100;
tform_power=1;
shift_power=1;

ET1=sqrt(errortform(:,:,1))./matchedcount;
ET1_p=prctile(ET1(:),pctile_thresh_tform);

ET2=sqrt(errortform(:,:,2))./matchedvcount;
ET2_p=prctile(ET2(:),pctile_thresh_tform);

ST1a=abs(squeeze(Tt(:,:,1,3,1))+(imsz*(1-overlap)));
ST1b=abs(squeeze(Tt(:,:,3,3,1))-(imsz*(1-overlap)));
ST1a_p=prctile(ST1a(:),pctile_thresh_shift);
ST1b_p=prctile(ST1b(:),pctile_thresh_shift);

ST2a=abs(squeeze(Tt(:,:,2,3,2))-(imsz*(1-overlap)));
ST2b=abs(squeeze(Tt(:,:,4,3,2))+(imsz*(1-overlap)));
ST2a_p=prctile(ST2a(:),pctile_thresh_shift);
ST2b_p=prctile(ST2b(:),pctile_thresh_shift);

ST_p=prctile([ST2a(:);ST2b(:);ST1a(:);ST1b(:)],pctile_thresh_shift);


isoutlier=zeros(colmax,rowmax);
%ff2a=[];
for qq1=cols_array
    for kk1=rowsforcol{colt==qq1}'
        
        if (ET1(qq1,kk1)>ET1_p || ST1a(qq1,kk1)>ST_p || ST1b(qq1,kk1)>ST_p) && isfinite(errortform(qq1,kk1,1))
            isoutlier(qq1,kk1)=1;
        end
        if ET1(qq1,kk1)<ET1_p && ST1a(qq1,kk1)<ST_p %ST1a_p %&& isfinite(Tt(qq1,kk1,3,1,1)) %kk1>1 &&
            in2a=EndNodes(:,1)==sub2ind([ny,nx],qq1,kk1-1) & EndNodes(:,2)==sub2ind([ny,nx],qq1,kk1);
            % sprintf('%d,%d : %d,%d',qq1,kk1-1,qq1,kk1)
            %ff2a=[ff2a,find(in2a)];
            G.Edges.Weight(in2a)=(ET1(qq1,kk1)./ET1_p).^tform_power+(ST1a(qq1,kk1)./ST_p).^shift_power; %(ET1(qq1,kk1)./ET1_p).^tform_power.*(ST1(qq1,kk1)./ST1_p).^shift_power; % errortform(qq1,kk1,1);
            G1.Edges.Weight(in2a)=1;
        end
        if ET1(qq1,kk1)<ET1_p && ST1b(qq1,kk1)<ST_p %ST1b_p
            in2a=EndNodes(:,1)==sub2ind([ny,nx],qq1,kk1) & EndNodes(:,2)==sub2ind([ny,nx],qq1,kk1-1);
            % sprintf('%d,%d : %d,%d',qq1,kk1,qq1,kk1-1)
            %ff2a=[ff2a,find(in2a)];
            G.Edges.Weight(in2a)=(ET1(qq1,kk1)./ET1_p).^tform_power+(ST1b(qq1,kk1)./ST_p).^shift_power; %(ET1(qq1,kk1)./ET1_p).^tform_power.*(ST1(qq1,kk1)./ST1_p).^shift_power; % errortform(qq1,kk1,1);
            G1.Edges.Weight(in2a)=1;
        end
        
        if (ET2(qq1,kk1)>ET2_p  || ST2a(qq1,kk1)>ST_p || ST2b(qq1,kk1)>ST_p) && isfinite(errortform(qq1,kk1,2))
            isoutlier(qq1,kk1)=1;
        end
        if ET2(qq1,kk1)<ET2_p && ST2a(qq1,kk1)<ST_p %ST2a_p %&& isfinite(Tt(qq1,kk1,4,1,1)) % qq1>1 &&
            in4a=EndNodes(:,1)==sub2ind([ny,nx],qq1-1,kk1) & EndNodes(:,2)==sub2ind([ny,nx],qq1,kk1);
            %ff2a=[ff2a,find(in4a)];
            %sprintf('%d,%d : %d,%d',qq1-1,kk1,qq1,kk1)
            G.Edges.Weight(in4a)=(ET2(qq1,kk1)./ET2_p).^tform_power+(ST2a(qq1,kk1)./ST_p).^shift_power; %((errortform(qq1,kk1,2)./matchedvcount_refine(qq1,kk1))./ET2_p).^tform_power.*(ST2(qq1,kk1)./ST2_p).^shift_power; % errortform(qq1,kk1,2);
            G1.Edges.Weight(in4a)=1;
        end
        if ET2(qq1,kk1)<ET2_p && ST2b(qq1,kk1)<ST_p %ST2b_p
            in4a=EndNodes(:,1)==sub2ind([ny,nx],qq1,kk1) & EndNodes(:,2)==sub2ind([ny,nx],qq1-1,kk1);
            %ff2a=[ff2a,find(in4a)];
            %sprintf('%d,%d : %d,%d',qq1,kk1,qq1-1,kk1)
            G.Edges.Weight(in4a)=(ET2(qq1,kk1)./ET2_p).^tform_power+(ST2b(qq1,kk1)./ST_p).^shift_power; %((errortform(qq1,kk1,2)./matchedvcount_refine(qq1,kk1))./ET2_p).^tform_power.*(ST2(qq1,kk1)./ST2_p).^shift_power; % errortform(qq1,kk1,2);
            G1.Edges.Weight(in4a)=1;
        end
        
    end
end

% [path1,d,edgepath]=shortestpath(G,sub2ind([ny,nx],111,56),sub2ind([ny,nx],114,58));  G.Edges.Weight(edgepath)

% FULL GRAPH:
[Tt2,isoutlier,dtotal] = calculateMosaicPositions(colmax,rowmax,cols_array,rowsforcol,G,ny,nx,qqo,kko,isoutlier,Tt,false);

% JUST ROWS:
% Tt2=nan(colmax,rowmax,3,3);
% dtotal=nan(colmax,rowmax);
%
% for qq=cols_array
%
%     [Tt2_temp,isoutlier_temp,dtotal_temp] = calculateMosaicPositions(colmax,rowmax,qq,rowsforcol,G,ny,nx,qq,maxc_index_t(qq),isoutlier,colt,Tt,true);
%
%     Tt2(qq,rowsforcol{colt==qq}',:,:)=Tt2_temp(qq,rowsforcol{colt==qq}',:,:);
%
%     isoutlier(qq,rowsforcol{colt==qq}',:,:)=isoutlier_temp(qq,rowsforcol{colt==qq}',:,:);
%
%     dtotal(qq,rowsforcol{colt==qq}',:,:)=dtotal_temp(qq,rowsforcol{colt==qq}',:,:);
%
%
% end






if dopause
    pause
end

%imcall=['echo ' converttexts converttext1 '| xargs ' converttext0]
imcall=[converttext0 converttexts converttext1]; % num2str(heightp)

%save(sprintf('mosaic_maker_%03d-%03d.mat',min(cols_array),max(cols_array)))
save(sprintf('mosaic_maker_%03d-%03d.mat',min(cols_array),max(cols_array)), '-regexp', '^(?!(gg)$).')

dzifilename=sprintf('deepzoom_%03d_%03d.dzi',min(cols_array),max(cols_array));
dzifilenametest='deepzoom.dzi';

% pyvipscall=['python3 ~/Documents/affinetransform4.py ~/Documents/MATLAB/' imageinfofilename ' ' dzifilenametest];
%
% if dovips
%     system(pyvipscall);
% else
%     disp(pyvipscall);
% end


if doimagemagick
    system(imcall);
    %imresizecall=sprintf('/usr/local/bin/vipsthumbnail strip%03d-%03d.png --size %d',min(cols_array),max(cols_array),pixelcountresizewidth);
    %system(imresizecall);
end

fid = fopen('imcall.sh','wt');
fprintf(fid, '%s\n\n%s','#!/bin/bash',imcall);
fclose(fid);

%toc

%system(['python ~/Documents/affinetransform.py ' samplename]);
%[Tt2, isoutlier, ~] = calculateMosaicPositions(colmax,rowmax,cols_array,rowsforcol{colt==qq1}',G,ny,nx,qqo,kko,isoutlier,Tt,false);
if dosavematfile
    save(sprintf('mosaic_maker_%03d-%03d.mat',min(cols_array),max(cols_array)))
end
% STRIPS:

% load(sprintf('mosaic_maker_%03d-%03d.mat',1,133)); %min(cols_array),max(cols_array)))
%%
% load(sprintf('mosaic_maker_%03d-%03d.mat',min(cols_array),max(cols_array))); %min(cols_array),max(cols_array)))
%load(sprintf('mosaic_maker_%03d-%03d.mat',90,110));
%load('mosaic_maker_100-105.mat');

fido=fopen('pyvipsstrips.sh','w');

imageinfofilename=sprintf('imageinformation_strips_%03d_%03d.txt',cols_array(1),cols_array(end));

fid=fopen(imageinfofilename,'w');

nrefineloops=4*2; %4*8;
etfT=zeros(nrefineloops,1);
nrefinepoints=zeros(nrefineloops,1);
% ROW ANALYSES:
etfTp=inf;
if dorefine
    for KK=1:nrefineloops
        doprintrowpoints=true;
        NN=0;
        ttrowT=nan(numel(cols_array),3,3);
        cols_array_=setdiff(cols_array,1);
        if mod(KK,4)==3 || mod(KK,4)==0
            cols_array_=fliplr(cols_array_);
        end
        
        for qq_=cols_array_
            mp_llpx=[];
            mp_llpy=[];
            llpx=[];
            llpy=[];
            rows_array=rowsforcol{colt==qq_}';
            if mod(KK,4)==2 || mod(KK,4)==0
                rows_array=fliplr(rows_array);
            end
            
            for kk_=rows_array
                %if (kk_-1)>=rowmin && (qq_-1)>=colmin && (kk_+1)<=rowmax && (qq_+1)<=colmax
                if isfinite(Tt2(qq_,kk_,1,1)) %&& isfinite(Tt2(qq_-1,kk_,1,1)) && isfinite(Tt2(qq_,kk_-1,1,1))  && isfinite(Tt2(qq_+1,kk_,1,1)) && isfinite(Tt2(qq_,kk_+1,1,1))
                    
                    etfp = errortformfun(qq_,kk_,Tt2,matchedPointsPrev_v,matchedPoints_v,matchedPointsPrev,matchedPoints,colmin,colmax,rowmin,rowmax);
                    
                    t00=projective2d;
                    t00.T=squeeze(Tt2(qq_,kk_,:,:));
                    if (qq_-1)>=colmin
                        if isfinite(Tt2(qq_-1,kk_,1,1))
                            t01=projective2d;
                            t01.T=squeeze(Tt2(qq_-1,kk_,:,:));
                            [llpx,llpy]        =transformPointsForward(t01,matchedPointsPrev_v{qq_,kk_}.Location(:,1),matchedPointsPrev_v{qq_,kk_}.Location(:,2));
                            [mp_llpx,mp_llpy] = transformPointsForward(t00,matchedPoints_v{qq_,kk_}.Location(:,1),matchedPoints_v{qq_,kk_}.Location(:,2));
                            errortformTp1=errortformT(qq_,kk_,1);
                            errortformT(qq_,kk_,1)=sum((llpx-mp_llpx).^2+(llpy-mp_llpy).^2);
                        else
                            llpx=[];
                            llpy=[];
                            mp_llpx=[];
                            mp_llpy=[];
                        end
                    end
                    
                    if (kk_-1)>=rowmin
                        if isfinite(Tt2(qq_,kk_-1,1,1))
                            t02=projective2d;
                            t02.T=squeeze(Tt2(qq_,kk_-1,:,:));
                            [llpx2,llpy2]        =transformPointsForward(t02,matchedPointsPrev{qq_,kk_}.Location(:,1),matchedPointsPrev{qq_,kk_}.Location(:,2));
                            [mp_llpx2,mp_llpy2]=transformPointsForward(t00,matchedPoints{qq_,kk_}.Location(:,1),matchedPoints{qq_,kk_}.Location(:,2));
                            errortformTp2=errortformT(qq_,kk_,2);
                            errortformT(qq_,kk_,2)=sum((llpx2-mp_llpx2).^2+(llpy2-mp_llpy2).^2);
                        else
                            llpx2=[];
                            llpy2=[];
                            mp_llpx2=[];
                            mp_llpy2=[];
                        end
                    end
                    
                    if (numel(llpx)+numel(llpx2))>mincountpoints
                        %                         if numel(llpx)>mincountpoints && numel(llpx2)>mincountpoints
                        %                             ttrf =estimateGeometricTransform([ [mp_llpx;mp_llpx2],[mp_llpy;mp_llpy2] ],[ [llpx;llpx2] ,[llpy;llpy2] ], 'affine', 'MaxNumTrials',1e5);
                        %                         else
                        [ttrf,iP1,iP2] =estimateGeometricTransform([ [mp_llpx;mp_llpx2],[mp_llpy;mp_llpy2] ],[ [llpx;llpx2] ,[llpy;llpy2] ], 'affine', 'MaxNumTrials',1e5);
                        %                         end
                        
                        
                        %[pfx,pfy]=transformPointsForward(ttrf,[mp_llpx;mp_llpx2],[mp_llpy;mp_llpy2]);
                        %plot([llpx;llpx2],[llpy;llpy2],'b.',[mp_llpx;mp_llpx2],[mp_llpy;mp_llpy2],'r.',pfx,pfy,'bo')
                        
                        %ttrf.T
                        Tt2qqkkold=squeeze(Tt2(qq_,kk_,:,:));
                        Tt2(qq_,kk_,:,:)=squeeze(Tt2(qq_,kk_,:,:))*ttrf.T;
                        %                         [AXX,AYY] = transformPointsForward(ttrf,matchedPoints_v{qq_,kk_}.Location(:,1),matchedPoints_v{qq_,kk_}.Location(:,2));
                        %                         isf=AXX>0 & AYY>0;
                        %                         matchedPoints_v{qq_,kk_}(~isf)=[];
                        %                         matchedPointsPrev_v{qq_,kk_}(~isf)=[];
                        %                         matchedPoints_v{qq_,kk_}.Location=[AXX(isf),AYY(isf)];
                        %                         [AXX,AYY]=transformPointsForward(ttrf,matchedPoints{qq_,kk_}.Location(:,1),matchedPoints{qq_,kk_}.Location(:,2));
                        %                         isf=AXX>0 & AYY>0;
                        %                         matchedPoints{qq_,kk_}(~isf)=[];
                        %                         matchedPointsPrev{qq_,kk_}(~isf)=[];
                        %                         matchedPoints{qq_,kk_}.Location=[AXX(isf),AYY(isf)];
                    end
                    % t22=projective2d;
                    % t22.T=squeeze(Tt(qq_,kk_,2,:,:));
                    % [llpx_,llpy_]=transformPointsForward(t22,llpx_,llpy_);
                    %                 llpx=[llpx;llpx_];
                    %                 llpy=[llpy;llpy_];
                    %                 mp_llpx=[mp_llpx;mp_llpx_];
                    %                 mp_llpy=[mp_llpy;mp_llpy_];
                    
                    %end
                    
                    
                    etf = errortformfun(qq_,kk_,Tt2,matchedPointsPrev_v,matchedPoints_v,matchedPointsPrev,matchedPoints,colmin,colmax,rowmin,rowmax);
                    if etf>etfp || isnan(etf) || isnan(etfp)
                        
                        Tt2(qq_,kk_,:,:)=Tt2qqkkold;
                        if (qq_-1)>=colmin
                            if isfinite(Tt2(qq_-1,kk_,1,1))
                                errortformT(qq_,kk_,1)=errortformTp1;
                            end
                        end
                        
                        if (kk_-1)>=rowmin
                            if isfinite(Tt2(qq_,kk_-1,1,1))
                                errortformT(qq_,kk_,2)=errortformTp2;
                            end
                        end
                        
                    else
                        
                        NN=NN+1;
                        
                    end
                    
                end
                
                %   end
                
            end
        end
        
        %         if KK==1
        %             keyboard
        %         end
        
        etfT(KK)=nansum(nansum(errortformT(:,:,2)))+nansum(nansum(errortformT(:,:,1)));
        fprintf('Iter = %d, Error = %.5e, Number Changed=%d',KK,etfT(KK),NN);
        if KK>1
            fracerror=(etfT(KK)-etfT(KK-1))./etfT(KK);
            fprintf(', Frac Error = %.5e\n',fracerror);
        else
            fprintf('\n');
        end
        nrefinepoints(KK)=NN;
        
        
    end
    
end

T2=Tt2;
for ii=1:3
    for jj=1:2
        if sum(sum(isnan(Tt2(colmin:colmax,rowmin:rowmax,ii,jj))))>0
            T2(colmin:colmax,rowmin:rowmax,ii,jj)=inpaint_nans(Tt2(colmin:colmax,rowmin:rowmax,ii,jj),0);
        end
    end
end

T2(colmin:colmax,rowmin:rowmax,1,3)=0;
T2(colmin:colmax,rowmin:rowmax,2,3)=0;
if sum(sum(isnan(Tt2(colmin:colmax,rowmin:rowmax,3,3))))>0
    T2(colmin:colmax,rowmin:rowmax,3,3)=double(single(inpaint_nans(Tt2(colmin:colmax,rowmin:rowmax,3,3),0)));
end


%         if numel(llpx)>mincountpoints
%
%             [ttrow, ~, ~, ~] =estimateGeometricTransform([mp_llpx,mp_llpy],[llpx,llpy], 'affine', 'MaxNumTrials',1e5);
%             for kk__=rows_array
%                 T2(qq_,kk__,:,:)=squeeze(T2(qq_,kk__,:,:))*ttrow.T;
%             end
%
%         else
%             ttrow=projective2d;
%
%         end
%
% %         yy=T2(qq_,rows_array,1,1);
% %         isf=~isfinite(yy);
% %         fisf=find(isf);
% %         for mm=1:sum(isf)
% %             T2(qq_,rows_array(fisf(mm)),:,:)=squeeze(T2(qq_,rows_array(fisf(mm)),:,:))*ttrow.T;
% %         end
%
%         %disp('TTROW:')
%         %disp(ttrow.T)
%         ttrowT(jjj,:,:)=ttrow.T;
%         %tformstage_(qq_,kk_,:,:)=ttrow.T*squeeze(tformstage_(qq_,kk_,:,:));
%         %ttrefine2 = images.geotrans.PolynomialTransformation2D([[mp_llpx,mp_llpy];[mp_llx,mp_lly]],[[llpx,llpy];[llx,lly]],3);
%         %ttrefine2i = images.geotrans.PolynomialTransformation2D([[llpx,llpy];[llx,lly]],[[mp_llpx,mp_llpy];[mp_llx,mp_lly]],3);
%         %II=imread(sprintf('/Users/ogliore/Data/deepzoom/strip%03d.png',qq_));
%         %outputView = imref2d(size(II));
%         %IIw  = imwarp(II,ttrow,'OutputView',outputView);
%
%         %(qq_,kk_)
%
%         %Tt2(qq_,kk_,:,:) = squeeze(Tt2(qq_,kk_,:,:))*ttrefine2.T;
%
%         if doprintrowpoints
%             f=figure('visible','off');
%             % [pfx,pfy]=transformPointsInverse(ttrefine2i,[mp_llpx;mp_llx],[mp_llpy;mp_lly]);
%             [pfx,pfy]=transformPointsForward(ttrow,mp_llpx,mp_llpy);
%             plot(llpx,llpy,'b.',mp_llpx,mp_llpy,'r.',pfx,pfy,'bo')
%             tformfname=sprintf('refine/rowpoints_%03d.png',qq_);
%             print('-dpng','-r300',tformfname);
%             close(f);
%         end
%
%
%     end
%
%     % COMMENT FOR NOW:
%     ttrow0index=round(median(1:size(ttrowT,1)));
%     for qq_=cols_array
%         for kk_=rows_array
%
%             T2(qq_,kk_,:,:)=squeeze(T2(qq_,kk_,:,:))/squeeze(ttrowT(ttrow0index,:,:));
%
%         end
%     end
%
%
%     % COLUMNS:
%     jjj=0;
%     ttcolT=nan(numel(cols_array),3,3);
%     rowmin_=max([2,rowmin]);
%     for kk_=rowmin_:rowmax
%         jjj=jjj+1;
%         for qq_=cols_array
%             [qq_,kk_]
%             mp_llpx=[];
%             mp_llpy=[];
%             llpx=[];
%             llpy=[];
%             if sum(ismember(rowsforcol{colt==qq_},kk_))
%
%                 if isfinite(Tt2(qq_,kk_,1,1)) && isfinite(Tt2(qq_,kk_-1,1,1))
%
%                     %disp(['COLS: ' num2str([qq_,kk_])]);
%
%                     t00=projective2d;
%                     t00.T=squeeze(T2(qq_,kk_,:,:));
%                     t01=projective2d;
%                     t01.T=squeeze(T2(qq_,kk_-1,:,:));
%
%                     %if isfinite(Tt(qq_,kk_,2,1,1)) % 4 or 2?
%
%                     % t01.T=squeeze(Tt(qq_,kk_,2,:,:))*t00.T;
%
%                     [mp_llpx_,mp_llpy_]=transformPointsForward(t00,matchedPoints{qq_,kk_}.Location(:,1),matchedPoints{qq_,kk_}.Location(:,2));
%                     [llpx_,llpy_]        =transformPointsForward(t01,matchedPointsPrev{qq_,kk_}.Location(:,1),matchedPointsPrev{qq_,kk_}.Location(:,2));
%                     % t22=projective2d;
%                     % t22.T=squeeze(Tt(qq_,kk_,2,:,:));
%                     % [llpx_,llpy_]=transformPointsForward(t22,llpx_,llpy_);
%                     llpx=[llpx;llpx_];
%                     llpy=[llpy;llpy_];
%                     mp_llpx=[mp_llpx;mp_llpx_];
%                     mp_llpy=[mp_llpy;mp_llpy_];
%
%                     %end
%
%                 end
%
%             end
%         end
%
%
%         if numel(llpx)>mincountpoints
%
%             [ttcol, ~, ~, ~] =estimateGeometricTransform([mp_llpx,mp_llpy],[llpx,llpy], 'affine', 'MaxNumTrials',1e5);
%             for qq__=cols_array
%                 T2(qq__,kk_,:,:)=squeeze(T2(qq__,kk_,:,:))*ttcol.T;
%             end
%
%         else
%             ttcol=projective2d;
%
%         end
%
% %         yy=Tt2(cols_array,kk_,1,1);
% %         isf=~isfinite(yy);
% %         fisf=find(isf);
% %         for mm=1:sum(isf)
% %             Tt2(cols_array(fisf(mm)),kk_,:,:)=squeeze(T2(cols_array(fisf(mm)),kk_,:,:))*ttcol.T;
% %         end
%
%         %disp('TTCOL:')
%         %disp(ttcol.T)
%         ttcolT(jjj,:,:)=ttcol.T;
%         %tformstage_(qq_,kk_,:,:)=ttrow.T*squeeze(tformstage_(qq_,kk_,:,:));
%         %ttrefine2 = images.geotrans.PolynomialTransformation2D([[mp_llpx,mp_llpy];[mp_llx,mp_lly]],[[llpx,llpy];[llx,lly]],3);
%         %ttrefine2i = images.geotrans.PolynomialTransformation2D([[llpx,llpy];[llx,lly]],[[mp_llpx,mp_llpy];[mp_llx,mp_lly]],3);
%         %II=imread(sprintf('/Users/ogliore/Data/deepzoom/strip%03d.png',qq_));
%         %outputView = imref2d(size(II));
%         %IIw  = imwarp(II,ttrow,'OutputView',outputView);
%
%         %(qq_,kk_)
%
%         %Tt2(qq_,kk_,:,:) = squeeze(Tt2(qq_,kk_,:,:))*ttrefine2.T;
%
%         if doprintrowpoints
%             f=figure('visible','off');
%             % [pfx,pfy]=transformPointsInverse(ttrefine2i,[mp_llpx;mp_llx],[mp_llpy;mp_lly]);
%             [pfx,pfy]=transformPointsForward(ttcol,mp_llpx,mp_llpy);
%             plot(llpx,llpy,'b.',mp_llpx,mp_llpy,'r.',pfx,pfy,'bo')
%             tformfname=sprintf('refine/colpoints_%03d.png',kk_);
%             print('-dpng','-r300',tformfname);
%             close(f);
%         end
%
%
%     end
%
%     % COMMENT FOR NOW:
%     ttcol0index=round(median(1:size(ttcolT,1)));
%     jjj=0;
%     for qq_=cols_array
%         jjj=jjj+1;
%         for kk_=rowmin_:rowmax
%             if sum(ismember(rowsforcol{colt==qq_},kk_))
%                 T2(qq_,kk_,:,:)=squeeze(T2(qq_,kk_,:,:))/squeeze(ttcolT(ttcol0index,:,:));
%             end
%         end
%     end




% else
%
%     for qq_=cols_array
%         rows_array=rowsforcol{colt==qq_}';
%         yy=Tt2(qq_,rows_array,1,1);
%         isf=~isfinite(yy);
%         fisf=find(isf);
%         for mm=1:sum(isf)
%             Tt2(qq_,rows_array(fisf(mm)),:,:)=squeeze(T2(qq_,rows_array(fisf(mm)),:,:));
%         end
%
%
%    end
% end

xMin=inf;
xMax=-inf;
yMin=inf;
yMax=-inf;
for qq=cols_array
    for kk=rowsforcol{colt==qq}'
        tt__=projective2d;
        tt__.T=squeeze(T2(qq,kk,:,:));
        [XLim, YLim] = outputLimits(tt__,  [1 imsz], [1 imsz]);
        xMin=min([XLim(:);xMin]);
        xMax=max([XLim(:);xMax]);
        yMin=min([YLim(:);yMin]);
        yMax=max([YLim(:);yMax]);
        
    end
    
end

fprintf(fid,'%.15f,%.15f,%.15f,%.15f\n',xMin,xMax,yMin,yMax);

XPos=nan(colmax,1);
for qq=cols_array
    
    rows_array=rowsforcol{colt==qq}';
    
    for kk=rows_array
        tt__=projective2d;
        tt__.T=squeeze(T2(qq,kk,:,:));
        %[XLim(kk,:), YL(qq,kk,:)] = outputLimits(tt__,  [1 imsz], [1 imsz]);
        [XL(qq,kk,:), YL(qq,kk,:)] = outputLimits(tt__,  [1 imsz], [1 imsz]);
    end
    
    XL_h(qq,rows_array(1))=round(XL(qq,rows_array(1),2));
    XL_l(qq,rows_array(1):(rows_array(end)-1))=round(mean([XL(qq,(rows_array(1)+1):rows_array(end),2);XL(qq,rows_array(1):rows_array(end)-1,1)],1));
    XL_h(qq,rows_array(2):rows_array(end))=XL_l(qq,(rows_array(2)-1):rows_array(end)-1);
    XL_l(qq,rows_array(end))=round(XL(qq,rows_array(end),1));
    XPos(qq)=round(min(min(XL(qq,:,:)))-xMin);
    
    % NaNs?
    XL_l(qq,isnan(XL_l(qq,:)))=round(XL(qq,isnan(XL_l(qq,:)),1));
    XL_h(qq,isnan(XL_h(qq,:)))=round(XL(qq,isnan(XL_h(qq,:)),2));
    
end


jjj=0;
cropY=nan(numel(cols_array)-1,1);
cropdiff=cell(numel(cols_array)-1,1);
CCYOVERLAP=nan(numel(cols_array)-1,1);
% maxYL1=nan(numel(cols_array)-1,1);
% minYL2=nan(numel(cols_array)-1,1);

for qq=cols_array(2:end)
    
    jjj=jjj+1;
    
    maxYL1(jjj)=max(YL(qq,:,1));
    minYL2(jjj)=min(YL(qq,:,2));
    
    if qq<cols_array(end)
        %cropdiff(jjj)=min(YL(qq-1,:,2))-max(YL(qq,:,1));
        cropdiff{jjj}=YL(qq-1,:,2)-YL(qq,:,1);
    end
    
    %     if max(YL(qq,:,1))<min(YL(qq-1,:,2)) %max(max(YL(qq-1,:,:)))-min(min(YL(qq,:,:)));
    %         %cropY(jjj)=round(mean([max(max(YL(qq-1,:,:))),min(min(YL(qq,:,:)))])); %round(mean([max(YL(qq,:,1)),min(YL(qq-1,:,2))]));
    %         cropY(jjj)=round(mean([max(YL(qq,:,1)),min(YL(qq-1,:,2))]));
    %     else
    %        %cropY(jjj)=ceil(max(YL(qq,:,1)));
    %        cropY(jjj)=floor(min(YL(qq-1,:,2)));
    %     end
    %     %end
    
    cropY(jjj)=round(nanmean(T2(qq,:,3,2)));
    
end

% fun0 = @(bb) rowimagecount(bb,colmax,cols_array,rowsforcol,colt,Tt2,imsz,xMin,xMax,XL,YL);
%
% nsteps=1e3;
% bb0=linspace(yMin,yMax,nsteps);
% nimmatch=zeros(nsteps,1);
%
% for ii=1:nsteps
%
%     nimmatch(ii)=fun0([bb0(1),bb0(ii)]);
%
% end
%
% % [~,nimmatch_ind]=findpeaks(diff(nimmatch));
% [~,nimmatch_ind]=findpeaks(diff(nimmatch),'NPeaks',numel(cols_array)-1);
%
% BB=bb0(nimmatch_ind);
% BB(end)=yMax;
%
% [~,rowindex]=rowimagecount(BB,colmax,cols_array,rowsforcol,colt,Tt2,imsz,xMin,xMax,XL,YL);

%BB=linspace(yMin,yMax,floor(numel(cols_array)/4)+1);
BB=linspace(yMin,yMax,17);
[~,rowindex]=rowimagecount(BB,colmax,cols_array,rowsforcol,colt,T2,imsz,xMin,xMax,XL,YL);



imageindices=cell(numel(BB)-1,1);

for jjj=1:numel(BB)-1
    
    index = false(size(rowindex));
    
    for qq=cols_array
        
        rows_array=rowsforcol{colt==qq}';
        
        for kk=rows_array
            
            if sum(rowindex{qq,kk}==jjj)>0
                
                imageindices{jjj}(end+1,:)=[qq,kk];
                
            end
            
        end
        
        
    end
    
end


%


% plot(1:numel(maxYL1),minYL2-maxYL1)
%
% ncompositetrials=1e6;
% ind0t=cell(ncompositetrials,1);
% compositest=zeros(ncompositetrials,1);
% for ii=1:ncompositetrials
% [compositest(ii),ind0t{ii}]=compositepath(minYL2,maxYL1);
% end
% min(compositest)

% squeeze(Tt2(105,59,:,:))
%
% squeeze(Tt2(105,58,:,:))*squeeze(Tt(105,59,1,:,:))

CCY2=[cropY;round(max(YL(cols_array(end),:,2)))];
CCY1=nan(size(CCY2));
CCY1(1)=round(min(YL(cols_array(1),:,1)));




% jjj=0;
compositenumber=zeros(size(CCY2));
% for qq=cols_array(2:end-1)
%
%     jjj=jjj+1;
%
%     kk=0;
%     if max(YL(qq+kk,:,1))>min(YL(qq-1,:,2)) && qq+kk<=size(YL,1)  %max(max(YL(qq-1,:,:)))-min(min(YL(qq,:,:)));
%
%     end
%     kk=kk-1;
%     qq+kk
%     if kk>0
%         CCY2(jjj)=round(mean([max(YL(qq+kk,:,1)),min(YL(qq,:,2))])); %round(mean([max(YL(qq,:,1)),min(YL(qq-1,:,2))]));
%         compositenumber(jjj:jjj+kk)=1;
%     end
% end


CCY1(2:end)=CCY2(1:end-1);

% ffstart=find(diff(compositenumber~=0)==1);
% ffend=find(diff(compositenumber~=0)==-1);
% for ii=1:numel(ffstart)
%     CCY1((ffstart(ii)+1):ffend(1))=CCY1((ffstart(ii)+1));
%     CCY2((ffstart(ii)+1):ffend(1))=CCY2(ffend(ii));
% end

%CCYHEIGHT=CCY2-CCY1;
% CCYHEIGHT(end)=imsz;
% CCYHEIGHT(1)=imsz;

%jjj=0;
dolimitsfigure=false;
if dolimitsfigure
    figure('visible','off');
    hold on
    for qq=cols_array
        jjj=jjj+1;
        rows_array=rowsforcol{colt==qq}';
        for kk=rows_array
            
            plot([XL(qq,kk,1),XL(qq,kk,2),XL(qq,kk,2),XL(qq,kk,1),XL(qq,kk,1)],[YL(qq,kk,1),YL(qq,kk,1),YL(qq,kk,2),YL(qq,kk,2),YL(qq,kk,1)],'k-')
            plot([XL_l(qq,kk),XL_h(qq,kk),XL_h(qq,kk),XL_l(qq,kk),XL_l(qq,kk)],[CCY1(jjj),CCY1(jjj),CCY2(jjj),CCY2(jjj),CCY1(jjj)],'r-')
            
        end
    end
    axis image
    hold off
    fname=sprintf('cropbounds_%03d_%03d.pdf',min(cols_array),max(cols_array));
    set(gca,'ColorScale','log')
    print('-dpdf',fname);
end


% jjj=0;
CCYHEIGHT=round(diff(BB));

fprintf(fid,'%d\n',numel(imageindices));


for jjj=1:numel(imageindices) %qq=cols_array
    
    %    jjj=jjj+1;
    
    fprintf(fid,'%d,%d,%d\n',size(imageindices{jjj},1),round(BB(jjj)),CCYHEIGHT(jjj)); %,xMin,xMax, YLimRmin(jjj),YLimRax(jjj));
    
    %rows_array=rowsforcol{colt==qq}';
    
    
    for mm=1:size(imageindices{jjj},1) %2:nf
        
        qq=imageindices{jjj}(mm,1);
        kk=imageindices{jjj}(mm,2);
        
        fname=[imgfolder sprintf(['%03d_%03d.' imgext],kk,qq)];
        
        tm=squeeze(T2(qq,kk,:,:));
        
        fprintf(fid,'%s,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f\n',fname,OFFSET(qq,kk),tm(1,1),tm(1,2),tm(2,1),tm(2,2),tm(3,1),tm(3,2));
        
    end
    
    
    
end

fclose(fid);




for qq=cols_array
    for kk=rowsforcol{colt==qq}'
        if qq>1
            if isfinite(Tt2(qq,kk,1,1)) && isfinite(Tt2(qq-1,kk,1,1))
                
                t2=projective2d;
                t2.T=squeeze(Tt2(qq,kk,:,:)); %*squeeze(Tt(qq,kk,4,:,:));
                
                t1=projective2d;
                t1.T=squeeze(Tt2(qq-1,kk,:,:));
                
                if numel(llp_t{qq,kk})>0
                    [llpx,llpy]=transformPointsForward(t1,llp_t{qq,kk}(:,1),llp_t{qq,kk}(:,2));
                    [mp_llpx,mp_llpy]=transformPointsForward(t2,mp_llp_t{qq,kk}(:,1),mp_llp_t{qq,kk}(:,2));
                    
                    errortformT(qq,kk,2)=sum((llpx-mp_llpx).^2+(llpy-mp_llpy).^2);
                    
                end
                
            end
        end
        
        if kk>1
            if isfinite(Tt2(qq,kk,1,1)) && isfinite(Tt2(qq,kk-1,1,1))
                
                t2=projective2d;
                t2.T=squeeze(Tt2(qq,kk,:,:)); %*squeeze(Tt(qq,kk,4,:,:));
                
                t1=projective2d;
                t1.T=squeeze(Tt2(qq,kk-1,:,:)); %squeeze(Tt(qq,kk,3,:,:))*t2.T;
                
                if numel(ll_t{qq,kk})>0
                    
                    
                    [llx,lly]=transformPointsForward(t1,ll_t{qq,kk}(:,1),ll_t{qq,kk}(:,2));
                    [mp_llx,mp_lly]=transformPointsForward(t2,mp_ll_t{qq,kk}(:,1),mp_ll_t{qq,kk}(:,2));
                    
                    errortformT(qq,kk,1)=sum((llx-mp_llx).^2+(lly-mp_lly).^2);
                end
                
                
            end
            
        end
    end
end

figure('visible','off');
imagesc(errortformT(colmin:colmax,rowmin:rowmax,1)); axis image; colorbar
fname=sprintf('etform1_%03d_%03d.png',min(cols_array),max(cols_array));
set(gca,'XDir','reverse');
%set(gca,'ColorScale','log')
print('-dpng','-r300',fname);
close(gcf);

figure('visible','off');
imagesc(errortformT(colmin:colmax,rowmin:rowmax,2)); axis image; colorbar
fname=sprintf('etform2_%03d_%03d.png',min(cols_array),max(cols_array));
set(gca,'XDir','reverse');
%set(gca,'ColorScale','log')
print('-dpng','-r300',fname);
close(gcf);


for qq=cols_array
    for kk=rowsforcol{colt==qq}' %2:nf
        
        t00=projective2d;
        t00.T=squeeze(T2(qq,kk,:,:));
        
        [imgshiftx(qq,kk),imgshifty(qq,kk)]=transformPointsForward(t00,0,0);
        
    end
end

cc=distinguishable_colors(numel(cols_array));
cc2=distinguishable_colors(rowmax);
figure('visible','off');
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
    text(max(imgshiftx(ii,rowsforcol{colt==ii}))+imsz,mean(imgshifty(ii,rowsforcol{colt==ii})),sprintf('%03d',ii))
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

figure('visible','off');
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

%imwrite(uint8(255*matchedcount./max(max(matchedcount))),sprintf('matchedcount_%03d_%03d.png',min(cols_array),max(cols_array)));
figure('visible','off');
imagesc(matchedvcount(colmin:colmax,rowmin:rowmax));
axis image
xticklabels = rowmin:rowmax;
xticks = 1:(rowmax-rowmin+1);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = colmin:colmax;
yticks = 1:(colmax-colmin+1);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
xlabel('Column');
ylabel('Row');
title('Vertical Matches')
colorbar
print('-dpng','-r300',sprintf('matchedvcount_%03d_%03d.png',min(cols_array),max(cols_array)));
close(gcf);

figure('visible','off');
imagesc(matchedcount(colmin:colmax,rowmin:rowmax));
axis image
xticklabels = rowmin:rowmax;
xticks = 1:(rowmax-rowmin+1);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = colmin:colmax;
yticks = 1:(colmax-colmin+1);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
xlabel('Column');
ylabel('Row');
title('Horizontal Matches')
colorbar
print('-dpng','-r300',sprintf('matchedcount_%03d_%03d.png',min(cols_array),max(cols_array)));
close(gcf);

[~,b]=min(abs(imgshifty(cols_array(1),:)));

origin=[qqo,kko];

[gridcoordsx,gridcoordsy]=meshgrid(linspace(0,imsz*(1-overlap)*rowmax,rowmax),linspace(0,imsz*(1-overlap)*colmax,colmax));
gridcoordsx=gridcoordsx(origin(1),origin(2))-gridcoordsx;
gridcoordsy=gridcoordsy-gridcoordsy(origin(1),origin(2));

grid_distance=sqrt((gridcoordsx-imgshiftx).^2 + (gridcoordsy-imgshifty).^2);

figure('visible','off');
imagesc(grid_distance(colmin:colmax,rowmin:rowmax)); axis image; colorbar
gridfname=sprintf('grid_error_%03d_%03d.png',min(cols_array),max(cols_array));
set(gca,'XDir','reverse');
print('-dpng','-r300',gridfname);
close(gcf);

figure('visible','off');
imagesc(sum(ET1(colmin:colmax,rowmin:rowmax,:),3)); axis image; colorbar
fname=sprintf('tform1_%03d_%03d.png',min(cols_array),max(cols_array));
set(gca,'XDir','reverse');
%set(gca,'ColorScale','log')
print('-dpng','-r300',fname);
close(gcf);

figure('visible','off');
imagesc(sum(ET2(colmin:colmax,rowmin:rowmax,:),3)); axis image; colorbar
fname=sprintf('tform2_%03d_%03d.png',min(cols_array),max(cols_array));
axis image;
set(gca,'XDir','reverse');
%set(gca,'ColorScale','log')
print('-dpng','-r300',fname);
close(gcf);

figure('visible','off');
imagesc(ST2a(colmin:colmax,rowmin:rowmax));
axis image;
set(gca,'XDir','reverse');
colorbar;
fname=sprintf('shifterry_%03d_%03d.png',min(cols_array),max(cols_array));
set(gca,'ColorScale','log')
print('-dpng','-r300',fname);
close(gcf);

figure('visible','off');
imagesc(ST1a(colmin:colmax,rowmin:rowmax));
axis image;
set(gca,'XDir','reverse');
colorbar;
fname=sprintf('shifterrx_%03d_%03d.png',min(cols_array),max(cols_array));
set(gca,'ColorScale','log')
print('-dpng','-r300',fname);
close(gcf);


pyvipscall=['python ~/Documents/affinetransform6.py ~/Documents/MATLAB/' imageinfofilename ' deepzoom.png'];

fprintf(fido,'%s\n',pyvipscall);

if dovips
    [~,~]=system(pyvipscall);
end
fclose(fido);

%%
dotesttransform=false;

if dotesttransform
    close all
    I1=imread('/Users/ogliore/Data/Tescan/Ryan/Acfer182/HighRes/055_113.png');
    I2=imread('/Users/ogliore/Data/Tescan/Ryan/Acfer182/HighRes/056_113.png');
    It={I1,I2};
    
    tt=projective2d;
    tt.T=eye(3); %squeeze(Tt(114,50,2,:,:))*squeeze(Tt2(qqo,kko,:,:));
    TTforms(1)=tt;
    %M=squeeze(Tt(110,59,1,:,:));
    tt=projective2d;
    tt.T=squeeze(Tt(113,56,1,:,:)); %squeeze(Tt(114,51,3,:,:))*squeeze(Tt2(qqo,kko,:,:));
    %tt.T=squeeze(Tt2(114,51,:,:));
    %tt.T=M;
    TTforms(2)=tt;
    
    % tt=projective2d;
    % TTforms(1)=squeeze(Tt2(110,54,:,:));
    % tt.T=squeeze(Tt2(112,55,:,:));
    % TTforms(2)=tt;
    
    [xlim(1,:), ylim(1,:)] = outputLimits(TTforms(1), [1 size(I1,1)], [1 size(I1,2)]);
    [xlim(2,:), ylim(2,:)] = outputLimits(TTforms(2), [1 size(I2,1)], [1 size(I2,2)]);
    
    
    maxImageSize = [size(I1,1),size(I1,2)];
    
    % Find the minimum and maximum output limits
    xMin = min([xlim(:)]);
    xMax = max([maxImageSize(2); xlim(:)]);
    
    yMin = min([ylim(:)]);
    yMax = max([maxImageSize(1); ylim(:)]);
    
    % Width and height of panorama.
    PANOwidth  = round(xMax - xMin);
    PANOheight = round(yMax - yMin);
    
    % Initialize the "empty" panorama.
    panorama = zeros([PANOheight PANOwidth 1], 'like', I1);
    blender = vision.AlphaBlender('Operation', 'Binary mask', 'MaskSource', 'Input port');
    
    % Create a 2-D spatial reference object defining the size of the panorama.
    xLimits = [xMin xMax];
    yLimits = [yMin yMax];
    panoramaView = imref2d([PANOheight PANOwidth], xLimits, yLimits);
    
    
    
    
    % Create the panorama.
    for ii = 1:2
        
        I = It{ii};
        
        % Transform I into the panorama.
        warpedImage = imwarp(I, TTforms(ii), 'OutputView', panoramaView);
        
        % Generate a binary mask.
        mask = imwarp(true(size(I,1),size(I,2)), TTforms(ii), 'OutputView', panoramaView);
        
        % Overlay the warpedImage onto the panorama.
        panorama = step(blender, panorama, warpedImage, mask);
    end
    
    figure
    imshow(panorama)
    squeeze(T2(113,56,:,:))/squeeze(T2(113,55,:,:))
    squeeze(Tt(113,56,1,:,:))
end

%%
% %%
% aa2=112;
% bb2=58;
%
% aa1=114;
% bb1=51;
%
% [A,~,~]=calculateMosaicPositions(colmax,rowmax,cols_array,rowsforcol{colt==qq1}',G,ny,nx,aa2,bb2,isoutlier,Tt,false);
% [B,~,~]=calculateMosaicPositions(colmax,rowmax,cols_array,rowsforcol{colt==qq1}',G,ny,nx,aa1,bb1,isoutlier,Tt,false);
% fprintf("\n\n\n\n");
% fprintf("%d,%d->%d,%d->%d,%d\n",aa2,bb2,aa1,bb1,aa2,bb2);
% squeeze(B(aa2,bb2,:,:))*squeeze(A(aa1,bb1,:,:))
%
% [A,~,~,TmA]=calculateMosaicPositions(colmax,rowmax,aa1,bb1,G,ny,nx,aa2,bb2,isoutlier,Tt,false);
% [B,~,~,TmB]=calculateMosaicPositions(colmax,rowmax,aa1,bb1,G1,ny,nx,aa2,bb2,isoutlier,Tt,false);
% fprintf("\n\n\n\n");
% fprintf("G:%d,%d->%d,%d,G1:%d,%d->%d,%d\n",aa2,bb2,aa1,bb1,aa2,bb2,aa1,bb1);
% squeeze(A(aa1,bb1,:,:))
% squeeze(B(aa1,bb1,:,:))
%
% % fprintf("\n\n\n\n");
% % inv(squeeze(A(110,57,:,:)))*squeeze(A(110,57,:,:))
