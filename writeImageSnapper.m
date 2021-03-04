function res=writeImageSnapper(varargin)

ptsxyz=varargin{1};
Zstage=varargin{2};
viewfield=varargin{3};
coordn=varargin{4};
imagepath=varargin{5};
imageformat=varargin{6};
dowritepythonfiles=varargin{7};

autofocus=0;
autogain=0;
stitching=0;
snapshotsize=4;
panoramasize=100;
panoramamaxw=1e4;
%Z=21.484375; % in mm (read off "Z" in Stage Control);
% WD=14.9333;  % in mm (read off "WD & Z" in Stage Control);
%viewfield=250;  % in *microns* (read off info panel, convert to microns)
imagebasename='Snap';

npts=size(ptsxyz,1);

fname='C:\Users\Tescan\Documents\ImageSnapper.xml';

fnamep='C:\Users\Tescan\Documents\PythonCoords.txt';

fnamepv='C:\Users\Tescan\Documents\PythonFOV.txt';

if dowritepythonfiles
    fid0=fopen(fnamepv,'w');
    fprintf(fid0,'%.15f',1e-6*viewfield);
    fclose(fid0);
end

fid=fopen(fname,'w');

if dowritepythonfiles
    fid2=fopen(fnamep,'w');
end

fprintf(fid,'<?xml version="1.0" encoding="utf-8" standalone="yes"?>\n');
fprintf(fid,'<ImageSnapperProject>\n');
fprintf(fid,'<Settings AutoFocus="%d" AutoGainBlack="%d" Stitching="%d" ShadingCorrection="0" AddInfobar="0" Path="%s" ImageFormat="%s" SnapshotSize="%d" PanoramaSize="%d" PanoramaMaxW="%d"/>\n',autofocus,autogain,stitching,imagepath,imageformat,snapshotsize,panoramasize,panoramamaxw);
fprintf(fid,'<Samples Count="%d">\n',npts);

formatcodex=sprintf('%%0%dd',max( [ceil(log10(max(coordn(:,1)))),3]) );
formatcodey=sprintf('%%0%dd',max( [ceil(log10(max(coordn(:,2)))),3]) );

if numel(varargin)==7
    mask=ones(npts,1);
else
    mask=varargin{8};
end

for ii=1:npts
    if mask(ii)
        samplenamei=[sprintf(formatcodex,coordn(ii,1)) '_' sprintf(formatcodey,coordn(ii,2))];
        fprintf(fid,['<PointSample Name="%s" Z="%s" WD="%.9f" ViewField="%.6f" Overlapping="0" ImageBaseName="%s" Pos="%.4f,%.4f"/>\n'],samplenamei,Zstage,ptsxyz(ii,3),viewfield*1e-6,imagebasename,ptsxyz(ii,1)*1e3,ptsxyz(ii,2)*1e3);
        if dowritepythonfiles
            fprintf(fid2,'%s,%.9f,%.4f,%.4f\n',samplenamei,ptsxyz(ii,3),ptsxyz(ii,1)*1e3,ptsxyz(ii,2)*1e3);
        end
    end
end


fprintf(fid,'</Samples>\n');
fprintf(fid,'</ImageSnapperProject>\n');

res=fclose(fid);

if dowritepythonfiles
    fclose(fid2);
end

fprintf('ImageSnapper XML written to: %s\n\n',fname);