function s = processEDAXmaps(mapfile)

mapnumber = str2double(extractBetween(mapfile,'Area ','/'));
mapname = extractBetween(mapfile,'EDAX_EDS/','/');
mapname=mapname{1};
dosaveplots=false;
shiftthresh=50;

%mapfile=mapfile{1};
% clearvars
% close all
%
% mapfile='~/Data/EDAX_EDS/ProgramData/EDAX/Ryan/Cole/Mapping/Lsm/Acfer094_EDS/New Sample/Area 20/map20200814034547577_0.spd';

s=readspc([mapfile(1:end-3) 'spc']);

fid=fopen(mapfile);
tag=fread(fid,16,'*char')';
version=fread(fid,1,'int32=>double');
nSpectra=fread(fid,1,'int32=>double');
nPoints=fread(fid,1,'int32=>double');
nLines=fread(fid,1,'int32=>double');
nChannels=fread(fid,1,'int32=>double');
countBytes=fread(fid,1,'int32=>double');
dataOffset=fread(fid,1,'int32=>double');
nFrames=fread(fid,1,'int32=>double');
fName=fread(fid,120,'*char')';
filler=fread(fid,900,'*char')';

dE=1e3*(s.energy(2)-s.energy(1));
energy=double((1e3*s.startEnergy):dE:dE*(nChannels-1));


if countBytes==1
    typename='uint8';
elseif countBytes==2
    typename='uint16';
elseif countBytes==4
    typename='uint32';
end

formatstr=['*' typename];

data=zeros(nSpectra,nChannels,typename);

sumspectrum=zeros(nChannels,1);
maxpixel=zeros(nChannels,1);

fseek(fid,dataOffset,'bof');
for  ii=1:nSpectra
    if mod(ii,1e5)==0
        fprintf('%d,...',ii);
    end
    dd=fread(fid,nChannels,formatstr);
    dd(dd>2^16-1e3)=0; % Error trapping
    if numel(dd)~=nChannels
        keyboard
    end
    data(ii,:)=dd;
    dd_double=double(dd);
    sumspectrum=sumspectrum+dd_double;
    maxpixel(dd>maxpixel)=dd_double(dd_double>maxpixel);
end
fclose(fid);
fprintf('%d\n',nSpectra);
data=reshape(data,nPoints,nLines,size(data,2));

% figure
% imagesc(sum(data,3))
% axis image
% title('Summed X-ray counts');


% Fe_Kalpha_energy=6408;
% energywindow=50;
%
% figure
% imagesc(sum(data(:,:,energy>Fe_Kalpha_energy-energywindow & energy<Fe_Kalpha_energy+energywindow),3))
% axis image
% title('Fe K Map');
%
% figure
% subplot(3,1,1)
% plot(energy,sumspectrum)
% xlabel('Energy (eV)');
% ylabel('Counts');
% title('Sum Spectrum')
%
% subplot(3,1,2)
% plot(energy,maxpixel);
% xlabel('Energy (eV)');
% ylabel('Counts');
% title('Max Pixel Spectrum')
%
% subplot(3,1,3)
% plot(energy,squeeze(data(200,200,:)));
% xlabel('Energy (eV)');
% ylabel('Counts');
% title('Spectrum of Pixel = 200,200')

C=277;
N=392.4;
O=524.9;
F=676.8;
Na=1040.98;
Mg=1253.6;
Al=1486.7;
Si=1739.98;
P=2013.7;
S=2307.84;
Cl=2622.39;
K=3313.8;
Ca=3691.68;
Sc=4090.6;
Ti=4510.84;
V=4952.2;
Cr=5414.72;
Mn=5898.75;
Fe=6403.84;
Co=6930.32;
Ni=7478.15;
Cu=8047.78;
Zn=8638.86;

% ZAF factors
C_zaf =0.24086028;
N_zaf =0.22780545;
O_zaf =0.3124552;
F_zaf =0.25247247;
Na_zaf =0.37853913552;
Mg_zaf =0.539219618;
Al_zaf =0.5861142472;
Si_zaf =0.70573243125;
P_zaf =0.699829600896;
S_zaf =0.787409296128;
Cl_zaf=0.7962099872;
K_zaf=0.86540901005;
Ca_zaf=0.903158142705;
Sc_zaf=0.857133971766;
Ti_zaf=0.85805826147;
V_zaf=0.860303081925;
Cr_zaf=0.897575196192;
Mn_zaf=0.87302754528;
Fe_zaf=0.823586657274;
Co_zaf=0.805451020617;
Ni_zaf=0.823216709047;
Cu_zaf=0.78968804969;
Zn_zaf=0.795078201294;

% Remove:
% N: sits on shoulder of large O peak
% F: Fe-L interference

allelementsnames={'C','N','O','F','Na','Mg','Al','Si','P','S','Cl','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn'};
allelements=[C,N,O,F,Na,Mg,Al,Si,P,S,Cl,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn];
allelementsZ=[6,7,8,9,11,12,13,14,15,16,17,19,20,21,22,23,24,25,26,27,28,29,30];
allelementszaf=[C_zaf,N_zaf,O_zaf,F_zaf,Na_zaf,Mg_zaf,Al_zaf,Si_zaf,P_zaf,S_zaf,Cl_zaf,K_zaf,Ca_zaf,Sc_zaf,Ti_zaf,V_zaf,Cr_zaf,Mn_zaf,Fe_zaf,Co_zaf,Ni_zaf,Cu_zaf,Zn_zaf];
nelements=numel(allelements);

[~,LOCS]=findpeaks(sumspectrum,'MinPeakProminence',max(sumspectrum)/500);

% plot(energy,sumspectrum,energy(LOCS),sumspectrum(LOCS),'*')
% yl=ylim;
% hold on
% plot([Fe,Fe],[yl(1),yl(2)],[O,O],[yl(1),yl(2)],[Mg,Mg],[yl(1),yl(2)],[S,S],[yl(1),yl(2)],[Si,Si],[yl(1),yl(2)],[Ca,Ca],[yl(1),yl(2)],[C,C],[yl(1),yl(2)])
% hold off

[~,b]=min(abs(energy(LOCS)-Fe));
fe_shift=energy(LOCS(b))-Fe;

[~,b]=min(abs(energy(LOCS)-O));
o_shift=energy(LOCS(b))-O;

[~,b]=min(abs(energy(LOCS)-Mg));
mg_shift=energy(LOCS(b))-Mg;

[~,b]=min(abs(energy(LOCS)-S));
s_shift=energy(LOCS(b))-S;

[~,b]=min(abs(energy(LOCS)-Si));
si_shift=energy(LOCS(b))-Si;

[~,b]=min(abs(energy(LOCS)-Ca));
ca_shift=energy(LOCS(b))-Ca;

[~,b]=min(abs(energy(LOCS)-C));
c_shift=energy(LOCS(b))-C;

[~,b]=min(abs(energy(LOCS)-Ni));
ni_shift=energy(LOCS(b))-Ni;

elshift=[C,O,Mg,Si,S,Ca,Fe,Ni];
shiftv=[c_shift,o_shift,mg_shift,si_shift,s_shift,ca_shift,fe_shift,ni_shift];

goodshifts=abs(shiftv)<shiftthresh;

shift_linear_fun= @(x) interp1(elshift(goodshifts),shiftv(goodshifts),x,'linear',0);
%plot([C,O,Mg,Si,S,Ca,Fe,Ni],[c_shift,o_shift,mg_shift,si_shift,s_shift,ca_shift,fe_shift,ni_shift],'O',energy,shift_linear_fun(energy))

if dosaveplots
    ll=linspace(0,10e3,1e3);
    ff=figure;
    plot(ll,shift_linear_fun(ll),allelements,shift_linear_fun(allelements),'s');
    xlabel('Energy (eV)');
    ylabel('Shift');
    print('shift.png','-dpng');
    close(ff)
end

F0=fit(energy(energy>Fe-200 & energy<Fe+200)',sumspectrum(energy>Fe-200 & energy<Fe+200),'gauss1');
energysigmairon=F0.c1;

F1=fit(energy(energy>C-50 & energy<C+50)',sumspectrum(energy>C-50 & energy<C+50),'gauss1');
energysigmacarbon=F1.c1;

energysigma_slope=(energysigmairon-energysigmacarbon)./(Fe-C);
energysigma_int=energysigmairon-energysigma_slope*Fe;


allcountsmap=sum(data,3);
pseudobse=zeros(size(allcountsmap));

allmaps=zeros(size(data,1),size(data,2),nelements);

nsigma=1.5;
bgrangesigma=0.2;

% bgv1t=zeros(nelements,1);
% bgv2t=zeros(nelements,1);
% bg1t=zeros(nelements,1);
% bg2t=zeros(nelements,1);
%
% for ii=1:nelements
%     energysigma=energysigma_slope*allelements(ii)+energysigma_int;
%     bg1=mean(sum(sum(data(:,:,energy>allelements(ii)+shift_linear_fun(allelements(ii))-(nsigmabg+bgenergygrange)*energysigma & energy<allelements(ii)+shift_linear_fun(allelements(ii))-(nsigmabg-bgenergygrange)*energysigma))));
%     bg2=mean(sum(sum(data(:,:,energy>allelements(ii)+shift_linear_fun(allelements(ii))+(nsigmabg-bgenergygrange)*energysigma & energy<allelements(ii)+shift_linear_fun(allelements(ii))+(nsigmabg+bgenergygrange)*energysigma))));
%     bgv1=allelements(ii)+shift_linear_fun(allelements(ii))-nsigmabg*energysigma;
%     bgv2=allelements(ii)+shift_linear_fun(allelements(ii))+nsigmabg*energysigma;
%     bgv1t(ii)=bgv1;
%     bgv2t(ii)=bgv2;
%     bg1t(ii)=bg1;
%     bg2t(ii)=bg2;
%
% end
%
% % manual:
% b0_e1=920;
% b0_e2=955;
% b0s=mean(sum(sum(data(:,:,energy>b0_e1 & energy<b0_e2))));
%
% xv=[bgv1t(1);mean([b0_e1,b0_e2]);bgv1t([8,10:17,19:nelements]);bgv2t([10:17,19:nelements])];
% yv=[bg1t(1);b0s;bg1t([8,10:17,19:nelements]);bg2t([10:17,19:nelements])];

%plot(xv,yv,'*',energy,interp1(xv,yv,energy,'pchip'),energy,sumspectrum); xlim([0,max(xv)])

bgcountsPerPixelTotal=zeros(nelements,1);
for ii=1:nelements
    energysigma=energysigma_slope*allelements(ii)+energysigma_int;
    bgelow=allelements(ii)+shift_linear_fun(allelements(ii))-nsigma*energysigma;
    bgehigh=bgelow+2*nsigma*energysigma;
    %bgelowm= repmat(bgelow,size(data,1),size(data,2));
    %bgehighm=repmat(bgehigh,size(data,1),size(data,2));
    bg1=mean(data(:,:,energy>(bgelow -bgrangesigma*energysigma) & energy<(bgelow +bgrangesigma*energysigma)),3);
    bg2=mean(data(:,:,energy>(bgehigh-bgrangesigma*energysigma) & energy<(bgehigh+bgrangesigma*energysigma)),3);
    %bg1=mean(sum(sum(data(:,:,energy>(bgelow -bgrangesigma*energysigma) & energy<(bgelow +bgrangesigma*energysigma)))));
    %bg2=mean(sum(sum(data(:,:,energy>(bgehigh-bgrangesigma*energysigma) & energy<(bgehigh+bgrangesigma*energysigma)))));
    bgt(1,:,:)=bg1;
    bgt(2,:,:)=bg2;
    bgcountsPerPixel=squeeze(sum(interp1([bgelow,bgehigh],bgt,energy(energy>bgelow & energy<bgehigh),'linear')));
    %     bg2=interp1(xv,yv,bgv2,'pchip');
    %     bgslope=(bg2-bg1)/(bgv2-bgv1);
    %     bgint=bg1-bgslope*bgv1;
    %     %bgcountstotal=(bgslope*bgehigh.^2/2 + bgint*bgehigh) - (bgslope*bgelow.^2/2 + bgint*bgelow);
    %     bgcountstotal=sum(bgslope*energy(energy>bgelow & energy<bgehigh)+bgint);
    
    % Average:
    % bgcountsPerPixel=bgcountstotal/(nPoints*nLines);
    % bgcountsPerPixelTotal(ii)=bgcountsPerPixel;
    
    % Nitrogen, no background:
    if allelements(ii)==N
        bgcountsPerPixel=0;
    end
    intensities_background_subtracted=(1/(1-2*(1-normcdf(nsigma)))).*( sum(data(:,:,energy>bgelow & energy<bgehigh),3) - bgcountsPerPixel ); %-bgcountsPerPixel;
    allmaps(:,:,ii)=intensities_background_subtracted; %.*allelementszaf(ii);
    pseudobse=pseudobse+allmaps(:,:,ii)*2^(-9./sqrt(allelementsZ(ii)));
    if dosaveplots
        ff=figure('visible','off');
        plot(energy(energy>bgelow*0.95 & energy<bgehigh*1.05),squeeze(sum(sum(data(:,:,energy>bgelow*0.95 & energy<bgehigh*1.05)))),...
            energy(energy>bgelow & energy<bgehigh),interp1([bgelow,bgehigh],[sum(bg1(:)),sum(bg2(:))],energy(energy>bgelow & energy<bgehigh),'linear'),'LineWidth',2);
        yl=ylim;
        hold on
        plot([bgelow-bgrangesigma*energysigma,bgelow-bgrangesigma*energysigma],[yl(1),yl(2)],'k:');
        plot([bgelow+bgrangesigma*energysigma,bgelow+bgrangesigma*energysigma],[yl(1),yl(2)],'k:');
        plot([bgehigh-bgrangesigma*energysigma,bgehigh-bgrangesigma*energysigma],[yl(1),yl(2)],'k:');
        plot([bgehigh+bgrangesigma*energysigma,bgehigh+bgrangesigma*energysigma],[yl(1),yl(2)],'k:');
        hold off
        xlabel('Energy (eV)');
        ylabel('Counts');
        set(gca,'LineWidth',2);
        title(allelementsnames(ii),'FontSize',16);
        saveas(ff,[mapname '_' num2str(mapnumber) '_bg_' allelementsnames{ii} '.png'],'png');
        close(ff)
    end
end


gamma=1.5;
pseudobse=pseudobse.^gamma;

s=struct();
s.allelementsZ=allelementsZ;
s.allmaps=allmaps;
s.maxpixel=maxpixel;
s.sumspectrum=sumspectrum;
s.data=data;
s.energy=energy;
s.allelementsnames=allelementsnames;


