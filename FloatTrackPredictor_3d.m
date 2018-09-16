

% This program computes the trajectories of floats launched on a given location 
%and date. 
% The velocity field used is the 3 dimensional climatological mean of the
%ECCO2 model. The horizontal component includes the seasonal variability
%while the vertical component is just the annual mean.
%ECCO2 model. But any velocity field product can be used just edting the corresponding
%variables. 
%The horizontal component includes the seasonal variability while the vertical component
%is just the annual mean.
% Artifitial eddy turbulence can be included based on the standard deviation
%of seasonal velocity components.Different turbulence can be added for
%the horizontal and vertical directions...


close all
FT=input('First time run? (y/n) ','s');
if strcmp(FT,'y');
    clear all
    FT='y';
end


%.........................................................................
%MAIN EXPERIMENT SETTINGS
loc0=[-56 -41 500]; %input('Enter longitude and latitude of float launching [lonE,latN,z]: ');
date0=[4 15]; %input('Enter date of launching [month day]: ');
Ndays=3300; %input('Enter number of days of the simulation: ');
SimuDirection=-1; %input('Select forward(+1) or backward(-1) simulation: ');
disptrack=input('Display simulation on live? (y/n): ','s');

%MINOR SETTINGS
Nfloats=1000; %Number of launched floats
tstep=3; %[days]
KturbXY=1; %Kturb introduces eddy variability folowing a normal distribution
KturbZ=0.5;
KmedXY=1.5; %Kmed amplified the mean velocity field. Kmed=1 no change
KmedZ=1.5;
Nnorm=100; %Number of element in normal distribution
lonlim=[-80 20]; %Region selection
latlim=[-50 10]; %Region selection
depthlim=[0 5000];
keepinside='n'; %Keep particles inside the regional domain. If not,'n', particles stop
file2save=input('Enter name of file to save: ','s');

%AVOIDING COAST
it_stuck=3; %times step stuck
Dstuck0=0.25; %[deg] Minimum radious of movement allowed within it_stuck timesteps
Njump=10; %jump to the Njump nearest grid point when stuck;

%HISTOGRAM PARAMETERS
dx=0.5; dy=0.5; %Horizontal resolution for the histogram
xbin=lonlim(1):dx:lonlim(2); %Histogram longitude grid
ybin=latlim(1):dy:latlim(2); %Histogram latitude grid
bin{1}=xbin; bin{2}=ybin;

%FIGURE PARAMETERS
coltrack=0.5*[1 1 1]; %track color [Red Green Blue]
wtrack=1; %width track
colstart=[1 0 0]; %color start point
sstart=30; %size start point
colend=[0 0 1]; %color end point
send=10; %size end point
%Bathimetry
contbat=[-100 -500 -1000 -2000]; %Bathymetry contours
colbat=0*[1 1 1]; %Color bathymetry
wbat=0.5; %width bathymetry
wcoast=3; %width coast



%.........................................................................
% NOTHING TO CHANGE BELOW THIS LINE.......................................
%.........................................................................

if strcmp(FT,'y');
    
    %Loading data....
    disp(' ')
    disp('Loading data...')
    disp(['Longitude= [' num2str(lonlim(1)) 'E , ' num2str(lonlim(2)) 'E]'])
    disp(['Latitude= [' num2str(latlim(1)) 'N , ' num2str(latlim(2)) 'N]'])
    disp(['Depth= [' num2str(depthlim(1)) 'm ' num2str(depthlim(2)) 'm]'])
    disp(['Seeding ' num2str(Nfloats) ' floats at:  [ '...
        num2str(loc0(1)) 'E , ' num2str(loc0(2)) 'N , ' num2str(loc0(3)) 'm ]'])
    disp(' ');
    cds('/ECCO2/OutputFiles');
    load Uvel_seasons19922011_ECCO2.mat; lonmod=lonmod-360;
    load Vvel_seasons19922011_ECCO2.mat vvel vvel_sd
    load Wvel_clim19922011_ECCO2.mat wvel wvel_sd
    
    %Selecting subregion
   
    
    ilon0=findprox(lonmod,lonlim(1),lonlim(2));
    ilat0=findprox(latmod,latlim(1),latlim(2));
    iz0=findprox(zmod,depthlim(1),depthlim(2));
    lonm=lonmod(ilon0);
    latm=latmod(ilat0);
    zm=zmod(iz0);
    uvel=(uvel(ilon0,ilat0,iz0,:));
    vvel=(vvel(ilon0,ilat0,iz0,:));
    wvel=(wvel(ilon0,ilat0,iz0));
    uvel_sd=(uvel_sd(ilon0,ilat0,iz0,:));
    vvel_sd=(vvel_sd(ilon0,ilat0,iz0,:));
    wvel_sd=(wvel_sd(ilon0,ilat0,iz0));
    
    load coastline
    
    zbat=bathimetry2d(lonm,latm);
end




%------------------------------------------
%Computing variables
lon0=loc0(1);
lat0=loc0(2);
z0=loc0(3);
dayseason=[0 1 2 3]*365/4 + 15;
day0=(date0(1)-1)*30+date0(2);
Nstep=ceil(Ndays/tstep);
tstepS=tstep*24*3600; %[seconds]
tstuck=tstep*it_stuck*24*3600; %[seconds]

if SimuDirection>0
    direction='forwards';
else
    direction='backwards';
end

%Ocean grid points to avoid coast
Llatm=length(latm); Llonm=length(lonm); Lzm=length(zm);
[LATm,LONm,Zm]=meshgrid(latm,lonm,zm);
for jz=1:Lzm
    
isea{jz}=find(~isnan(uvel(:,:,iz0(jz),1)) & ~isnan(wvel(:,:,iz0(jz))))+Llonm*Llatm*(jz-1);
Usea{jz}=uvel(isea{jz}); Usea_sd{jz}=uvel_sd(isea{jz}); 
Vsea{jz}=vvel(isea{jz}); Vsea_sd{jz}=vvel_sd(isea{jz}); 
Wsea{jz}=wvel(isea{jz}); Wsea_sd{jz}=wvel_sd(isea{jz}); 
LATsea{jz}=LATm(isea{jz});
LONsea{jz}=LONm(isea{jz});
end

clear uF vF wF LATm LONm lonJ latJ floatJ


%COMPUTING TRACK
%Initial step
day=day0;
lonF=single(lon0*ones(1,Nfloats));
latF=single(lat0*ones(1,Nfloats));
zF=single(z0*ones(1,Nfloats));
J=0;

disp('Computing....')
disp(' ')
for it=2:Nstep
    
    %day of the year
    day=day + tstep*SimuDirection;
    if day<0
        day=365+day;
    end
    day=rem(day,365);
    
    %season of day
    iseason=findprox(dayseason,day);
    
    %Keep inside domain
    if keepinside=='y'
        in=1:Nfloats;
    else
        out=(latF(it-1,:)<latlim(1) | latF(it-1,:)>latlim(2) |...
            lonF(it-1,:)<lonlim(1) | lonF(it-1,:)>lonlim(2) );%|...
            %zF(it-1,:)<depthlim(1) | zF(it-1,:)>depthlim(2));
        iout=find(out);
        if ~isempty(iout)
            lonF(it,iout)=lonF(it-1,iout);
            latF(it,iout)=latF(it-1,iout);
            zF(it,iout)=zF(it-1,iout);
        end
        in=find(~out);
    end
    Nin=length(in);
    
    if Nin
        
        %Mean velocity of the day
        umF=interp3(latm,lonm,zm,uvel(:,:,:,iseason),latF(it-1,in),lonF(it-1,in),zF(it-1,in));
        vmF=interp3(latm,lonm,zm,vvel(:,:,:,iseason),latF(it-1,in),lonF(it-1,in),zF(it-1,in));
        wmF=interp3(latm,lonm,zm,wvel(:,:,:),latF(it-1,in),lonF(it-1,in),zF(it-1,in));
        
        %Eddy variability
        ueF=interp3(latm,lonm,zm,uvel_sd(:,:,:,iseason),latF(it-1,in),lonF(it-1,in),zF(it-1,in));
        veF=interp3(latm,lonm,zm,vvel_sd(:,:,:,iseason),latF(it-1,in),lonF(it-1,in),zF(it-1,in));
        weF=interp3(latm,lonm,zm,wvel_sd(:,:,:),latF(it-1,in),lonF(it-1,in),zF(it-1,in));
         
        %Avoiding coast
        icoast=find(isnan(umF) | isnan(wmF));
        icoastin=in(icoast);
        Lcoast=length(icoast);
        if Lcoast
            for i=1:Lcoast
                izcoast=findprox(zm,zF(it-1,icoastin(i)));
                
                D=abs(distance([latF(it-1,icoastin(i)) lonF(it-1,icoastin(i))],[LATsea{izcoast} LONsea{izcoast}]));
                
                %Reavoiding coast
                if it>it_stuck
                    mcoast_lon=nanmean(lonF((it-it_stuck):(it-1),icoastin(i)));
                    mcoast_lat=nanmean(latF((it-it_stuck):(it-1),icoastin(i)));
                    Dstuck_lonit=abs(lonF((it-it_stuck):(it-1),icoastin(i))-mcoast_lon);
                    Dstuck_latit=abs(latF((it-it_stuck):(it-1),icoastin(i))-mcoast_lat);
                    Dstuck_lon=mean(Dstuck_lonit);
                    Dstuck_lat=mean(Dstuck_latit);
                    Dstuck=hypot(Dstuck_lon,Dstuck_lat);
                    
                    if Dstuck<Dstuck0
                        [Ds,imin]=sort(D,1,'ascend');
                        imin=imin(Njump);
                        
%                         J=J+1;
%                         lonJ(J,1)=lonF(it-1,icoastin(i));
%                         latJ(J,1)=latF(it-1,icoastin(i));
%                         lonJ(J,2)=LONsea{izcoast}(imin);
%                         latJ(J,2)=LATsea{izcoast}(imin);
%                         floatJ(J)=icoastin(i);
                        
                    elseif Dstuck>=Dstuck0
                        [Dmin,imin]=min(D);
                    end
                else
                    [Dmin,imin]=min(D);
                end
                
                lonF(it-1,icoastin(i))=LONsea{izcoast}(imin);
                latF(it-1,icoastin(i))=LATsea{izcoast}(imin);
                %Avoiding particles flying or digging
                if zF(it-1,icoastin(i))<0 || zF(it-1,icoastin(i))>6000
                    zF(it-1,icoastin(i))=zm(izcoast);
                end
                
                %Recomputing mean velocities within the sea
                umF(icoast(i))=Usea{izcoast}(imin);
                vmF(icoast(i))=Vsea{izcoast}(imin);
                wmF(icoast(i))=Wsea{izcoast}(imin);
                
                %Recomputing eddy variability within the sea
                ueF(icoast(i))=Usea_sd{izcoast}(imin);
                veF(icoast(i))=Vsea_sd{izcoast}(imin);
                weF(icoast(i))=Wsea_sd{izcoast}(imin);
            end
        end
        
        %Mean velocity amplification
        umF=umF*KmedXY;
        vmF=vmF*KmedXY;
        wmF=wmF*KmedZ;
        
            
        %Instantaneous velocity
        clear uF vF wF
        for inF=1:Nin
            %Normal distribution centerd on the mean and with eddy std
            ueFnorm=normrnd(umF(inF),ueF(inF),[Nnorm,1]) - umF(inF);
            veFnorm=normrnd(vmF(inF),veF(inF),[Nnorm,1]) - vmF(inF);
            weFnorm=normrnd(wmF(inF),weF(inF),[Nnorm,1]) - wmF(inF);
            %Random index
            irnd=ceil((Nnorm-1)*rand(1,1));
            uF(inF)=umF(inF) + KturbXY*ueFnorm(irnd);
            vF(inF)=vmF(inF) + KturbXY*veFnorm(irnd);
            wF(inF)=wmF(inF) + KturbZ*weFnorm(irnd);
        end
        
        
        %Time step advection
        %Zonal component
        uFlon=km2deg(uF/1000,6371*cosd(lat0)); %[deg/seconds];
        dlonF=tstepS*uFlon;
        lonF(it,in)=lonF(it-1,in) + dlonF*SimuDirection;
        %Meridional component
        vFlat=km2deg(vF/1000,6371); %[deg/seconds];
        dlatF=tstepS*vFlat;
        latF(it,in)=latF(it-1,in) + dlatF*SimuDirection;
        %Vertical component
        wFz=-wF; %[deg/seconds];
        dzF=tstepS*wFz;
        zF(it,in)=zF(it-1,in) + dzF*SimuDirection;
        
        
        
        disp(['Day ' num2str((it)*tstep) ' of ' num2str(Ndays) ' computed.    '...    
            'Tracking ' num2str(Nfloats) ' ' direction ' floats from:  [ '...
            num2str(loc0(1)) 'E, ' num2str(loc0(2)) 'N, ' num2str(loc0(3)) 'm ]'])
            
   
        
    %No more particles to compute tracks
    else 
        disp(' ')
        disp(['Day ' num2str((it)*tstep) '. No more particles inside domain.'])
        
        Nstep=it;
        Ndays=it*tstep;  
    end
    
     if disptrack=='y'
        %PLOTTING.......
        
        %Floats tracks
        figure(1)
        if it==2
            contour(lonm,latm,zbat,contbat,'color',colbat,'linewidth',wbat); hold on
            line(loncoast,latcoast,'color','k','linewidth',wcoast);
            p1=plot(lonF,latF,'color',coltrack,'linewidth',wtrack);
            p2=plot(lonF(end,:),latF(end,:),'.','color',colend,'markersize',send);
            p3=plot(lonF(1),latF(1),'.','color',colstart,'markersize',sstart);
        else
            delete(p1,p2,p3);
            p1=plot(lonF,latF,'color',coltrack,'linewidth',wtrack);
            p2=plot(lonF(end,:),latF(end,:),'.','color',colend,'markersize',send);
            p3=plot(lonF(1),latF(1),'.','color',colstart,'markersize',sstart);
        end
        axis equal;xlim(lonlim); ylim(latlim);
        xlabel('Longitude [ºE]'); ylabel('Latitude [ºN]')
        title(['Day ' num2str((it)*tstep)]);
        FormatFig
    end
end



%PLOTTING.......

%Floats tracks
figure
contour(lonm,latm,zbat,contbat,'color',colbat,'linewidth',wbat); hold on
line(loncoast,latcoast,'color','k','linewidth',wcoast)
plot(lonF,latF,'color',coltrack,'linewidth',wtrack)
plot(lonF(end,:),latF(end,:),'.','color',colend,'markersize',send);
plot(lonF(1),latF(1),'.','color',colstart,'markersize',sstart);
axis equal;xlim(lonlim); ylim(latlim);
xlabel('Longitude [ºE]'); ylabel('Latitude [ºN]')
title('Float tracks');
FormatFig



%Density of all locations
figure
F=repmat(1:Nfloats,[Nstep,1]);
Nf=hist3Sfloats(lonF(:),latF(:),F(:),xbin,ybin); 
Nf=Nf/Nfloats; Nf(Nf==0)=nan;
imagescnan(xbin,ybin,Nf); axis xy; caxis([0 1]);
hold on
contour(lonm,latm,zbat,contbat,'color',colbat,'linewidth',wbat);
line(loncoast,latcoast,'color','k','linewidth',wcoast);
plot(lonF(1),latF(1),'.','color',colstart,'markersize',sstart);
axis equal;xlim(lonlim); ylim(latlim); colorbar; colormap(jet(20))
xlabel('Longitude [ºE]'); ylabel('Latitude [ºN]')
title('Density of all tracks');
FormatFig

% %Floats tracks 3d
% figure
% surf(lonm,latm,zbat,'linestyle','none'); hold on
% plot3(lonF(:,1:10),latF(:,1:10),-zF(:,1:10),'color',coltrack,'linewidth',wtrack)
% plot3(lonF(end,:),latF(end,:),-zF(end,:),'*','color',colend,'markersize',send);
% plot3(lonF(1),latF(1),-zF(1),'.','color',colstart,'markersize',sstart);
% xlim(lonlim); ylim(latlim); %axis equal;
% xlabel('Longitude [ºE]'); ylabel('Latitude [ºN]')
% title('Float tracks');
% FormatFig
% 



% %Density of final locations
figure(3)
Nend=hist3([lonF(end,:)',latF(end,:)'],'ctrs',bin);
Nend=Nend/Nfloats; Nend(Nend==0)=nan;
imagescnan(xbin,ybin,Nend'); axis xy; caxis([0 max(Nend(:))]);
hold on
contour(lonm,latm,zbat,contbat,'color',colbat,'linewidth',wbat);
line(loncoast,latcoast,'color','k','linewidth',wcoast)
plot(lonF(1),latF(1),'.','color',colstart,'markersize',sstart);
axis equal;xlim(lonlim); ylim(latlim); colorbar; colormap(jet(20))
xlabel('Longitude [ºE]'); ylabel('Latitude [ºN]')
title('Density of final locations');
FormatFig



cds('/Lagrangian/OutputFiles');
save(file2save,'lonm','latm','zm','lonF','latF','zF','loc0','date0',...
    'Ndays','Nstep','SimuDirection','Nfloats','tstep','KturbXY','KmedXY',...
    'KturbZ','KmedZ','lonlim','latlim','depthlim','ilon0','ilat0','iz0',...
    'xbin','ybin','Nin','Nf','Nend');



