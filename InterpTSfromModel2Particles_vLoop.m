
%Este programa interpola los campos de temperatura y salinidad a todas las
%posiciones de todas las particulas.Las salidas son vectores T y S de la
%misma forma que lat, lon, depth en la salida de CMS.

%---------------------------------------------
% CMS POSTPROCESSING....
% 0. CMS2mat.m
% 1. FilterRipParticles.m
% 2. InterpTSfromModel2Particles_v2.m
% 3. FillingNans.m
% 4. TS2TcSaSIG.m
%-----------------------------------------------


FT=input('First time? ');
if FT==1
    
    clear all
    close all
    
    %cd /srv/ccrc/data15/z3439946/Matlab/MOW/
    cd /srv/ccrc/data15/z3439946/Matlab/CMS/OutputFiles/
    filetraj=input('File to be loaded? ','s');
    %'traj_OMZ.mat';%'trajs_ecco_alive.mat';%cropped_alive.mat';
    file2save=input('File to be saved?: ','s');
    FB=input('Forward or Backward? [f/b] ','s');
    dateLOOP=input('Select initial and final loop dates, [anoI mesI diaI; anoF mesF diaF]: ');
    XYZlim=0; %input('Select geographical subregion big enough to include all location, [lonlim,latlim,zlim] or ''0=auto'' : ');
    XYZlim_model=[-98 21; -50 74; 5 5900]; %Total region of the model
    disp('Loading data...')
    load(filetraj); %,'date','lat','lon','depth','exitcode')
    disp('Data loaded.')
    disp('')
else
    file2save=input('File to be saved?: ','s');
    FB=input('Forward or Backward? [f/b] ','s');
    dateLOOP=input('Select initial and final loop dates, [anoI mesI diaI; anoF mesF diaF]: '); %I always before F, even in backward simulation    
    XYZlim=0; %input('Select geographical subregion big enough to include all location, [lonlim,latlim,zlim] or ''0=auto'' : ');
    XYZlim_model=[-98 21; -50 74; 5 5900]; %Total region of the model
end


disp('Sorting data...')
if FB=='f'
    [date_s,i_s]=sort(date,'ascend'); 
elseif FB=='b'
    [date_s,i_s]=sort(date,'descend'); 
end
lat_s=lat(i_s);
lon_s=lon(i_s);
z_s=z(i_s);
ntrac_s=ntrac(i_s);
%clear date lat lon depth
days(:,1)=unique(date_s);
if FB=='b'
    days(:,1)=flipud(days);
end
Ldays=length(days);


cd /srv/ccrc/data15/z3396816/RawData/ECCO2/Cube92_NAtl/perday/raw/
 
filesT=dir('*t.nc'); 
Lt=length(filesT);
filesS=dir('*s.nc');
Ls=length(filesS);

%-------------------------------------
%Dates sorting and grouping.....
if Lt~=Ls
    disp('Error, missing Temp o Salt files!')
end

if Lt<Ldays
    disp('Less days of model than lagrangian tracks???!!!');
end


%Preparing loop
%-------------------------------
dateLi=datenum(dateLOOP(1,:));
dateLf=datenum(dateLOOP(2,:));
anoLi=dateLOOP(1,1);
anoLf=dateLOOP(2,1);
danoL=anoLf-anoLi+1;

% for im=1:Lt
%     ano_m(im)=str2num(filesT(im).name(8:11));
%     mes_m(im)=str2num(filesT(im).name(12:13));
%     dia_m(im)=str2num(filesT(im).name(14:15));
%     date_m(im)=datenum([ano_m(im) mes_m(im) dia_m(im) 0 0 0]);
% end
% 
% idateLi=find(date_m==dateLi);
% idateLf=find(date_m==dateLf);
% 
% Dloop=date_m(idateLi:idateLf);
% if FB=='b' %in descending order
%     Dloop=fliplr(Dloop);
% end
%---------------------------------------
    


%----------------------------------
%Region limits
if XYZlim==0
    XYZlim=[min(lon)-1 max(lon)+1; min(lat)-1 max(lat)+1; min(z) max(z)];
    iout=(XYZlim>XYZlim_model | XYZlim<XYZlim_model);
    if ~isempty(iout)
        XYZlim(iout)=XYZlim_model(iout);
    end
end
%Reading the subregion grid coordinates
lonm=ncread(filesT(1).name,'Longitude')-360;
latm=ncread(filesT(1).name,'Latitude');
zm=ncread(filesT(1).name,'Depth');

ilonlim=[findprox(lonm,XYZlim(1,1)) findprox(lonm,XYZlim(1,2))];
ilatlim=[findprox(latm,XYZlim(2,1)) findprox(latm,XYZlim(2,2))];
izlim=[findprox(zm,XYZlim(3,1)) findprox(zm,XYZlim(3,2))];

lonm1=lonm(ilonlim(1):ilonlim(2));
latm1=latm(ilatlim(1):ilatlim(2));
zm1=zm(izlim(1):izlim(2));
%------------------------------------------
disp('Done.')
disp('')

T_s=nan(length(lat_s),1); S_s=nan(length(lat_s),1);
ntrac_s2=ones(length(lat_s),1);
jnan=0;
hoy=days(1);
wb=waitbar(0,'Please wait...');
kday=0;
disp('Interpolating...')
for iday=1:Ldays
    
    hoyvec=datevec(hoy);
    ano_hoy=num2str(hoyvec(1),'%4.0f');
    mes_hoy=num2str(hoyvec(2),'%02.0f');
    dia_hoy=num2str(hoyvec(3),'%02.0f');
    
    
    
    
    %Loop
    %--------------------------
    ano_hoy_num=str2num(ano_hoy);
    if FB=='f'
        dano_hoy_loop=ano_hoy_num-anoLf;
        nLoop=fix(dano_hoy_loop/danoL);    
        ano_hoy_loop=ano_hoy_num-nLoop*danoL;
    elseif FB=='b'
        dano_hoy_loop=anoLi-ano_hoy_num;
        nLoop=fix(dano_hoy_loop/danoL);    
        ano_hoy_loop=ano_hoy_num+nLoop*danoL;
    end
    ano_hoy_loop=num2str(ano_hoy_loop);
    
    
    
    filename=['nest_1_' ano_hoy_loop mes_hoy dia_hoy '000000'];
    %-------------------------
    filenameT=[filename 't.nc'];
    filenameS=[filename 's.nc'];  
   
    
    jday=find(date_s==hoy);
    
    lat_phoy=lat_s(jday);
    lon_phoy=lon_s(jday);
    z_phoy=z_s(jday);
    ntrac_phoy=ntrac_s(jday);
    
    Lkday=length(lat_phoy);
    
    %if nLoop
    
    if isempty(dir(filenameT)) || isempty(dir(filenameS))
        
        kday=kday+Lkday;
        
        T_s(kday-Lkday+1:kday)=nan;
        S_s(kday-Lkday+1:kday)=nan; 
        ntrac_s2(kday-Lkday+1:kday)=ntrac_phoy;
        
        disp(['No file found for this date: ' num2str(hoyvec(1:3))])
        disp(['NaN''s will be use instead.'])
        
    else 
    
    
        
        Tgridhoy=ncread(filenameT,'temp', ...%);
            [ilonlim(1) ilatlim(1) izlim(1) 1],diff([ilonlim' ilatlim' izlim' [1 1]'])+1);
        Sgridhoy=ncread(filenameS,'salt', ...%);
            [ilonlim(1) ilatlim(1) izlim(1) 1],diff([ilonlim' ilatlim' izlim' [1 1]'])+1);
           
        T_phoy=interp3(latm1,lonm1,zm1,Tgridhoy,lat_phoy,lon_phoy,z_phoy);
        S_phoy=interp3(latm1,lonm1,zm1,Sgridhoy,lat_phoy,lon_phoy,z_phoy);  
        

        
        kday=kday+Lkday;
        T_s(kday-Lkday+1:kday)=single(T_phoy);
        S_s(kday-Lkday+1:kday)=single(S_phoy); 
        ntrac_s2(kday-Lkday+1:kday)=ntrac_phoy;
        
    end
    
    %end
    
    waitbar(iday/Ldays,wb,['Please wait... ' datestr(hoy)]);
    
    if iday<Ldays
        hoy=date_s(jday(end)+1);     
    end      
    
end
close(wb);
disp('T and S field interpolated to all particles locations.')
disp('')

[i_s_s,ii_s]=sort(i_s);


ntrac2=ntrac_s2(ii_s);
disp(['This number should be zero... ' num2str(sum(ntrac2-single(ntrac)))])
disp('')

T=single(T_s(ii_s)); clear T_s
S=single(S_s(ii_s)); clear S_s


disp('Saving...')
%cd /srv/ccrc/data15/z3439946/Matlab/MOW/
cd /srv/ccrc/data15/z3439946/Matlab/CMS/OutputFiles/
save(file2save,'T','S','date','lat','lon','z','ntrac','transp','starts','ends','-v7.3');
disp('Done! ;)')

