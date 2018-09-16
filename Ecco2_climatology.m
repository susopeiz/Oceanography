%This Program compute the climatological mean of the u,v,w, T, S and SIG 
%fields of ECCO2 model.

FT=input('First time run? ');

%ECCO2 model dataset
dirdata='/Volumes/TOSH1TB/ECCO2_data_CMS/Atlantic/';

cd(dirdata) 
if FT==1
    clear all
    close all  
    
    dirout='/Volumes/Suso/oCeaNoS/Matlab/ECCO2/OutputFiles/';

    
    ivarnam=input('Select variable to compute, (u,v,w,t,s,g): ','s');
    ivariability=input('Select time scale, (1=tot,2=anual,3=seas,4=month): ');
    anos=input('Select years for the climatology, [anoI anoF]: '); %ECCO2 = [1992 2011]
    anoI=anos(1); anoF=anos(2); Lanos=anoF-anoI+1;
    file2save=input('File to be saved? ','s');
    
    
    
    if ivarnam=='g'
        filesT=dir(['*t.nc']); filesT=filesT(2:end); files=filesT;
        filesS=dir(['*s.nc']); filesS=filesS(2:end);
        Ldays=length(filesT);  
    else
        files=dir(['*' ivarnam '.nc']); %files=files(2:end); %because the first day is not used in the experiment
        Ldays=length(files);  
    end
     
    
    seasonRef=[12 1 2; 3 4 5; 6 7 8; 9 10 11];
    
    lonmod=ncread(files(1).name,'Longitude'); Llon=length(lonmod);
    latmod=ncread(files(1).name,'Latitude'); Llat=length(latmod);
    zmod=ncread(files(1).name,'Depth'); Lz=length(zmod);
       
    
    if ivarnam=='u'
        var='uvel';
    elseif ivarnam=='v'
        var='vvel';
    elseif ivarnam=='w'
        var='wvel';
    elseif ivarnam=='t'
        var='temp';
    elseif ivarnam=='s';
        var='salt';
    elseif ivarnam=='g'
        var='sigm';
    end
        
        
        
    
end



km=zeros(12,1);
ks=zeros(4,1);
ka=zeros(Lanos,1);
kt=0;
j=0;
wb=waitbar(0,'Please wait...');
tic
for iday=1:Ldays

    
    ano=str2num(files(iday).name(8:11));
    mes=str2num(files(iday).name(12:13));
    dia=str2num(files(iday).name(14:15));
    [season,messeason]=find(seasonRef==mes);
    
    if ano>=anoI && ano<=anoF
        j=j+1;
        
        if ivarnam=='g'
            mt=ncread(filesT(iday).name,'temp');
            ms=ncread(filesS(iday).name,'salt');
            
            lonmod2=lonmod-360;
            m=TS2sigma(lonmod2,latmod,zmod,mt,ms);
        else
            m=ncread(files(iday).name,var);
        end
        
        %Monthly
        %-------------------------------------
        if ivariability==4
            if j==1
                M=zeros(Llon,Llat,Lz,12);
                M_sd=zeros(Llon,Llat,Lz,12);
            end
            
            
            %counter
            km(mes)=km(mes)+1;
            
            %Standard Deviation (always first than the mean)
            M_sd(:,:,:,mes)=M_sd(:,:,:,mes) + ...
                ((km(mes)-1)/km(mes)) * (m - M(:,:,:,mes)).^2;
            %Mean
            M(:,:,:,mes)=M(:,:,:,mes) + ...
                (m - M(:,:,:,mes)) / km(mes);
        end
        
        
        %------------------------------------------------
        
        %Seasonal
        %-------------------------------------
        if ivariability==3
            if j==1
                M=zeros(Llon,Llat,Lz,4);
                M_sd=zeros(Llon,Llat,Lz,4);
            end
            
            %counter
            ks(season)=ks(season)+1;
            
            %Standard Deviation (always first than the mean)
            M_sd(:,:,:,season)=M_sd(:,:,:,season) + ...
                ((ks(season)-1)/ks(season)) * (m - M(:,:,:,season)).^2;
            %Mean
            M(:,:,:,season)=M(:,:,:,season) + ...
                (m - M(:,:,:,season)) / ks(season);
        end
        %------------------------------------------------
        
        %Interanual
        %-------------------------------------
        if ivariability==2
            if j==1
                M=zeros(Llon,Llat,Lz,Lanos);
                M_sd=zeros(Llon,Llat,Lz,Lanos);
            end
            
            %counter
            a=ano-anoI+1;
            ka(a)=ka(a)+1;
            
            %Standard Deviation (always first than the mean)
            M_sd(:,:,:,a)=M_sd(:,:,:,a) + ...
                ((ka(a)-1)/ka(a)) * (m - M(:,:,:,a)).^2;
            %Mean
            M(:,:,:,a)=M(:,:,:,a) + ...
                (m - M(:,:,:,a)) / ka(a);
        end
        %------------------------------------------------
        
        
        %Total
        %-------------------------------------
        if ivariability==1
            if j==1
                M=zeros(Llon,Llat,Lz);
                M_sd=zeros(Llon,Llat,Lz);
            end
            
            %counter
            kt=kt+1;
            
            %Standard Deviation (always first than the mean)
            M_sd(:,:,:)=M_sd(:,:,:) + ...
                ((kt-1)/kt) * (m - M(:,:,:)).^2;
            %Mean
            M(:,:,:)=M(:,:,:) + ...
                (m - M(:,:,:)) / kt;     
            
        end
        %------------------------------------------------
        
        %Checking nans in fields
        reg(j,:)=[ano mes dia sum(isnan(m(:)))];        
    end
    
    t=toc;
    tr=t*(Ldays/iday-1);
    waitbar(iday/Ldays,wb,['Please wait... ' num2str(tr/60,'%1.1f') 'minutes left.']);
end
close(wb)

if length(unique(reg(:,end)))>1
    disp('NaN problems in the grid!!!')    
else
    disp('Number of nans ok in all files!')    
end
disp( [num2str(unique(reg(:,end)')) ' nans'])


%Final computation of std
K=sum(km(1)+ks(1)+ka(1)+kt);
M_sd=sqrt(M_sd/K);

eval([var '=M;']); clear M
eval([var '_sd=M_sd;']); clear M_sd

cds('/ECCO2/OutputFiles')
save(file2save,var,[var '_sd'],'lonmod','latmod','zmod','reg');

    
    
    
   
   