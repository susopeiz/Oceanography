clear all
close all


%Este Programa selecciona y guarda los perfiles que cumplen ciertos requisitos
% especificados en el programa. Se trata ahora de ser poco restrictivo para 
% trabajar luego con los perfiles. Las condicones:
%Localizaci�n geogr�fica, posQC, timeQC, Plim, perfiles con T S o O2 data,
%QC medio del perfil de P o de T....VER PROGRAMA

% [SE RECOMIENDA USAR ANTES ProfilePreSelector.m PARA CONOCER
% CUALES Y CUANTOS PERFILES Y PERFILADORES SON BUENOS Y MALOS Y PORQU� LO
% SON Y CON ELLOS AUMENTAR O DISMINUIR LAS REQUISITOS A LA HORA DE 
% SELECCIONAR LOS PERFILES EN ESTE PROGRAMA]

%Estos flotadores han sido previamente seleccionado con SelectorS.m.
%La salida son vectores con las caracteristicas de cada uno de los perfiles,
% todos ellos ya puestos en un mismo vector por lo que ya no ordenados por perfiladores. 

%DIMENSIONES DE LA REGION Y LIMITES VERTICALES DE LOS PERFILES
lonmin=-80; %-80; %Pos este Neg oeste
lonmax=60; %20; %20;
latmin=-50; %-15; %-15; %Pos norte Neg Sur
latmax=30; %30;
PlimUp=100; %en [db]
PlimDown=1800; %en [db]
%Calidad seleccionada de los datos de T y P del perfil: A y B. Ver m�s
%abajo para cambiar

%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREAMOS PERFIL EST�NDAR

%Esto se uso para calcular inicialmente una primera matriz ordenada en
%press
% PrefMin=10; %input('Minimum pressure of the standard profile? ');
% PrefMax=2000; %input('Maximum pressure of the standard profile? ');
% dPref=10; %input('Vertical resolution of the standard profile? ');
% Pref(:,1)=PrefMin:dPref:PrefMax;

%Con eso, calculamos con "Creando_SigmaRef.m" un perfil estandar en sigma0
cd /export/suso/oCeaNoS/Matlab/ARGO/OutputFiles/
load Sigma0Ref.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Nombre del fichero .mat
SaveFileName=input('Name of file to be save? ','s');



%cd S:\oCeaNoS\Matlab\ARGO\SelectedRegion\prof
cd /export/suso/oCeaNoS/Matlab/ARGO/SelectedRegion_50S30N80W60E/
files=dir('*.nc');
Lf=length(files);

isi=0; %Contador de perfiladores seleccionados
ino=0; %contador de perfiladores desechados
kR=0; kDA=0; %Contador de perfilers Real and Delayed/Adjusted modes
kmal_R=0; kmal_DA=0;
Li_R=0; Li_DA=0; %Contador perfiles Real y Delayed/Adjusted
wb=waitbar(0, ['Please wait...0 de ' num2str(Lf)]);
seasons=[1 2 3; 4 5 6; 7 8 9; 10 11 12];

%Cargamos la costa ahora para no pisar nuestras variables
load coast
latcoast=lat;
loncoast=long;
clear lat long

for f=1:Lf %Lf:-1:1 %1:Lf
    waitbar(f/Lf,wb,['Please wait... ' num2str(f) ' de ' num2str(Lf)]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LEEMOS LAS VARIABLES DE TODOS LOS PERFILES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Tomamos todos los perfiles
    infoperf=nc_info(files(f).name);
    varnameperf=strvcat(infoperf.Dataset.Name); %Cada fila un nombre de una variable
    lonperf=nc_varget(files(f).name,'LONGITUDE');
    latperf=nc_varget(files(f).name,'LATITUDE');
    posperfQC=nc_varget(files(f).name,'POSITION_QC'); inodata=find(posperfQC==' '); posperfQC(inodata)='0';
    posperfQC=str2num(posperfQC);
    timeperf=floor(nc_varget(files(f).name,'JULD'));
    timeperfQC=str2num(nc_varget(files(f).name,'JULD_QC'));

    Ls0(f)=length(lonperf); %Numero total de perfiles en el flotador

    modeperf=nc_varget(files(f).name,'DATA_MODE'); %Real, Delayed or Adjusted
    DCrefperf=nc_varget(files(f).name,'DC_REFERENCE'); %Ref del profile
    IDperf=repmat(str2num(files(f).name(1:end-8)),length(modeperf),1);
    PperfQC=nc_varget(files(f).name,'PROFILE_PRES_QC'); %Profile Param_QC
    TperfQC=nc_varget(files(f).name,'PROFILE_TEMP_QC');
    %SperfQC=nc_varget(files(f).name,'PROFILE_PSAL_QC');

    Pperf=nc_varget(files(f).name,'PRES');
    PperfMin=min(Pperf')'; PperfMax=max(Pperf')';
    Pperf_QC=nc_varget(files(f).name,'PRES_QC');
    Pperf_adj=nc_varget(files(f).name,'PRES_ADJUSTED');
    Pperf_adj_QC=nc_varget(files(f).name,'PRES_ADJUSTED_QC');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %BUSCAMOS QU� VARIABLES TIENE PERFILADOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    e=varnameperf'; varnameperf=e(:)'; %Y ahora tenemos todos nombres de variables en una misma fila
    if findstr(varnameperf,'TEMP')
        Tperf=nc_varget(files(f).name,'TEMP');
        Tperfsi=sum(~isnan(Tperf'))'; %numero de no nans por perfil
        Tperf_QC=nc_varget(files(f).name,'TEMP_QC');
        Tperf_adj=nc_varget(files(f).name,'TEMP_ADJUSTED');
        Tperf_adj_QC=nc_varget(files(f).name,'TEMP_ADJUSTED_QC');
    else
        Tperfsi=zeros(Ls0(f),1);
    end
    if findstr(varnameperf,'PSAL')
        Sperf=nc_varget(files(f).name,'PSAL');
        Sperfsi=sum(~isnan(Sperf'))'; %numero de no nans por perfil
        SperfQC=nc_varget(files(f).name,'PROFILE_PSAL_QC'); %QC promedio del perfil (lo metemos aqui por sino hay PSAL
        Sperf_QC=nc_varget(files(f).name,'PSAL_QC');
        Sperf_adj=nc_varget(files(f).name,'PSAL_ADJUSTED');
        Sperf_adj_QC=nc_varget(files(f).name,'PSAL_ADJUSTED_QC');
    else
        Sperfsi=zeros(Ls0(f),1);
        SperfQC=repmat('N',Ls0(f),1);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SELECCIONAMOS LOS PERFILES QUE NOS INTERESAN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%   CAMBIAR CRITERIOS DE SELECCION AKI %%%%%%%%%%%%%
    %Seleccionamos s�lo los perfiles que cumplen las condiciones siguientes
    is_s{f}=find(lonperf<=lonmax & lonperf>=lonmin & latperf>=latmin & latperf<=latmax ... %Limites geogr�ficos
        & (posperfQC==1 | posperfQC==2 )...%| posperfQC==5 | posperfQC==8 | posperfQC==0) ... % | posperfQC==3) ... %Calidad posici�n
        & (timeperfQC==1 | timeperfQC==2 )...%| timeperfQC==5 | timeperfQC==8 | timeperfQC==0) ... % | timeperfQC==3) ... %Calidad tiempo
        & PperfMin<PlimUp & PperfMax>PlimDown ...%Limites verticales
        & Tperfsi>0 ... %Perfiles con valores no nans de Temp
        & Sperfsi>0 ... %Perfiles con valores no nans de Sal  %& O2perfsi>0 ... %Perfiles con valores no nans de O2
        & (PperfQC=='A' )...%| PperfQC=='B' | PperfQC=='C' | PperfQC=='D' | PperfQC=='E') ... %Calidad de media del perfil de P
        & (TperfQC=='A' )...%| TperfQC=='B' | TperfQC=='C' | TperfQC=='D' | TperfQC=='E') ... %Calidad de T
        & (SperfQC=='A' ));%| SperfQC=='B' | SperfQC=='C' | SperfQC=='D' | SperfQC=='E')); %Calidad de S

    %%%%%%%%%%%% INFO %%%%%%%%%%%%%%%%%%%%%%%
    % 0=No QC, 1=Good QC, 2=probably goo, 3=Bad but potentially correctable,
    % 4=Bad data, 5=Value Change, 8=Interpolated value
    % QC = 1,2,5,8 are considered good data in the mean profile QC
    % Mean profile QC ==>> A.=100%good, B.>75%,C.>50%,D.>25%,E.>0%,F.=0%)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    Ls(f)=length(is_s{f}); %N�mero de perfiles buenos para el perfilador f

    if Ls(f)==0
        ino=ino+1; % N�mero de perfiladores sin perfiles buenos
    else

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TRABAJAMOS CON LOS PERFILADORES QUE TIENEN PERFILES BUENOS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Perfiladores con al menos un perfil bueno!!!
        isi=isi+1; %Numero de perfiladores con al menos un perfil bueno
        filename_si{isi}=files(f).name; %Nombre del fichero del perfilador con al menos un perfil bueno
        iperf_si(isi)=is_s(f); %Indices de los perfiles buenos en el perfilador bueno
        
        %Cada columna un perfil
        Pperf=Pperf'; Pperf_QC=Pperf_QC'; Tperf=Tperf'; Tperf_QC=Tperf_QC'; Sperf=Sperf'; Sperf_QC=Sperf_QC';
        %Seleccionamos los perfiles
        Pperf_s=Pperf(:,is_s{f}); Pperf_QC_s=Pperf_QC(:,is_s{f});
        Tperf_s=Tperf(:,is_s{f}); Tperf_QC_s=Tperf_QC(:,is_s{f});
        Sperf_s=Sperf(:,is_s{f}); Sperf_QC_s=Sperf_QC(:,is_s{f}); 
        
        %Quitamos los valores con P<0
        iPSmal=find(Pperf_s<0 | Sperf_s>100 | Tperf_s<=0);
        Pperf_s(iPSmal)=nan;
        Tperf_s(iPSmal)=nan;
        Sperf_s(iPSmal)=nan;

        
        % Trabajamos con Salinidad Absoluta, Temperatura Conservativa y
        % Densidad potencial (TEOS-10)
        Sperf2_s=gsw_SA_from_SP(Sperf_s,Pperf_s,lonperf(is_s{f}),latperf(is_s{f})); %Ahora Absolute Salinity
        Tperf2_s=gsw_CT_from_t(Sperf2_s,Tperf_s,Pperf_s); %Ahora Conservative Temperature
        SIGperf2_s=gsw_sigma0(Sperf2_s,Tperf2_s);  % Densidad potencial con nivel ref z=0
        
        
                    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%Guardamos los perfiles de P, T y S interpolados a un perfil estandar%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Separamos entre Real y Delayed-Adjusted modes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        modeperf_s=modeperf(is_s{f});
        
        %Real Time
        iR=find(modeperf_s=='R'); %Dentro de los seleccionados, los reales
        if ~isempty(iR)
            Pperf_s_R=Pperf_s(:,iR); 
            Pperf_QC_s_R=Pperf_QC_s(:,iR); Pperf_QC_s_R(Pperf_QC_s_R==' ')='7';
            Tperf_s_R=Tperf2_s(:,iR); 
            Tperf_QC_s_R=Tperf_QC_s(:,iR); Tperf_QC_s_R(Tperf_QC_s_R==' ')='7';
            Sperf_s_R=Sperf2_s(:,iR); 
            Sperf_QC_s_R=Sperf_QC_s(:,iR); Sperf_QC_s_R(Sperf_QC_s_R==' ')='7';  
            SIGperf_s_R=SIGperf2_s(:,iR);
            
            kR2=0;
            for ip=1:length(iR)
                Pperf_QC_s_R_n(:,ip)=str2num(Pperf_QC_s_R(:,ip)); Pperf_QC_s_R_n(Pperf_QC_s_R_n(:,ip)==7)=nan;
                Tperf_QC_s_R_n(:,ip)=str2num(Tperf_QC_s_R(:,ip)); Tperf_QC_s_R_n(Tperf_QC_s_R_n(:,ip)==7)=nan;
                Sperf_QC_s_R_n(:,ip)=str2num(Sperf_QC_s_R(:,ip)); Sperf_QC_s_R_n(Sperf_QC_s_R_n(:,ip)==7)=nan;
                iPTS=find(~isnan(Pperf_s_R(:,ip)) & ~isnan(Tperf_s_R(:,ip)) & ~isnan(Sperf_s_R(:,ip)) ... %niveles no todos buenos P y T y S
                    & (Pperf_QC_s_R_n(:,ip)==1 | Pperf_QC_s_R_n(:,ip)==2 | Pperf_QC_s_R_n(:,ip)==5 | Pperf_QC_s_R_n(:,ip)==8) ... %buen P QC
                    & (Tperf_QC_s_R_n(:,ip)==1 | Tperf_QC_s_R_n(:,ip)==2 | Tperf_QC_s_R_n(:,ip)==5 | Tperf_QC_s_R_n(:,ip)==8) ... %Buen T QC
                    & (Sperf_QC_s_R_n(:,ip)==1 | Sperf_QC_s_R_n(:,ip)==2 | Sperf_QC_s_R_n(:,ip)==5 | Sperf_QC_s_R_n(:,ip)==8)); %buen S QC
                %Interpolamos s�lo pares de puntos nonans y con buen QC
                if length(iPTS)<2
                    kmal_R=kmal_R+1;
                else
                    %Comprobamos que no hay malos datos
                    if nansum(nansum([(Tperf_s_R(iPTS,ip)>35) , (Tperf_s_R(iPTS,ip)<0) , (Sperf_s_R(iPTS,ip)<0) , (Sperf_s_R(iPTS,ip)>42) ],2) )
                        %3
                    end
                    
                    kR=kR+1;
                    kR2=kR2+1;
                    %Ordenamos SIGperf para que no haya lios en la
                    %interpolacion
                    [SIGsort,iSIGsort]=sort(SIGperf_s_R(iPTS,ip));
                    %Ordenamos de la misma forma Tperf y Sperf
                    Tsort=Tperf_s_R(iPTS,ip); Ssort=Sperf_s_R(iPTS,ip);
                    Tsort=Tsort(iSIGsort); Ssort=Ssort(iSIGsort);
                    %Interpolamos
                    T_R(:,kR)=interp1q(SIGsort,Tsort,SIGref);                    
                    S_R(:,kR)=interp1q(SIGsort,Ssort,SIGref);
                    %T_QC_R(:,kR)=interp1q(SIGperf_s_R(iPTS,ip),Tperf_QC_s_R_n(iPTS,ip),SIGref);
                    %S_QC_R(:,kR)=interp1q(SIGperf_s_R(iPTS,ip),Sperf_QC_s_R_n(iPTS,ip),SIGref);
                    %Indice perfiles
                    ip_s_R{f}(kR2)=ip;
                end
            end
            if kR2==0 %La data de ningun perfil cumple los requisitos. O sea A,B,C... est� mal.
                ip_s_R{f}=[];
                i_R{f}=[];
            else
                Lip_s_R=length(ip_s_R{f});
                i_R{f}=is_s{f}(iR(ip_s_R{f}));
                Li_R=Li_R+Lip_s_R;
                %i_R(Li_R-Lip_s_R+1:Li_R)=is_s{f}(iR(ip_s_R));
            end
        else
            i_R{f}=[];
        end

        %Delayed time and Adjusted
        iDA=find(modeperf_s=='D' | modeperf_s=='A');
        if ~isempty(iDA)
            Pperf_s_DA=Pperf_s(:,iDA); 
            Pperf_QC_s_DA=Pperf_QC_s(:,iDA); Pperf_QC_s_DA(Pperf_QC_s_DA==' ')='7';
            Tperf_s_DA=Tperf2_s(:,iDA); 
            Tperf_QC_s_DA=Tperf_QC_s(:,iDA); Tperf_QC_s_DA(Tperf_QC_s_DA==' ')='7';
            Sperf_s_DA=Sperf2_s(:,iDA); 
            Sperf_QC_s_DA=Sperf_QC_s(:,iDA); Sperf_QC_s_DA(Sperf_QC_s_DA==' ')='7';  
            SIGperf_s_DA=SIGperf2_s(:,iDA);
            
            kDA2=0;
            for ip=1:length(iDA)                                
                Pperf_QC_s_DA_n(:,ip)=str2num(Pperf_QC_s_DA(:,ip)); Pperf_QC_s_DA_n(Pperf_QC_s_DA_n(:,ip)==7,ip)=nan;
                Tperf_QC_s_DA_n(:,ip)=str2num(Tperf_QC_s_DA(:,ip)); Tperf_QC_s_DA_n(Tperf_QC_s_DA_n(:,ip)==7,ip)=nan;
                Sperf_QC_s_DA_n(:,ip)=str2num(Sperf_QC_s_DA(:,ip)); Sperf_QC_s_DA_n(Sperf_QC_s_DA_n(:,ip)==7,ip)=nan;
                iPTS=find(~isnan(Pperf_s_DA(:,ip)) & ~isnan(Tperf_s_DA(:,ip)) & ~isnan(Sperf_s_DA(:,ip)) ... %niveles con: no nans en P y T
                    & (Pperf_QC_s_DA_n(:,ip)==1 | Pperf_QC_s_DA_n(:,ip)==2 | Pperf_QC_s_DA_n(:,ip)==5 | Pperf_QC_s_DA_n(:,ip)==8) ... %buen P QC
                    & (Tperf_QC_s_DA_n(:,ip)==1 | Tperf_QC_s_DA_n(:,ip)==2 | Tperf_QC_s_DA_n(:,ip)==5 | Tperf_QC_s_DA_n(:,ip)==8) ... %Buen T QC
                    & (Sperf_QC_s_DA_n(:,ip)==1 | Sperf_QC_s_DA_n(:,ip)==2 | Sperf_QC_s_DA_n(:,ip)==5 | Sperf_QC_s_DA_n(:,ip)==8)); %buen S QC
                %Interpolamos s�lo pares de puntos nonans y con buen QC
                if length(iPTS)<2
                    kmal_DA=kmal_DA+1;
                else
                    %Comprobamos que no hay malos datos
                    if nansum(nansum([(Tperf_s_DA(iPTS,ip)>35) , (Tperf_s_DA(iPTS,ip)<0) , (Sperf_s_DA(iPTS,ip)<0) , (Sperf_s_DA(iPTS,ip)>42) ],2) )
                        [Tperf_s(iPTS,ip) str2num(Tperf_QC_s(iPTS,ip)) Tperf2_s(iPTS,ip) Sperf_s(iPTS,ip) str2num(Sperf_QC_s(iPTS,ip)) Sperf2_s(iPTS,ip)]
                        ip
                    end
                    
                    kDA=kDA+1;
                    kDA2=kDA2+1;
                    %Ordenamos SIGperf para que no haya lios en la
                    %interpolacion
                    [SIGsort,iSIGsort]=sort(SIGperf_s_DA(iPTS,ip));
                    %Ordenamos de la misma forma Tperf y Sperf
                    Tsort=Tperf_s_DA(iPTS,ip); Ssort=Sperf_s_DA(iPTS,ip);
                    Tsort=Tsort(iSIGsort); Ssort=Ssort(iSIGsort);
                    %Interpolamos
                    T_DA(:,kDA)=interp1q(SIGsort,Tsort,SIGref);                    
                    S_DA(:,kDA)=interp1q(SIGsort,Ssort,SIGref);
                    %T_QC_DA(:,kDA)=interp1q(SIGperf_s_DA(iPTS,ip),Tperf_QC_s_DA_n(iPTS,ip),SIGref);
                    %S_QC_DA(:,kDA)=interp1q(SIGperf_s_DA(iPTS,ip),Sperf_QC_s_DA_n(iPTS,ip),SIGref);
                    %Indice perfiles
                    ip_s_DA{f}(kDA2)=ip;
                    
                    if nansum(nansum([(T_DA>35) , (T_DA<0) , (S_DA<0) , (S_DA>42) ],2) )
                        [Tperf_s(iPTS,ip) str2num(Tperf_QC_s(iPTS,ip)) Tperf2_s(iPTS,ip) Sperf_s(iPTS,ip) str2num(Sperf_QC_s(iPTS,ip)) Sperf2_s(iPTS,ip) SIGperf2_s(iPTS,ip)]
                        [T_DA(:,kDA) S_DA(:,kDA) SIGref]
                        ip
                    end
                end
            end
            if kDA2==0 %La data de ningun perfil cumple los requisitos. O sea A,B,C... est� mal.
                ip_s_DA{f}=[];
                i_DA{f}=[];
            else
                Lip_s_DA=length(ip_s_DA{f});
                i_DA{f}=is_s{f}(iDA(ip_s_DA{f}));
                Li_DA=Li_DA+Lip_s_DA;
                %i_DA(Li_DA-Lip_s_DA+1:Li_DA)=is_s{f}(iDA(ip_s_DA));
            end
        else
            i_DA{f}=[];
        end
        disp(['kmalR/LiR ' num2str(100*kmal_R/sum(Ls),4) '%  kmalDA/LiDA ' num2str(100*kmal_DA/sum(Ls),4) '%      '...
            num2str([min(min(T_DA)) max(max(T_DA)) min(min(S_DA)) max(max(S_DA))]) ])
        clear Pperf_QC_s_R_n Pperf_QC_s_DA_n Tperf_QC_s_R_n Tperf_QC_s_DA_n Sperf_QC_s_R_n Sperf_QC_s_DA_n



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TRABAJAMOS YA SOLO CON LOS PERFILES BUENOS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Usamos ya s�lo estos perfiles. Generamos vectores de 1dim con info
        % de cada uno de los perfiles seleccionados

        %tiempo de referencia
        timeperf0=nc_varget(files(f).name,'REFERENCE_DATE_TIME'); timeperf0=timeperf0';
        anoperf0=(timeperf0(1:4));
        mesperf0=(timeperf0(5:6));
        diaperf0=(timeperf0(7:8));
        timeperf_s=timeperf(is_s{f});
        clear mesperf anoperf seasonperf
        for s=1:Ls(f)
            time=datestr(addtodate(datenum([mesperf0 '/' diaperf0 '/' anoperf0 ' 00:00']),floor(timeperf_s(s)),'day'),2);
            mesperf(s,1)=str2num(time(1:2));
            anoperf(s,1)=2000+str2num(time(7:8));
            [seasonperf(s,1),r]=find(seasons==mesperf(s));
        end

        for jmode=1:2
            if (jmode==1 && ~isempty(i_DA{f})) || (jmode==2 && ~isempty(i_R{f}))
                if jmode==1 && max(i_DA{f})
                    j=i_DA{f}; %perfiles seleccionados definitivamente de is_s{f}
                    jt=iDA(ip_s_DA{f}); %perfiles seleccionados (tiempo)
                    ij=Li_DA-Lip_s_DA+1:Li_DA;
                elseif jmode==2 && max(i_R{f})
                    j=i_R{f}; %perfiles seleccionados definitivamente de is_s{f}
                    jt=iR(ip_s_R{f}); %perfiles seleccionados (tiempo)
                    ij=Li_R-Lip_s_R+1:Li_R;
                end
                lon{jmode}(ij)=lonperf(j);
                lat{jmode}(ij)=latperf(j);
                posQC{jmode}(ij)=posperfQC(j);
                mes{jmode}(ij)=mesperf(jt);
                ano{jmode}(ij)=anoperf(jt);
                season{jmode}(ij)=seasonperf(jt);
                timeQC{jmode}(ij)=timeperfQC(j);
                mode{jmode}(ij)=modeperf(j);
                DCref{jmode}(ij,:)=DCrefperf(j,:);
                ID{jmode}(ij)=IDperf(j);
                PQC{jmode}(ij)=PperfQC(j);
                TQC{jmode}(ij)=TperfQC(j);
                SQC{jmode}(ij)=SperfQC(j);
                PMin{jmode}(ij)=PperfMin(j);
                PMax{jmode}(ij)=PperfMax(j);
            end
        end
    end
end
Lsi=isi; %numero total de perfiladores buenos

close(wb)

genprogname=mfilename;
daterunned=date;

cd /export/suso/oCeaNoS/Matlab/ARGO/OutputFiles
save(SaveFileName) %TropAtlanProf_20N15S70W20E_TS.mat
