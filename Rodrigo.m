clear all
close all

fitType='poly2'; % Función de fiteo
alturas_fit=[-1.5 0 1.5]; % Alturas que se usan para el fiteo

load('sujetos.mat')
color=['r','b','g','k','m'];
mech_sizes=[33 40 48 57 65 75 100];
asyn=-30; % Asincronía buscada

for s=1:length(sujeto)
    figure(s) % Una figura por sujeto
    c_data=1;
    color(s)='b';
    %% Grafico las asyn de todos los trials

    for t=1:length(sujeto(s).exp) % En todos los trials
        if sujeto(s).exp(t).out==0 % Me quedo con aquellos que no son outliers
        
        angulo=sujeto(s).exp(t).mech_size; % Ángulo del Servo      
        altura=ang2alt(angulo); % Convierto el ángulo en una altura
        
        % Busco la asyn del bip perturbado
        pert_indx=find(sujeto(s).exp(t).asyn(2,:)==sujeto(s).exp(t).mech_bip+1);
        % Asigno una condicion de acuerdo a la perturbacion
        cond_indx=find(mech_sizes==sujeto(s).exp(t).mech_size);  
        
        % Grafico todos las asyn forzadas por escalón
        y=sujeto(s).exp(t).asyn(1,pert_indx)-sujeto(s).condicion(cond_indx).mean_prebaseline;
        
        subplot(1,2,1) % Datos de todos los trials
        plot(altura,y,'Color',color(s))
        hold on
        subplot(1,2,2) % Datos de todos los trials invertidos
        plot(y,altura,'Color',color(s))
        hold on
        end
    end

%% Grafico las asyn medias por condicion
    for c=1:length(sujeto(s).condicion) % En todas las condiciones

        angulo=sujeto(s).condicion(c).mech_size; % Ángulo del Servo
        altura=ang2alt(angulo);  % Convierto el ángulo en una altura
        % Busco la asyn del bip perturbado
        idx=find(sujeto(s).condicion(c).serie_prom(2,:)==1);
        % Grafico la asyn media forzadas por escalón
        subplot(1,2,1) % Medias
        plot(altura,sujeto(s).condicion(c).serie_prom_bs(1,idx),'*','Color',color(s))
        hold on   
        subplot(1,2,2) % Medias invertidas
        plot(sujeto(s).condicion(c).serie_prom_bs(1,idx),altura,'*','Color',color(s))
        hold on
        
        data(c,1,s)=altura; % x=altura
        data(c,2,s)=sujeto(s).condicion(c).serie_prom_bs(1,idx); % y=Asincronia
        
        if ismember(altura,alturas_fit) % Guardo los x e y de interés
        data_min(c_data,1,s)=altura; % x=altura
        data_min(c_data,2,s)=sujeto(s).condicion(c).serie_prom_bs(1,idx); % y=Asincronia
        c_data=c_data+1;
        end
    end
     
    %% Calibro una cuadrática utilizando solo los puntos -M,0,+M
    [fiteo,gof_d] = fit(data_min(:,1,s),data_min(:,2,s),fitType);
    subplot(1,2,1)
    plot(fiteo) 
    hold on
    
    coef(s,:) = [fiteo.p1 fiteo.p2 fiteo.p3];
    p=coef;
    p(s,3) = coef(s,3)-asyn;
    r = roots(p(s,:));   
    if isreal(r)==1 % Si no tiene parte imaginaria, graficar el menor escalón que satisface
        subplot(1,2,1)
        plot(min(r),asyn)
        hold on
    end    
    
    % Matriz de coeficientes y escalones predichos
    Pred(s,:,1)=[fiteo.p1 fiteo.p2 fiteo.p3 r(1) r(2)];
    
    %% Calibro una cuadrática invierto x e y con -M,0,+M y con todos los datos
    [fiteo2m,gof_d2m] = fit(data_min(:,2,s),data_min(:,1,s),fitType);
    [fiteo2,gof_d2] = fit(data(:,2,s),data(:,1,s),fitType);

    subplot(1,2,2)
    plot(fiteo2m,data_min(:,2,s),data_min(:,1,s))
    hold on
    plot(fiteo2,data(:,2,s),data(:,1,s))
    
    coef2(s,:) = [fiteo2.p1 fiteo2.p2 fiteo2.p3];
    p2=coef2;
    p2(s,3) = coef2(s,3)-asyn;
    r2 = roots(p2(s,:));   
    coef2m(s,:) = [fiteo2m.p1 fiteo2m.p2 fiteo2m.p3];
    p2m=coef2m;
    p2m(s,3) = coef2m(s,3)-asyn;
    r2m = roots(p2m(s,:)); 
    
    Pred(s,:,2)=[fiteo2m.p1 fiteo2m.p2 fiteo2m.p3 r2m(1) r2m(2)];
    Pred(s,:,3)=[fiteo2.p1 fiteo2.p2 fiteo2.p3 r2(1) r2(2)];
    
    if isreal(r2)==1 % Si no tiene parte imaginaria, graficar el menor escalón que satisface
        subplot(1,2,1)
        plot(min(r2),asyn,'o','Color','g')
        hold on
    end
    if isreal(r2m)==1 % Si no tiene parte imaginaria, graficar el menor escalón que satisface
        plot(min(r2m),asyn,'o','Color','k')
        hold on
    end       
    
    
    %% Invierto x e y, (despeje Rodrigo) con -M,0,+M
    
    M=1.5; % Tamaño del escalón
    t_plus=data_min(find(data_min(:,1,s)==M),2,s); % Asincronía de +M
    t_minus=data_min(find(data_min(:,1,s)==-M),2,s); % Asincronía de -M
    a(s)=(2*M)/(t_minus*t_plus);
    v(s)=((t_minus-t_plus)/(t_minus*t_plus))*M;
    
    escalon(s)=0.5*a(s)*(asyn)^2+v(s)*asyn;
    
    coef3(s,:) = [a(s) v(s) escalon(s)];
    p3=coef3;
    p3(s,3) = coef3(s,3)-asyn;
    r3 = roots(p3(s,:));   
     
    
    Pred(s,:,4)=[fiteo2.p1 fiteo2.p2 fiteo2.p3 r3(1) r3(2)];
    
    
    subplot(1,2,1)
    plot(escalon(s),asyn,'+','Color','m')
    ylabel('asyn [ms]')
    xlabel('altura [cm]')
    xlim([-2 2]) % Limite de alturas
    ylim([-150 150])
    title('Fiteo')
    hold on
    subplot(1,2,2)
    plot(asyn,escalon(s),'+','Color','m')
    xlabel('asyn [ms]')
    ylabel('altura [cm]')
    ylim([-2 2]) % Limite de alturas
    xlim([-150 150])
    title('Fiteo invertido')
    hold on
    
end
