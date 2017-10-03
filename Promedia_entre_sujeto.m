clear all
close all
load('sujetos_alturas.mat')

%% Promedia las series temporales entre sujetos

s=length(sujeto);
c=length(sujeto(s).condicion);
max_trial=[-12:10]; % coordenadas del trial
Cond_Mat=nan(s,length(max_trial),c); % Matriz de series promedios de condición de cada sujeto

for s=1:length(sujeto)
    
    for c=1:length(sujeto(s).condicion);   
        
        idx=sujeto(s).condicion(c).serie_prom(2,1)+13;
     
     Cond_Mat(s,idx:length(sujeto(s).condicion(c).serie_prom_bs)+idx-1,c)=sujeto(s).condicion(c).serie_prom_bs;
       
    end
    
end

Cond_Prom=nanmean(Cond_Mat);
Cond_Std=nanstd(Cond_Mat);

for c=1:length(sujeto(s).condicion);

    condicion(c).serie_prom(1,:)=Cond_Prom(:,:,c);
%     condicion(c).serie_prom(2,:)=[-7:7];
 condicion(c).serie_prom(2,:)=max_trial;
    condicion(c).serie_prom_std(1,:)=Cond_Std(:,:,c);
     condicion(c).serie_prom_std(2,:)=max_trial;
%     condicion(c).serie_prom_std(2,:)=[-7:7];
    condicion(c).n=length(sujeto);
    
end

%% Grafica las series temporales promedio por condicion
figure(1)
for i=1:7
    plot(condicion(i).serie_prom(2,:),condicion(i).serie_prom(1,:))
%         mseb(condicion(i).serie_prom(2,:),condicion(i).serie_prom(1,:),condicion(i).serie_prom_std(1,:))

    hold all
    xlim([-6 8])
end
legend('-1.5','-1.0','-0.5','0','+0.5','+1.0','+1.5')

% xlswrite('cond_1.xls',Cond_Mat(:,:,1))
save('promedios.mat','condicion')