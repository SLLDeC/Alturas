clear all
close all
load('sujetos.mat')

%% Promedia las series temporales entre sujetos

for s=1:length(sujeto)
    
    for c=1:length(sujeto(s).condicion);   
        
        idx=find(sujeto(s).condicion(c).serie_prom(2,:)==-7);
     
     Cond_Mat(s,:,c)=sujeto(s).condicion(c).serie_prom_bs(1,idx:end);
       
    end
    
end

Cond_Prom=mean(Cond_Mat);
Cond_Std=std(Cond_Mat);

for c=1:length(sujeto(s).condicion);

    condicion(c).serie_prom(1,:)=Cond_Prom(:,:,c);
    condicion(c).serie_prom(2,:)=[-7:7];
    condicion(c).serie_prom_std(1,:)=Cond_Std(:,:,c);
    condicion(c).serie_prom_std(2,:)=[-7:7];
    condicion(c).n=length(sujeto);
    
end

%% Grafica las series temporales promedio por condicion
figure(1)
for i=1:7
    plot(condicion(i).serie_prom(2,:),condicion(i).serie_prom(1,:))
    hold all
end
legend('-1.5','-1.0','-0.5','0','+0.5','+1.0','+1.5')

save('promedios.mat','condicion')