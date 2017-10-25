%% Grafica todos los trials por sujeto, para cada condición
clear all
close all

load('sujetos_alturas.mat')
s = input('Numero de sujeto a graficar ');

perts = unique(horzcat( sujeto(s).exp.mech_size ));

for t=1:length(sujeto(s).exp)
    
    for m=1:length(perts)
        if sujeto(s).exp(t).mech_size == perts(m)
            if sujeto(s).exp(t).out == 0
        figure(m)
        title(['Condicion= ' num2str(m)])
        plot(sujeto(s).exp(t).asyn_pro(2,:),sujeto(s).exp(t).asyn_pro(1,:)-sujeto(s).condicion(m).mean_prebaseline)
        hold on
        plot(sujeto(s).condicion(m).serie_prom(2,:),sujeto(s).condicion(m).serie_prom_bs(1,:),'r')
        hold on
            end
        end
    end
end

