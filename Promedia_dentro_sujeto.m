clear all
close all
load('sujetos.mat')

%% Par�metros----------------------------------------------
% s=2;
pre_baseline_bips = 5;        % N° de bips usados para calcular el prebaseline
pos_baseline_bips = 5;        % N° de bips usados para calcular el posbaseline
out_limit = 150;              % Asincronia maxima permitida
%%----------------------------------------------------------


for s=1:length(sujeto)
    
    N=length(sujeto(s).exp);
A=[];B=[];C=[];D=[];F=[];E=[];G=[];c1=1;c2=1;c3=1;c4=1;c5=1;c6=1;c7=1;

% Arma los vectores para cada trial
for m=1:N
    mech_sizes(m)=sujeto(s).exp(m).mech_size; % vector de pert. mech
    pert_bip(m)=sujeto(s).exp(m).mech_bip;
    pert_indx(m)=find(sujeto(s).exp(m).asyn(2,:)==pert_bip(m));
    min_bip(m,sujeto(s).exp(m).mech_size)=min(sujeto(s).exp(m).asyn(2,:)-pert_bip(m));
    max_bip(m,mech_sizes(m))=max(sujeto(s).exp(m).asyn(2,:)-pert_bip(m));
    pre_baseline(m)=mean(sujeto(s).exp(m).asyn(1,pert_indx(m)-pre_baseline_bips-1:pert_indx(m)-1));
    pos_baseline(m)=mean(sujeto(s).exp(m).asyn(1,end-pos_baseline_bips:end));    
    
    if isempty(find(abs(sujeto(s).exp(m).asyn(1,:)-pre_baseline(m))>=out_limit))==1
        out(m)=0;
    else
        out(m)=1;
    end
end

steps=unique(mech_sizes);

for n=1:length(steps)
    min_bip_cond(n)=max(nonzeros(min_bip(:,steps(n))));
    max_bip_cond(n)=min(nonzeros(max_bip(:,steps(n))));
end

for t=1:N
    
    for i=1:length(sujeto(s).exp(t).resp)-1
        resp_n=sujeto(s).exp(t).resp(1,i);
        resp_n1=sujeto(s).exp(t).resp(1,i+1);
        if i<pert_bip(t)
            preITI(i)=resp_n1-resp_n;
        elseif i>=pert_bip(t)
            posITI(i)=resp_n1-resp_n;
        end
    end
    
    meanpreITI(t)=mean(preITI);
    meanposITI(t)=mean(posITI);
    
    if sujeto(s).exp(t).mech_size==steps(1)
        lim_inf_indx=find(sujeto(s).exp(t).asyn(2,:)-pert_bip(t)==min_bip_cond(1));
        lim_sup_indx=find(sujeto(s).exp(t).asyn(2,:)-pert_bip(t)==max_bip_cond(1));
        sujeto(s).exp(t).asyn_pro(1,:)=sujeto(s).exp(t).asyn(1,lim_inf_indx:lim_sup_indx);
        sujeto(s).exp(t).asyn_pro(2,:)=sujeto(s).exp(t).asyn(2,lim_inf_indx:lim_sup_indx)-pert_bip(t);
        sujeto(s).condicion(1).serie_prom(2,:)=sujeto(s).exp(t).asyn_pro(2,:);
        sujeto(s).condicion(1).mech_size=steps(1);
          if out(t)==0
            A(c1,:,1)=sujeto(s).exp(t).asyn_pro(1,:);
            A(c1,:,2)=sujeto(s).exp(t).asyn_pro(1,:)-pre_baseline(t);
            A_prebaseline(c1)=pre_baseline(t);
            A_posbaseline(c1)=pos_baseline(t);
            A_preITI(c1)=meanpreITI(t);
            A_posITI(c1)=meanposITI(t);
            c1=c1+1;

        else
        end
                    
        
    elseif sujeto(s).exp(t).mech_size==steps(2)
        lim_inf_indx=find(sujeto(s).exp(t).asyn(2,:)-pert_bip(t)==min_bip_cond(2));
        lim_sup_indx=find(sujeto(s).exp(t).asyn(2,:)-pert_bip(t)==max_bip_cond(2));
        sujeto(s).exp(t).asyn_pro(1,:)=sujeto(s).exp(t).asyn(1,lim_inf_indx:lim_sup_indx);
        sujeto(s).exp(t).asyn_pro(2,:)=sujeto(s).exp(t).asyn(2,lim_inf_indx:lim_sup_indx)-pert_bip(t);
        sujeto(s).condicion(2).serie_prom(2,:)=sujeto(s).exp(t).asyn_pro(2,:);
        sujeto(s).condicion(2).mech_size=steps(2);
        if out(t)==0
            B(c2,:,1)=sujeto(s).exp(t).asyn_pro(1,:);
            B(c2,:,2)=sujeto(s).exp(t).asyn_pro(1,:)-pre_baseline(t);
            B_prebaseline(c2)=pre_baseline(t);
            B_posbaseline(c2)=pos_baseline(t);
            B_preITI(c2)=meanpreITI(t);
            B_posITI(c2)=meanposITI(t);
            c2=c2+1;
        else
        end
        
    elseif sujeto(s).exp(t).mech_size==steps(3)
        lim_inf_indx=find(sujeto(s).exp(t).asyn(2,:)-pert_bip(t)==min_bip_cond(3));
        lim_sup_indx=find(sujeto(s).exp(t).asyn(2,:)-pert_bip(t)==max_bip_cond(3));
        sujeto(s).exp(t).asyn_pro(1,:)=sujeto(s).exp(t).asyn(1,lim_inf_indx:lim_sup_indx);
        sujeto(s).exp(t).asyn_pro(2,:)=sujeto(s).exp(t).asyn(2,lim_inf_indx:lim_sup_indx)-pert_bip(t);
        sujeto(s).condicion(3).serie_prom(2,:)=sujeto(s).exp(t).asyn_pro(2,:);
        sujeto(s).condicion(3).mech_size=steps(3);
        if out(t)==0
            C(c3,:,1)=sujeto(s).exp(t).asyn_pro(1,:);
            C(c3,:,2)=sujeto(s).exp(t).asyn_pro(1,:)-pre_baseline(t);
            C_prebaseline(c3)=pre_baseline(t);
            C_posbaseline(c3)=pos_baseline(t);
            C_preITI(c3)=meanpreITI(t);
            C_posITI(c3)=meanposITI(t);
            c3=c3+1;
        else
        end
        
    elseif sujeto(s).exp(t).mech_size==steps(4)
        lim_inf_indx=find(sujeto(s).exp(t).asyn(2,:)-pert_bip(t)==min_bip_cond(4));
        lim_sup_indx=find(sujeto(s).exp(t).asyn(2,:)-pert_bip(t)==max_bip_cond(4));
        sujeto(s).exp(t).asyn_pro(1,:)=sujeto(s).exp(t).asyn(1,lim_inf_indx:lim_sup_indx);
        sujeto(s).exp(t).asyn_pro(2,:)=sujeto(s).exp(t).asyn(2,lim_inf_indx:lim_sup_indx)-pert_bip(t);
        sujeto(s).condicion(4).serie_prom(2,:)=sujeto(s).exp(t).asyn_pro(2,:);
          sujeto(s).condicion(4).mech_size=steps(4);
        if out(t)==0
            D(c4,:,1)=sujeto(s).exp(t).asyn_pro(1,:);
            D(c4,:,2)=sujeto(s).exp(t).asyn_pro(1,:)-pre_baseline(t);
            D_prebaseline(c4)=pre_baseline(t);
            D_posbaseline(c4)=pos_baseline(t);
            D_preITI(c4)=meanpreITI(t);
            D_posITI(c4)=meanposITI(t);
            c4=c4+1;
        else
        end
        
    elseif sujeto(s).exp(t).mech_size==steps(5)
        lim_inf_indx=find(sujeto(s).exp(t).asyn(2,:)-pert_bip(t)==min_bip_cond(5));
        lim_sup_indx=find(sujeto(s).exp(t).asyn(2,:)-pert_bip(t)==max_bip_cond(5));
        sujeto(s).exp(t).asyn_pro(1,:)=sujeto(s).exp(t).asyn(1,lim_inf_indx:lim_sup_indx);
        sujeto(s).exp(t).asyn_pro(2,:)=sujeto(s).exp(t).asyn(2,lim_inf_indx:lim_sup_indx)-pert_bip(t);
        sujeto(s).condicion(5).serie_prom(2,:)=sujeto(s).exp(t).asyn_pro(2,:);
          sujeto(s).condicion(5).mech_size=steps(5);
          if out(t)==0
            E(c5,:,1)=sujeto(s).exp(t).asyn_pro(1,:);
            E(c5,:,2)=sujeto(s).exp(t).asyn_pro(1,:)-pre_baseline(t);
            E_prebaseline(c5)=pre_baseline(t);
            E_posbaseline(c5)=pos_baseline(t);
            E_preITI(c5)=meanpreITI(t);
            E_posITI(c5)=meanposITI(t);
            c5=c5+1;
        else
        end
        
    elseif sujeto(s).exp(t).mech_size==steps(6)
        lim_inf_indx=find(sujeto(s).exp(t).asyn(2,:)-pert_bip(t)==min_bip_cond(6));
        lim_sup_indx=find(sujeto(s).exp(t).asyn(2,:)-pert_bip(t)==max_bip_cond(6));
        sujeto(s).exp(t).asyn_pro(1,:)=sujeto(s).exp(t).asyn(1,lim_inf_indx:lim_sup_indx);
        sujeto(s).exp(t).asyn_pro(2,:)=sujeto(s).exp(t).asyn(2,lim_inf_indx:lim_sup_indx)-pert_bip(t);
        sujeto(s).condicion(6).serie_prom(2,:)=sujeto(s).exp(t).asyn_pro(2,:);
          sujeto(s).condicion(6).mech_size=steps(6);
          if out(t)==0
            F(c6,:,1)=sujeto(s).exp(t).asyn_pro(1,:);
            F(c6,:,2)=sujeto(s).exp(t).asyn_pro(1,:)-pre_baseline(t);
            F_prebaseline(c6)=pre_baseline(t);
            F_posbaseline(c6)=pos_baseline(t);
            F_preITI(c6)=meanpreITI(t);
            F_posITI(c6)=meanposITI(t);
            c6=c6+1;
        else
        end
        
    elseif sujeto(s).exp(t).mech_size==steps(7)
        lim_inf_indx=find(sujeto(s).exp(t).asyn(2,:)-pert_bip(t)==min_bip_cond(7));
        lim_sup_indx=find(sujeto(s).exp(t).asyn(2,:)-pert_bip(t)==max_bip_cond(7));
        sujeto(s).exp(t).asyn_pro(1,:)=sujeto(s).exp(t).asyn(1,lim_inf_indx:lim_sup_indx);
        sujeto(s).exp(t).asyn_pro(2,:)=sujeto(s).exp(t).asyn(2,lim_inf_indx:lim_sup_indx)-pert_bip(t);
        sujeto(s).condicion(7).serie_prom(2,:)=sujeto(s).exp(t).asyn_pro(2,:);
        sujeto(s).condicion(7).mech_size=steps(7);
        if out(t)==0
            G(c7,:,1)=sujeto(s).exp(t).asyn_pro(1,:);
            G(c7,:,2)=sujeto(s).exp(t).asyn_pro(1,:)-pre_baseline(t);
            G_prebaseline(c7)=pre_baseline(t);
            G_posbaseline(c7)=pos_baseline(t);
            G_preITI(c7)=meanpreITI(t);
            G_posITI(c7)=meanposITI(t);
            c7=c7+1;
        else
        end
    end
    sujeto(s).exp(t).out=out(t);
end

sujeto(s).condicion(1).serie_prom(1,:)=mean(A(:,:,1));
sujeto(s).condicion(1).serie_prom_bs(1,:)=mean(A(:,:,2));
sujeto(s).condicion(1).mean_prebaseline=mean(A_prebaseline);
sujeto(s).condicion(1).std_prebaseline=std(A_prebaseline);
sujeto(s).condicion(1).mean_posbaseline=mean(A_posbaseline);
sujeto(s).condicion(1).std_posbaseline=std(A_posbaseline);
sujeto(s).condicion(1).mean_preITI=mean(A_preITI);
sujeto(s).condicion(1).std_preITI=std(A_preITI);
sujeto(s).condicion(1).mean_posITI=mean(A_posITI);
sujeto(s).condicion(1).std_posITI=std(A_posITI);
sujeto(s).condicion(1).n=c1;

sujeto(s).condicion(2).serie_prom(1,:)=mean(B(:,:,1));
sujeto(s).condicion(2).serie_prom_bs(1,:)=mean(B(:,:,2));
sujeto(s).condicion(2).mean_prebaseline=mean(B_prebaseline);
sujeto(s).condicion(2).std_prebaseline=std(B_prebaseline);
sujeto(s).condicion(2).mean_posbaseline=mean(B_posbaseline);
sujeto(s).condicion(2).std_posbaseline=std(B_posbaseline);
sujeto(s).condicion(2).mean_preITI=mean(B_preITI);
sujeto(s).condicion(2).std_preITI=std(B_preITI);
sujeto(s).condicion(2).mean_posITI=mean(B_posITI);
sujeto(s).condicion(2).std_posITI=std(B_posITI);
sujeto(s).condicion(2).n=c2;

sujeto(s).condicion(3).serie_prom(1,:)=mean(C(:,:,1));
sujeto(s).condicion(3).serie_prom_bs(1,:)=mean(C(:,:,2));
sujeto(s).condicion(3).mean_prebaseline=mean(C_prebaseline);
sujeto(s).condicion(3).std_prebaseline=std(C_prebaseline);
sujeto(s).condicion(3).mean_posbaseline=mean(C_posbaseline);
sujeto(s).condicion(3).std_posbaseline=mean(C_posbaseline);
sujeto(s).condicion(3).mean_preITI=mean(C_preITI);
sujeto(s).condicion(3).std_preITI=std(C_preITI);
sujeto(s).condicion(3).mean_posITI=mean(C_posITI);
sujeto(s).condicion(3).std_posITI=std(C_posITI);
sujeto(s).condicion(3).n=c3;

sujeto(s).condicion(4).serie_prom(1,:)=mean(D(:,:,1));
sujeto(s).condicion(4).serie_prom_bs(1,:)=mean(D(:,:,2));
sujeto(s).condicion(4).mean_prebaseline=mean(D_prebaseline);
sujeto(s).condicion(4).std_prebaseline=std(D_prebaseline);
sujeto(s).condicion(4).mean_posbaseline=mean(D_posbaseline);
sujeto(s).condicion(4).std_posbaseline=std(D_posbaseline);
sujeto(s).condicion(4).mean_preITI=mean(D_preITI);
sujeto(s).condicion(4).std_preITI=std(D_preITI);
sujeto(s).condicion(4).mean_posITI=mean(D_posITI);
sujeto(s).condicion(4).std_posITI=std(D_posITI);
sujeto(s).condicion(4).n=c4;

sujeto(s).condicion(5).serie_prom(1,:)=mean(E(:,:,1));
sujeto(s).condicion(5).serie_prom_bs(1,:)=mean(E(:,:,2));
sujeto(s).condicion(5).mean_prebaseline=mean(E_prebaseline);
sujeto(s).condicion(5).std_prebaseline=std(E_prebaseline);
sujeto(s).condicion(5).mean_posbaseline=mean(E_posbaseline);
sujeto(s).condicion(5).std_posbaseline=std(E_posbaseline);
sujeto(s).condicion(5).mean_preITI=mean(E_preITI);
sujeto(s).condicion(5).std_preITI=std(E_preITI);
sujeto(s).condicion(5).mean_posITI=mean(E_posITI);
sujeto(s).condicion(5).std_posITI=std(E_posITI);
sujeto(s).condicion(5).n=c5;

sujeto(s).condicion(6).serie_prom(1,:)=mean(F(:,:,1));
sujeto(s).condicion(6).serie_prom_bs(1,:)=mean(F(:,:,2));
sujeto(s).condicion(6).mean_prebaseline=mean(F_prebaseline);
sujeto(s).condicion(6).std_prebaseline=std(F_prebaseline);
sujeto(s).condicion(6).mean_posbaseline=mean(F_posbaseline);
sujeto(s).condicion(6).std_posbaseline=std(F_posbaseline);
sujeto(s).condicion(6).mean_preITI=mean(F_preITI);
sujeto(s).condicion(6).std_preITI=std(F_preITI);
sujeto(s).condicion(6).mean_posITI=mean(F_posITI);
sujeto(s).condicion(6).std_posITI=std(F_posITI);
sujeto(s).condicion(6).n=c6;

sujeto(s).condicion(7).serie_prom(1,:)=mean(G(:,:,1));
sujeto(s).condicion(7).serie_prom_bs(1,:)=mean(G(:,:,2));
sujeto(s).condicion(7).mean_prebaseline=mean(G_prebaseline);
sujeto(s).condicion(7).std_prebaseline=std(G_prebaseline);
sujeto(s).condicion(7).mean_posbaseline=mean(G_posbaseline);
sujeto(s).condicion(7).std_posbaseline=std(G_posbaseline);
sujeto(s).condicion(7).mean_preITI=mean(G_preITI);
sujeto(s).condicion(7).std_preITI=std(G_preITI);
sujeto(s).condicion(7).mean_posITI=mean(G_posITI);
sujeto(s).condicion(7).std_posITI=std(G_posITI);
sujeto(s).condicion(7).n=c7;

%% Series temporales de cada sujeto
colores=['g' 'r' 'k' 'm' 'c' 'b' 'y'];
figure(s)
    for p=1:7
        title(['sujeto ' num2str(s)]) 
        color=colores(p);
        plot(sujeto(s).condicion(p).serie_prom(2,:),sujeto(s).condicion(p).serie_prom(1,:)-sujeto(s).condicion(p).mean_prebaseline,'.-','Color',color)
        hold on
        plot(sujeto(s).condicion(p).serie_prom(2,:),sujeto(s).condicion(p).serie_prom_bs(1,:),'--','Color',color)
        hold all
    end
% legend('-1.5','-1.0','-0.5','0','+0.5','+1.0','+1.5')

save('sujetos.mat','sujeto')

end

%% Efecto de la perturbacion espacial en NMA
figure()
x=[-1.5 -1.0 -0.5 0 0.5 1.0 1.5];

for s=1:2
    for p=1:7
        
        y=sujeto(s).condicion(p).serie_prom_bs(1,(sujeto(s).condicion(p).serie_prom(2,:)==1));

        plot(x(p),y,'LineWidth',5,'MarkerSize',15)
        hold all
        xlim([-2 2])
    end
end

legend('s1', 's2') 
