%% Controle de Qualidade de Dados MOVAR
%% Introdução
% Dados coletados com instrumentos muitas vezes podem apresentar falhas. O objetivo 
% deste documento é trabalhar uma forma automatica lendo os dados de temperatura 
% da linha de monitoramento AX97, identificando e retirando falhas nos dados. 
% Esse documento está divido em:
% 
% * Importação dos dados brutos em .txt.
% * Importação da Climatologia WOA2018.
% * Pré-processamento básico.
% * Identificação e tratamento de perfis de temperatura irregulares.
% * Verificação final e gravação dos dados.
% 
% Nota: o software utilizado na coleta de dados é o AMVERSEAS V9, outras 
% versões ou outros software podem gerar arquivos de saida diferentes, logo será 
% necessário realizar ajustes para utilizar esse Script.
%% Primeira parte
%% Importação e organização dos dados brutos

clearvars
% Caminho onde estão localizados os .txt:
location="/media/samantha/OCEAN/Desktop/NOTE-MOVAR/ENVIAR/MOVAR_atualizado_201709/2 - cruzeiros/POIT_2016_10/txt/";
files=dir(location+"*.txt"); % pega todos os arquivos no fomato txt

% importa o arquivo XBTLog, se houver
xbtlog= string(importdata('/media/samantha/OCEAN/Desktop/NOTE-MOVAR/ENVIAR/MOVAR_atualizado_201709/2 - cruzeiros/POIT_2016_10/XBTLog_PWAZ.txt', '\t'));
xbtlog(1)=[];xbtlog=split(xbtlog,',');


% inicializando a tabela:
Poit(1,:)=table(1,{0.5},{0.5},0.5,0.5,datetime('now'),categorical("UNVALUED"),'VariableNames'...
    ,{'Number','Depth','Temp','Latitude','Longitude','Date','QC'}); 


disp("Importando "+length(files)+" perfis de temperatura")
for i=1:length(files)
    [Header,XBT] = import_txt(location+files(i).name);
    
    lat=char(Header.Value(2));lat([3,9,10])=[];lat=-str2double(lat)/100;
    lon=char(Header.Value(3));lon([4,10,11])=[];lon=-str2double(lon)/100;
    
    Poit(i,:)=table(i,{XBT.Depth},{XBT.Temperature},lat,lon,datetime(strjoin(Header.Value(7:11))...
        ,'InputFormat','yyyy MM dd HH mm'),categorical("UNVALUED"));
    if(contains(xbtlog(i,8),'test','IgnoreCase',true))
        Poit.QC(i)="DROPTEST";
    end
end
clear i location files Header XBT lat lon
%% Importação da Climatologia WOA2018
% A climatologia mensal é carregada aqui (mude o caminho para a localização 
% do arquivo no seu computador), de acordo com o mês em que houve o cruzeiro.
% 
% Faço aqui um corte nas latitude e longitudes próximas da área onde são 
% coletados os dados do MOVAR, assim fica mais leve para trabalharmos, porém a 
% cada perfil de temperatura será interpolado um perfil de temperatura da climatologia 
% para comparar.
% 
% *Atenção: o desvio padão da climatologia muitas vezes não existe em vários 
% pontos, ficando com valores NaNs.*

mes=Poit.Date.Month(1); if mes<10; mes="0"+mes;end
woa="/media/samantha/OCEAN/Pesquisa/Climatology/WOA18/temp/woa18_decav_t"+mes+"_04.nc";
disp("Using: "+ncreadatt(woa,'/',"title"))
woaLat=double(ncread(woa,'lat',264,25)); %end 289
woaLon=double(ncread(woa,'lon',544,65)); %end 609
woaProf=double(ncread(woa,'depth'));
woaTemp=ncread(woa,'t_an',[544,264,1,1],[65,25,57,1]);
woaStd=ncread(woa,'t_sd',[544,264,1,1],[65,25,57,1]); clear woa;
%% Segunda Parte
%% Pré-processamento básico
% Nesta etapa são feitas várias verificações nos perfis de temperatura para 
% identificar possíveis erros e identificar cada um. Em cada um dos perfis, cortar 
% aquele "traço" característico do xbt quando perde contato com o computador, 
% e identificar os perfir em categorias:
% 
% * UNVALUED: quando não foi possível avaliar o perfil.
% * GOOD: quando o perfil está comprovadamente dentro dos padrões.
% * BAD: quando mais de um indicativo que o perfil está ruim.
% * QUEST: quando algum ídice pode indicar que o perfil está ruim.
% * DROPTEST: os identificados no XBTLog.

for i=1:size(Poit,1)
    
    depth=cell2mat(Poit.Depth(i));
    temp=cell2mat(Poit.Temp(i));
    
    % Limpa o corte de profundidade
    if(~isempty(find(diff(temp)>0.1 & depth(2:end)>800, 1)))
        e=find(diff(temp)>0.1 & depth(2:end)>800);
        temp(e(1):end)=[]; depth(e(1):end)=[];
    end
    
    % Gerando um perfil a partir do WOA
    Twoa=squeeze(interp3(woaLat,woaLon,woaProf,woaTemp,Poit.Latitude(i),Poit.Longitude(i),depth));
    Stdwoa=squeeze(interp3(woaLat,woaLon,woaProf,woaStd,Poit.Latitude(i),Poit.Longitude(i),depth));

    % Atualizando variaveis
    Poit.Twoa(i)={Twoa};
    Poit.Stdwoa(i)={Stdwoa};
    Poit.RMSEwoa(i)=rmse(Twoa',temp');
    Poit.Temp(i)={temp};
    Poit.Depth(i)={depth};
    Poit.isnan(i)=sum(isnan(Twoa))/length(Twoa);
    
    % Cálculo do erro quadrático médio com o perfil anterior e posterior
    if i>1
        temp0=cell2mat(Poit.Temp(i-1));l=min(size(temp,1),size(temp0,1));
        Poit.RMSEanterior(i)=rmse(temp0(1:l),temp(1:l));
    end
    if i<size(Poit,1)
        temp0=cell2mat(Poit.Temp(i+1));l=min(size(temp,1),size(temp0,1));
        Poit.RMSEposterior(i)=rmse(temp0(1:l),temp(1:l));
    end
    
    
    if(Poit.QC(i)==categorical("DROPTEST"))
    else
        if((Poit.RMSEwoa(i)<=2 && Poit.isnan(i)==0) && (Poit.RMSEposterior(i)<=1 && Poit.RMSEanterior(i)<=1))
            Poit.QC(i)=categorical("GOOD");
            plot(temp,-depth);hold on; title('Perfis já aprovados')
            xlabel('Temperature');ylabel('Depth')
        elseif(Poit.RMSEwoa(i)>2 || Poit.isnan(i)>0.50)
            Poit.QC(i)=categorical("QUEST");
        end
        if(Poit.RMSEwoa(i)>2 && Poit.isnan(i)==0)
            Poit.QC(i)=categorical("BAD");
        end
    end
    
end
a=summary(Poit); a1=a.QC.Categories; a2=a.QC.Counts;
drop=a2(strcmp(a1,'DROPTEST'));if isempty(drop); drop=0; end
unv=a2(strcmp(a1,'UNVALUED'));if isempty(unv); unv=0; end
good=a2(strcmp(a1,'GOOD'));if isempty(good); good=0; end
bad=a2(strcmp(a1,'BAD')); if isempty(bad); bad=0; end
quest=a2(strcmp(a1,'QUEST'));if isempty(quest); quest=0; end
disp(good+" Perfis aprovados")
disp(bad+quest+" Perfis ruins ou questionáveis")
disp(unv+" Perfis que não foram avaliados")
disp(drop+" Testes de canistes encontrados")
hold off
Poit = movevars(Poit, 'Twoa', 'Before', 'Latitude');
Poit = movevars(Poit, 'Stdwoa', 'Before', 'Latitude');
Poit = movevars(Poit, 'RMSEwoa', 'Before', 'QC');
Poit.Properties.VariableDescriptions{9} = 'Root mean square error from Climatology.';
clear e i temp depth Stdwoa woaLat woaLon woaProf woaTemp woaStd a drop a1 a2 quest unv good bad
%% Pré-processamento avançado
%% Identificação e tratamento de perfis de temperatura irregulares.
% Todos os perfis já estão indicados com uma Flag, agora o usuário vai apenas 
% verificar se houve algum erro de análise e o próprio usuário poderá alterar 
% a flag de um perfil específico.

% Separando os perfins bons, dos ruins, dos drop test
PoitGood = Poit(Poit.QC == 'GOOD',:);

PoitBad = Poit(Poit.QC ~= 'GOOD' & Poit.QC ~= 'DROPTEST',:)

PoitDrop=Poit(Poit.QC == 'DROPTEST',:)
for i=1:size(PoitDrop,1);plot(cell2mat(PoitDrop.Temp(i)),-cell2mat(PoitDrop.Depth(i)));hold on; title('Drop Tests encontrados')
        xlabel('Temperature');ylabel('Depth');end

for i=1:size(PoitBad,1)
    depth=cell2mat(PoitBad.Depth(i));
    temp=cell2mat(PoitBad.Temp(i));
    Twoa=cell2mat(PoitBad.Twoa(i));
    figure
    for k=1:size(PoitGood,1)
        plot(cell2mat(PoitGood.Temp(k)),-cell2mat(PoitGood.Depth(k)),'Color',[0.790 0.790 0.790],'HandleVisibility','off');hold on
    end
    plot(Twoa,-depth,'r','LineWidth',1.5);hold on
    plot(temp,-depth,'k','LineWidth',1.5);hold on;
    legend('show','Location','northwest','String',{'WOA','XBT'})

    axis([-5 30 -900 0])
    grid on
    hold off
    
    title('Launch # '+string(PoitBad.Number(i))+" - "+string(PoitBad.QC(i)))
    xlabel('Temperature')
    ylabel('Depth')
    set(gca,'FontSize',12)
    if(~isnan(PoitBad.RMSEwoa(i)))
        text(-4,-200,"RMSE="+PoitBad.RMSEwoa(i),'FontSize',16)
    end
    
end
%% Correção de avaliação errada do computador
%%
% Aqui você pode mudar a bandeira de algum perfil ruim para bom.
Poit.QC(5)='GOOD';

% Atualizando a lista de perfins bons
PoitGood = Poit(Poit.QC == 'GOOD',:);
PoitBad = Poit(Poit.QC ~= 'GOOD' & Poit.QC ~= 'DROPTEST',:)
%% Terceira parte
%% Verificação final
% Nessa etapa o script apenas vai sobrepor todos os perfis bons em um gráfico 
% e todos os ruins em outro grafico, além de mostrar o balanço final e savar os 
% dados tanto em .mat quanto em .csv

hold off
for k=1:size(PoitGood,1)
    plot(cell2mat(PoitGood.Temp(k)),-cell2mat(PoitGood.Depth(k)));hold on; title('Perfis já aprovados')
    xlabel('Temperature');ylabel('Depth')
end
hold off
for k=1:size(PoitBad,1)
    plot(cell2mat(PoitBad.Temp(k)),-cell2mat(PoitBad.Depth(k)));hold on; title('Perfis reprovados')
    xlabel('Temperature');ylabel('Depth')
end

a=summary(Poit); a1=a.QC.Categories; a2=a.QC.Counts;
drop=a2(strcmp(a1,'DROPTEST'));if isempty(drop); drop=0; end
unv=a2(strcmp(a1,'UNVALUED'));if isempty(unv); unv=0; end
good=a2(strcmp(a1,'GOOD'));if isempty(good); good=0; end
bad=a2(strcmp(a1,'BAD')); if isempty(bad); bad=0; end
quest=a2(strcmp(a1,'QUEST'));if isempty(quest); quest=0; end
disp(good+" Perfis aprovados")
disp(bad+quest+" Perfis ruins ou questionáveis")
disp(unv+" Perfis que não foram avaliados")
disp(drop+" Testes de canistes encontrados")
%% Gravação dos dados
%%
cruzeiro=mes+string(Poit.Date.Year(1));
% retirando informações desnecessárias da tabela
PoitGood.isnan=[];PoitGood.RMSEanterior=[];PoitGood.Twoa=[];PoitGood.RMSEposterior=[];PoitGood.RMSEwoa=[];
save("poit"+cruzeiro+".mat",'PoitGood')
writetable(PoitGood,'poit'+cruzeiro+'.csv')