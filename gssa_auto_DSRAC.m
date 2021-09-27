%% Algoritmo para automatizão do processo de GSSA
% CONVERSOR DSRAC
%%
close all;
clear global; clear variables; clc;
%% Define equações generalizadas
syms Lm Rm Cc Rs Ls Cx R Cb Ca n sf sl q1 q2 Vin ws real;
Ag = [-Rm/Lm -sl/Lm -(n*Rm*q1/Lm + n*Rm*q2/Lm) 0 0;...
      sl/Cc 0 n*q2/Cc 0 0;...
      -(n*Rm*q1/Ls + n*Rm*q2/Ls) -n*q2/Ls -(n*n*Rm*q1/Ls + n*n*Rm*q2/Ls) (q1/Ls + q2/Ls) -q1/Ls;...
      0 0 (q1/Cb - q2/Ca) 0 -1/(R*Cx);...
      0 0 (q1/Cx - q2/Cx) 0 -2/(R*Cx)];
Bg = [sf/Lm 0 n*q1/Ls 0 0];
%% Definições gerais do algoritmo
nh = 1; % Número de harmonicas
nsv = length(Ag); % Número de variáveis de estado das eqs generalizadas
niv = 1; % Número de entradas (Fontes do circuito, não o duty)
%% Inicializa matrizes do GAM
% SISTEMA FINAL LINEARIZADO NA FORMA
% dot{x}= A*x + Bd*deltaD + BvdeltaV
h = (2*nh + 1)*nsv;    % Número de estados do sistema médio
A = sym(zeros(h,h));   % Matriz A do sistema final linearizado
B = sym(zeros(h,niv)); % Vetor B do sistema final linearizado
x_i = sym('x', [1 h]); % Variáveis do sistema final
dx = sym(zeros(size(x_i)))'; % Eqs. Dif. N-Lin Finais do Modelo.
dx_i = sym(zeros(1,(nsv*(1+nh))))'; %Eqs de cada coeficiente complexo
assume(x_i,'real');

j = i; % Variável complexa j
% Cria vetor x auxiliar -> vetor de coeficientes k para as variáveis de
% estado
display('Comeco. Definindo vetores e variáeis de estado...');
for k = 1:nsv
    l = (2*nh+1)*(k-1)+1;
    ki = (nh+1)*(k-1)+1;
    x(ki) = x_i(l);
    for kh = 1:nh
        kk = ki + kh;
        p_real = sym(['x' num2str(l+2*kh-1)]); 
        assume(p_real,'real');
        p_img = sym(['x' num2str(l+2*kh)]); 
        assume(p_img,'real');
        x(kk) = p_real + j*p_img;
     end
 end
x = conj(x'); % Vetor linha para os estados
display('Definindo vetores com as EDFs...');
%% Define vetor com EDFs
% edf(1) = s(t); edf(2) = (1-s(t)); edf(3) = q1(t); edf(4) = q2(t);
syms Do wsr kedf Pi real;
% Define vetor de "nomes" das EDFs
edf_nome = [sf;sl;q1;q2];
% Define vetor com os termos das EDFs
edfs = sym(zeros(nh+1,4));
    for k = 2:nh+1
        kn = k-1;
        % #1 EDF s(t)
        edfs(1,1) = (Do);
        s_edf_R = (1/2)*sin(2*Do*kn*Pi)/(kn*Pi);
        s_edf_I = (1/2)*(-1+cos(2*Do*kn*Pi))/(Pi*kn);
        edfs(k,1) = s_edf_R + j*s_edf_I;
        % #2 EDF 1-s(t)
        edfs(1,2) = (1-Do);
        s_edf_R = -(1/2)*(sin(2*Do*kn*Pi)-sin(2*kn*Pi))/(Pi*kn);
        s_edf_I = -(1/2)*(cos(2*Do*kn*Pi)-cos(2*kn*Pi))/(kn*Pi);
        edfs(k,2) = s_edf_R + j*s_edf_I;
        % #3 EDF q1(t)
        edfs(1,3) = wsr/2;
        s_edf_R = (1/2)*sin(kn*Pi*wsr)/(Pi*kn);
        s_edf_I = (1/2)*(-1+cos(kn*Pi*wsr))/(Pi*kn);
        edfs(k,3) = s_edf_R + j*s_edf_I;
        % #4 EDF q2(t)
        edfs(1,4) = wsr/2;
        s_edf_R = -(1/2)*(sin(2*Do*kn*Pi)-sin((2*Do+wsr)*kn*Pi))/(kn*Pi);
        s_edf_I = -(1/2)*(cos(2*Do*kn*Pi)-cos((2*Do+wsr)*kn*Pi))/(kn*Pi);
        edfs(k,4) = s_edf_R + j*s_edf_I;
    end
% return
%%
% %%%%% Fim das definições %%%%%%%%
%  %%%%% COMEÇO DO ALGORITMO %%%%%%
%%
display('Começo do algoritmo...');
% para extração dos termos
% q1 = 1; q2 = 1; sf = 1;sl = 1;
%Itera sobre as equação generalizadas de cada estado
for ki = 1:nsv
    % Adicion termos referentes a Bg nas eqs
    % Como para esse conversor em B só aparece s(t), nao é
    % feito a verificação se o termo é uma soma. 
    Kb = Bg(ki);
    if Kb ~=0
        % identifica EDF do termo
        iEDF = idEDF(Kb,edf_nome);
        if iEDF == 0 % Multiplica apenas por Vin -> Apenas DC
           i = (nh+1)*(ki-1)+1;
           dx_i(i) = dx_i(i) + Kb*Vin;
        else % Multiplica por função de chaveamento. Add Kb*<x>k*Vin
            s = edfs(:,iEDF);
            Kb = Kb/edf_nome(iEDF);
            for kh=0:nh
                i = (nh+1)*(ki-1)+1+kh;
                dx_i(i) = dx_i(i) + Kb*s(kh+1)*Vin;
            end
        end
    end
   % itera termo a termo para adicionar termos de Ag
        for kj = 1:nsv
%             if ki == 1 && kj == 2
%                 display('pause')
%             end
           Kij = Ag(ki,kj);
           if Kij ~= 0
              % identifica EDF do termo
              [childs, issum, isprod, ispow, isfun] = children2(Kij);
              if issum == 1
                somas = length(childs);
                Ki = childs;
              else
                somas = 1;
                Ki = Kij;
              end
              for kt = 1:somas
                    Kaux = Ki(kt);
                    iEDF = idEDF(Kaux,edf_nome);
                for kh = 0:nh
                    l = (nh+1)*(ki-1)+1+kh;
                    i = (nh+1)*(kj-1)+1+kh;
                    if iEDF == 0 % Apenas multiplica pelo estado K*x_i
                        K = Kaux;
                        dx_i(l) = dx_i(l)+ K*x(i);
                    else % Função d existência * estado => K*s*x_i
                      s = edfs(:,iEDF);
                      K = simplify(Kaux/edf_nome(iEDF));
                      cov = gen_cov(kh,nh);
                      [si vo_pos] = size(cov);
                      for c = 1:si
                          es = cov(c,1); %harmonica do estado
                          qh = cov(c,2); %harmonica da função de existencia
                          index = (nh+1)*(kj-1)+1+abs(es);
                          if es < 0
                              es_add = conj(x(index));
                          else
                              es_add = x(index);
                          end
                          i_qh = abs(qh)+1;
                          if qh < 0
                              q_add = conj(s(i_qh));
                          else
                              q_add = s(i_qh);
                          end   
                         dx_i(l) = dx_i(l) + K*es_add*q_add;
                      end
                    end
                end
              end
           end
         % Adiciona os termo +- kjws<x>k
        if ki == kj
            for kh = 1:nh
              l = (nh+1)*(ki-1)+1+kh;
              dx_i(l) = dx_i(l) - j*ws*kh*x(l);
            end
        end
        end
end
% return
%% SEPARA PARTES REAIS E IMAGINÁRIAS PARA CHEGAR NAS EQS FINAIS.
% dx -> contém as equações FINAIS do modelo não linear
display('Separando partes reais e imaginárias....');
for ki = 1:nsv
    i = (2*nh+1)*(ki-1) + 1; % Inicio das variáveis correspondentes ao estado
    for kh = 0:nh
        eq_i = (nh+1)*(ki-1)+1+kh; % ACESSO AO dx_i correto
        eq_j = i+2*(kh-1) + 1; % Acesso ao dx correto para harmonicas
        x_real = simplify(real(dx_i(eq_i)));
        x_imag = simplify(imag(dx_i(eq_i)));
        if kh == 0 %% Insere termo DC do estado
            dx(i) = x_real;
            % Confere se a parte imag é zero
            if x_imag ~= 0 
                dx(i)
                msbox = 'Parte imaginária de um coeficiente DC diferente de zero';
                error(msbox)
            end
        else %% Insere demais termos
            dx(eq_j) = x_real;
            dx(eq_j+1) = x_imag;
        end
    end
end
%% 
% Linearização do modelo não linear das equaçãos em dx -> A e B são as
% matrizes finais para o modelo: dot{x}= A*x_i + Bu => u = [deltaVin deltaD]
display('Linearizando o modelo...');
A = jacobian(dx,x_i);
u = [Vin;Do];
B = jacobian(dx,u);

% %% OBTENÇÃO MODELO NUMÉRICO

% Op Lab
R = 250
Vin = 15.5;
Do = 0.3; n = 4;
Rm = 0.25; 
Lm = 15e-6;  Ls = 1e-6;
Cc = 150e-6; Cr = 1e-6; 
Co = 150e-6;
fs = 50e3;
% 
Cx = 2*Co + Cr;
wr = 1/sqrt(Ls*2*Cr);
fr = wr/(2*pi);
Tr = 1/fr;

ws = 2*pi*fs;
Pi = pi;
wsr = ws/wr;
Ca =1/((Cr+Co)/(Cx*Cr));    
Cb = 1/(1/Ca - 1/Cr);

display('Verifica se o modelo obtido é estável...');
An = eval(subs(A));
autoV = eig(An);
if (sum(autoV > 0) > 0)
    msbox = 'Sistema Instável.';
    error(msbox)
end
% Bss = subs(B(:,1));
% xss = eval(-(An^-1*Bss)*Vin);
display('Modelo estável!');

display('Calculando ponto de operação...')
% remove is0
Ads = A;
Bds = B;
r = (1+2*nh)*2 + 1; % posição de is0 dependendo do número de harmonicas
Ads(r,:) = []; Ads(:,r) = []; Bds(r,:) = [];
A_linear = (subs(Ads));
A_linear = eval(A_linear);
B_mimo = (subs(Bds));
% Calcula ponto de operação. Considerando que Do = constante, sistema é
% linear com u = Vin.
BssN = B_mimo(:,1); % Pega apenas coluna da entrada Vin
BssN = eval(BssN);
% xss = -(A_linear^-1*BssN)*Vin; 
xss = -(A_linear\BssN)*Vin;

% Insere 0 no lugar do antigo is0 para ponto de operação
rowToInsert = r;
rowVectorToInsert = 0;
xssB = [xss(1:rowToInsert-1,:); rowVectorToInsert; xss(rowToInsert:end,:)];
% Monta sistema no espaço de estados
B_linear = subs(Bds,x_i,xssB');
B_mimo = eval(B_linear);
B_linear = eval(B_linear(:,2));
C = zeros(1,length(A_linear));
vo_pos = length(A) - (1+2*nh); % Posição do Vo q é a saída
C(vo_pos) = 1;

% Sistema final em espaço de estados
sist = ss(A_linear,B_linear,C,0);
sis_mimo = ss(A_linear,B_mimo,C,0);

display('Fim do algoritmo =)')


