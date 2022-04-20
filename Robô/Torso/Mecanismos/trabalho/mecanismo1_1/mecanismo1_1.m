%função x2 = l2 + ((l1 + l2 - sqrt(l1^2 + l2^2))/2)*cos(t) e y2 = - l1 + ((l1 + l2 - sqrt(l1^2 + l2^2))/2)*sin(t)
%função escolhida para não haver singularidades
%%Simulador de mecanismos planos
clear
%%
%Construção do modelo baseado em variáveis simbólicas
%Definição das variáveis
syms x1 y1 phi1 x2 y2 phi2 real %posição
syms x1p y1p phi1p x2p y2p phi2p real %velocidade
syms l l1 l2 %constantes do modelo
syms t %tempo

%Construção do vetor de coordenadas generalizadas e a sua derivada
q = [x1;y1;phi1;x2;y2;phi2];

qp = [x1p;y1p;phi1p;x2p;y2p;phi2p];

%Construção do vetor de restrições cinemáticas e diretoras
PhiK = [x1;
        y1;
        x1 + l1*cos(phi1) - x2 + l2*cos(phi2)
        y1 + l1*sin(phi1) - y2 + l2*sin(phi2)];

PhiD=[x2 - l2 - ((l1 + l2 - sqrt(l1^2 + l2^2))/2)*cos(t);
      y2 + l1 - ((l1 + l2 - sqrt(l1^2 + l2^2))/2)*sin(t)];

Phi=[PhiK;PhiD];

%Construção da jacobiana, nu e gamma

[Phiq, nu, Gamma] = construtor (q, qp, t, Phi);

%Os resultados tem que ser passados para as respectivas funções e nos nomes
%das variáveis convertidos para x(i).
%%
clear
%Entrada da estimativa da posição inicial e constantes do modelo
l1 = 5;
l2 = 5;
l = 10;
x0 = [0;0;0.1;0;0;0];
t0 = 0;

%Início da simulação
x = x0;
t = t0;
V=[];
Y=[];
i=2;

while t<2*pi
    %Análise de posição
    %Definir as funções funcaoseqnln e jacobiana
    A=jacobiana(x,t);
    if det(A)==0
    warning('jacobiana singular');
    end
    x=newtonraphson(x,t);

    %Análise de velocidade
    b=nu(x,t);
    A=jacobiana(x,t);
    xp = eliminacaogauss(A,b);

    %Análise de aceleração
    b = Gamma(x,xp,t);
    xpp = eliminacaogauss(A,b);
    
    %velocidade do cursor perpendicular à barra
    vc = -xp(4)*sin(x(6))+xp(5)*cos(x(6))-l*xp(6);
    V(:,i) = vc;

    %Atualização do tempo e estimativa da nova posição
    dt=0.0001;
    x = x + xp*dt + xpp*dt^2/2;
    t = t + dt;
    X(:,i)=[x;t];
    i=i+1;
end

%área da função
fprintf('A área da curva é %.4f\n',polyarea(X(4,2:end),X(5,2:end)))

%lrtheta
fprintf('l*r*theta = %.4f\n',l2*trapz(X(7,2:end),V(2:end)))
