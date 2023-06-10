% DynamicSeismicAnalysis_2DFrames_Damping_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To compute the Dynamic Seismic analysis for a Reinforced Concrete
%    Plane Frame with/without damping effects.
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-06-07
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------

clc 
clear all

nnodes=12;
nbars=15;

%% Materials
fpc=[300;
     300;
     300;
     300;
     300;
     300;
     300;
     300;
     300;
     300;
     300;
     300;
     300;
     300;
     300];

% Modulus of Elasticity of each element
E=14000*(fpc).^0.5;

%% Geometry
dimensions=[20 50;
            20 50;
            20 50;
            20 50;
            20 50;
            20 50;
            20 50;
            20 50;
            20 50;
            20 50;
            20 50;
            20 50;
            20 50;
            20 50;
            20 50];

A=dimensions(:,1).*dimensions(:,2);
I=1/12.*dimensions(:,1).*dimensions(:,2).^3;

% Coordinates of each node for each bar
coordxy=[0 0;
         0 300;
         0 600;
         0 900;
         600 0;
         600 300;
         600 600;
         600 900;
         1200 0;
         1200 300;
         1200 600;
         1200 900]; 
                 
%% Topology
% Node conectivity
ni=[1;2;3;4;3;2;5;6;7;8; 7; 6; 9; 10;11];
nf=[2;3;4;8;7;6;6;7;8;12;11;10;10;11;12];

% Length of each element
L=sqrt((coordxy(nf,1)-coordxy(ni,1)).^2+(coordxy(nf,2)-coordxy(ni,2)).^2);

% Topology matrix
Edof=zeros(nbars,7);
for i=1:nbars
    Edof(i,1)=i;
    Edof(i,2)=ni(i)*3-2;
    Edof(i,3)=ni(i)*3-1;
    Edof(i,4)=ni(i)*3;
    
    Edof(i,5)=nf(i)*3-2;
    Edof(i,6)=nf(i)*3-1;
    Edof(i,7)=nf(i)*3;
end

%% Boundary conditions
% Prescribed DOF - [N-DOF, displacement]
bc=[1 0;
    2 0;
    3 0;
    13 0;
    14 0;
    15 0;
    25 0;
    26 0;
    27 0];

%% Loads

type_elem=[1 "Col";
           2 "Col";
           3 "Col";
           4 "Beam";
           5 "Beam";
           6 "Beam";
           7 "Col";
           8 "Col";
           9 "Col";
           10 "Beam";
           11 "Beam";
           12 "Beam";
           13 "Col";
           14 "Col";
           15 "Col"];
       
% Distributed Loads on beams
beamsLoads=[1 -50;
          2 -50;
          3 -50;
          4 -50;
          5 -50;
          6 -50];

elem_cols=[];
elem_beams=[];
beams=0;
cols=0;
for j=1:nbars % To identify which elements are beams and which are columns
    if type_elem(j,2)=="Beam"
        beams=beams+1;
        elem_beams=[elem_beams,j];
    elseif type_elem(j,2)=="Col"
        cols=cols+1;
        elem_cols=[elem_cols,j];
    end
end
qbarxy=zeros(nbars,2);
for i=1:beams
    qbarxy(elem_beams(i),2)=1.1*(beamsLoads(i,2));
end

modal=1;

%% Damping matrix (for the damped case)
omega1=5;
omega2=100;

zeta1=0.25;
zeta2=0.2;

D=[1/(2*omega1) omega1/2;
    1/(2*omega2) omega2/2];

theta=[zeta1;
    zeta2];

AlfaBeta=D\theta; % Rayleigh coefficients

%% Max Seismic response from the CFE-15 spectrum

g=981; % gravity acceleration (cm/sec^2)
Fsit=2.4; FRes=3.8; % Factores de sitio y de respuesta
a0_tau=200; % cm/sec^2

ro=0.8; % Redundance factor
alf=0.9; % Irregularity factor
Q=4; % Seismic behaviour factor

Ta=0.1;
Tb=0.6;
Te=0.5; % Structure's period
k=1.5; % Design spectrum slope
Qp=1+(Q-1)*sqrt(Te/(k*Tb)); % Ductility factor

Ro=2.5; % Over-resistance index
R=Ro+1-sqrt(Te/Ta); % Over-resistance factor

sa=a0_tau*Fsit*FRes/(R*Qp*alf*ro); % Reduced pseudo-acceleration (cm/sec^2)

%% Modal analysis
pvconc=0.0024; % unit weight of concrete
unitWeightElm=zeros(nbars,1)+pvconc;

% Consistent mass method
[Cgl,Mgl,Kgl]=SeismicModalMDOF2DFrames2...
(coordxy,A,unitWeightElm,qbarxy,Edof,E,I,ni,nf,AlfaBeta,g);

Dof=zeros(nnodes,3);
for i=1:nnodes
    Dof(i,1)=3*i-2;
    Dof(i,2)=3*i-1;
    Dof(i,3)=3*i;
end
[Ex,Ey]=coordxtr(Edof,coordxy,Dof,2);

%% Dynamic analysis
% Time discretization
dt=0.05;
ttotal=10;
t=0:dt:ttotal;
npoints=length(t);

% Ground acceleration history
tload=1.5; % duration of external excitation

g=sa*cos(5*t); % Acceleration in time
for i=1:length(g)
    if t(i)>tload
        g(i)=20*cos(30*t(i));
    end
end

figure(1)
grid on
plot(t,g,'b -','LineWidth',1.8)
hold on
xlabel('Time (sec)')
ylabel('Acceleration (Kg/cm^2)')
title('Ground acceleration in time')

dofhist=[16 19 22]; % dof to evaluate

% Forces history
f(:,1)=zeros(3*nnodes,1);
for i=1:npoints
    % Modal analysis without Damping
    [f(:,i+1),T(:,i),La(:,i),Egv]=ModalsMDOF2DFrames2(Mgl,Kgl,...
        bc,g(i),modal);
    
end
% Analysis in time without viscous damping
beta=0.25;
gamma=0.5;

d0=zeros(3*nnodes,1);
v0=zeros(3*nnodes,1);

% Solving the motion equation with the Newmark-Beta method
[Dsnap,D,V,A]=NewmarkBetaMDOF2(Kgl,Cgl,Mgl,d0,v0,dt,beta,gamma,t,f,...
                              dofhist,bc);

%% Dynamic displacement analysis per DOF
figure(2)
grid on
plot(t,D(1,:),'b -','LineWidth',1.8)
legend(strcat('DOF-',num2str(dofhist(1))))
hold on
for i=2:length(dofhist)
    plot(t,D(i,:),'LineWidth',1.8,'DisplayName',...
        strcat('DOF-',num2str(dofhist(i))))
end
hold on
xlabel('Time (sec)')
ylabel('Displacements (cm)')
title('Displacements in time per DOF')

%% Deformation history of structures
dtstep=5;

WidthStruc=max(coordxy(:,1));
HeightStruc=max(coordxy(:,2));
figure(3)
axis('equal')
axis off
sfac=10;
title(strcat('Deformed structures in time. Scale x ',num2str(sfac)))
for i=1:5
    Ext=Ex+(i-1)*(WidthStruc+300);
    eldraw2(Ext,Ey,[2 3 0]);
    Edb=extract(Edof,Dsnap(:,dtstep*i-(dtstep-1)));
    eldisp2(Ext,Ey,Edb,[1 2 2],sfac);
    Time=num2str(t(5*i-4));
    NotaTime=strcat('Time(seg)= ',Time);
    text((WidthStruc+150)*(i-1)+50,150,NotaTime);
end

Eyt=Ey-(HeightStruc+200);
for i=6:10
    Ext=Ex+(i-6)*(WidthStruc+300);
    eldraw2(Ext,Eyt,[2 3 0]);
    Edb=extract(Edof,Dsnap(:,dtstep*i-(dtstep-1)));
    eldisp2(Ext,Eyt,Edb,[1 2 2],sfac);
    Time=num2str(t(5*i-4));
    NotaTime=strcat('Time(seg)= ',Time);
    text((WidthStruc+150)*(i-6)+50,-250,NotaTime)
    
end
% ------------------------------------------------------------------