function [Cgl,Mgl,Kgl]=SeismicModalMDOF2DFrames2...
    (coordxy,A,unitWeightEl,qbarxy,Edof,E,I,ni,nf,AlfaBeta,g)
% SYNTAX : 
% [Cgl,Mgl,Kgl]=SeismicModalMDOF2DFrames2...
%  (coordxy,A,unitWeightEl,qbarxy,Edof,E,I,ni,nf,AlfaBeta)
%---------------------------------------------------------------------
%    PURPOSE
%     To compute the global stiffness matrix of a plane frame as well as
%     the global mass matrix and the global damping matrix.
% 
%    INPUT:  coordxy:           Node coordinates of the structure [x,y]
%
%            A:                 Cross-sectional elements' area
%
%            E,I:               Modulus of Elasticity and Cross-sectional 
%                               inertia of the frame's elements
%
%            unitWeightEl:      unit weight material of each element:
%                               Size: nbars x 1
%
%            qbarxy:            uniformly distributed loads. 
%                               Size: nbars x 2.
%                               The first column corresponds to the
%                               distributed loads in the local X' direction
%                               and the second column the the loads
%                               distributed in the local Y' direction.
%
%            Edof:              Topology matrix. Size: nbars x 7
%                               Format: [bar-i, dof-ni, dof-nf,
%                                        ...]
%
%            ni, nf:            Vectors containing the initial node and 
%                               final node of each element
%
%            AlfaBeta:          Rayleigh's coefficients for the computation
%                               of the damping matrix: [alpha, beta]
%
%    OUTPUT: Cgl:               Global damping matrix
%
%            Mgl:               Global Mass matrix
%            Kgl:               Global Stiffness matrix
%
%--------------------------------------------------------------------

% LAST MODIFIED: L.Verduzco    2023-06-07
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro
%--------------------------------------------------------------------
nnodes=length(coordxy(:,1)); nbars=length(E);

%% Stiffness and Mass matrices
Kgl=zeros(3*nnodes);
Mgl=zeros(3*nnodes);
Cgl=zeros(3*nnodes);
for i=1:nbars      
    ex=[coordxy(ni(i),1) coordxy(nf(i),1)];
    ey=[coordxy(ni(i),2) coordxy(nf(i),2)];
    PVe=unitWeightEl(i,1)*A(i)/g-qbarxy(i,2)/g; % unit weight of each element
                                          % The distributed downward loads
                                          % on the BEAMS are considered.
    ep=[E(i) A(i) I(i) PVe AlfaBeta'];
    
    %% Stiffness matrix, Mass matrix and Damping matrix
    
    [Kebar,Mebar,Cebar]=beam2d(ex,ey,ep); % CALFEM function
                                          % To download the CALFEM package
                                          % visit its repository at:
                                          % https://github.com/CALFEM/calfem-matlab
    
    [Kgl]=assem(Edof(i,:),Kgl,Kebar);
    [Mgl]=assem(Edof(i,:),Mgl,Mebar);
    [Cgl]=assem(Edof(i,:),Cgl,Cebar);
end 

