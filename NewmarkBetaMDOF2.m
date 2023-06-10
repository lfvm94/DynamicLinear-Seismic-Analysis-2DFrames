function [Dsnap,D,V,A]=NewmarkBetaMDOF2(K,C,M,d0,v0,dt,beta,gamma,...
            t,f,dofhist,bc)
% SYNTAX : 
% [Dsnap,D,V,A]=NewmarkBetaMDOF2(K,C,M,d0,v0,dt,beta,gamma,...
%               t,f,dofhist,bc)
%---------------------------------------------------------------------
%    PURPOSE
%     To solve a dynamic system of 2nd order with the numerical
%     integration method "Beta-Newmark".
% 
%    INPUT:  K:                 Global stiffness matrix
%            C:                 Global damping matrix (if any).
%                               Set C=[] for vibration free systems.
%
%            M:                 Global mass matrix
%
%            d0,v0:             Initial displacements and initial
%                               velocities. Vectors of size: n-dof x 1
%
%            beta, gamma:       Chosen parameters for the Beta-Newmark   
%                               method
%
%            t:                 time vector: t0,t1,t2,t3,....tn
%
%            f:                 forces history f(t). Vector of size:
%                               n-dof x n
%
%            dofhist:           degrees of freedom in question
%
%            bc:                boundary condition array. 
%                               Size: n-prescribed-dof x 2
%
%    OUTPUT: D,V,A:             Displacement, velocity and acceleration
%                               history for each DOF in question (taken 
%                               from the vector dofhist)
%
%            Dsnap:             Displacement history for all DOF at each
%                               time step
%
%--------------------------------------------------------------------

% LAST MODIFIED: L.Verduzco    2023-06-09
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro
%--------------------------------------------------------------------

[nd,nd]=size(K);
if isempty(C)==1
    C=zeros(nd,nd);  % In case damping is not considered (vibration free)
end
  
nstep=length(t);
b1 = dt^2*0.5*(1-2*beta);  
b2 = (1-gamma)*dt;
b3 = gamma*dt;              
b4 = beta*dt^2;     

ef=(f(:,1)-C*v0-K*d0);
a0=solveq(M,ef,bc); % -> This is a CALFEM function, to download
                    % CALFEM visit its repository at:
                    % https://github.com/CALFEM/calfem-matlab

D(:,1) = d0(dofhist');  % Initial values of solution 
V(:,1) = v0(dofhist');   
A(:,1) = a0(dofhist');

Keff = M+b3*C+b4*K;         

dnew=d0;    
vnew=v0;    
anew=a0; 
for isnap = 1:nstep
    dpred=dnew+dt*vnew+b1*anew;     
    vpred=vnew+b2*anew;
     
    eff=f(:,isnap+1)-C*vpred-K*dpred;
    
    % Update a,v,d for the next iteration
    anew=solveq(Keff,eff,bc); % -> This is a CALFEM function, to download
                              % CALFEM visit its repository at:
                              % https://github.com/CALFEM/calfem-matlab
    dnew=dpred+b4*anew;  
    vnew=vpred+b3*anew;    

    % Store a,v,d in respective arrays
    D(:,isnap+1) = dnew(dofhist);  % Time-History for each DOF in question
    V(:,isnap+1) = vnew(dofhist);  
    A(:,isnap+1) = anew(dofhist); 

    Dsnap(:,isnap) = dnew; % Time-History displacement for all DOF
end
D=D(:,2:nstep+1);
V=V(:,2:nstep+1);
A=A(:,2:nstep+1);

%--------------------------End--------------------------------