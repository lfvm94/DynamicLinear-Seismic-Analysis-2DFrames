function [FmaxMDOF,Te,lambda,fi,Ma]=ModalsMDOF2DFrames2(M,K,bc,Sa,mode)
% SYNTAX : [FmaxMDOF,Te,lambda,fi,Ma]=ModalsMDOF2DFrames2(M,K,bc,mode)
%---------------------------------------------------------------------
%    PURPOSE
%     To compute the equivalent inertial forces at the DOF's of a 
%     plane frame due to an acceleration at its base.
% 
%    INPUT:  M:                 Global Mass matrix
%            K:                 Global Stiffness matrix
%
%            bc:                Boundary condition array
%
%            Sa:                Design Spectrum 
%                               (1) from the Eurocode EC8
%
%            mode:              Mode of vibration of interest:
%                               [mode-1,mode-2,...] -> The equivalent 
%                                    inertial forces are computed with the
%                                    contribution of the given modes of 
%                                    vibration
%                               
%                               mode-i -> The equivalent inertial forces are
%                                       computed with the mode of
%                                       vibration inserted
%
%    OUTPUT: lambda :           Modal of vibration for each DOF. 
%                               Size: Nmodals x 1
%
%            fi :               DOF's eigenvalues: NDOF x Nmodals
%
%            Te :               Structure's periods for each modal      
%
%            FmaxMDOF :         Equivalent DOF's forces for the modal
%                               in question
%
%            Ma:                Modal mass contribution for every DOF
%
%--------------------------------------------------------------------

% LAST MODIFIED: L.Verduzco    2023-06-07
% Copyright (c)  School of Engineering
%                The Hong Kong University of Science and Technology
%--------------------------------------------------------------------

%% Solving eigenvalues (frequencies) and eigenvectors (modals)
[lambda,fi]=eigen(K,M,bc(:,1)); % Eigenvalue/eigenvectors

[ndof,nmodes]=size(fi);

for i=1:nmodes
    [factor,ifactor]=max(abs(fi(:,i))); % Eigenvectors - vibration modals
    factor=factor*sign(fi(ifactor,i));
    fi(:,i)=fi(:,i)./factor;
end

% Circular frequencies
omega=sqrt(lambda);

% Frequencies
freq=omega/(2*pi);

% Periods
Te=1./freq;

%% Lateral equivalent inertial loads caused by the soil acceleration
fmax=zeros(ndof,nmodes);
for i=1:nmodes
    Mn=fi(:,i)'*M*fi(:,i);
    fmaxn=fi(:,i)'*M;
    vector1=abs(fi(:,i)'); % influence vector
    
    vector1(1,3:3:ndof)=0; % the roation DOF are not considered influential
    vector1=vector1/max(vector1);
    
    Ln=fmaxn*vector1'; 
    rn=Ln/Mn; % Modal participation factor
    
    fmax(:,i)=-fmaxn.*vector1*Sa; % fmax(:,i)=rn*Sd;
    
    Ma(i)=rn^2*Mn; % Effective modal mass
end

% Lateral equivalent inertial loads considering the constribution of all 
% required modes 
npmodes=length(mode);
FmaxMDOF=zeros(ndof,1);
for i=1:ndof
    for j=1:npmodes
        FmaxMDOF(i,1)=FmaxMDOF(i,1)+fmax(i,mode(j))^2;
    end
end

if npmodes>1 % To consider the contribution of all modals
    FmaxMDOF=sqrt(FmaxMDOF);
else
    FmaxMDOF=fmax(:,mode); % When only a certain modal is taken
end