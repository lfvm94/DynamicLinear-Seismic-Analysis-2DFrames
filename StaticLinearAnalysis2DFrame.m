function [displacements,reactions,Ex,Ey,esbarsnormal,esbarsshear,esbarsmoment]=...
    StaticLinearAnalysis2DFrame(Eelem,areaElem,inertia,bc,...
    fglobal,ni,nf,qbarxy,np,coordxy,plotAnalysisResults,eRefN,eRefV,eRefM)

%------------------------------------------------------------------------
% Syntax:
% Option 1:
% [displacements,reactions,Ex,Ey,esbarsnormal,esbarsshear,esbarsmoment]=...
%   StaticLinearAnalysis2DFrame(Eelem,areaElem,inertia,bc,...
%   fglobal,ni,nf,qbarxy,np,coordxy,plotAnalysisResults,eRefN,eRefV,eRefM)
%
% Option 2:
% [displacements,reactions,Ex,Ey,esbarsnormal,esbarsshear,esbarsmoment]=...
%   StaticLinearAnalysis2DFrame(Eelem,areaElem,inertia,bc,...
%   fglobal,ni,nf,qbarxy,np,coordxy,plotAnalysisResults)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%------------------------------------------------------------------------
% PURPOSE: To execute a linear static analysis of plane frame, and compute
% the mechanic element diagrams for each  of its elements.
% 
% OUTPUT: displacements:        is a vector containing the EDOF displacements.
%                               Size: [nEDOF,1]
%
%         reactions:            is the EDOF reaction forces. 
%                               Size = [nEDOF,1]
%
%         Ex,Ey:                are the element end nodes' coordinates
%
%         esbarsnormal:         are the normal mechanic elements for each
%                               bar. Size: [np,nbars]
%
%         esbarsshear:          are the shear mechanic elements for each
%                               bar. Size: [np,nbars]
%
%         esbarsmoment:         are the bending moment mechanic elements 
%                               for each bar. Size: [np,nbars]
%
% INPUT:  nnodes:               are the number of nodes of the structure
%
%         nbars:                are the number of structural elements of 
%                               the structural frame (neglecting the 
%                               footings)
%
%         Eelem:                is the vector containing the elasticity
%                               modulus for each bar
%
%         areaElem:             is the vector containing the cross-section
%                               area of each element
%
%         inertia:              is the vector containing the cross-section
%                               inertia momentum of each element
%
%         bc:                   is the array containing the boundary 
%                               conditions for the respective prescribed 
%                               (or restricted) DOF. Size=[nRestrictedDOF,2]
%                               in format [DOF,prescribed-displacement]
%
%         fglobal:              is the global external force vector for
%                               each DOF. Size: [NDOF,1]
%
%         ni,nf:                are the vectors containing the initial
%                               and final nodes for each element. Size: 
%                               [nbars,1] for each
%
%         qbarxy:               are the distributed loads on the
%                               elements. Size: nbars x 2. The 1st column 
%                               is for the loads in the local x
%                               direction and 2nd column for the loads in
%                               the local y direction of each element.
%
%         Edof:                 is the topology matrix, consisting of the
%                               DOF for each element. The DOF of the initial
%                               node are placed first as: 
%           ________________________________________________________
%           [nbar,DOFxi,DOFyi,DOF(theta-i),DOFxf,DOFyf,DOF(theta-f)]
%           --------------------------------------------------------
%
%         np:                   are the number of points for the evaluation
%                               of the mechanic elements for each structural 
%                               element
%
%         coordxy:              is the array containing the node coordinates.
%                               Size = [nNodes,2] in format [xi,yi]
%
%         plotAnalysisResults:  is the parameter that indicates if the
%                               structural analysis results should be 
%                               plotted or not. Options are: (1) the plots
%                               are required, (2) the plots are not required
%
%         eRefN,eRefV,eRefM;    These are the structural elements of 
%                               reference to scale the diagrams, 
%                               respectively for axial loads (N), shear 
%                               forces (V) and bending moments (M) 
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------
nnodes=length(coordxy(:,1)); nbars=length(Eelem);
% Topology matrix 
% Format:
% ______________________________________________________________
% [nbar, DOFxi, DOFyi, DOF(theta-i), DOFxf, DOFyf, DOF(theta-f)]
% --------------------------------------------------------------
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

Kglobal=zeros(3*nnodes);
 ep_bars=zeros(nbars,3); 
 eq_bars=zeros(nbars,2);
 for i=1:nbars      
     
     ex=[coordxy(ni(i),1) coordxy(nf(i),1)];
     ey=[coordxy(ni(i),2) coordxy(nf(i),2)];
     ep=[Eelem(i) areaElem(i) inertia(i)];
     eq=qbarxy(i,:);
     
     ep_bars(i,:)=ep;
     eq_bars(i,:)=eq;
     [Kebar,febar]=beam2e(ex,ey,ep,eq); % CALFEM function
                                          % To download the CALFEM package
                                          % visit its repository at:
                                          % https://github.com/CALFEM/calfem-matlab
    
     [Kglobal,fglobal]=assem(Edof(i,:),Kglobal,Kebar,fglobal,febar);
     
 end     

[displacements,reactions]=solveq(Kglobal,fglobal,bc);
Ed=extract(Edof,displacements);
ex=coordxy(:,1);
ey=coordxy(:,2);

Ex=zeros(nbars,2);
Ey=zeros(nbars,2);

for j=1:nbars
    Ex(j,1)=ex(Edof(j,4)/3);
    Ex(j,2)=ex(Edof(j,7)/3);

    Ey(j,1)=ey(Edof(j,4)/3);
    Ey(j,2)=ey(Edof(j,7)/3);

end

%% Forces diagrams
esbarsnormal=zeros(np,nbars);
esbarsshear=zeros(np,nbars);
esbarsmoment=zeros(np,nbars);
for i=1:nbars
    es_bar=beam2s(Ex(i,:),Ey(i,:),ep_bars(i,:),Ed(i,:),eq_bars(i,:),np);
    esbarsnormal(:,i)=es_bar(:,1);
    esbarsshear(:,i)=es_bar(:,2);
    esbarsmoment(:,i)=es_bar(:,3);
end

if plotAnalysisResults==1
     scstruc=50;
     scstructxt=num2str(scstruc);
     %-----Undeformed mesh-----%
     figure(1)
     xlabel('Width')
     ylabel('Height')
     title(strcat('Deformed structure. Scale: 1-',scstructxt));
     plotpar=[2 1 0];
     eldraw2(Ex,Ey,plotpar);

       %----Deformed mesh----%

     plotpar=[1 2 1];
     eldisp2(Ex,Ey,Ed,plotpar,100);  
     sfac=scalfact2(Ex(eRefN,:),Ey(eRefN,:),esbarsnormal(:,eRefN),0.2);
     for i=1:nbars
         figure(2)
         plotpar=[2 1];
         eldia2(Ex(i,:),Ey(i,:),esbarsnormal(:,i),plotpar,sfac);
         title('Axial force')
     end
     sfac=scalfact2(Ex(eRefV,:),Ey(eRefV,:),esbarsshear(:,eRefV),0.2);
     for i=1:nbars
         figure(3)
         plotpar=[2 1];
         eldia2(Ex(i,:),Ey(i,:),esbarsshear(:,i),plotpar,sfac);
         title('Shear force') 
    end

    sfac=scalfact2(Ex(eRefM,:),Ey(eRefM,:),esbarsmoment(:,eRefM),0.2);

    for i=1:nbars
         figure(4)
         plotpar=[2 1];
         eldia2(Ex(i,:),Ey(i,:),esbarsmoment(:,i),plotpar,sfac);
         title('Bending moment')
    end
end