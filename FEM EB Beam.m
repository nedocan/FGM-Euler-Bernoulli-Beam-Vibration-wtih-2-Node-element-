clc
clear
tic
syms A11 B11 D11 I0 I1 I2 x L E A z w real 

%%% Shape Functions %%%
N1=1-x/L;
N2=x/L;
N3=1-3*x^2/L^2+2*x^3/L^3;
N4=x-2*x^2/L+x^3/L^2;
N5=3*x^2/L^2-2*x^3/L^3;
N6=-x^2/L+x^3/L^2;

N_u=[N1 N2];
N_w=[N3 N4 N5 N6];

M11=I0*int(N_u'*N_u,x,[0,L]);
M12=-I1*int(N_u'*diff(N_w,x),x,[0,L]);
M21=-I1*int(diff(N_w',x)*N_u,x,[0,L]);
M22=I0*int(N_w'*N_w,[0,L])+I2*int(diff(N_w',x)*diff(N_w,x),x,[0,L]);
M=[M11 M12;M21 M22]; %% Local mass matrix (symbolic)

K11=A11*int(diff(N_u',x)*diff(N_u,x),x,[0,L]);
K12=-B11*int(diff(N_u',x)*diff(N_w,x,x),x,[0,L]);
K21=-B11*int(diff(N_w',x,x)*diff(N_u,x),x,[0,L]);
K22=D11*int(diff(N_w',x,x)*diff(N_w,x,x),x,[0,L]);
K=[K11 K12;K21 K22]; %% Local stiffness matrix (symbolic)


K([1 2 3 4 5 6],[1 2 3 4 5 6])=K([1 3 4 2 5 6],[1 3 4 2 5 6]);
M([1 2 3 4 5 6],[1 2 3 4 5 6])=M([1 3 4 2 5 6],[1 3 4 2 5 6]);
format long
%%%% Geometry %%%%%
h=  0.1; %height of beam
b=   0.05; %width of beam
l=   2;%length of beam  

%%% Shear Considerations %%%
f=0 ;   %Euler Bernoulli Beam

%%%% Material Properties %%%%%
k=2;  %volume fraction

E_metallic=210*10^9;
rho_metallic=7800;
nu_metallic=0.3;

E_ceramic=380*10^9;
rho_ceramic=3800;
nu_ceramic=0.14;


%%%%%  Case1; upper face is ceramic; lower face is metallic %%%%%
E_effective=(E_ceramic-E_metallic)*(z/h+0.5)^k+E_metallic;
rho_effective=(rho_ceramic-rho_metallic)*(z/h+0.5)^k+rho_metallic;
nu_effective=(nu_ceramic-nu_metallic)*(z/h+0.5)^k+nu_metallic;


%{
%%%%%  Case2; upper face is metallic; lower face is ceramic %%%%%
E_effective=(-E_ceramic+E_metallic)*(z/h+0.5)^k+E_ceramic;
rho_effective=(-rho_ceramic+rho_metallic)*(z/h+0.5)^k+rho_ceramic;
nu_effective=(-nu_ceramic+nu_metallic)*(z/h+0.5)^k+nu_ceramic;
%}


%%% Inertial Properties %%%
I0_numeric=double(int(rho_effective,[-h/2 h/2])*b);
I1_numeric=double(int(z*rho_effective,[-h/2 h/2])*b);
I2_numeric=double(int(z^2*rho_effective,[-h/2 h/2])*b);

I=[I0_numeric,I1_numeric,I2_numeric] ;%% Ineartias

%%%Coefficients of Stress Resultants%%%
A11_numeric=double(int(E_effective,[-h/2 h/2])*b);
B11_numeric=double(int(z*E_effective,[-h/2 h/2])*b);
D11_numeric=double(int(z^2*E_effective,[-h/2 h/2])*b);

R=[A11_numeric,B11_numeric,D11_numeric]; %% Coefficinet of Stress Resultants

nel=40;  %%Number of Divided Elements

mm=double(subs(M,{A11,B11,D11,I0,I1,I2,L},{R(1),R(2),R(3),I(1),I(2),I(3),l/nel})); %% Local mass matrix
kk=double(subs(K,{A11,B11,D11,I0,I1,I2,L},{R(1),R(2),R(3),I(1),I(2),I(3),l/nel})); %% Local stiffness matrix

n_nodes=nel+1;  %%Node Number
Dof=3;  %%Degree of Freedom of a Node

%%% Boundary Conditions %%%
BCs=[1,2,3,n_nodes*Dof-2,n_nodes*Dof-1,n_nodes*Dof]; % Clamped-Clamped
%BCs=[1,2,n_nodes*Dof-2,n_nodes*Dof-1]; % Simply Supported-Simply Supported
%BCs=[1,2,3,n_nodes*Dof-2,n_nodes*Dof-1]; %Clamped-Supported
%BCs=[1,2,3]; % Clamped-Free



           
Kglob=zeros(Dof*n_nodes,Dof*n_nodes)  ;
Mglob=zeros(Dof*n_nodes,Dof*n_nodes)  ;

for i=1:nel
    top(i,:)=[3*i-2 3*i-1 3*i 3*i+1 3*i+2 3*i+3]; %% Topology Matrix
end


%%% Assembly Process %%%

for n=1:nel
    for i=1:6
        for j=1:6
            Kglob(top(n,i),top(n,j))=...
                Kglob(top(n,i),top(n,j))+kk(i,j);
        end
    end
end

for n=1:nel
    for i=1:6
        for j=1:6
            Mglob(top(n,i),top(n,j))=...
                Mglob(top(n,i),top(n,j))+mm(i,j);
        end
    end
end

activeDof=setdiff(1:Dof*n_nodes',BCs); %% Rows and columns are omitted which are boundary conditions

Kglobactive=Kglob(activeDof,activeDof);
Mglobactive=Mglob(activeDof,activeDof);

[vectorc,freqc]=eig(Kglobactive,Mglobactive); %% Eigenvalue problem
nat_freq=freqc;
freqc=diag(freqc);
freqc=sort(sqrt(freqc)) %% Natural Frequncies

disp=linspace(0,l,n_nodes);  %% Displacement along x-direction

%%% Mass normalisation of mass and stiffness matrices %%%

mod=[];
for i=1:length(freqc)
    v=vectorc(:,i)/norm(vectorc(:,i));
    v=v/sqrt(v'*Mglobactive*v);
    mod=[mod,v];
end


%%% Mode Shapes %%%

number_of_mode_shapes=15; %% Number of modes which are requested


%%%%%%%   Clamped-Clamped Boundary Condition   %%%%%%%%%

mod_shapes=[zeros(3,length(freqc));mod;zeros(3,length(freqc))];

for i=1:n_nodes
    for j=1:length(freqc);
      u_mod(i,j)=mod_shapes(3*i-2,j);
    end
end

for i=1:n_nodes
    for j=1:length(freqc);
      w_mod(i,j)=mod_shapes(3*i-1,j);
    end
end

for i=1:number_of_mode_shapes
figure
subplot(2,1,1)
plot(disp,u_mod(:,i))
grid on
title('Clamped-Clamped')
xlabel('x')
ylabel('u')
subplot(2,1,2)
plot(disp,w_mod(:,i))
grid on
xlabel('x')
ylabel('w')
end




%{



%%%%%%%   Simply Supported-Simply Supported Boundary Condition   %%%%%%%%%

mod_shapes=[zeros(2,length(freqc));mod(1:length(freqc)-1,:);zeros(2,length(freqc));mod(length(freqc),:)];

for i=1:n_nodes
    for j=1:length(freqc);
      u_mod(i,j)=mod_shapes(3*i-2,j);
    end
end

for i=1:n_nodes
    for j=1:length(freqc);
      w_mod(i,j)=mod_shapes(3*i-1,j);
    end
end


for i=1:number_of_mode_shapes
figure
subplot(2,1,1)
plot(disp,u_mod(:,i))
grid on
title('Simply Supported-Simply Supported')
xlabel('x')
ylabel('u')
subplot(2,1,2)
plot(disp,w_mod(:,i))
grid on
xlabel('x')
ylabel('w')
end








%%%%%%%   Clamped-Free Boundary Condition   %%%%%%%%%

mod_shapes=[zeros(3,length(freqc));mod];

for i=1:n_nodes
    for j=1:length(freqc);
      u_mod(i,j)=mod_shapes(3*i-2,j);
    end
end

for i=1:n_nodes
    for j=1:length(freqc);
      w_mod(i,j)=mod_shapes(3*i-1,j);
    end
end


for i=1:number_of_mode_shapes
figure
subplot(2,1,1)
plot(disp,u_mod(:,i))
grid on
title('Clamped-Free')
xlabel('x')
ylabel('u')
subplot(2,1,2)
plot(disp,w_mod(:,i))
grid on
xlabel('x')
ylabel('w')
end




%%%%%%%   Clamped-Simply Supported Boundary Condition   %%%%%%%%%

mod_shapes=[zeros(3,length(freqc));mod(1:length(freqc)-1,:);zeros(2,length(freqc));mod(length(freqc),:)];
number_of_mode_shapes=15
for i=1:n_nodes
    for j=1:length(freqc);
      u_mod(i,j)=mod_shapes(3*i-2,j);
    end
end

for i=1:n_nodes
    for j=1:length(freqc);
      w_mod(i,j)=mod_shapes(3*i-1,j);
    end
end

for i=1:number_of_mode_shapes
figure
subplot(2,1,1)
plot(disp,u_mod(:,i))
grid on
title('Clamped-Simpþy Supported')
xlabel('x')
ylabel('u')
subplot(2,1,2)
plot(disp,w_mod(:,i))
grid on
xlabel('x')
ylabel('w')
end
%}



%%% Frequency Response Function %%%

i=140; %% Measured point
j=400; %% Excitation point
frf_modal=0;
for k=1:length(freqc)
frf_modal=mod(i,k)*mod(j,k)/(-w^2+nat_freq(k,k))+frf_modal;
end

frequency=0:5:freqc(10)+10000;
figure
plot(frequency,20*log(double(subs(frf_modal,frequency))))
title('Receptance')
xlabel('Frequency')
ylabel('Amplitude (dB)')
grid on

toc
