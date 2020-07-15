function main
global LIGHTSPEED
global sin_theta
global m_polar

LIGHTSPEED=2.997925458e8;   % m/s

%%%%%%%%  HELMHOLTZ PARAMS %%%%%%%%%%%%%%
m_polar=1;        % s-pol==1, p-pol==2
lambda=800e-9;    % in meters
theta=68.4;       % angle of incidence in deg
EDAMP=20;         % damping param for divergence-protection
nelements=1000;   % nr of 1-dimensional cells
delta=1.0;         % Cell width in nm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
format long;

theta=theta*pi/180; %deg to rad
sin2theta=complex(sin(theta)*sin(theta),0.0);
sin_theta=sin(theta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k0=1e-9*(2.0*pi/lambda); % wave-number in 1/nm
ECUT=nelements;


%1st cell needs to be vacuum
fd(1).eps=1+i*0;
fd(1).k=sqrt(fd(1).eps)*k0;

%Semi-inf vacuum for 1st and last cell
fd(nelements).delta=1e10;  
fd(1).delta=1e10;
fd(1).Te=0;       % electron temperature
fd(1).Ti=0;       % lattice/ion temperature


% prepare cells (all vacuum)
for n=2:nelements
    fd(n).rho=0;
    fd(n).eps=1+i*0;
    fd(n).k=sqrt(fd(n).eps)*k0;
    fd(n).delta=delta; 
    fd(n).Te=0;       % electron temperature
    fd(n).Ti=0;       % lattice/ion temperature
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Example for custom density-profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho0=1000; % in kg/m^3
dprof=ones(nelements,1)*0.0;
dprof(100:120)=rho0- [0:20]*0.9*rho0/20;   
dprof(121:140)=0.1*rho0 + [1:20]*0.9*rho0/20;
dprof(181:200)=rho0- [1:20]*0.9*rho0/20;   
dprof(201:220)=0.1*rho0 + [1:20]*0.9*rho0/20;
dprof(250:nelements)=500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n=1:nelements
    fd(n).rho=dprof(n);       
    fd(n).Ti=300; % in Kelvin
    fd(n).Te=300;
    fd(n).eps=getEpsilon(lambda,fd(n).Te,fd(n).Ti,fd(n).rho);
    fd(n).k=sqrt(fd(n).eps)*k0;    
    fd(n).delta=delta;
end
Q_vs_n=zeros(nelements,1); %power density vs cell-nr

%PREPARE SOLVER MATRICES
E=[1 0;0 1];
B0=1+i*0;
BR=0;

for n=1:nelements
    epsl=fd(n).eps;
    if n<nelements
        epsr=fd(n+1).eps;
    else
        epsr=1.0+i*0; %vacuum
    end
    
    kl=sqrt(epsl-sin2theta)*k0;
    kr=sqrt(epsr-sin2theta)*k0;
    phi1=real(fd(n).delta*kl*(0+i*1));
    phi2=imag(fd(n).delta*kl*(0+i*1));
    eiphi=exp(phi1)*cos(phi2)+i*exp(phi1)*sin(phi2);
    
    fd(n).P=[eiphi 1/eiphi;eiphi -1/eiphi];    
    
    if m_polar==1 %s-pol
        fd(n).Cinv=[0.5 0.5*kl/kr;0.5 -0.5*kl/kr];       
    elseif m_polar==2 %p-pol
        fd(n).Cinv=[0.5 0.5*kl/kr*epsr/epsl;0.5 -0.5*kl/kr*epsr/epsl];               
    end
    
    E=fd(n).Cinv*fd(n).P*E;      
    
    BR=-E(2,1)*B0/E(2,2);
    BT=E(1,1)*B0+E(1,2)*BR;
        
    if abs(BT)^2<exp(-EDAMP)
        ECUT=n;
        break;
    end
end    

F=[B0;BR];

for n=1:ECUT   
    F=(fd(n).Cinv*fd(n).P)*F;
    fd(n+1).Bplus=F(1);
    fd(n+1).Bminus=F(2);                
end

for n=2:ECUT
    kl=sqrt(fd(n).eps-sin2theta)*k0;    
    Q=k0*fd(n).delta*imag(fd(n).eps)*Runge5(fd(n).delta,kl,fd(n).eps,fd(n).Bplus,fd(n).Bminus);           
    Q_vs_n(n)=Q;
end

R=real(BR/B0)^2+imag(BR/B0)^2 % integral reflection
T=real(BT)^2+imag(BT)^2       % integral transmission
A=1-R-T                       % integral absorption

x=linspace(0,nelements,nelements).*delta; %absorbed power-density profile
yyaxis left
semilogy(x,Q_vs_n);
title('Absorbed power-density profile');
xlabel('x in nm');
ylabel('Relative power-density/Intensity');

yyaxis right
plot(x,[fd.rho],'--')
ylabel('Density in kg/m^3');

txt={strcat('A: ', num2str(A)), strcat('R: ', num2str(R)),strcat('T: ', num2str(T))};
text(nelements*delta*0.7, max([fd.rho])*0.85,txt);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SOLVER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval=EE(z,delta,kl,epsl,Bplus,Bminus)
global sin_theta
global m_polar

    if z<0
        disp('ERROR:z<0')
        return
    end
    if z>1
        disp('ERROR:z>1')
        return
    end
        
    im=0+i*1;
    phi1=real(delta*z*kl*im);
    phi2=imag(delta*z*kl*im);
    eiphi=exp(phi1)*cos(phi2)+i*exp(phi1)*sin(phi2);
    if m_polar==2 %p-polar
        Ez=sin_theta/epsl*(Bplus*eiphi+Bminus/eiphi);
        Ey=-sqrt(epsl-sin_theta*sin_theta)/epsl*(Bplus*eiphi-Bminus/eiphi);
    elseif m_polar==1 %s-polar
        Ez=Bplus*eiphi+Bminus/eiphi;        
        Ey=0;
    end    
    retval=abs(Ez)^2+abs(Ey)^2;    
end

% Adaptive runge-kutta algorithm
function retval=Runge5(delta,kl,epsl,Bplus,Bminus)
    dz=1;
    cur_zpos=0;
    result=0;
    errval=1e-5;
    k1=0;
    k2=0;
    k3=0;
    k4=0;
    k5=0;
    while (cur_zpos<1 & dz> 1e-5)
        
        k1 = dz / 3.0 * EE(cur_zpos, delta, kl, epsl, Bplus, Bminus);
        k3 = dz / 3.0 * EE(cur_zpos + dz / 3.0, delta, kl, epsl, Bplus, Bminus);
        k4 = dz / 3.0 * EE(cur_zpos + 0.5 * dz, delta, kl, epsl, Bplus, Bminus);
        k5 = dz / 3.0 * EE(cur_zpos + dz, delta, kl, epsl, Bplus, Bminus);
        ERR = k1 - 4.5 * k3 + 4 * k4 - 0.5 * k5;
        
        if ERR<5/32*errval
            result=result+0.5*(k1+4*k4+k5);
            cur_zpos=cur_zpos+dz;
            dz=dz*1.1;
        elseif ERR>5*errval
            dz=dz*0.5;
            continue;
        else
            result=result+0.5*(k1+4*k4+k5);
            cur_zpos=cur_zpos+dz;
        end
        
        if cur_zpos+dz>1
            dz=1-cur_zpos;
        end
    end
    retval=result;            
end


% Function for relative permittivity calculation
% Customize this for density-,temperature-, and wavelength-dependent permittivity
function retval=getEpsilon(lambda,Te,Ti,rho)
    if rho>0
        retval=-111 + i*52.6;
    else
        retval=1.0;
    end
end

