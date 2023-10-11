close all
clear all
clear all 


%Initial Parameters 
s=150; %Maximum stock price 
K=50; %Strike 
sigma=0.2; %Volatility
r=0.04; %Interest rate
T=1; %Maturity time 
 

for Nx=[4,5,6,7,8,10] %Nx is the number of qubits 

    %Definition of operators 

    %Position 
    xmax=2*log(s);
    deltax=4*log(s)/(2^Nx-1);
    x=[-2*log(s)];  
    for j=1:2^Nx-1    
        x(j+1)=-2*log(s)+j*deltax;
    end
    sv=exp(x+xmax/2);

    %Momentum
    P=diag(ones(2^Nx-1,1),1)+diag(-ones(2^Nx-1,1),-1);
    P(1,2^Nx)=-1; %Periodic Boundary condition
    P(2^Nx,1)=1; %Periodic Boundary condition
    P=-1i/(2*deltax)*P;

    %Initial condition
    C0=zeros(2^(Nx-1),1);
    for j=1:length(C0)
        if exp(log(1/s)+(j-1)*deltax)<K
           C0(j)=K-exp(log(1/s)+(j-1)*deltax);
        end
    end  
    C0=[C0;flip(C0)];


    %Evolution operators
    Unit=expm(-i*T*(sigma^2/2-r)*P);
    O=expm(-T*((sigma^2)/2*P*P+r*eye(2^Nx))) ;  

    %Evolved state
    Phif=Unit*O*C0;
    Phif=Phif(1:length(Phif)/2);
    sv=sv(1:length(Phif));
    x=x(1:length(Phif));
    C0=C0(1:length(Phif));

    %Analytical Solution
    [Call,Put] = blsprice(sv(1:length(Phif)),K,r,T,sigma);
    
    if Nx==10

        [Call,Put] = blsprice(sv,K,r,T,sigma);
        h1=plot(sv,Put,'-');
        axis ([0 135 0 40]);
        hold on
        set([h1],'LineWidth',2);
        set(gca, 'fontsize', 12);
       
        xlabel('Stock Value (S)','FontSize',17,'Interpreter','latex');
        ylabel('Option Price','FontSize',17,'Interpreter','latex');
      
    
    end
   
    
    
    if Nx~=10
    
        h2=plot(sv,Phif,'--');
        hold on
        axis ([0 135 0 40]);
        set([h2],'LineWidth',2);
        set(gca, 'fontsize', 12);
        xlabel('Log Stock Value (x=Log(S))','FontSize',17,'Interpreter','latex')
        ylabel('Initial Option Price','FontSize',17,'Interpreter','latex')
  
    end


end
legend({'4 Qubits','5 Qubits','6 Qubits','7 Qubits','8 Qubits','Analytical'},'Interpreter','latex')



















