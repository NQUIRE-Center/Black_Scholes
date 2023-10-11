close all
clear all

%Initial Parameters 
s=150; %Maximum stock price 
K=50; %Strike 
sigma=0.2; %Volatility
r=0.04; %Interest rate
T=1; %Maturity time 
  

for Nx=8

    %Definition of operators 

    %Position 
    x=[-log(s)];  
    deltax=2*log(s)/(2^Nx-1);
    for j=1:2^Nx-1
        x(j+1)=-log(s)+j*deltax;
    end
    sv=exp(x);

    %Momentum
    P=diag(ones(2^Nx-1,1),1)+diag(-ones(2^Nx-1,1),-1);
%     P(1,2^Nx)=-1;
%     P(2^Nx,1)=1;
    P=-1*i/(2*deltax)*P;
    P2=P*P;
    
    %Initial condition WITHOUT DUPLICATION
    C0=zeros(2^Nx,1);
    for j=1:length(C0)
    
            if exp(log(1/s)+(j-1)*deltax)<K
            
                   C0(j)=K-exp(log(1/s)+(j-1)*deltax);
                   
            end
    end  

    %Evolution operators
    Unit=expm(-i*T*(sigma^2/2-r)*P);
    O=expm(-T*((sigma^2)/2*P2+r*eye(2^Nx)))   ;   
    Phif=Unit*O*C0;      
    [Call,Put] = blsprice(sv,K,r,T,sigma);


    figure
    
    h2=plot(sv(1:end),Phif(1:end));
    set([h2],'Color', '#c7c5c5','LineWidth',2)
    hold on
    
    h1=plot(sv(1:end),Put(1:end),':');
    
    set([h1],'LineWidth',2,'Color', '#6e55a5') 
    
    xlabel('Stock Value (S)','FontSize',17,'Interpreter','latex')
    ylabel('Option Price (V)','FontSize',17,'Interpreter','latex')
    legend({'8 Qubits without duplication','Analytical'},'FontSize',17,'Interpreter','latex')
    
    legend boxoff
    set(gcf, 'PaperPosition', [0 0 20 15]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [20 15]); %Set the paper to have width 5 and height 5.


end






