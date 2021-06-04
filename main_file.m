close all;
clc;
lower_bound=[0 0 1 1 40];
upper_bound=[0.8 0.5 2 2 105];
Np=90;
T=300;
t=1;
Pc=0.3;
F=0.42;
f=NaN(Np,1);
fu=NaN(Np,1);
D=length(lower_bound);
U=NaN(Np,D);
P= repmat(lower_bound,Np,1)+repmat((lower_bound-upper_bound),Np,1).*rand(Np,D);
P(1)=0.6;
prob=@fitnessfunction;
for p=1:Np
    f(p)= prob(P(p,:));
end
for T = 1:T
        for i = 1:Np
        %% Mutation
        Candidates = [1:i-1 i+1:Np];
        idx = Candidates(randperm(Np-1,3));      %selection of random parameters
        X1 =  P(idx(1),:);
        X2 = P(idx(2),:);
        X3 = P(idx(3),:);
        V = X1 + F*(X2-X3);                 %mutuant vector generated
        
        %% Cross-Over
        del = randi(D,1);                    %generate random variable data i
        for j = 1:D
            
            if( rand <=Pc || del==j)
                U(i,j)=V(j);                      %Mutuant vector assign
            else
                U(i,j)=P(i,j);                     % Target vector assign
            end
        end
        end
      %% Greedy Selection
      for j = 1:Np
          U(j,:) == min(upper_bound,U(j,:));
          U(j,:) == max(lower_bound,U(j,:));
          fu(j) = prob(U(j,:));                               %fitness function
          if fu(j) < f(j)
              P(j,:) = U(j,:);
              f(j) = fu(j);
          end
      end;
      best(t,1)=min(f);
      t=t+1;
end
[bestfitness,ind] = min(f);
bestac1=P(ind,:);
bestac
t=1:T;
plot(t,best,'linewidth',2);
title('Fitness Convergence');
figure;
Rs=bestac(1);
R=bestac(2);
Xs=bestac(3);
X=bestac(4);
Xm=bestac(5);
n=1;
for s=0:0.01:1
   T(n,1)=s*(400^2)*R*((R^2+(s*X)^2)^-1)*3*((2*pi*25)^-1);
   T1(n,1)=s*(400^2)*0.23*((0.23^2+(s*1.27)^2)^-1)*3*((2*pi*25)^-1);
   n=n+1;
end
s=1:-0.01:0;
plot(s,T,'--.r');
hold on
plot(s,T1,'linewidth',2);
xticks=[1 0.8 0.6 0.4 0.2 0];
xlabel("Slip");
ylabel("Torque of the induction machine");
title("Torque-Slip Charecteristics of induction motor");
legend('predicted-value','actual-value');
grid on;
figure;
s=0;
x=1;
while s<0.0451
    t(x,1)=(s*400*400*3*bestac(2))/((2*pi*25)*((bestac(2)*bestac(2))+(s*s*bestac(4)*bestac(4))));
    w(x,1)=(4*pi*50*(1-s)/4);
    sl(x,1)=s;
    prctsl(x,1)=s*100/0.0451;
    
  
  
    z(x,1)=0.72+2.155i+(((0.37/s)+2.155i)*(96.805i)/((0.37/s)+96.805i));
   
    current(x,1)=400/(sqrt(3)*sqrt(((real(z(x,1)))^2)+(imag(z(x,1)))^2));
    prctcurrent(x,1)=current(x,1)*100/(13.26);
    
    pf(x,1)=cos(atan(imag(z(x,1))/(real(z(x,1)))));
    prctpf(x,1)=pf(x,1)*100/(0.8676);
   
    power(x,1)=(t(x,1)*w(x,1));
    prctpower(x,1)=power(x,1)*100/(93877);
    
    
    op(x,1)=current(x,1)*400*pf(x,1)*sqrt(3);
    prctop(x,1)=op(x,1)/80; 
    inputpower(x,1)=op(x,1)+(3*current(x,1)*current(x,1)*(2.155+0.72));
    eff(x,1)=op(x,1)*100/inputpower(x,1);
    prcteff(x,1)=eff(x,1)*100/75.16;
    
    
    x=x+1;
    s=s+0.0001;
end
plot(prctop,prctsl,prctop,prctcurrent,prctop,prcteff,prctop,prctpf,'linewidth',1.5)
legend('Slip','Current','Efficiency','Power Factor')
title('Performance charecteristics of induction motor')
hold on;
s=0;
x=1;
while s<0.0451
    t(x,1)=(s*400*400*3*0.23)/((2*pi*25)*((0.23*0.23)+(s*s*1.27*1.27)));
    w(x,1)=(4*pi*50*(1-s)/4);
    sl(x,1)=s;
    prctsl(x,1)=s*100/0.0451;
    
  
  
    z(x,1)=0.72+2.155i+(((0.37/s)+2.155i)*(96.805i)/((0.37/s)+96.805i));
   
    current(x,1)=400/(sqrt(3)*sqrt(((real(z(x,1)))^2)+(imag(z(x,1)))^2));
    prctcurrent(x,1)=current(x,1)*100/(13.26);
    
    pf(x,1)=cos(atan(imag(z(x,1))/(real(z(x,1)))));
    prctpf(x,1)=pf(x,1)*100/(0.8676);
   
    power(x,1)=(t(x,1)*w(x,1));
    prctpower(x,1)=power(x,1)*100/(93877);
    
    
    op(x,1)=current(x,1)*400*pf(x,1)*sqrt(3);
    prctop(x,1)=op(x,1)/80; 
    inputpower(x,1)=op(x,1)+(3*current(x,1)*current(x,1)*(2.155+0.72));
    eff(x,1)=op(x,1)*100/inputpower(x,1);
    prcteff(x,1)=eff(x,1)*100/75.16;
    
    
    x=x+1;
    s=s+0.0001;
end
plot(prctop,prctsl,prctop,prctcurrent,prctop,prcteff,prctop,prctpf,'linewidth',1)
grid



          
            
        
        
        
