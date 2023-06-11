clc
clear all

% Important Note:
% This code Tessted on MAtlab R2021a successfully


x0=[0.629610834937804;0.239207893767559;0.00199999999999997;-2.68855075384503e-21;0.00853182278410227;145495.000013079;-13991.3912640677;0.0100000000008989;-1.54606158016281e-16;15107.6087385481;85.0721740273458;1510760.87371901;8.84523538757923e-06;7.41913341880677e-06;0.00100000000234248;5.56849995061363e-06;1.77453877292112e-06;1.00000000234248;0.199999999999999;-2.43777028428438e-22;634205.429840439;871.715546812033;2641247.42473264;9000;350;0;0;0.01]
xd0=zeros(28,1)
xd0(7)=[1392.95080573910];
xd0(26)=[0.000999952490545232];
xd0(27)=[0.000999984859949650];
xd0(28)=[-0.00853182280408788];

 Tspan=[0 60]
F_15i(0,x0,xd0)
options= odeset('abstol',5e-6)
[t,x] = ode15i (@F_15i,Tspan,x0, xd0,options)
 
figure;
plot(t,x(:,13),'LineWidth',2)
hold on
plot(t,x(:,14),'LineWidth',2)
plot(t,x(:,15))
title('X Mol fraction')
legend('x1','x2','x3')

figure;
plot(t,x(:,16),'LineWidth',2)
hold on
plot(t,x(:,17),'LineWidth',2)
plot(t,x(:,18),'LineWidth',2)
title('Y Mol fraction')
legend('y1','y2','y3')

figure;
yyaxis right
plot(t,x(:,25),'LineWidth',2)
yyaxis left
plot(t,x(:,6),'LineWidth',2)
legend('P','T')

figure;
yyaxis right
plot(t,x(:,5),'LineWidth',2)
yyaxis left
plot(t,x(:,4),'LineWidth',2)
legend('FL','Fv')

figure;
yyaxis right
plot(t,x(:,8),'LineWidth',2)
yyaxis left
plot(t,x(:,9),'LineWidth',2)
legend('ML','Mv')


%%%%%%%%%%%%% Phase Detector

Mv=x(:,8);
ML=x(:,9);
x1=x(:,13);
x2=x(:,14);
x3=x(:,15);
y1=x(:,16);
y2=x(:,17);
y3=x(:,18);
n=length(t);
State=[];
for i=1:1:n
    temp=median([Mv(i)/(ML(i)+Mv(i)),(x1(i)+x2(i)+x3(i))-(y1(i)+y2(i)+y3(i)),Mv(i)/(ML(i)+Mv(i))-1]);
    switch temp
        case Mv(i)/(ML(i)+Mv(i))
            State(i)=1;
        case (x1(i)+x2(i)+x3(i))-(y1(i)+y2(i)+y3(i))
            State(i)=2;
        case Mv(i)/(ML(i)+Mv(i))-1
             State(i)=3;
        otherwise
            display('Error happend')
            
    end
            
    
end



