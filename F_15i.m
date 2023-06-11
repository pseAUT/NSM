function F= F_15i(t,z,zp)
VT=0.2;
Fin0=0.002;
Q0=9000;
tsd= 20;
cv= 0.01;
cL=0.141;
Vvmin=0.02*0.2;
VLmin=0.01*0.2;
Tin=325;
A=0.1;
P0=1e5;
Rho1=11.22;
Rho2=9.41;
R= 8314;
Mw1=78.11;
Mw2=92.14;
z1=0.5;
z2= 0.5;
z3=0;
g=9.81;
Tref=298.15;
epsi=1e-6;
k3=1000;
dHform=[82.88e3,50.17e3,0 ];	    
alpha=[3.551,3.866,3.539;-6.184e-3,3.558e-3,-0.261e-3;14.365e-5,13.356e-5,0.007e-5;-19.807e-8,-18.659e-8,0.157e-8;8.234e-11, 7.690e-11,-0.099e-11];
Tc=[562.16,591.80];		   
a=[-7.01433 ,-7.31600];          
b=[1.55256 ,1.59425];             
c=[-1.8479 ,-1.93165];            
d=[-3.7130 , -3.72220];      	   
Pc=[48.98 , 41.06];	

k1=z(1);
k2=z(2);
Fin=z(3);
FL=z(4);
Fv=z(5);
P=z(6);
U=z(7);
Mv=z(8);
ML=z(9);
H=z(10);
hL=z(11);
hv=z(12);
x1=z(13);
x2=z(14);
x3=z(15);
y1=z(16);
y2=z(17);
y3=z(18);
Vv=z(19);
VL=z(20);
RhoL=z(21);
RhoPL=z(22);
hin=z(23);
Q=z(24);
T=z(25);
M1=z(26);
M2=z(27);
M3=z(28);

du=zp(7);
dm1=zp(26);
dm2=zp(27);
dm3=zp(28);



AA=[Mv/(ML+Mv), (x1+x2+x3)-(y1+y2+y3), Mv/(ML+Mv)-1];
tau = (Tc-T)./Tc;
Pvp=1e5*exp((log(Pc) + Tc ./ T .* (a .* tau + b .* tau .^1.5 + c .* tau .^2.5 + d .* tau .^5) ) );
dHvap(1) = 1e3*47.41*exp(-0.1231*(1-tau(1)))*(tau(1))^0.3602  ;

dHvap(2) = 1e3*53.09*exp(-0.2774*(1-tau(2)))*(tau(2))^0.2774  ;	
dHTref = R* ( alpha(1,1:3).*(T-Tref) + alpha(2,1:3)./2.*(T^2-Tref^2) + alpha(3,1:3)./3.*(T^3-Tref^3) +alpha(4,1:3)./4*(T^4-Tref^4) + alpha(5,1:3)./5.*(T^5-Tref^5) ) ;

h_Vcmp = dHform + dHTref;
h_Lcmp = h_Vcmp(1:2) - dHvap ;


dHTref_in = R * ( alpha(1,1:2)*(Tin-Tref) + alpha(2,1:2)/2*(Tin^2-Tref^2) + alpha(3,1:2)/3*(Tin^3-Tref^3) +alpha(4,1:2)/4*(Tin^4-Tref^4) + alpha(5,1:2)/5*(Tin^5-Tref^5) );		
dHvap_in(1) = 1e3*47.41*exp(-0.1231*(Tin/Tc(1)))*(1 - Tin/Tc(1))^0.3602  ;
dHvap_in(2) = 1e3*53.09*exp(-0.2774*(Tin/Tc(2)))*(1 - Tin/Tc(2))^0.2774  ;	
h_Vcmp_in = dHform(1:2) + dHTref_in;
h_Lcmp_in = h_Vcmp_in - dHvap_in;

F=[x1*k1-y1;
x2*k2-y2;
x3*k3-y3;
x1*ML+y1*Mv-M1;
x2*ML+y2*Mv-M2;
x3*ML+y3*Mv-M3;
P*Vv-Mv*R*T;
Mv/(ML+Mv)+(x1+x2+x3)-(y1+y2+y3)+Mv/(ML+Mv)-1-min(AA)-max(AA);
k1*P-Pvp(1);
k2*P-Pvp(2);
M1+M2+M3-Mv-ML;
VL+Vv-VT;
VL*RhoL-ML;
(x1/Rho1+x2/Rho2)*RhoL-1;
RhoPL-RhoL*(x1*Mw1+x2*Mw2);
Fin0+Fin0-1e-4*(t-tsd)+0-min([Fin0,Fin0-1e-4*(t-tsd),0])-max([Fin0,Fin0-1e-4*(t-tsd),0])-Fin;
cv*min(Vvmin,Vv)*max(0,(P-P0)/(abs(P-P0)+epsi)^0.5)-Fv;
cL*min(VLmin,VL)*max(0,(g*VL/A+(P-P0)/(RhoPL))/(abs(g*VL/A+(P-P0)/(RhoPL))+epsi)^0.5)-FL;
Q0+Q0-500*(t-tsd)+0-min([Q0,Q0-500*(t-tsd),0])-max([Q0,Q0-500*(t-tsd),0])-Q;
z1*h_Lcmp_in(1)+z2*h_Lcmp_in(2)-hin;
x1*h_Lcmp(1)+x2*h_Lcmp(2)-hL;
y1*h_Vcmp(1)+y2*h_Vcmp(2)+y3*h_Vcmp(3)-hv;
ML*hL+ Mv*hv-H;
U+P*(VT)-H;
Fin*z1-FL*x1-Fv*y1-dm1;
Fin*z2-FL*x2-Fv*y2-dm2;
Fin*z3-FL*x3-Fv*y3-dm3;
Fin*hin-FL*hL-Fv*hv+Q-du;
]
end

