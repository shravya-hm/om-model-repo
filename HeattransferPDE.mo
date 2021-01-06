model HeattransferPDE

parameter Real rho = 1000;
parameter Real Cp = 1;
parameter Real L =1;
parameter Real k =1;
parameter Real Tc = 0;
parameter Real Th = 100;
parameter Real Tinitial = 30;
parameter Integer N = 10;
Real T[N-1]; //state variables representing the rod temperature
Real dX = L/N; //Discretization length
Real i =1;

initial equation

for i in 1:N-1 loop
 T[i] = Tinitial;
end for;

equation

rho*Cp*der(T[1]) = k*(T[2]-(2*T[1])+Th)/(dX^2);
rho*Cp*der(T[N-1]) = k*(Tinitial-(2*T[N-1])+T[N-2])/(dX^2);

for i in 2:N-2 loop
 rho*Cp*der(T[i]) = k*(T[i+1]-(2*T[i])*T[i-1])/(dX^2);
 end for;
 

equation

end HeattransferPDE;
