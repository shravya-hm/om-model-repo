package Nanofiltration

 connector Inlet
  flow Real Q(quantity="volumetric flowrate", unit="L/h") "Flowrate"; 
  input Real c(quantity="concentration", unit="mg/L") "concentration of contaminant";
  Real P(quantity="pressure", unit="Pa") "Potential variable";
  
 end Inlet;

 connector Outlet
  flow Real Q(quantity="volumetric flowrate", unit="L/h")"Flowrate"; 
  output Real c(quantity="concentration", unit="mg/L") "concentration of  contaminant";
  output Real P(quantity="pressure", unit="Pa") "Potential variable";
 
 end Outlet;
 
 record membrane_properties

 parameter Real n = 7"Number of capillaries per fiber";
 parameter Real m = 1"Number of Membrane Fibers";
 parameter Modelica.SIunits.Length d_cap = 0.0009 "Capillary diameter";
 parameter Modelica.SIunits.Length l_m = 1.5 "Capillary membrane length";
 parameter Real sn = 4 "number of segments";

end membrane_properties;

record set_flow_parameters
  
parameter Modelica.SIunits.Velocity u_cf = 0.4; //Crossflow velocity at the membrane inlet
//parameter Real J_p = 35; //Average permeate flux as the permeate leaves the membrane L/hm^2
//parameter Real Y = 0.40; //Water conversion factor (Q_p/Q_raw)
//parameter Real R = 0.81; //rejection rate of the membrane
//parameter Real R_int = 0.90; //intrinsic retention of membrane
//parameter Real k_w = 18.3;   //water permeability through the membrane in L/hm^2
//parameter Real c_raw = 100 "Concentration of contaminant / mg/L";

end set_flow_parameters;

class calculated_flows
  
  constant Real pi = 3.141592;
  //Real Q_p "Permeate flowrate / m3/s";
  Real d_h "hydraulic diameter";
  extends membrane_properties;
  extends set_flow_parameters;
  
  protected 
  output Real A_cf "Crosssectional flow area / m²";
  output Real A_act "Active Membrane area / m²";
  
  equation
  d_h = d_cap * m * n;
  A_cf = pi * (d_cap / 2) ^ 2 * n * m;
  A_act = pi * d_h * l_m;
  //Q_p = J_p * A_act; //permeate flowrate calulated based on J_p
end calculated_flows;

 
model Raw_water

Inlet In_MO;
Outlet Out_F;

Real V(unit="m3") "Volume of water in the raw water tank varying with progression of the filtration process with time";
Real V_initial(unit="m3")= 50*10^(-3) "Volume of water present in the tank";
parameter Real P = 1;

Real c_raw(unit="mg/L") "concentration of salt in the raw water tank varying with progression of the filtration process with time";

Real c_raw_initial(unit="mg/L") = 100 "initial concentration of salt in the raw water tank";

initial equation
V = V_initial;
c_raw = c_raw_initial;

equation 

Out_F.P = P;
In_MO.P = Out_F.P;

if V > abs(Out_F.Q) then

der(V) = In_MO.Q + Out_F.Q ;

else 
Out_F.Q = 0; // I am unsure what the steady state equation would be.

end if;

end Raw_water;

model Feed

 Inlet In_RW, In_cyc;
 Outlet Out_M;
 
 extends calculated_flows;
 Real Q_f = u_cf * A_cf;   //feed flowrate calculated based on fixed feed velocity in m^3/s
    
  equation 
  Out_M.Q = -abs(Q_f);
  In_cyc.Q + Out_M.Q + In_RW.Q = 0; 
  In_cyc.P = Out_M.P;
  In_RW.P = Out_M.P;

end Feed;

model Membrane
  
  Inlet In_F;
  Outlet Out_P, Out_MO;

  
  equation

  In_F.Q + Out_P.Q + Out_MO.Q = 0;
  Out_P.Q = -abs(0.15 * In_F.Q);  
  In_F.P = Out_P.P;
  Out_MO.P = Out_P.P;
  
  end Membrane;
  
  model Permeate
  
  Inlet In_M;
  output Real Q_permeate;
  output Real P_permeate;

  
  equation
  Q_permeate = abs(In_M.Q);
  P_permeate = In_M.P;

  
  end Permeate;

  model Membrane_Outlet
  Inlet In_M;
  Outlet Out_RW, Out_cyc;
  
  equation
  Out_cyc.Q + Out_RW.Q + In_M.Q = 0;
  Out_cyc.Q = -abs(0.7 * In_M.Q);
  
  end Membrane_Outlet;

  class Run
  
  Raw_water RW;
  Feed F;
  Membrane M;
  Permeate P;
  Membrane_Outlet MO;
  
  
  equation
  connect (RW.Out_F,F.In_RW);
  connect (F.Out_M, M.In_F);
  connect (M.Out_P, P.In_M);
  connect (M.Out_MO, MO.In_M);
  connect (MO.Out_RW, RW.In_MO);
  connect (MO.Out_cyc, F.In_cyc); 
  

end Run;

end Nanofiltration;
