package Continuous_Nanofiltration_Concentrations

 connector Inlet
  flow Real Q(quantity="volumetric flowrate", unit="m3/s") "Flowrate"; 
  input Real c(quantity="concentration", unit="mg/L") "concentration of contaminant";
 end Inlet;

 connector Outlet
 flow Real Q(quantity="volumetric flowrate", unit="m3/s")"Flowrate"; 
 output Real c(quantity="concentration", unit="mg/L") "concentration of  contaminant";
 end Outlet;
 
record membrane_properties

 parameter Real n = 7"Number of capillaries per fiber";
 parameter Real m = 1"Number of Membrane Fibers";
 parameter Modelica.SIunits.Length d_cap = 0.0009 "Capillary diameter";
 parameter Modelica.SIunits.Length l_m = 1.5 "Capillary membrane length";
 parameter Real sn = 4 "number of segments";

end membrane_properties;

record set_flow_parameters
  
parameter Modelica.SIunits.Velocity u_cf = 2; //Crossflow velocity at the membrane inlet
parameter Real J_p = 35; //Average permeate flux as the permeate leaves the membrane L/hm^2
parameter Real Y = 0.40; //Water conversion factor (Q_p/Q_raw)
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


model Raw_water_tank

Inlet In_MO;
Outlet Out_F;
extends calculated_flows;

Real V(unit="m3") "Volume of water in the raw water tank varying with progression of the filtration process with time";
Real V_initial(unit="m3")= 5*10^(-3) "Volume of water present in the tank";
Real c_raw(unit="mg/L") "concentration of salt in the raw water tank varying with progression of the filtration process with time";
Real c_raw_initial(unit="mg/L") = 100 "initial concentration of salt in the raw water tank";
Real Q_p(unit="m3/s") "Permeate flow rate based on constant flux value";

initial equation

V = V_initial;
c_raw = c_raw_initial;

equation 

Out_F.c = c_raw;
Q_p = J_p * A_act * 10 ^ (-3)/3600; //permeate flowrate based on J_p
Out_F.Q = -abs(Q_p/Y);

if V > abs(Out_F.Q) then

der(V) - In_MO.Q + Out_F.Q = 0 ;
der(c_raw * V) + Out_F.c * Out_F.Q - In_MO.Q * In_MO.c = 0;

else 

Out_F.Q + Out_F.Q = 0; // I am unsure what the steady state equation would be.
Out_F.c = 0;

end if;

end Raw_water_tank;

model Feed

 Inlet In_RW, In_cyc;
 Outlet Out_M;
 
 extends calculated_flows;
 Real Q_f = u_cf * A_cf;   //feed flowrate calculated based on fixed feed velocity in m^3/s
 
  equation 
  Out_M.Q = -abs(Q_f);
  In_cyc.Q + Out_M.Q + In_RW.Q = 0; 
  In_cyc.Q * In_cyc.c + Out_M.Q * Out_M.c + In_RW.Q * In_RW.c = 0;

end Feed;

  model Membrane
  
  Inlet In_F;
  Outlet Out_P, Out_MO;
  parameter Real R = 0.81; //rejection rate of the membrane;
  
  equation

  In_F.Q + Out_P.Q + Out_MO.Q = 0;
  Out_P.Q = -abs(0.15 * In_F.Q);
  R = 1 -(Out_P.c/In_F.c);
  In_F.Q * In_F.c + Out_P.Q * Out_P.c + Out_MO.Q * Out_MO.c = 0;
  
  
  end Membrane;

  model Permeate
  
  Inlet In_M;
  output Real Q_permeate;
  output Real c_permeate;
  
  equation
  Q_permeate = abs(In_M.Q);
  c_permeate = In_M.c;
  
  end Permeate;

  model Membrane_Outlet
  Inlet In_M;
  Outlet Out_RW, Out_cyc;
  
  equation
  Out_cyc.Q + Out_RW.Q + In_M.Q = 0;
  //Out_cyc.Q = -abs(0.7 * In_M.Q);
Out_cyc.c = In_M.c;
  Out_RW.c = In_M.c;
  
  
  
  end Membrane_Outlet;

  class Run
  
  Raw_water_tank RW;
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


end Continuous_Nanofiltration_Concentrations;
