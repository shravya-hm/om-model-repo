package Theoretical_nanofiltration
 
 connector Inlet
  flow Real Q(quantity="volumetric flowrate", unit="L/h") "Flowrate"; 
  input Real c(quantity="concentration", unit="mg/L") "concentration of contaminant";
 end Inlet;

 connector Outlet
 flow Real Q(quantity="volumetric flowrate", unit="L/h")"Flowrate"; 
 output Real c(quantity="concentration", unit="mg/L") "concentration of  contaminant";
 end Outlet;
 
 record membrane_properties

  parameter Real n = 7"Number of capillaries per fiber";
  parameter Real m = 1"Number of Membrane Fibers";
  parameter Modelica.SIunits.Length d_cap = 0.0009 "Capillary diameter";
  parameter Modelica.SIunits.Length l_m =1.3 "Capillary membrane length";
 
  end membrane_properties;
  
  record set_flow_parameters
    
  parameter Modelica.SIunits.Velocity u_cf = 8;     //Crossflow velocity at the membrane inlet
  parameter Real J_p = 35;     //Average permeate flux as the permeate leaves the membrane L/hm^2
  parameter Real Y = 0.601;   //Water conversion factor (Q_p/Q_raw)
    //parameter Real R = 0.81; //rejection rate of the membrane
  parameter Real R_int = 0.78;     //intrinsic retention of membrane
  parameter Real k_w = 10;       //water permeability through the membrane in L/hm^2
  parameter Real c_f = 100 "Concentration of contaminant / mg/L";
  end set_flow_parameters;
  
  record Conc_pol_parameters
  //output
  Real Re "Reynoldsnumber";
  Real Sc "Schmidtnumber";
  Real Sh "Sherwoodnumber";
  Real CP "concentration polarization";
  Real delta "d_h/Sh";
  Real k "D_i.Sh/d_h";   //mass transfer coefficient
    //input
  constant Real D_i = 1.07*10^(-9) "Diffusion coefficient / m²/s";
  
  end Conc_pol_parameters;
  
  record osmopressure_solute
 constant Real R_gasconstant = Modelica.Constants.R;
    parameter Real alpha = 1 "disscociation constant between 0 and 1";
    parameter Real T = 293 "temperature /K";
    parameter Real v = 2 "stoichiometric coefficient";
    //c_f is the only varaiation from one membrane module to the next
    Real P_osmotic;
    Real B;
 end osmopressure_solute;
    
    class calculated_flows
 constant Real pi = 3.141592;
    Real Q_p "Permeate flowrate / L/h";
    Real d_h "hydraulic diameter";
    extends membrane_properties;
    extends set_flow_parameters;
  protected
    Real A_cf "Crosssectional flow area / m²";
    Real A_act "Active Membrane area / m²";
  equation
    d_h = d_cap * m * n;
    A_cf = pi * (d_cap / 2) ^ 2 * n * m;
    A_act = pi * d_h * l_m;
    Q_p = J_p * A_act;
//permeate flowrate calulated based on J_p
  end calculated_flows;
    
   model Feed
  
   Outlet Out;
   extends calculated_flows;
   Real Q_f = u_cf * A_cf * 3600 * 1000;       //feed flowrate calculated based on fixed feed velocity
   equation 
    Out.Q = -abs(Q_f);
    Out.c = c_f;

  end Feed;
  
   model Membrane
  
  Inlet In;
  Outlet Out_P, Out_R;
  Real R_obs "rejection considering CP";
  Real c_m "concentration along the membrane / mg/L";
  Real t "thickness of boundary layer / m";
  Real TMP "transmembrane pressure /bar";
  Real WCF "water conversion factor";
  Real P_loss "pressure loss from inlet to outlet of the membrane due to friction, only accounted for by the concentration polarization in this model";
  Real R_raw "Rejecion of conductivity with respect to raw water";
  parameter Real rho = 998.2 "density of water kg/m^3";
  parameter Real mu = 1 * 10^ (-3) "dynamic viscosity / Pa.s" ;
  Real P_loss_HP "Pressure loss along the membrane length using Hagen-poiseuille equation";
  
  extends calculated_flows;
  extends Conc_pol_parameters;
  extends osmopressure_solute;
  
  protected
    Real f_D = 64 / Re "friction factor";
    
  equation
//Out_P.c = In.c * (1-R); // permeate from the membrane calculated based on assumed R value
    Out_P.Q = -abs(Q_p);
  In.Q + Out_P.Q + Out_R.Q = 0;
  (Out_P.c * Out_P.Q) + (Out_R.c * Out_R.Q) +  (In.c * In.Q) = 0;
  WCF = abs(Out_P.Q/In.Q);
  
  (CP, Re, Sc, Sh, delta, k)= Concentration_polarization(d_cap, D_i, l_m, c_f, u_cf, J_p, R_int);
  c_m = CP * c_f; 
  R_int = 1 - (Out_P.c/c_m); // calculates permeate concentration
  R_obs = 1 - (Out_P.c/c_f); // calculates Rejection
  t = D_i/k;
  
  P_loss = f_D * rho * l_m * u_cf ^ 2 / (2 * d_cap) * 10^(-5);
  P_loss_HP = (In.Q * (8/(m*n)) * mu * l_m) * 10 ^(-5)/ ((d_cap/2)^4 * pi * 3.6*10^(6)) ;
// convert from pa to bar - Q from l/h to m3/s, flow is considered to be equally distributed through the m capillaries in the n fibres
   
   (B, P_osmotic) = osmotic_pressure(R_gasconstant, alpha, T, c_f, v);  
  
   k_w = J_p/(TMP-P_osmotic);//This equation calculates the TMP
  R_raw = 1 - (Out_P.c/c_f);
  
  end Membrane;
  
  model Permeate

  Inlet In; 
  output Real Q_permeate;
  output Real c_permeate;
  equation
  Q_permeate = In.Q;
  c_permeate = In.c;
  
  end Permeate;
  
   model Retentate

  output Real Q_retentate;
  output Real c_retentate;
  Inlet In;
  equation 
  Q_retentate = In.Q;
  c_retentate = In.c;

  end Retentate;
  
  class Run
  
  Feed F;
  Membrane M;
  Permeate P;
  Retentate R;

  equation
  
  connect (F.Out, M.In);
  connect (M.Out_P, P.In);
  connect (M.Out_R, R.In);

  end Run;

  function Concentration_polarization
  
  input Real d_cap, D_i, l_m, c_f, u_cf, J_p, R_int;
  output Real CP "Concentration polarization / c_0/c_f", Re "Reynold's number", Sc "Schmidt number", Sh "Sherwood number", delta "thickness of boundary layer laminar flow / m", k "Mass transfer coefficient";
  
   protected
   
  constant Real e = Modelica.Constants.e;
  Real mu = 1 * 10^ (-3) "dynamic viscosity / Pa.s" ;
  Real rho = 998.2 "density of water kg/m^3";
  
  algorithm
  
  Re := d_cap * u_cf * rho / mu;
  Sc:= mu /(rho*D_i);
  
  if Re<2300 
  then Sh:=1.62*(Re*Sc*d_cap/l_m)^(1/3);
  else
  Sh:=0.04*Re^(3/4)*Sc;
  end if;
  delta:=d_cap/Sh;
  k:=D_i/delta;
  CP:=e^(J_p/1000/3600/k)/(R_int+(1-R_int)*e^(J_p/1000/3600/k));
  
  end Concentration_polarization;

  function osmotic_pressure
  
  input Real R_gasconstant, alpha, T, c_f, v "c_f is the only varaiation from one membrane module to the next";
  output Real B, P_osmotic;

  algorithm
  B := 1 + alpha*(v-1);
  P_osmotic := B * R_gasconstant *10^(-2) * T * c_f / (1000 * 120.3676); // Use magnesium sulphate in moles
  end osmotic_pressure;
    
end Theoretical_nanofiltration;
