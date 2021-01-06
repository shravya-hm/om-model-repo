package Nanofiltration_split_continuous
  
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
  parameter Real m = 15"Number of Membrane Fibers";
  parameter Modelica.SIunits.Length d_cap = 0.0009 "Capillary diameter";
  parameter Modelica.SIunits.Length l_m = 1.5 "Capillary membrane length";
 
 end membrane_properties;
  
 record set_flow_parameters
    
  parameter Modelica.SIunits.Velocity u_cf =0.4; //Crossflow velocity at the membrane inlet
  //parameter Real J_p = 35; //Average permeate flux as the permeate leaves the membrane L/hm^2
  parameter Real Y = 0.412; //Water conversion factor (Q_p/Q_raw)
  parameter Real k_w = 18.3;   //water permeability through the membrane in L/hm^2
  parameter Real P_feed = 3.31; // used to calculate potential presure terms
  
 end set_flow_parameters;
  
  record Conc_pol_parameters
  //output
  Real Re "Reynoldsnumber";
  Real Sc "Schmidtnumber";
  Real Sh "Sherwoodnumber";
  Real CP "concentration polarization";
  Real delta "d_h/Sh";
  Real k "D_i.Sh/d_h";
  //input
  constant Real D_i = 1.07*10^(-9) "Diffusion coefficient / m²/s";
  
  end Conc_pol_parameters;

  record osmopressure_solute
  
  constant Real R_gasconstant = Modelica.Constants.R;
  parameter Real alpha = 1 "disscociation constant
  between 0 and 1";
  parameter Real T = 293 "temperature /K";
  //why error
  parameter Real v = 2 "stoichiometric coefficient";
 //c_f is the only variation from one membrane module to the next
  Real P_osmotic;
  Real B;
  
  end osmopressure_solute;
  
    class calculated_flows
      
      constant Real pi = 3.141592;
      //Real Q_p "Permeate flowrate / L/h";
      Real d_h "hydraulic diameter";
      extends membrane_properties;
      extends set_flow_parameters;
      
      output Real A_cf "Crosssectional flow area / m²";
      output Real A_act "Active Membrane area / m²";
      
      equation
      d_h = d_cap * m * n;
      A_cf = pi * (d_cap / 2) ^ 2 * n * m;
      A_act = pi * d_h * l_m;
      
    
      
    end calculated_flows;
    
    model Raw_water

    Inlet In_MO;
    Outlet Out_F; 
    extends calculated_flows;
    parameter Real c_raw_init = 100 "Concentration of contaminant / mg/m3";
    Real Q_raw_f "raw water flowrate calculated based on Y / m3/s " ;
    Real Q_p;
    Real J_p;
    Real c_raw;
    
    Real V(unit="m3") "Volume of water in the raw water tank varying with progression of the filtration process with time";
    Real V_initial(unit="m3")= 500*10^(-3) "Volume of water present in the tank";
    
    initial equation
    c_raw = c_raw_init;
    V = V_initial;
    
    equation
    
    Out_F.c = c_raw;
    Q_p = J_p * A_act; //permeate flowrate calulated based on J_p
    Q_raw_f = Q_p/Y;
    Out_F.Q = -abs(Q_raw_f);
    
    if V > Q_raw_f then
    
    der(V) = In_MO.Q + Out_F.Q ;
    der(c_raw * V) = In_MO.Q * In_MO.c ;
    In_MO.P = 1;
    Out_F.P = 1;
    
   else
   In_MO.Q - Out_F.Q = 0; // I am unsure what the steady state equation would be.
    end if;
    
    end Raw_water;

  model Feed
  
   Inlet In_raw, In_cyc;
   Outlet Out;
   extends calculated_flows;
   Real Q_f = u_cf * A_cf * 3600 * 1000;   //feed flowrate calculated based on fixed feed velocity
   //Real P_head;
   
   equation 
    (In_raw.c * In_raw.Q) + (In_cyc.c * In_cyc.Q) + (Out.c * Out.Q) = 0;
    Out.Q = -abs(Q_f);
    In_cyc.Q + Out.Q + In_raw.Q = 0; 
    In_raw.P =1;
    In_cyc.P =1;
    Out.P = 1;

  end Feed;

  model Membrane
  
  Inlet In;
  Outlet Out_P, Out_M;
  
  Real sn "segment number";
  
  Real R_obssn "rejection considering CP";
  Real c_msn "concentration along the membrane / mg/L";
  Real c_fsn "concentration of feed or bulk / mg/L";
  Real tsn "thickness of boundary layer / m";
  Real TMP "transmembrane pressure /bar";
  Real WCF "water conversion factor";
  
  Real P_loss "Pressure loss";
  Real P_loss_HP "Pressure loss along the membrane length using Hagen-poiseuille equation";
  Real P_loss_sn "Pressure loss in each segment";
  Real e_sn "energy used in each segment unit";
  Real TMP_test;
  
  parameter Real rho = 998.2 "density of water kg/m^3";
  parameter Real mu = 1 * 10^ (-3) "dynamic viscosity / Pa.s" ;
  
  Real Q_psn "permeate leaving segment n";
  Real J_psn "permeate flux through segment sn";
  Real u_cfsn;
  //Parameter assumption to test the model-
    //Pressure varies with velocity and a relationship must be attained experimentally, or J_p(TMP) must be calculated by another method and not Darcy's law
    //Intrinsic retetntion can be obtained from experimental calculations
  Real P_in = P_feed;
  Real P_out;
  
  parameter Real R_int = 0.9; //intrinsic retention of membrane
  
  extends calculated_flows;
  extends Conc_pol_parameters;
  extends osmopressure_solute;
  
  
  protected
  Real f_D = 64 / Re "friction factor";
  
  equation
  
  abs(In.c) = c_fsn;
  //Out_P.c = In.c * (1-R); // permeate from the membrane calculated based on assumed R value
  
  Q_psn = J_psn * (A_act / 4);
  u_cfsn = In.Q / (A_cf * 3600 * 1000);
  
  Out_P.Q = -abs(Q_psn);
  In.Q + Out_P.Q + Out_M.Q = 0;
  (Out_P.c * Out_P.Q) + (Out_M.c * Out_M.Q) + (In.c * In.Q) = 0;
  WCF = abs(Out_P.Q/In.Q);
  
  (CP, Re, Sc, Sh, delta, k) = Concentration_polarization(d_cap, D_i, l_m*sn/4, c_fsn, u_cfsn, J_psn, R_int);
  
  
  (B, P_osmotic) = osmotic_pressure(R_gasconstant, alpha, T, c_msn, v);
  
  c_msn = CP * c_fsn; 
  R_int = 1 - (Out_P.c/c_fsn);
  R_obssn = 1 - (Out_P.c/c_msn);
  tsn = D_i/k;

  In.P - Out_M.P = P_loss_sn;
  Out_P.P = 0;
  
  P_loss_sn =  f_D * rho * l_m/4 * u_cfsn ^ 2 / (2 * d_cap) * 10^(-5);
  P_loss = f_D * rho * l_m * sn/4 * u_cfsn ^ 2 / (2 * d_cap) * 10^(-5);
  P_loss_HP = (In.Q * (8/(m*n)) * mu * l_m * sn/4) * 10 ^(-5)/ ((d_cap/2)^4 * pi * 3.6*10^(6)) ;
// convert from pa to bar - Q from l/h to m3/s, flow is considered to be equally distributed through the m capillaries in the n fibres
  P_in - P_out = P_loss_HP;
  TMP = (P_in + P_out)/2;
  TMP_test = (In.P + Out_M.P)/2;
  
  k_w = J_psn/(TMP-P_osmotic);
  e_sn = P_loss_sn * u_cfsn * A_cf * 3600 * 1000;
  
    
  end Membrane;
  
  model Permeate

  Inlet In_P1, In_P2, In_P3, In_P4; 
  output Real Q_permeate;
  output Real c_permeate;
  
  equation
  
  Q_permeate = In_P1.Q + In_P2.Q + In_P3.Q + In_P4.Q;
  c_permeate = (In_P1.c * In_P1.Q + In_P2.c * In_P2.Q + In_P3.c * In_P3.Q + In_P4.c * In_P4.Q)/Q_permeate;
  
  
  end Permeate;
  
  model Membrane_outlet

  Inlet In;
  Outlet Out_RW, Out_cyc;
  //parameter Real rec = -1 "recirculation";
  equation
  Out_RW.c = In.c;
  Out_cyc.c = In.c;
//Out_cyc.Q = -abs(rec);
  Out_cyc.Q + Out_RW.Q + In.Q = 0;
  Out_cyc.P + Out_RW.P = In.P;
  
  end Membrane_outlet;
  
  class Run
  
  Real WCF_overall;
  Real R_overall;
  Real P_loss_overall;
  Real R_raw_overall;
  Real J_p_overall;
  
  
  Raw_water RW (J_p = 35);
  Feed F;
  Membrane M1(sn=1);
  Membrane M2(sn=2);
  Membrane M3(sn=3);
  Membrane M4(sn=4);
  Permeate Pn;
  Membrane_outlet MO;
  
  
  equation
  connect (RW.Out_F,F.In_raw);
  connect (F.Out, M1.In);
  
  connect (M1.Out_P, Pn.In_P1);
  connect (M2.Out_P, Pn.In_P2);
  connect (M3.Out_P, Pn.In_P3);
  connect (M4.Out_P, Pn.In_P4);
    
  //final permeate stream is made up of the 4 permeate streams
  
  connect (M1.Out_M, M2.In);
  connect (M2.Out_M, M3.In);
  connect (M3.Out_M, M4.In);
  connect (M4.Out_M, MO.In);
  connect (MO.Out_RW, RW.In_MO);
  connect (MO.Out_cyc, F.In_cyc); 
  
  WCF_overall = abs(Pn.Q_permeate/F.Out.Q) ;
  R_overall = 1 - (Pn.c_permeate/F.Out.c);
  R_raw_overall = 1 - (Pn.c_permeate/RW.c_raw_init);
  P_loss_overall = M1.P_loss_sn + M2.P_loss_sn + M3.P_loss_sn + M4.P_loss_sn;
  J_p_overall = M1.J_psn + M2.J_psn + M3.J_psn + M4.J_psn;
  
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
  delta:=d_cap/Sh "thickness of boundary layer laminar flow / m";
  k:=D_i/delta "Mass transfer coefficient";
  CP:=e^(J_p/1000/3600/k)/(R_int+(1-R_int)*e^(J_p/1000/3600/k)) "Concentration polarization / c_0/c_f";
  
  end Concentration_polarization;

  function osmotic_pressure
  
  input Real R_gasconstant, alpha, T, c_fsn, v;//c_f is the only varaiation from one membrane module to the next
  output Real B, P_osmotic;

  algorithm
  B := 1 + alpha*(v-1);
  P_osmotic := B * R_gasconstant *10^(-2) * T * c_fsn / (1000 * 120.3676); // Use magnesium sulphate in moles
  
  end osmotic_pressure;

  
end Nanofiltration_split_continuous;
