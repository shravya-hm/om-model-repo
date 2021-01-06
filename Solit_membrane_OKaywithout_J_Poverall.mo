package Nanofiltration_CP  
  connector Inlet  
    flow Real Q(quantity = "volumetric flowrate", unit = "L/h") "Flowrate";
    input Real c(quantity = "concentration", unit = "mg/L") "concentration of contaminant";
    Real P;
  end Inlet;

  connector Outlet  
    flow Real Q(quantity = "volumetric flowrate", unit = "L/h") "Flowrate";
    output Real c(quantity = "concentration", unit = "mg/L") "concentration of  contaminant";
    output Real P;
  end Outlet;

  record membrane_properties  
    parameter Real n = 7 "Number of capillaries per fiber";
    parameter Real m = 15 "Number of Membrane Fibers";
    parameter Modelica.SIunits.Length d_cap = 0.0009 "Capillary diameter";
    parameter Modelica.SIunits.Length l_m = 1.5 "Capillary membrane length";
  end membrane_properties;

  record set_flow_parameters  
    parameter Modelica.SIunits.Velocity u_cf = 0.4;
    parameter Real Y = 0.412;
    parameter Real k_w = 18.3;
    parameter Real P_feed = 3.31;
  end set_flow_parameters;

  record Conc_pol_parameters  
    Real Re "Reynoldsnumber";
    Real Sc "Schmidtnumber";
    Real Sh "Sherwoodnumber";
    Real CP "concentration polarization";
    Real delta "d_h/Sh";
    Real k "D_i.Sh/d_h";
    constant Real D_i = 1.07 * 10 ^ (-9) "Diffusion coefficient / m²/s";
  end Conc_pol_parameters;

  record osmopressure_solute  
    constant Real R_gasconstant = Modelica.Constants.R;
    parameter Real alpha = 1 "disscociation constant
      between 0 and 1";
    parameter Real T = 293 "temperature /K";
    parameter Real v = 2 "stoichiometric coefficient";
    Real P_osmotic;
    Real B;
  end osmopressure_solute;

  class calculated_flows  
    constant Real pi = 3.141592;
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

  model Membrane  
    Inlet In;
    Outlet Out_P;
    Outlet Out_M;
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
    parameter Real mu = 1 * 10 ^ (-3) "dynamic viscosity / Pa.s";
    Real Q_psn "permeate leaving segment n";
    Real J_psn "permeate flux through segment sn";
    Real u_cfsn;
    Real P_in = P_feed;
    Real P_out;
    parameter Real R_int = 0.9;
    extends calculated_flows;
    extends Conc_pol_parameters;
    extends osmopressure_solute;
  protected
    Real f_D = 64 / Re "friction factor";
  equation
    abs(In.c) = c_fsn;
    Q_psn = J_psn * (A_act / 4);
    u_cfsn = In.Q / (A_cf * 3600 * 1000);
    Out_P.Q = -abs(Q_psn);
    In.Q + Out_P.Q + Out_M.Q = 0;
    Out_P.c * Out_P.Q + Out_M.c * Out_M.Q + In.c * In.Q = 0;
    WCF = abs(Out_P.Q / In.Q);
    (CP, Re, Sc, Sh, delta, k) = Concentration_polarization(d_cap, D_i, l_m * sn / 4, c_fsn, u_cfsn, J_psn, R_int);
    (B, P_osmotic) = osmotic_pressure(R_gasconstant, alpha, T, c_msn, v);
    c_msn = CP * c_fsn;
    R_int = 1 - Out_P.c / c_fsn;
    R_obssn = 1 - Out_P.c / c_msn;
    tsn = D_i / k;
    In.P - Out_M.P = P_loss_sn;
    Out_P.P = 0;
    P_loss_sn = f_D * rho * l_m / 4 * u_cfsn ^ 2 / (2 * d_cap) * 10 ^ (-5);
    P_loss = f_D * rho * l_m * sn / 4 * u_cfsn ^ 2 / (2 * d_cap) * 10 ^ (-5);
    P_loss_HP = In.Q * (8 / (m * n)) * mu * l_m * sn / 4 * 10 ^ (-5) / ((d_cap / 2) ^ 4 * pi * 3.6 * 10 ^ 6);
    P_in - P_out = P_loss_HP;
    TMP_test = (P_in + P_out) / 2;
    TMP = (In.P + Out_M.P) / 2;
    k_w = J_psn / (TMP - P_osmotic);
    e_sn = P_loss_sn * u_cfsn * A_cf * 3600 * 1000;
  end Membrane;

  function Concentration_polarization  
    input Real d_cap;
    input Real D_i;
    input Real l_m;
    input Real c_f;
    input Real u_cf;
    input Real J_p;
    input Real R_int;
    output Real CP "Concentration polarization / c_0/c_f";
    output Real Re "Reynold's number";
    output Real Sc "Schmidt number";
    output Real Sh "Sherwood number";
    output Real delta "thickness of boundary layer laminar flow / m";
    output Real k "Mass transfer coefficient";
  protected
    constant Real e = Modelica.Constants.e;
    Real mu = 1 * 10 ^ (-3) "dynamic viscosity / Pa.s";
    Real rho = 998.2 "density of water kg/m^3";
  algorithm
    Re := d_cap * u_cf * rho / mu;
    Sc := mu / (rho * D_i);
    if Re < 2300 then
      Sh := 1.62 * (Re * Sc * d_cap / l_m) ^ (1 / 3);
    else
      Sh := 0.04 * Re ^ (3 / 4) * Sc;
    end if;
    delta := d_cap / Sh "thickness of boundary layer laminar flow / m";
    k := D_i / delta "Mass transfer coefficient";
    CP := e ^ (J_p / 1000 / 3600 / k) / (R_int + (1 - R_int) * e ^ (J_p / 1000 / 3600 / k)) "Concentration polarization / c_0/c_f";
  end Concentration_polarization;

  function osmotic_pressure  
    input Real R_gasconstant;
    input Real alpha;
    input Real T;
    input Real c_fsn;
    input Real v;
    output Real B;
    output Real P_osmotic;
  algorithm
    B := 1 + alpha * (v - 1);
    P_osmotic := B * R_gasconstant * 10 ^ (-2) * T * c_fsn / (1000 * 120.3676);
  end osmotic_pressure;
end Nanofiltration_CP;

package ModelicaServices  "ModelicaServices (OpenModelica implementation) - Models and functions used in the Modelica Standard Library requiring a tool specific implementation" 
  extends Modelica.Icons.Package;

  package Machine  "Machine dependent constants" 
    extends Modelica.Icons.Package;
    final constant Real eps = 1e-15 "Biggest number such that 1.0 + eps = 1.0";
    final constant Real small = 1e-60 "Smallest number such that small and -small are representable on the machine";
    final constant Real inf = 1e60 "Biggest Real number such that inf and -inf are representable on the machine";
    final constant Integer Integer_inf = OpenModelica.Internal.Architecture.integerMax() "Biggest Integer number such that Integer_inf and -Integer_inf are representable on the machine";
  end Machine;
  annotation(version = "3.2.3", versionBuild = 3, versionDate = "2019-01-23", dateModified = "2019-09-21 12:00:00Z"); 
end ModelicaServices;

package Modelica  "Modelica Standard Library - Version 3.2.3" 
  extends Modelica.Icons.Package;

  package Math  "Library of mathematical functions (e.g., sin, cos) and of functions operating on vectors and matrices" 
    import SI = Modelica.SIunits;
    extends Modelica.Icons.Package;

    package Icons  "Icons for Math" 
      extends Modelica.Icons.IconsPackage;

      partial function AxisCenter  "Basic icon for mathematical function with y-axis in the center" end AxisCenter;
    end Icons;

    function asin  "Inverse sine (-1 <= u <= 1)" 
      extends Modelica.Math.Icons.AxisCenter;
      input Real u;
      output SI.Angle y;
      external "builtin" y = asin(u);
    end asin;

    function exp  "Exponential, base e" 
      extends Modelica.Math.Icons.AxisCenter;
      input Real u;
      output Real y;
      external "builtin" y = exp(u);
    end exp;
  end Math;

  package Constants  "Library of mathematical constants and constants of nature (e.g., pi, eps, R, sigma)" 
    import SI = Modelica.SIunits;
    import NonSI = Modelica.SIunits.Conversions.NonSIunits;
    extends Modelica.Icons.Package;
    final constant Real e = Modelica.Math.exp(1.0);
    final constant Real pi = 2 * Modelica.Math.asin(1.0);
    final constant SI.Velocity c = 299792458 "Speed of light in vacuum";
    final constant SI.FaradayConstant F = 9.648533289e4 "Faraday constant, C/mol (previous value: 9.64853399e4)";
    final constant Real R(final unit = "J/(mol.K)") = 8.3144598 "Molar gas constant (previous value: 8.314472)";
    final constant Real N_A(final unit = "1/mol") = 6.022140857e23 "Avogadro constant (previous value: 6.0221415e23)";
    final constant Real mue_0(final unit = "N/A2") = 4 * pi * 1.e-7 "Magnetic constant";
  end Constants;

  package Icons  "Library of icons" 
    extends Icons.Package;

    partial package Package  "Icon for standard packages" end Package;

    partial package IconsPackage  "Icon for packages containing icons" 
      extends Modelica.Icons.Package;
    end IconsPackage;
  end Icons;

  package SIunits  "Library of type and unit definitions based on SI units according to ISO 31-1992" 
    extends Modelica.Icons.Package;

    package Conversions  "Conversion functions to/from non SI units and type definitions of non SI units" 
      extends Modelica.Icons.Package;

      package NonSIunits  "Type definitions of non SI units" 
        extends Modelica.Icons.Package;
        type Temperature_degC = Real(final quantity = "ThermodynamicTemperature", final unit = "degC") "Absolute temperature in degree Celsius (for relative temperature use SIunits.TemperatureDifference)" annotation(absoluteValue = true);
      end NonSIunits;
    end Conversions;

    type Angle = Real(final quantity = "Angle", final unit = "rad", displayUnit = "deg");
    type Length = Real(final quantity = "Length", final unit = "m");
    type Velocity = Real(final quantity = "Velocity", final unit = "m/s");
    type Acceleration = Real(final quantity = "Acceleration", final unit = "m/s2");
    type ElectricCharge = Real(final quantity = "ElectricCharge", final unit = "C");
    type FaradayConstant = Real(final quantity = "FaradayConstant", final unit = "C/mol");
  end SIunits;
  annotation(version = "3.2.3", versionBuild = 3, versionDate = "2019-01-23", dateModified = "2019-09-21 12:00:00Z"); 
end Modelica;

model Membrane_total
  extends Nanofiltration_CP.Membrane;
end Membrane_total;
