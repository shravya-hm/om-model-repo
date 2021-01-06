model Calculations
parameter Real Q_f[4] = {9.918, 10.0242, 10.0146, 10.7382};
Real Q_raw[4];
Real Q_mo[4];
Real Q_rec[4];
parameter Real Q_r[4] = {2.244, 2.4324, 1.725, 1.605};
parameter Real Q_p[4] = {0.48, 0.9528, 1.548, 0.7432};


equation

Q_f = Q_mo + Q_p;
Q_f = Q_raw + Q_rec;
Q_mo= Q_rec + Q_r;


annotation(
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,newInst");
end Calculations;
