TITLE potassium current
 
COMMENT
Based on https://neuron.yale.edu/neuron/docs/hodgkin-huxley-using-rxd

ENDCOMMENT
 
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)	
}
 
NEURON {
 	SUFFIX potassium
	USEION k WRITE ik
	RANGE gkbar, gk, ninf, n, ik, ntau
	: RANGE ek
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
	gkbar = 0.036 (mho/cm2)
	ek	  = -77	(mV)

	: HH params
    
    An = 0.01
    v0_An = 55
    k_An = 10
    Bn = 0.125
    v0_Bn = 65
    k_Bn = 80

}
 
STATE {
    n
}
 
ASSIGNED {
    v           (mV)
    ik          (mA/cm2)
    celsius	    (degC)
 	ninf
    ntau
    gk
    q10
}
 
BREAKPOINT {
    SOLVE states METHOD cnexp
    gk = gkbar * n^4
    ik = gk * (v - ek)
}
 
UNITSOFF
 
INITIAL {
    rates(v)
    n = ninf
    q10 = 3.0^((celsius - 6.3)/10.0)
}

DERIVATIVE states {  
    :Computes state variable n at the current v and t.
    rates(v)
    n' = (ninf - n) / ntau
}
 
PROCEDURE rates(v) {  
    :Computes rate and other constants at current v.
    :Call once from HOC to initialize inf at resting v.
    LOCAL alpha_n, beta_n
    : n
    
    alpha_n = An * vtrap(-(v + v0_An), k_An)
    beta_n = Bn * exp(-(v + v0_Bn)/k_Bn)
    ntau = 1 / ( q10 * ( alpha_n + beta_n ))
    ninf = alpha_n / ( alpha_n + beta_n )    
    
}
 
FUNCTION vtrap(x (mV),y (mV)) (mV) {
    : vtrap(x,y) is 1/(exp(x/y)-1) if |x/y|>=1e-6 or y*(1.0 - x/y/2.0) otherwise.
    if (fabs(x/y) < 1e-6) {
            vtrap = y*(1 - x/y/2)
    }else{
            vtrap = x/(exp(x/y) - 1)
    }
}
 
UNITSON
