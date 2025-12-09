TITLE sodium current
 
COMMENT
Based on https://neuron.yale.edu/neuron/docs/hodgkin-huxley-using-rxd


q10 = 3^((celsius - 6.3)/10)
:"m" sodium activation system
alpha = .1 * vtrap(-(v+40),10)
beta =  4 * exp(-(v+65)/18)
sum = alpha + beta
mtau = 1/(q10*sum)
minf = alpha/sum
:"h" sodium inactivation system
alpha = .07 * exp(-(v+65)/20)
beta = 1 / (exp(-(v+35)/10) + 1)
sum = alpha + beta
htau = 1/(q10*sum)
hinf = alpha/sum


:"n" potassium activation system
alpha = .01*vtrap(-(v+55),10) 
beta = .125*exp(-(v+65)/80)
sum = alpha + beta
ntau = 1/(q10*sum)
ninf = alpha/sum



ENDCOMMENT
 
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)	
}
 
NEURON {
 	SUFFIX sodium
	USEION na WRITE ina
	RANGE gnabar, gna, minf, hinf, mtau, htau, m, h, ina, Am, v0_Am, k_Am, Bm, v0_Bm, k_Bm, Ah, v0_Ah, k_Ah, Bh, v0_Bh, k_Bh
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
  	
	ena	= 55	(mV)    : CHECK
	gnabar	= 0.12 (mho/cm2) : OK
	
	: HH params
	Am = 0.1
    v0_Am = 40
    k_Am = 10
    Bm = 4
    v0_Bm = 65
    k_Bm = 18 
    
	Ah = 0.07
    v0_Ah = 65
    k_Ah = 20
    Bh = 1
    v0_Bh = 35
    k_Bh = 10
     
}
 
STATE {
    m 
    h
}
 
ASSIGNED {
    v           (mV)
    ina         (mA/cm2)
    celsius		(degC)
 	minf
    hinf
    mtau
    htau
    gna
    q10
}
 
BREAKPOINT {
    SOLVE states METHOD cnexp
    gna = gnabar * m^3 * h
    ina = gna * (v - ena)
}
 
UNITSOFF
 
INITIAL {
    rates(v)
    m = minf
    h = hinf
    q10 = 3.0^((celsius - 6.3)/10.0)
}

DERIVATIVE states {  
    :Computes states variable m and h at the current v and t.
    rates(v)
    m' = (minf - m) / mtau
    h' = (hinf - h) / htau
}
 
PROCEDURE rates(v) {  
    :Computes rate and other constants at current v.
    :Call once from HOC to initialize inf at resting v.
    LOCAL alpha_m, beta_m, alpha_h, beta_h
    : m
    alpha_m = Am * vtrap(-(v + v0_Am), k_Am)
    beta_m = Bm * exp(-(v + v0_Bm) / k_Bm)
    mtau = 1.0/(q10 * (alpha_m + beta_m))
    minf = alpha_m / (alpha_m + beta_m)
    : h
    alpha_h = Ah * exp(-(v + v0_Ah)/k_Ah)
    beta_h = Bh / ( 1 + exp(-(v + v0_Bh)/k_Bh))
    htau = 1 / ( q10 * ( alpha_h + beta_h ) )
    hinf = alpha_h / ( alpha_h + beta_h )    

  
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
