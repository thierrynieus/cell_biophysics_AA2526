TITLE leak current

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}
 
NEURON {
    SUFFIX leak
    NONSPECIFIC_CURRENT il
    RANGE gl, el
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {       
    gl = 0.0003 (mho/cm2) 
    el = -54.3 (mV)
}
  
ASSIGNED {
    v (mV)
    il (mA/cm2)
}
 
BREAKPOINT {
    il = gl*(v - el)
}
