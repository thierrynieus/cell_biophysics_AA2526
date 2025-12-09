TITLE Electrode for ramp current clamping

NEURON {
	POINT_PROCESS ramp
	RANGE delay, dur, tpeak, amp0, amp1, i
	ELECTRODE_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	delay = 100 (ms)
	dur = 800 (ms)
	tpeak = 200 (ms)	
	amp0 = 0.006 (nA)
	amp1 = 0.012 (nA)
}
ASSIGNED { i (nA) }

INITIAL {
	i = 0
}

BREAKPOINT {
	at_time(delay)
	at_time(delay+dur)
	if (t < delay+dur && t > delay) {
	    if (t< (delay+tpeak)) {
    	    i = amp0 + (amp1-amp0)/ tpeak * (t-delay)    
	    }else{
    	    i = amp1
	    }
		
	}else{
		i = 0
	}
}
