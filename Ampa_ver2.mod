TITLE AMPA

COMMENT
ENDCOMMENT

NEURON {
	POINT_PROCESS AMPA_ver2

	NONSPECIFIC_CURRENT i
	RANGE Cdur, Erev, g, gmax, kB
	RANGE r1FIX, r6FIX, r1, r2, r5, r6
	:RANGE tau_1, tau_rec, tau_facil, U	 
	RANGE T, Tmax, Trelease		
	:RANGE NTdiffusion 	
	:RANGE xview,yview,zview,Pview
}

UNITS {
	(nA) 	= (nanoamp)
	(mV) 	= (millivolt)
	(umho) 	= (micromho)
	(mM) 	= (milli/liter)
	(pS)	= (picosiemens)
	(nS) 	= (nanosiemens)
	(um) 	= (micrometer)
	PI 	= (pi)		(1)
}

PARAMETER {
	: postsynaptic parameters
	gmax		= 1200  (pS)		
	Cdur		= 0.3	(ms)	
	Erev		= 0	(mV)
	kB		    = 0.44	(mM)
		 
	r1FIX		= 5.4	(/ms/mM) 
	r6FIX		= 1.12	(/ms/mM)		
	r2		    = 0.82	(/ms)		
	r5		    = 0.013	(/ms)		 
	
	: presynaptic parameters
	:tau_1 		= 3 (ms) 	< 1e-9, 1e9 >
	:tau_rec 	= 35.1 (ms) 	< 1e-9, 1e9 > 	
	:tau_facil 	= 10.8 (ms) 	< 0, 1e9 > 	

	:U 		    = 0.416 (1) 	< 0, 1 >
	:u0 		    = 0 (1) 	< 0, 1 >	 
	Tmax		= 1  (mM)
	
	 
}


ASSIGNED {
	v		(mV)		 
	i 		(nA)		 
	g 		(pS)		 
	T		(mM)
	
	r1		(/ms)
	r6		(/ms)

	Trelease	(mM)
 
	:x 
	:tsyn		(ms)
	
	:xview   
	:yview  
	:zview
	:Pview   
}

STATE {	
	C
	O
	D
}

INITIAL {
	C=1
	O=0
	D=0
	
	T=0 (mM)
	Trelease=0 (mM)
	
	:xview = 1 
	:yview = 0 
	:zview = 0 
	:Pview = 0	
}

BREAKPOINT {	
	Trelease = T  
	SOLVE kstates METHOD sparse
	g = gmax * O
	i = (1e-6) * g * (v-Erev) 
}

KINETIC kstates {
	r1 = r1FIX * Trelease^2 / (Trelease + kB)^2
	r6 = r6FIX * Trelease^2 / (Trelease + kB)^2
	~ C  <-> O (r1,r2)
	~ D  <-> C (r5,r6)
	CONSERVE C+O+D = 1
}

:NET_RECEIVE(weight, on, nspike, t0 (ms),y, z, u, tsyn (ms)) {
NET_RECEIVE(weight, on, nspike, t0 (ms)) { 
	INITIAL {
		:y = 0
		:z = 0
		:u = u0
		:tsyn = t
		nspike = 1
	}
  	if (flag == 0) { 
      	: printf("Event at time %g\n",t)
		: presynaptic modulation
		nspike = nspike + 1
		if (!on) {
			t0 = t
			on = 1		
					
			:z = z*exp( - (t - tsyn) / tau_rec )	
			:z = z + ( y*(exp(-(t - tsyn)/tau_1) - exp(-(t - tsyn)/tau_rec))/((tau_1/tau_rec)-1) )  
			:y = y*exp(-(t - tsyn)/tau_1)			
			:x = 1-y-z
			
			:if (tau_facil > 0) { 
			:	u = u*exp(-(t - tsyn)/tau_facil)
			:	u = u + U * ( 1 - u )							
			:} else { u = U }
			:y = y + x * u

			:xview = x	 
			:yview = y  
			:Pview = u

			T = Tmax :* y			
			 
			:tsyn = t
						
		}
		net_send(Cdur, nspike)	 
   	}
	if (flag == nspike) { 
			t0  = t
			T   = 0
			on  = 0
	}
}

