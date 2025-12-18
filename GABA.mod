TITLE Gaba mod file to test the potential role of singly versus doubly bound receptors spatially distributed like on a synaptic disk

COMMENT
    GABA model adapted from Petrini et al. 2011 J.Neuroscience
    
	This represents a cleaned version aimed for demo purpose. Kinetic rates come from v4 -> v6 with cooperativity (kon1,kon2 with kon2>>kon1)
		
	Grc[0].synGa1_slice[0].konC0C1 = 2.5
	Grc[0].synGa1_slice[0].konC1C2 = 20
	Grc[0].synGa1_slice[0].kon = 5
	Grc[0].synGa1_slice[0].koff = 0.3
	Grc[0].synGa1_slice[0].a1 = 0.35
	Grc[0].synGa1_slice[0].b1 = 0.1
	Grc[0].synGa1_slice[0].a2 = 0.375
	Grc[0].synGa1_slice[0].b2 = 8
	Grc[0].synGa1_slice[0].r1 = 0.00013
	Grc[0].synGa1_slice[0].r2 = 0.04
	Grc[0].synGa1_slice[0].d1 = 0.013
	Grc[0].synGa1_slice[0].d2 = 1.5

ENDCOMMENT



NEURON {
	POINT_PROCESS GABA
	 
	NONSPECIFIC_CURRENT i
	RANGE g,Open 
	RANGE gmax,Cdur,Tmax,T,Erev,U		 
	RANGE kon,koff,d3,r3,d1d2,r1r2,a1,b1,a2,b2,r1,r2,d1,d2
	RANGE konDEFAULT,konC0C1,konC1C2
	RANGE Prel,xview,yview,zview,Pview
	RANGE tau_rec,tau_facil 
}

UNITS {
	(nA) 	= (nanoamp)
	(mV) 	= (millivolt)
	(umho)  = (micromho)
	(mM) 	= (milli/liter)
	(pS) 	= (picosiemens)
	PI   	= (pi)(1)
}

PARAMETER {

	gmax	    = 700	(pS)	
	Cdur	    = 0.3		(ms)	
    Tmax        =   1   (mM)
    Prel	    =   0.364	

	konDEFAULT	= 0	: choose cooperativity instead of the usual 2*kon followed by kon
	konC0C1 	= 2.5 	(/ms/mM) 
	konC1C2 	= 20 	(/ms/mM)
	kon		    = 5	    (/ms/mM)	 
	koff		= 0.3	(/ms) 		 
	a1		    = 0.35	(/ms)	
	b1		    = 0.1	(/ms)
	a2		    = 0.375	(/ms)	
	b2		    = 8	    (/ms)	
	d1		    = 0.013	(/ms)
	r1		    = 1.3e-4(/ms)
	d2		    = 1.5	(/ms)	
	r2		    = 0.04	(/ms)	 
	
	: In the Westbrook kinetic scheme these rates are zero
	d3		    = 0	(/ms) 		 
	r3		    = 0	(/ms) 		 
	d1d2		= 0	(/ms/mM) 
	r1r2		= 0	(/ms)	

	: presynaptic parameters
	tau_1 		= 3 (ms) 	    < 1e-9, 1e9 >
	tau_rec 	= 35.1 (ms) 	< 1e-9, 1e9 > 	
	tau_facil 	= 10.8 (ms) 	< 0, 1e9 > 	

	U 		    = 0.416 (1) 	< 0, 1 >
	u0 		    = 0 (1) 	    < 0, 1 >	 
		
	Erev		= -65		    (mV)	
}


ASSIGNED {
	v		(mV)	: postsynaptic voltage
	i 		(nA)	: current = g*(v - Erev)
	g 		(pS)	: conductance	
	Open
	T		(mM)		 
	kon1
	kon2
	x
	xview
	yview
	zview
	Pview
}

STATE {	
	C
	CA1
	CA2
	DA1
	DA2
	DA2f
	OA1
	OA2	
}

INITIAL {
	: kinetic states
	C	    =	1
	CA1	    =	0
	CA2	    =	0
	DA1	    =	0
	DA2	    =	0
	DA2f	=	0
	OA1	    =	0  	
	OA2	    =	0
	CA1	    =	0
	CA2	    =	0
	Open	=	0
	  
	T	    =	0 	(mM)

	if (konDEFAULT) {
		kon1 = 2 * kon
		kon2 = kon
	} else {
		kon1 = konC0C1
		kon2 = konC1C2
	}
}

BREAKPOINT {
	SOLVE kstates METHOD sparse
	Open = OA1 + OA2
	g = gmax * Open
	i = (1e-6) * g * ( v - Erev )
}

KINETIC kstates {
	: second row
	~	C  	    <-> 	CA1	(kon1*T,koff)
	~	CA1 	<-> 	CA2	(kon2*T,2*koff) 
	~	CA2	    <->	    DA2f (d3,r3)
	: third row
	~ 	DA1  	<-> 	DA2	(d1d2*T,r1r2)
	: first <=> second row
	~ 	OA1  	<-> 	CA1	(a1,b1)
	~ 	OA2  	<-> 	CA2	(a2,b2)
	: third <=> second row
	~	DA1	    <->	CA1	(r1,d1)
	~	DA2	    <->	CA2	(r2,d2)	
	CONSERVE C+CA1+CA2+DA1+DA2+DA2f+OA1+OA2=1
}

NET_RECEIVE(weight, on, nspike, t0 (ms),y, z, u, tsyn (ms)) {
	INITIAL {
		y = 0
		z = 0
		u = u0
		tsyn = t
		nspike = 1
	}
  	if (flag == 0) { 
		: presynaptic modulation
		nspike = nspike + 1
		if (!on) {
			t0 = t
			on = 1		
			z = z*exp( - (t - tsyn) / tau_rec )	
			z = z + ( y*(exp(-(t - tsyn)/tau_1) - exp(-(t - tsyn)/tau_rec))/((tau_1/tau_rec)-1) )  
			y = y*exp(-(t - tsyn)/tau_1)			
			x = 1-y-z
			
			if (tau_facil > 0) { 
				u = u*exp(-(t - tsyn)/tau_facil)
				u = u + U * ( 1 - u )							
			} else { u = U }
			y = y + x * u
			xview = x	 
			yview = y  
			Pview = u
			T = Tmax * y			
			tsyn = t
						
		}
		net_send(Cdur, nspike)	 
   	}
	if (flag == nspike) { 
			t0 = t
			T = 0
			on = 0
	}
}

