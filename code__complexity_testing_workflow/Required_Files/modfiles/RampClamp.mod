: Ramp current clamp

NEURON {
	POINT_PROCESS RampClamp
	RANGE del, dur, alpha, beta, i
	ELECTRODE_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	PI = (pi) (1)
}


PARAMETER {
	del (ms)
	dur (ms) < 0, 1e9 >
	alpha
	beta
}

ASSIGNED {
	i (nA)
}

INITIAL {
	i = 0
}

BREAKPOINT {
	at_time(del)
	at_time(del+dur)
	if (t < del + dur && t > del) {
		i = alpha*(t-del) + beta
	} 
	else {
		i = 0
	}
}