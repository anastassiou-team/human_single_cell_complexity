: Sinewave current clamp

NEURON {
	POINT_PROCESS Chirp
	RANGE del, dur, amp, freq, fdt, i
	ELECTRODE_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	PI = (pi) (1)
}


PARAMETER {
	del (ms)
	dur (ms) < 0, 1e9 >
	amp (nA)
	freq (Hz)
	fdt (1)
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
:		i = amp*sin(2*PI*freq/1000 * (fdt^(t-del) -1 ) / log(fdt)) *  fdt^((t-del)*0.66)
		i = amp*sin(2*PI*freq/1000 * (fdt^(t-del) -1 ) / log(fdt))
	} 
	else {
		i = 0
	}
}