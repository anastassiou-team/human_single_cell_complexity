: Sinewave current clamp

NEURON {
	POINT_PROCESS SineClamp
	RANGE del, dur, amp, freq, off, i
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
	off (nA)
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
		i = amp*sin(2*PI*t*freq/1000)+off
	} 
	else {
		i = 0
	}
}