package pnormGO

import (
	"math"
)

const (
	M1Sqrt2PI  = 0.398942280401432677939946059934
	MSqrt32    = 5.656854249492380195206754896838
	DBLMin     = math.SmallestNonzeroFloat64
	DBLEpsilon = 2.2204460492503131e-16
)

func _R_D__0(logP bool) float64 {
	if logP {
		return math.Inf(-1)
	} else {
		return 0
	}
}

func _R_D__1(logP bool) float64 {
	if logP {
		return 0
	} else {
		return 1
	}
}

func _R_DT_0(lowerTail bool, logP bool) float64 {
	if lowerTail {
		return _R_D__0(logP)
	} else {
		return _R_D__1(logP)
	}
}

func _R_DT_1(lowerTail bool, logP bool) float64 {
	if lowerTail {
		return _R_D__1(logP)
	} else {
		return _R_D__0(logP)
	}
}

func _DoDel(ccum *float64, cum *float64, logP bool, X float64, temp float64, upper bool, lower bool, x float64) {
	xsq := math.Ldexp(math.Trunc(math.Ldexp(X, 4)), -4)
	del := (X - xsq) * (X + xsq)
	if logP {
		*cum = -xsq*math.Ldexp(xsq, -1) - math.Ldexp(del, -1) + math.Log(temp)
		if (lower && x > 0) || (upper && x <= 0) {
			*ccum = math.Log1p(-math.Exp(-xsq*math.Ldexp(xsq, -1)) * math.Exp(-math.Ldexp(del, -1)) * temp)
		}
	} else {
		*cum = math.Exp(-xsq*math.Ldexp(xsq, -1)) * math.Exp(-math.Ldexp(del, -1)) * temp
		*ccum = 1.0 - *cum
	}
}

func pnormBoth(x float64, cum *float64, ccum *float64, iTail bool, logP bool) {
	a := [...]float64{2.2352520354606839287, 161.02823106855587881, 1067.6894854603709582, 18154.981253343561249, 0.065682337918207449113}
	b := [...]float64{47.20258190468824187, 976.09855173777669322, 10260.932208618978205, 45507.789335026729956}
	c := [...]float64{0.39894151208813466764, 8.8831497943883759412, 93.506656132177855979, 597.27027639480026226, 2494.5375852903726711, 6848.1904505362823326, 11602.651437647350124, 9842.7148383839780218, 1.0765576773720192317e-8}
	d := [...]float64{22.266688044328115691, 235.38790178262499861, 1519.377599407554805, 6485.558298266760755, 18615.571640885098091, 34900.952721145977266, 38912.003286093271411, 19685.429676859990727}
	p := [...]float64{0.21589853405795699, 0.1274011611602473639, 0.022235277870649807, 0.001421619193227893466, 2.9112874951168792e-5, 0.02307344176494017303}
	q := [...]float64{1.28426009614491121, 0.468238212480865118, 0.0659881378689285515, 0.00378239633202758244, 7.29751555083966205e-5}

	var xden, xnum, temp, xsq float64
	const min = DBLMin
	const eps = DBLEpsilon * 0.5
	lower := !iTail
	upper := iTail
	y := math.Abs(x)

	if y <= 0.67448975 {
		if y > eps {
			xsq = x * x
			xnum = a[4] * xsq
			xden = xsq
			for i := 0; i < 3; i++ {
				xnum = (xnum + a[i]) * xsq
				xden = (xden + b[i]) * xsq
			}
		} else {
			xnum = 0.0
			xden = 0.0
		}
		temp = x * (xnum + a[3]) / (xden + b[3])
		if lower {
			*cum = 0.5 + temp
		}
		if upper {
			*ccum = 0.5 - temp
		}
		if logP {
			if lower {
				*cum = math.Log(*cum)
			}
			if upper {
				*ccum = math.Log(*ccum)
			}
		}
	} else if y <= MSqrt32 {
		xnum = c[8] * y
		xden = y
		for i := 0; i < 7; i++ {
			xnum = (xnum + c[i]) * y
			xden = (xden + d[i]) * y
		}
		temp = (xnum + c[7]) / (xden + d[7])
		_DoDel(ccum, cum, logP, y, temp, upper, lower, x)
		if x > 0 {
			temp = *cum
			if lower {
				*cum = *ccum
			}
			*ccum = temp
		}
	} else if (logP && y < 1e170) || (lower && -37.5193 < x && x < 8.2924) || (upper && -8.2924 < x && x < 37.5193) {
		xsq = 1 / (x * x)
		xnum = p[5] * xsq
		xden = xsq
		for i := 0; i < 4; i++ {
			xnum = (xnum + p[i]) * xsq
			xden = (xden + q[i]) * xsq
		}
		temp = xsq * (xnum + p[4]) / (xden + q[4])
		temp = (M1Sqrt2PI - temp) / y
		_DoDel(ccum, cum, logP, y, temp, upper, lower, x)
		if x > 0 {
			temp = *cum
			if lower {
				*cum = *ccum
			}
			*ccum = temp
		}
	} else {
		if x > 0 {
			*cum = _R_D__1(logP)
			*ccum = _R_D__0(logP)
		} else {
			*cum = _R_D__0(logP)
			*ccum = _R_D__1(logP)
		}
	}
	if logP {
		if *cum > -min {
			*cum = 0
		}
		if *ccum > -min {
			*ccum = 0
		}
	} else {
		if *cum < min {
			*cum = 0
		}
		if *ccum < min {
			*ccum = 0
		}
	}
}

func pnorm5(x float64, mu float64, sigma float64, lowerTail bool, logP bool) float64 {
	if math.IsNaN(x) || math.IsNaN(mu) || math.IsNaN(sigma) {
		return x + mu + sigma
	}
	if math.IsInf(x, 0) && mu == x {
		return math.NaN()
	}
	if sigma <= 0 {
		if sigma < 0 {
			return math.NaN()
		}
		if x < mu {
			return _R_DT_0(lowerTail, logP)
		} else {
			return _R_DT_1(lowerTail, logP)
		}
	}
	p := (x - mu) / sigma
	var cp float64
	if math.IsInf(p, 0) {
		if x < mu {
			return _R_DT_0(lowerTail, logP)
		} else {
			return _R_DT_1(lowerTail, logP)
		}
	}
	x = p
	pnormBoth(x, &p, &cp, !lowerTail, logP)
	if lowerTail {
		return p
	} else {
		return cp
	}
}

func Pnorm(q float64, mean float64, sd float64, lowerTail bool, logP bool) float64 {
	return pnorm5(q, mean, sd, lowerTail, logP)
}
