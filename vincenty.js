const wgs84 = {
	earth: {a: 6378137.0, f: 1/298.257223563, b: 6356752.314245}//🜃
};

// https://en.wikipedia.org/wiki/Vincenty%27s_formulae#Inverse_problem
function distVincenty(lat1, lon1, lat2, lon2, params=wgs84.earth) {
	const a = params.a;
	const b = params.b;
	const f = params.f;

	//U1, U2: Reduced latitudes (i.e. latitudes on the auxiliary sphere) 
	const U1 = Math.atan((1 - f) * Math.tan(lat1 * Math.PI / 180));
	const U2 = Math.atan((1 - f) * Math.tan(lat2 * Math.PI / 180));

	//(precompute their sine and cosines)
	const sin_U1 = Math.sin(U1);
	const sin_U2 = Math.sin(U2);
	const cos_U1 = Math.cos(U1);
	const cos_U2 = Math.cos(U2);

	//L: Difference (in radians) in longitude of two points
	const L = (lon2 - lon1) * Math.PI / 180;

	//λ: Difference (in radians) in longitude of the points on the auxiliary sphere
	let la = L;
	let sin_la, cos_la;//Precomputation of sin(λ) and cos(λ)

	//α: Forward azimuth of the geodesic at the equator, if it were extended that far
	let aa;
	let sin_aa, cos2_aa;//Precomputation of sin(α) and cos²(α)

	//α_1, α_2: forward azimuths at the points
	let aa_1, aa_2;

	//s: Ellipsoidal distance between the two points
	let s = Infinity;

	//σ: Angular separation between points
	let sa = Infinity;
	let sin_sa, cos_sa;//Precomputation of sin(σ) and cos(σ)
	let d_sa;//Δσ

	//σ_1: angular separation between the point and the equator
	let s1 = Infinity;

	//σ_m: angular separation between the midpoint of the line and the equator
	let sam = Infinity;
	let cos_2sam, cos2_2sam;//Precomputation of cos(2 * σ_m) and cos²(2 * σ_m)

	let laPrev = -Infinity;
	let A, B, C, u2;
	while( (la - laPrev) > 1e-12 ) {

		//Precompute sin and cos of λ
		cos_la = Math.cos(la);
		sin_la = Math.sin(la);

		// σ is not evaluated directly from sin σ or cos σ
		// to preserve numerical accuracy near the poles and equator 
		sin_sa = Math.sqrt(
			( cos_U2 * sin_la ) ** 2 +
			( cos_U1 * sin_U2  -  sin_U1 * cos_U2 * cos_la ) ** 2
		);
		cos_sa = sin_U1 * sin_U2  +  cos_U1 * cos_U2 * cos_la;
		sa = Math.atan2(sin_sa, cos_sa);

		// If sin(σ) == 0 the value of sin α is indeterminate:
		// it represents an end point coincident with, or
		// diametrically opposed to, the start point. 
		sin_aa = cos_U1 * cos_U2 * sin_la / sin_sa;

		//Precompute this intermediate value
		cos2_aa = (1 - sin_aa ** 2);

		cos_2sam = cos_sa - 2 * sin_U1 * sin_U2 / cos2_aa;
		cos2_2sam = cos_2sam ** 2;

		C = (f/16) * (cos2_aa) * (4 + f * (4 - 3 * cos2_aa));

		laPrev = la;
		la = L + (1 - C) * f * sin_aa *\
                     (sa + C * sin_sa *\
		      (cos_2sam + C * cos_sa *\
		       (-1 + 2 * cos2_2sam)));
	}

	u2 = cos2_aa * ( (a ** 2 - b ** 2) / b ** 2);

	A = 1 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)));

	B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)));

	d_sa = B * sin_sa *\
	       (cos_2sam + (B/4) *\
	        (cos_sa * (-1 + 2 * cos2_2sam) -\
	         (B/6) *\
	         (cos_2sam) *\
	         (-3 + 4 * sin_sa ** 2) *\
	         (-3 + 4 * cos2_2sam)));

	s = b * A * (sa - d_sa);

	//aa_1 = Math.atan2(cos_U2 * sin_la, cos_U1 * sin_U2 - sin_U1 * cos_U2 * cos_la);

	//aa_2 = Math.atan2(cos_U1 * sin_la, -sin_U1 * cos_U2);

	return s;
}

module.exports = {distVincenty};
