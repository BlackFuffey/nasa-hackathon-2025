import { Vec3 } from './vector'

class PlanetCordinates {
    lat!: f64; lon!: f64; alt!: f64
}

// Bowring itration for ECEF -> Geodetic conversion
class EcefToGeodeticParams {
    pos!: Vec3; axis_m!: f64; ecc!: f64; flat!: f64;
}
export function ecefToGeodetic(params: EcefToGeodeticParams): PlanetCordinates {
    const ecc2 = params.ecc ** 2;

    const x = params.pos.x, y = params.pos.y, z = params.pos.z;
    const lon = Math.atan2(y, x);
    const r = Math.sqrt(x*x + y*y);
    const b = params.axis_m * (1.0 - params.flat);
    const ep2 = (params.axis_m*params.axis_m - b*b)/(b*b);
    let lat = Math.atan2(z, r * (1 - ecc2));
    let N: f64 = 0;
    for (let i=0; i<5; i++) {
        N = params.axis_m / Math.sqrt(1 - ecc2 * Math.sin(lat)*Math.sin(lat));
        lat = Math.atan2(z + ep2*N*Math.sin(lat), r);
    }
    const h = r / Math.cos(lat) - N;
    return { lat, lon, alt: h };
}

class GravityAtParams {
    pos!: Vec3; eqtgrav_mS2!: f64; axis_m!: f64; ecc!: f64; flat!: f64;
}

export function gravityAt(params: GravityAtParams): f64 {
    const lat = ecefToGeodetic({
        pos: params.pos, axis_m: params.axis_m,
        ecc: params.ecc, flat: params.flat
    }).lat;

    return params.eqtgrav_mS2 * 
        (1 + 0.0053024*Math.sin(lat)**2 - 0.0000058*Math.sin(2*lat)**2
        );
}

