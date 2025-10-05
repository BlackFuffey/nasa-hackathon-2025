import { vec3, Vec3 } from './vector'

export class KeplerOrbit {
    axis_m!: f64;
    ecc!: f64;
    incl_rad!: f64;
    peri_rad!: f64;
    lon_rad!: f64;
    mean_rad!: f64
}

export class CartesianMechState {
    pos!: Vec3; vel!: Vec3
}

// Convert kepler orbit data to cartesian position + velocity
class KeplerToCartesianParams {
    mu!: f64;
    orbital!: KeplerOrbit;
}

// Rotation: Rz(Ω) * Rx(i) * Rz(ω)
class Rotations {
    R11!: f64;
    R12!: f64;
    R13!: f64;
    R21!: f64;
    R22!: f64;
    R23!: f64;
    R31!: f64;
    R32!: f64;
    R33!: f64;
}
export function keplerToCartesian(params: KeplerToCartesianParams): CartesianMechState {
    const a = params.orbital.axis_m;
    const e = params.orbital.ecc;
    const i = params.orbital.incl_rad;
    const omega = params.orbital.lon_rad;   
    const w = params.orbital.peri_rad;     
    const mu = params.mu;

    const rp = a * (1.0 - e);
    const vp = Math.sqrt(mu * (1.0 + e) / (a * (1.0 - e)));

    const r_orb = vec3(rp, 0.0, 0.0);
    const v_orb = vec3(0.0, vp, 0.0);

    const cosO = Math.cos(omega), sinO = Math.sin(omega);
    const cosi = Math.cos(i), sini = Math.sin(i);
    const cosw = Math.cos(w), sinw = Math.sin(w);


    const rotVals: Rotations = {

        R11: cosO*cosw - sinO*sinw*cosi,
        R12: -cosO*sinw - sinO*cosw*cosi,
        R13: sinO*sini,
        R21: sinO*cosw + cosO*sinw*cosi,
        R22: -sinO*sinw + cosO*cosw*cosi,
        R23: -cosO*sini,
        R31: sinw*sini,
        R32: cosw*sini,
        R33: cosi,
    }

    function rot(v: Vec3, rot: Rotations): Vec3 {
        return vec3(
            rot.R11*v.x + rot.R12*v.y + rot.R13*v.z,
            rot.R21*v.x + rot.R22*v.y + rot.R23*v.z,
            rot.R31*v.x + rot.R32*v.y + rot.R33*v.z
        );
    }

    return { pos: rot(r_orb, rotVals), vel: rot(v_orb, rotVals) };
}
