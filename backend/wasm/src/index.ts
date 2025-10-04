// Types
type Vec3 = { x: f64, y: f64, z: f64 };

type Mesh3D = {
    vertices: StaticArray<Vec3>,
    faces: StaticArray<{ 0: i32, 1: i32, 2: i32 }>
}

type State = {
    id: i32,
    origin: Vec3,
    rot: f64,
    spinAxis: Vec3,
    shape: Mesh3D
}[]

// Vector Math Helpers
function vec3(x: f64, y: f64, z: f64): Vec3 { return { x, y, z }; }

function add(a: Vec3, b: Vec3): Vec3 { return { x: a.x + b.x, y: a.y + b.y, z: a.z + b.z }; }

function sub(a: Vec3, b: Vec3): Vec3 { return { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z }; }

function scale(a: Vec3, s: f64): Vec3 { return { x: a.x * s, y: a.y * s, z: a.z * s }; }

function dot(a: Vec3, b: Vec3): f64 { return a.x*b.x + a.y*b.y + a.z*b.z; }

function norm(a: Vec3): f64 { return Math.sqrt(dot(a,a)); }

function normalize(a: Vec3): Vec3 { let n = norm(a); return scale(a, 1.0/n); }


// Placeholder exponential atmosphere (replace with NRLMSISE-00)
function airDensity(alt_km: f64): f64 {
    if (alt_km < 0.0) return 1.225; // sea level
    return 1.225 * Math.exp(-alt_km / 7.5); // scale height ~7.5 km
}

// Drag coefficient (replace with simulated drag)
const CD: f64 = 1.3;

function keplerToCartesian(params: {
    mu: f64,
    orbital: {
        axis_m: f64,
        ecc: f64,
        incl_rad: f64,
        peri_rad: f64,
        lon_rad: f64,
        mean_rad: f64
    }
}): { pos: Vec3, vel: Vec3 } {
    const a = params.orbital.axis_m;
    const e = params.orbital.ecc;
    const i = params.orbital.incl_rad;
    const omega = params.orbital.lon_rad;   // Longitude of ascending node
    const w = params.orbital.peri_rad;      // Argument of periapsis
    const mu = params.mu;

    // Periapsis distance
    const rp = a * (1.0 - e);
    // Velocity magnitude at periapsis
    const vp = Math.sqrt(mu * (1.0 + e) / (a * (1.0 - e)));

    // In orbital plane (2D)
    const r_orb = vec3(rp, 0.0, 0.0);
    const v_orb = vec3(0.0, vp, 0.0);

    // Rotation matrices
    const cosO = Math.cos(omega), sinO = Math.sin(omega);
    const cosi = Math.cos(i), sini = Math.sin(i);
    const cosw = Math.cos(w), sinw = Math.sin(w);

    // Rotation: Rz(Ω) * Rx(i) * Rz(ω)
    const R11 = cosO*cosw - sinO*sinw*cosi;
    const R12 = -cosO*sinw - sinO*cosw*cosi;
    const R13 = sinO*sini;
    const R21 = sinO*cosw + cosO*sinw*cosi;
    const R22 = -sinO*sinw + cosO*cosw*cosi;
    const R23 = -cosO*sini;
    const R31 = sinw*sini;
    const R32 = cosw*sini;
    const R33 = cosi;

    function rot(v: Vec3): Vec3 {
        return vec3(
            R11*v.x + R12*v.y + R13*v.z,
            R21*v.x + R22*v.y + R23*v.z,
            R31*v.x + R32*v.y + R33*v.z
        );
    }

    return { pos: rot(r_orb), vel: rot(v_orb) };
}

// Main simulation
export function meteorsim(params: {
    planet: { gravity: f64, radius: f64 },
    orbital: { axis_m: f64, ecc: f64, incl_rad: f64, peri_rad: f64, lon_rad: f64, mean_rad: f64 },
    asteroid: { mass_kg: f64, rot_radS: f64, spinAxis: Vec3, shape: Mesh3D },
    resol_ms: i32
}): State[] {

    const mu = params.planet.gravity * Math.pow(params.planet.radius, 2);
    const dt = params.resol_ms as f64 / 1000.0; // timestep in seconds
    const planetR = params.planet.radius;
    const mass0 = params.asteroid.mass_kg;

    let { pos, vel } = keplerToCartesian({ mu, orbital: params.orbital });

    let m = mass0;
    let t: f64 = 0.0;
    let rot = 0.0;

    let results = new Array<State>();

    // Integrate until hit ground or burn-up
    for (let step = 0; step < 200000; step++) {
        let rmag = norm(pos);
        let alt = (rmag - planetR) / 1000.0; // km altitude

        if (alt < 0.0) break; // ground impact
        if (m < 1.0) break;   // burned up

        // Atmosphere
        let rho = airDensity(alt);
        let vmag = norm(vel);

        // Forces
        let Fd = 0.5 * CD * rho * vmag * vmag * Math.pow(m / 2000.0, 2/3); // area scaling
        let drag = scale(normalize(vel), -Fd / m);
        let gravity = scale(normalize(pos), -params.planet.gravity);

        // Ablation (mass loss)
        let heatFlux = 0.5 * rho * vmag * vmag * vmag;
        let dm = -heatFlux * 1e-8 * dt; // tuning coefficient

        m = Math.max(0.0, m + dm);

        // Integrate
        let acc = add(gravity, drag);
        vel = add(vel, scale(acc, dt));
        pos = add(pos, scale(vel, dt));

        // Rotation update
        rot += params.asteroid.rot_radS * dt;
        if (rot > 2.0*Math.PI) rot -= 2.0*Math.PI;

        // Record every Nth state for memory
        if (step % 10 == 0) {
            let s: State = [{
                id: 0,
                origin: pos,
                rot,
                spinAxis: params.asteroid.spinAxis,
                shape: params.asteroid.shape
            }];
            results.push(s);
        }

        t += dt;
    }

    return results;
}
