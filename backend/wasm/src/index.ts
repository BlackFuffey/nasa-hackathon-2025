// Types

@unmanaged
class Vec3 { 
    x!: f64;
    y!: f64;
    z!: f64;
}

type Mesh3D = {
    vertices: StaticArray<Vec3>,
    faces: StaticArray<{ 0: i32, 1: i32, 2: i32 }>
}

type Simulation = {
    t: i32,
    state: {
        id: i32,
        origin: Vec3,
        rot_rad: f64,
        spinAxis: Vec3,
        shape: Mesh3D
    }[]
}[]

type Asteroid = {
    id: i32,
    mass_kg: f64, shape: Mesh3D,
    rotr_radS: f64, rot_rad: f64, spinAxis: Vec3, 
    pos: Vec3, vel: Vec3
}


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

// Convert kepler orbit data to cartesian position + velocity
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

// Bowring itration for ECEF -> Geodetic conversion
function ecefToGeodetic(params: {
    pos: Vec3, axis_m: f64, ecc: f64, flat: f64
}): { lat: f64, lon: f64, alt: f64 } {

    const { pos, axis_m, ecc, flat } = params;
    const ecc2 = ecc ** 2;

    const x = pos.x, y = pos.y, z = pos.z;
    const lon = Math.atan2(y, x);
    const r = Math.sqrt(x*x + y*y);
    const b = axis_m * (1.0 - flat);
    const ep2 = (axis_m*axis_m - b*b)/(b*b);
    let lat = Math.atan2(z, r * (1 - ecc2));
    let N: f64 = 0;
    for (let i=0; i<5; i++) {
        N = axis_m / Math.sqrt(1 - ecc2 * Math.sin(lat)*Math.sin(lat));
        const h = r / Math.cos(lat) - N;
        lat = Math.atan2(z + ep2*N*Math.sin(lat), r);
    }
    const h = r / Math.cos(lat) - N;
    return { lat, lon, alt: h };
}

function gravityAt(params: {
    pos: Vec3, eqtgrav_mS2: f64, axis_m: f64, ecc: f64, flat: f64
}): f64 {
    const { lat } = ecefToGeodetic({
        pos: params.pos, axis_m: params.axis_m,
        ecc: params.ecc, flat: params.flat
    });

    return params.eqtgrav_mS2 * 
        (1 + 0.0053024*Math.sin(lat)**2 - 0.0000058*Math.sin(2*lat)**2
    );
}

// Main simulation
export function meteorsim(params: {
    planet: {
        mu: f64, radius: f64, axis_m: f64,
        ecc: f64, rotr_radS: f64,
        eqtgrav_mS2: f64, flat: f64
    },
    orbital: { 
        axis_m: f64, ecc: f64, incl_rad: f64, 
        peri_rad: f64, lon_rad: f64, mean_rad: f64
    },
    asteroid: {
        mass_kg: f64, rotr_radS: f64, 
        spinAxis: Vec3, shape: Mesh3D
    },
    resol_ms: i32
}): Simulation {

    const dt = params.resol_ms as f64 / 1000.0; // timestep in seconds
    const planetR = params.planet.radius;

    let { pos, vel } = keplerToCartesian({ mu: params.planet.mu, orbital: params.orbital });

    let fragId = 0;
    function nextId() { return fragId++; }

    let frags: Asteroid[] = [{
        rotr_radS: params.asteroid.rotr_radS,
        spinAxis: params.asteroid.spinAxis,
        rot_rad: 0.0,

        mass_kg: params.asteroid.mass_kg,
        shape: params.asteroid.shape,

        pos, vel,

        id: nextId()
    }]

    let t: f64 = 0.0;

    let results: Simulation = [];

    // runge-kutta 4th order integration
    function rk4Step(pos: Vec3, vel: Vec3, m: f64, dt: f64): { pos: Vec3, vel: Vec3 } {
        
        // Local acceleration
        function accel(p: Vec3, v: Vec3, m: f64): Vec3 {
            let rmag = norm(p);
            let alt = (rmag - planetR) / 1000.0;
            let rho = airDensity(alt);
            let vmag = norm(v);

            // Drag and gravity
            let Fd = 0.5 * CD * rho * vmag * vmag * Math.pow(m / 2000.0, 2/3);
            let drag = scale(normalize(v), -Fd / m);
            let gravity = scale(normalize(p), -gravityAt({
                pos: p, eqtgrav_mS2: params.planet.eqtgrav_mS2,
                axis_m: params.planet.axis_m, 
                ecc: params.planet.ecc, flat: params.planet.flat
            }))

            return add(gravity, drag);
        }

        let k1v = accel(pos, vel, m);
        let k1r = vel;

        let k2v = accel(add(pos, scale(k1r, dt * 0.5)), add(vel, scale(k1v, dt * 0.5)), m);
        let k2r = add(vel, scale(k1v, dt * 0.5));

        let k3v = accel(add(pos, scale(k2r, dt * 0.5)), add(vel, scale(k2v, dt * 0.5)), m);
        let k3r = add(vel, scale(k2v, dt * 0.5));

        let k4v = accel(add(pos, scale(k3r, dt)), add(vel, scale(k3v, dt)), m);
        let k4r = add(vel, scale(k3v, dt));

        let newVel = add(vel, scale(
            add(add(k1v, scale(add(k2v, k3v), 2.0)), k4v), dt / 6.0
        ));
        let newPos = add(pos, scale(
            add(add(k1r, scale(add(k2r, k3r), 2.0)), k4r), dt / 6.0
        ));
        return { pos: newPos, vel: newVel };
    }

    function integrate({ id, rot_rad, mass_kg, pos, vel, shape, spinAxis, rotr_radS }: Asteroid): Asteroid[] {
        let rmag = norm(pos);
        let alt = (rmag - planetR) / 1000.0; // km altitude

        if (alt < 0.0) return []; // ground impact
        if (mass_kg < 1.0) return [];   // burned up

        // Atmosphere
        let rho = airDensity(alt);
        let vmag = norm(vel);

        // Mass loss
        let massPrev = mass_kg;

        let heatFlux = 0.5 * rho * vmag * vmag * vmag;
        let dm = -heatFlux * 1e-8 * dt; // tuning coefficient

        mass_kg = Math.max(0.0, mass_kg + dm);

        // Adjust shape
        const dir = normalize(vel);

        for (let i = 0; i < shape.vertices.length; i++) {
            const v = shape.vertices[i];
            // bias shrinkage slightly toward direction of flight (erosion asymmetry)
            const bias = 1.0 - 0.05 * dot(normalize(v), dir); // 5% faster erosion on leading side
            shape.vertices[i] = scale(v, Math.pow(mass_kg/massPrev, 1.0/3.0) * bias);
        }

        // Integrate
        const integ = rk4Step(pos, vel, mass_kg, dt);
        vel = integ.vel;
        pos = integ.pos;

        // Rotation update
        rot_rad += rotr_radS * dt;
        if (rot_rad > 2.0*Math.PI) rot_rad -= 2.0*Math.PI;

        return [{ id, rot_rad, mass_kg, pos, vel, shape, spinAxis, rotr_radS }];
    }

    // Integrate until hit ground or burn-up
    for (let step = 0; step < 200000; step++) {

        for (let i = 0; i < frags.length; i++) {
            integrate(frags[i])
        }
        
        results.push({
            t: (t * 1000) as i32,
            state: frags.map(frag => ({
                id: frag.id,
                origin: frag.pos,
                rot_rad: frag.rot_rad,
                spinAxis: frag.spinAxis,
                shape: frag.shape
            }))
        })

        t += dt;
    }

    return results;
}
