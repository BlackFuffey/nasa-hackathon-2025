// Types

class Vec3 { 
    x!: f64;
    y!: f64;
    z!: f64;
}

class i32Vec3 { 
    x!: f64;
    y!: f64;
    z!: f64;
}

class Mat3x3 {
    a!: Vec3;
    b!: Vec3;
    c!: Vec3;
}

class Mesh3D {
    vertices!: StaticArray<Vec3>;
    faces!: StaticArray<i32Vec3>;
}

class FragmentState {
    id!: i32;
    origin!: Vec3;
    rot_rad!: f64;
    spinAxis!: Vec3;
    shape!: Mesh3D;
}

class SimulationFrame {
    t_ms!: i32;
    state!: StaticArray<FragmentState>
}

class AsteroidInitial {
    mass_kg!: f64; shape!: Mesh3D; strength_MPa!: f64
    rotr_radS!: f64; spinAxis!: Vec3; 
}

class Asteroid extends AsteroidInitial {
    id!: i32;
    rot_rad!: f64;
    pos!: Vec3; vel!: Vec3
}

class PlanetConsts {
    mu!: f64; radius!: f64; axis_m!: f64;
    ecc!: f64; rotr_radS!: f64;
    eqtgrav_mS2!: f64; flat!: f64
}

class KeplerOrbit {
    axis_m!: f64;
    ecc!: f64;
    incl_rad!: f64;
    peri_rad!: f64;
    lon_rad!: f64;
    mean_rad!: f64
}

class PlanetCordinates {
    lat!: f64; lon!: f64; alt!: f64
}

class CartesianParams {
    pos!: Vec3; vel!: Vec3
}

// Vector Math Helpers
function vec3(x: f64, y: f64, z: f64): Vec3 { return { x, y, z }; }

function add(a: Vec3, b: Vec3): Vec3 { return { x: a.x + b.x, y: a.y + b.y, z: a.z + b.z }; }

function sub(a: Vec3, b: Vec3): Vec3 { return { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z }; }

function scale(a: Vec3, s: f64): Vec3 { return { x: a.x * s, y: a.y * s, z: a.z * s }; }

function dot(a: Vec3, b: Vec3): f64 { return a.x*b.x + a.y*b.y + a.z*b.z; }

function norm(a: Vec3): f64 { return Math.sqrt(dot(a,a)); }

function normalize(a: Vec3): Vec3 { let n = norm(a); return scale(a, 1.0/n); }

// Matrix math helpers
function rotationMat(axis: Vec3, ang_rad: f64): Mat3x3 {
    axis = normalize(axis);
    const c = Math.cos(ang_rad);
    const s = Math.sin(ang_rad);
    const t = 1 - c;

    return {
        a: {
            x: t*axis.x*axis.x + c,
            y: t*axis.x*axis.y - s*axis.z,
            z: t*axis.x*axis.z + s*axis.y
        },
        b: {
            x: t*axis.x*axis.y + s*axis.z,  
            y: t*axis.y*axis.y + c,         
            z: t*axis.y*axis.z - s*axis.x 
        },
        c: {
            x: t*axis.x*axis.z - s*axis.y, 
            y: t*axis.y*axis.z + s*axis.x, 
            z: t*axis.z*axis.z + c  
        },
    };
}

function rotVec(M: Mat3x3, v: Vec3): Vec3 {
    return vec3(
        M.a.x*v.x + M.a.y*v.y + M.a.z*v.z,
        M.b.x*v.x + M.b.y*v.y + M.b.z*v.z,
        M.c.x*v.x + M.c.y*v.y + M.c.z*v.z
    );
}

// Placeholder exponential atmosphere (replace with NRLMSISE-00)
function airDensity(alt_km: f64): f64 {
    if (alt_km < 0.0) return 1.225; // sea level
    return 1.225 * Math.exp(-alt_km / 7.5); // scale height ~7.5 km
}

// Drag coefficient (replace with simulated drag)
const CD: f64 = 1.3;

// Convert kepler orbit data to cartesian position + velocity
class KeplerToCartesian {
    mu!: f64;
    orbital!: KeplerOrbit;
}

let R11: f64, R12: f64, R13: f64;
let R21: f64, R22: f64, R23: f64;
let R31: f64, R32: f64, R33: f64;

function keplerToCartesian(params: KeplerToCartesian): CartesianParams {
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
    R11 = cosO*cosw - sinO*sinw*cosi;
    R12 = -cosO*sinw - sinO*cosw*cosi;
    R13 = sinO*sini;
    R21 = sinO*cosw + cosO*sinw*cosi;
    R22 = -sinO*sinw + cosO*cosw*cosi;
    R23 = -cosO*sini;
    R31 = sinw*sini;
    R32 = cosw*sini;
    R33 = cosi;

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
class EcefToGeodeticParams {
    pos!: Vec3; axis_m!: f64; ecc!: f64; flat!: f64;
}
function ecefToGeodetic(params: EcefToGeodeticParams): PlanetCordinates {
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
function gravityAt(params: GravityAtParams): f64 {
    const lat = ecefToGeodetic({
        pos: params.pos, axis_m: params.axis_m,
        ecc: params.ecc, flat: params.flat
    }).lat;

    return params.eqtgrav_mS2 * 
        (1 + 0.0053024*Math.sin(lat)**2 - 0.0000058*Math.sin(2*lat)**2
        );
}

// Main simulation
class MeteorsimParams {
    planet!: PlanetConsts;
    orbital!: KeplerOrbit;
    asteroid!: AsteroidInitial;
    resol_ms!: i32;
}

let fragId = 0;
export function meteorsim(params: MeteorsimParams): Array<SimulationFrame> {

    const dt = params.resol_ms as f64 / 1000.0; // timestep in seconds
    const planetR = params.planet.radius;

    let movement = keplerToCartesian({ mu: params.planet.mu, orbital: params.orbital });

    function nextId(): i32 { return fragId++; }

    const frags = new Map<i32, Asteroid>();
    frags.set(0, {
        rotr_radS: params.asteroid.rotr_radS,
        spinAxis: params.asteroid.spinAxis,
        rot_rad: 0.0,

        mass_kg: params.asteroid.mass_kg,
        shape: params.asteroid.shape,
        strength_MPa: params.asteroid.strength_MPa,

        pos: movement.pos, vel: movement.vel,

        id: 0
    })

    let t: f64 = 0.0;

    let results = new Array<SimulationFrame>();

    function compute(params: Asteroid, dt: f64): StaticArray<Asteroid> {
        let rmag = norm(params.pos);
        let alt = (rmag - planetR) / 1000.0; // km altitude

        if (alt < 0.0) return []; // ground impact
        if (params.mass_kg < 1.0) return [];   // burned up

        // Atmosphere
        let rho = airDensity(alt);
        let vmag = norm(params.vel);

        // Mass loss
        let massPrev = params.mass_kg;

        let heatFlux = 0.5 * rho * vmag * vmag * vmag;
        let dm = -heatFlux * 1e-8 * dt; // tuning coefficient

        params.mass_kg = Math.max(0.0, params.mass_kg + dm);

        // Fragmentation 
        // O(sobbing terrible)
        if (0.5 * rho * vmag * vmag > params.strength_MPa * 1e6) {

            // --- 1. Compute aerodynamic stress direction (shock front) ---
            let stressDir = vec3(0.0, 0.0, 0.0);
            let maxStress = 0.0;
            const vdir = normalize(params.vel);

            for (let i = 0; i < params.shape.faces.length; i++) {
                const f = params.shape.faces[i];
                const a = params.shape.vertices[f.x];
                const b = params.shape.vertices[f.y];
                const c = params.shape.vertices[f.z];

                // Face normal
                const ab = sub(b, a);
                const ac = sub(c, a);
                const n = normalize(vec3(
                    ab.y * ac.z - ab.z * ac.y,
                    ab.z * ac.x - ab.x * ac.z,
                    ab.x * ac.y - ab.y * ac.x
                ));

                // Pressure-induced stress along local normal
                const stress = 0.5 * rho * vmag * vmag * Math.max(0.0, dot(n, vdir));
                if (stress > maxStress) {
                    maxStress = stress;
                    stressDir = n;
                }
            }

            // --- 2. Define fracture plane through centroid ---
            let cx = 0.0, cy = 0.0, cz = 0.0;
            for (let i = 0; i < params.shape.vertices.length; i++) {
                const v = params.shape.vertices[i];
                cx += v.x; cy += v.y; cz += v.z;
            }
            const center = vec3(cx / params.shape.vertices.length, cy / params.shape.vertices.length, cz / params.shape.vertices.length);
            const nfract = normalize(stressDir);

            // --- 3. Split mesh into two fragments along the fracture plane ---
            const vA: Vec3[] = [], vB: Vec3[] = [];
            const fA: i32Vec3[] = [], fB: i32Vec3[] = [];

            for (let i = 0; i < params.shape.vertices.length; i++) {
                const v = params.shape.vertices[i];
                const side = dot(sub(v, center), nfract);
                if (side >= 0.0) vA.push(v);
                    else vB.push(v);
            }

            for (let i = 0; i < params.shape.faces.length; i++) {
                const f = params.shape.faces[i];
                let countA = 0;
                if (vA.includes(params.shape.vertices[f.x])) countA++;
                if (vA.includes(params.shape.vertices[f.y])) countA++;
                if (vA.includes(params.shape.vertices[f.z])) countA++;
                if (countA >= 2) fA.push(f);
                    else fB.push(f);
            }

            const mesh1: Mesh3D = {
                vertices: StaticArray.fromArray(vA),
                faces: StaticArray.fromArray(fA)
            };
            const mesh2: Mesh3D = {
                vertices: StaticArray.fromArray(vB),
                faces: StaticArray.fromArray(fB)
            };

            // --- 4. Rotation-aware fragmentation ---
            const rotMatrix = rotationMat(params.spinAxis, params.rot_rad);

            // Rotate fragment meshes into world orientation
            for (let i = 0; i < mesh1.vertices.length; i++) { 
                const v = mesh1.vertices[i];
                const vr = rotVec(rotMatrix, v); 
                v.x = vr.x; v.y = vr.y; v.z = vr.z; 
            }
            for (let i = 0; i < mesh2.vertices.length; i++) { 
                const v = mesh2.vertices[i];
                const vr = rotVec(rotMatrix, v); 
                v.x = vr.x; v.y = vr.y; v.z = vr.z; 
            }

            const newSpinAxis = normalize(params.spinAxis);
            const newRotRadS = params.rotr_radS;

            // --- 5. Momentum-conserving separation ---
            const mass1 = params.mass_kg * 0.5;
            const mass2 = params.mass_kg * 0.5;
            const sepVel = scale(nfract, vmag * 0.02);

            // --- 6. Fill fracture surface ---
            const eps = 1e-6;
            const fractureVerts: Vec3[] = [];
            for (let i = 0; i < params.shape.vertices.length; i++) {
                const v = params.shape.vertices[i];
                const dist = dot(sub(v, center), nfract);
                if (Math.abs(dist) < eps * 100.0) {
                    fractureVerts.push(v);
                }
            }

            let fx = 0.0, fy = 0.0, fz = 0.0;
            for (let i = 0; i < fractureVerts.length; i++) {
                const v = fractureVerts[i];
                fx += v.x; fy += v.y; fz += v.z;
            }
            const fcenter = vec3(fx / fractureVerts.length, fy / fractureVerts.length, fz / fractureVerts.length);

            function hashNoise(v: Vec3): f64 {
                const s = Math.sin(dot(v, vec3(12.9898, 78.233, 45.164)) * 43758.5453);
                return s - Math.floor(s);
            }

            const offsetMag = norm(sub(fcenter, center)) * 0.01;
            const fractureInnerA: Vec3[] = [];
            const fractureInnerB: Vec3[] = [];
            for (let i = 0; i < fractureVerts.length; i++) {
                const v = fractureVerts[i];
                const n = hashNoise(v);
                const offset = scale(nfract, (n - 0.5) * offsetMag);
                fractureInnerA.push(add(v, offset));
                fractureInnerB.push(sub(v, offset));
            }

            const fTriA: i32Vec3[] = [];
            const fTriB: i32Vec3[] = [];
            for (let i = 0; i < fractureVerts.length; i++) {
                const i1 = i;
                const i2 = (i + 1) % fractureVerts.length;
                fTriA.push({ x: i1, y: i2, z: fractureVerts.length + i1 });
                fTriB.push({ x: i1, y: i2, z: fractureVerts.length + i1 });
            }

            for (let i = 0; i < fractureInnerA.length; i++) vA.push(fractureInnerA[i]);
            for (let i = 0; i < fractureInnerB.length; i++) vB.push(fractureInnerB[i]);
            for (let i = 0; i < fTriA.length; i++) fA.push(fTriA[i]);
            for (let i = 0; i < fTriB.length; i++) fB.push(fTriB[i]);

            const frag1 = compute({
                id: nextId(),
                rot_rad: 0.0,
                mass_kg: mass1,
                pos: params.pos,
                vel: add(params.vel, sepVel),
                shape: mesh1,
                spinAxis: newSpinAxis,
                rotr_radS: newRotRadS,
                strength_MPa: params.strength_MPa
            }, dt);
            
            const frag2 = compute({
                id: nextId(),
                rot_rad: 0.0,
                mass_kg: mass2,
                pos: params.pos,
                vel: sub(params.vel, sepVel),
                shape: mesh2,
                spinAxis: newSpinAxis,
                rotr_radS: newRotRadS,
                strength_MPa: params.strength_MPa
            }, dt);

            const result = new StaticArray<Asteroid>(frag1.length + frag2.length);

            for (let i = 0; i < frag1.length; i++) {
                result[i] = frag1[i];
            }

            for (let i = 0; i < frag2.length; i++) {
                result[frag1.length + 1] = frag2[i];
            }

            return result;
        }

        // Adjust params.shape (account for rotation)
        const dir = normalize(params.vel);

        // Build rotation matrix around spin axis
        const R = rotationMat(params.spinAxis, params.rot_rad);
        const Rinv = rotationMat(params.spinAxis, -params.rot_rad);

        for (let i = 0; i < params.shape.vertices.length; i++) {
            let v = params.shape.vertices[i];

            // Transform vertex to world orientation
            v = rotVec(R, v);

            // Apply erosion bias in world frame
            const bias = 1.0 - 0.05 * dot(normalize(v), dir); // leading side erodes faster
            const scaleFactor = Math.pow(params.mass_kg / massPrev, 1.0 / 3.0) * bias;

            v = scale(v, scaleFactor);

            // Transform back to local params.shape frame
            params.shape.vertices[i] = rotVec(Rinv, v);
        }

        // runge-kutta 4th order integration
        
        // Local acceleration
        function accel(planet: PlanetConsts, pos: Vec3, p: Vec3, v: Vec3, m: f64): Vec3 {
            let alt = ecefToGeodetic({
                pos, 
                axis_m: planet.axis_m, 
                ecc: planet.ecc, 
                flat: planet.flat
            }).alt;
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

        let k1v = accel(params.pos, params.pos, params.vel, params.mass_kg);
        let k1r = params.vel;

        let k2v = accel(params.pos, add(params.pos, scale(k1r, dt * 0.5)), add(params.vel, scale(k1v, dt * 0.5)), params.mass_kg);
        let k2r = add(params.vel, scale(k1v, dt * 0.5));

        let k3v = accel(params.pos, add(params.pos, scale(k2r, dt * 0.5)), add(params.vel, scale(k2v, dt * 0.5)), params.mass_kg);
        let k3r = add(params.vel, scale(k2v, dt * 0.5));

        let k4v = accel(params.pos, add(params.pos, scale(k3r, dt)), add(params.vel, scale(k3v, dt)), params.mass_kg);
        let k4r = add(params.vel, scale(k3v, dt));

        params.vel = add(params.vel, scale(
            add(add(k1v, scale(add(k2v, k3v), 2.0)), k4v), dt / 6.0
        ));
        params.pos = add(params.pos, scale(
            add(add(k1r, scale(add(k2r, k3r), 2.0)), k4r), dt / 6.0
        ));

        params.vel = integ.vel;
        params.pos = integ.pos;

        // Rotation update
        params.rot_rad += params.rotr_radS * dt;
        if (params.rot_rad > 2.0*Math.PI) params.rot_rad -= 2.0*Math.PI;

        return [params];
    }

    // Integrate until hit ground or burn-up
    for (let step = 0; step < 200000; step++) {

        let keys = frags.keys();
        let state = new StaticArray<FragmentState>(keys.length);

        for (let i = 0; i < keys.length; i++) {
            let result = compute(frags.get(keys[i]), dt)
            frags.delete(keys[i]);
            result.forEach(frag => frags.set(frag.id, frag));
        }

        keys = frags.keys();

        for (let i = 0; i < keys.length; i++) {
            const frag = frags.get(keys[i])
            state[i] = {
                id: frag.id,
                origin: frag.pos,
                rot_rad: frag.rot_rad,
                spinAxis: frag.spinAxis,
                shape: frag.shape
            }
        }

        results.push({
            t_ms: (t * 1000) as i32,
            state
        })

        t += dt;
    }

    return results;
}
