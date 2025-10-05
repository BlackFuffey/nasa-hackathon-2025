// Types

@unmanaged
class Vec3 { 
    x!: f64;
    y!: f64;
    z!: f64;
}

@unmanaged
class i32Vec3 { 
    x!: f64;
    y!: f64;
    z!: f64;
}

@unmanaged
class Mat3x3 {
    a!: Vec3;
    b!: Vec3;
    c!: Vec3;
}

interface Mesh3D {
    vertices: StaticArray<Vec3>;
    faces: StaticArray<i32Vec3>;
}

interface FragmentState {
    id: i32;
    origin: Vec3;
    rot_rad: f64;
    spinAxis: Vec3;
    shape: Mesh3D;
}

interface SimulationFrame {
    t_ms: i32;
    state: FragmentState[]
}

interface AsteroidInitial {
    mass_kg: f64, shape: Mesh3D, strength_MPa: f64
    rotr_radS: f64, spinAxis: Vec3, 
}

interface Asteroid extends AsteroidInitial {
    id: i32,
    rot_rad: f64,
    pos: Vec3, vel: Vec3
}

interface PlanetConsts {
    mu: f64, radius: f64, axis_m: f64,
    ecc: f64, rotr_radS: f64,
    eqtgrav_mS2: f64, flat: f64
}

interface KeplerOrbit {
    axis_m: f64,
    ecc: f64,
    incl_rad: f64,
    peri_rad: f64,
    lon_rad: f64,
    mean_rad: f64
}

interface PlanetCordinates {
    lat: f64, lon: f64, alt: f64
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
interface KeplerToCartesian {
    mu: f64,
    orbital: KeplerOrbit
}
function keplerToCartesian(params: KeplerToCartesian): { pos: Vec3, vel: Vec3 } {
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
interface EcefToGeodeticParams {
    pos: Vec3, axis_m: f64, ecc: f64, flat: f64
}
function ecefToGeodetic(params: EcefToGeodeticParams): PlanetCordinates {

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

interface GravityAtParams {
    pos: Vec3, eqtgrav_mS2: f64, axis_m: f64, ecc: f64, flat: f64
}
function gravityAt(params: GravityAtParams): f64 {
    const { lat } = ecefToGeodetic({
        pos: params.pos, axis_m: params.axis_m,
        ecc: params.ecc, flat: params.flat
    });

    return params.eqtgrav_mS2 * 
        (1 + 0.0053024*Math.sin(lat)**2 - 0.0000058*Math.sin(2*lat)**2
        );
}

// Main simulation
interface MeteorsimParams {
    planet: PlanetConsts,
    orbital: KeplerOrbit,
    asteroid: AsteroidInitial,
    resol_ms: i32
}
export function meteorsim(params: MeteorsimParams): Array<SimulationFrame> {

    const dt = params.resol_ms as f64 / 1000.0; // timestep in seconds
    const planetR = params.planet.radius;

    let { pos, vel } = keplerToCartesian({ mu: params.planet.mu, orbital: params.orbital });

    let fragId = 0;
    function nextId() { return fragId++; }

    const frags = new Map<i32, Asteroid>();
    frags.set(0, {
        rotr_radS: params.asteroid.rotr_radS,
        spinAxis: params.asteroid.spinAxis,
        rot_rad: 0.0,

        mass_kg: params.asteroid.mass_kg,
        shape: params.asteroid.shape,
        strength_MPa: params.asteroid.strength_MPa,

        pos, vel,

        id: 0
    })

    let t: f64 = 0.0;

    let results = new Array<SimulationFrame>;

    // runge-kutta 4th order integration
    function rk4Step(pos: Vec3, vel: Vec3, m: f64, dt: f64): { pos: Vec3, vel: Vec3 } {

        // Local acceleration
        function accel(p: Vec3, v: Vec3, m: f64): Vec3 {
            let rmag = norm(p);
            let { lat, lon, alt } = ecefToGeodetic({
                pos, 
                axis_m: params.planet.axis_m, 
                ecc: params.planet.ecc, 
                flat: params.planet.flat
            });
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

    function compute({ id, rot_rad, mass_kg, pos, vel, shape, spinAxis, rotr_radS }: Asteroid, dt: f64): Asteroid[] {
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

        // Fragmentation 
        // O(sobbing terrible)
        if (0.5 * rho * vmag * vmag > params.asteroid.strength_MPa * 1e6) {

            // --- 1. Compute aerodynamic stress direction (shock front) ---
            let stressDir = vec3(0.0, 0.0, 0.0);
            let maxStress = 0.0;
            const vdir = normalize(vel);

            for (let i = 0; i < shape.faces.length; i++) {
                const f = shape.faces[i];
                const a = shape.vertices[f.x];
                const b = shape.vertices[f.y];
                const c = shape.vertices[f.z];

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
            for (let i = 0; i < shape.vertices.length; i++) {
                const v = shape.vertices[i];
                cx += v.x; cy += v.y; cz += v.z;
            }
            const center = vec3(cx / shape.vertices.length, cy / shape.vertices.length, cz / shape.vertices.length);
            const nfract = normalize(stressDir);

            // --- 3. Split mesh into two fragments along the fracture plane ---
            const vA: Vec3[] = [], vB: Vec3[] = [];
            const fA: i32Vec3[] = [], fB: i32Vec3[] = [];

            for (let i = 0; i < shape.vertices.length; i++) {
                const v = shape.vertices[i];
                const side = dot(sub(v, center), nfract);
                if (side >= 0.0) vA.push(v);
                    else vB.push(v);
            }

            for (let i = 0; i < shape.faces.length; i++) {
                const f = shape.faces[i];
                let countA = 0;
                if (vA.includes(shape.vertices[f.x])) countA++;
                if (vA.includes(shape.vertices[f.y])) countA++;
                if (vA.includes(shape.vertices[f.z])) countA++;
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
            const rotMatrix = rotationMat(spinAxis, rot_rad);

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

            const newSpinAxis = normalize(spinAxis);
            const newRotRadS = rotr_radS;

            // --- 5. Momentum-conserving separation ---
            const mass1 = mass_kg * 0.5;
            const mass2 = mass_kg * 0.5;
            const sepVel = scale(nfract, vmag * 0.02);

            // --- 6. Fill fracture surface ---
            const eps = 1e-6;
            const fractureVerts: Vec3[] = [];
            for (let i = 0; i < shape.vertices.length; i++) {
                const v = shape.vertices[i];
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

            return [
                ...compute({
                    id: nextId(),
                    rot_rad: 0.0,
                    mass_kg: mass1,
                    pos,
                    vel: add(vel, sepVel),
                    shape: mesh1,
                    spinAxis: newSpinAxis,
                    rotr_radS: newRotRadS,
                }, dt),
                ...compute({
                    id: nextId(),
                    rot_rad: 0.0,
                    mass_kg: mass2,
                    pos,
                    vel: sub(vel, sepVel),
                    shape: mesh2,
                    spinAxis: newSpinAxis,
                    rotr_radS: newRotRadS,
                }, dt)
            ];
        }

        // Adjust shape (account for rotation)
        const dir = normalize(vel);

        // Build rotation matrix around spin axis
        const R = rotationMat(spinAxis, rot_rad);
        const Rinv = rotationMat(spinAxis, -rot_rad);

        for (let i = 0; i < shape.vertices.length; i++) {
            let v = shape.vertices[i];

            // Transform vertex to world orientation
            v = rotVec(R, v);

            // Apply erosion bias in world frame
            const bias = 1.0 - 0.05 * dot(normalize(v), dir); // leading side erodes faster
            const scaleFactor = Math.pow(mass_kg / massPrev, 1.0 / 3.0) * bias;

            v = scale(v, scaleFactor);

            // Transform back to local shape frame
            shape.vertices[i] = rotVec(Rinv, v);
        }

        // Integrate
        const integ = rk4Step(pos, vel, mass_kg, dt);
        vel = integ.vel;
        pos = integ.pos;

        // Rotation update
        rot_rad += rotr_radS * dt;
        if (rot_rad > 2.0*Math.PI) rot_rad -= 2.0*Math.PI;

        return [{ id, rot_rad, mass_kg, pos, vel, shape, spinAxis, rotr_radS, strength_MPa }];
    }

    // Integrate until hit ground or burn-up
    for (let step = 0; step < 200000; step++) {

        frags.keys().forEach(key => {
            let result = compute(frags.get(key), dt)
            frags.delete(key);
            result.forEach(frag => frags.set(frag.id, frag));
        })

        results.push({
            t_ms: (t * 1000) as i32,
            state: frags.keys().map(key => {
                const frag = frags.get(key);
                return {
                    id: frag.id,
                    origin: frag.pos,
                    rot_rad: frag.rot_rad,
                    spinAxis: frag.spinAxis,
                    shape: frag.shape
                }
            })
        })

        t += dt;
    }

    return results;
}
