import { 
    Vec3,
    vec3, add, sub, scale, dot, norm, normalize,
    rotationMat, rotVec, cross
} from './vector'

import {
    ecefToGeodetic, gravityAt
} from './geography'

import {
    KeplerOrbit, keplerToCartesian
} from './mechanics'

// Types
class Mesh3D {
    vertices!: StaticArray<Vec3>;
    faces!: StaticArray<i32Vec3>;
}

class i32Vec3 { 
    x!: i32;
    y!: i32;
    z!: i32;
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


// Placeholder exponential atmosphere (replace with NRLMSISE-00)
function airDensity(alt_km: f64): f64 {
    if (alt_km < 0.0) return 1.225; // sea level
    return 1.225 * Math.exp(-alt_km / 7.5); // scale height ~7.5 km
}

// Drag coefficient (replace with simulated drag)
const CD: f64 = 1.3;

// Main simulation
let fragId = 0;
function nextId(): i32 { return fragId++; }

function compute(params: Asteroid, planet: PlanetConsts, dt: f64): StaticArray<Asteroid> {
    // Rotation matrices
    const R_bodyToWorld = rotationMat(params.spinAxis, params.rot_rad);
    const R_worldToBody = rotationMat(params.spinAxis, -params.rot_rad);

    // ---- 1. Basic kinematics ----
    const rmag = norm(params.pos);
    const alt = (rmag - planet.radius) / 1000.0; // km altitude

    if (alt < 0.0) return []; // ground impact
    if (params.mass_kg < 1.0) return []; // burned up

    // ---- 2. Atmospheric and heating ----
    const rho = airDensity(alt);
    const vmag = norm(params.vel);
    const vdir_world = normalize(params.vel);
    const vdir_body = rotVec(R_worldToBody, vdir_world); // body-frame airflow dir

    const massPrev = params.mass_kg;
    const heatFlux = 0.5 * rho * vmag * vmag * vmag;
    const dm = -heatFlux * 1e-8 * dt; // tuning coefficient
    params.mass_kg = Math.max(0.0, params.mass_kg + dm);

    // ---- 3. Fragmentation ----
    const dynPressure = 0.5 * rho * vmag * vmag;
    if (dynPressure > params.strength_MPa * 1e6) {
        // --- 1. Compute aerodynamic stress direction in body frame ---
        let stressDir_body = vec3(0.0, 0.0, 0.0);
        let maxStress = 0.0;

        for (let i = 0; i < params.shape.faces.length; i++) {
            const f = params.shape.faces[i];
            const a = params.shape.vertices[f.x];
            const b = params.shape.vertices[f.y];
            const c = params.shape.vertices[f.z];

            const ab = sub(b, a);
            const ac = sub(c, a);
            const n = normalize(cross(ab, ac)); // local face normal

            const stress = dynPressure * Math.max(0.0, dot(n, vdir_body));
            if (stress > maxStress) {
                maxStress = stress;
                stressDir_body = n;
            }
        }

        // --- 2. Define fracture plane through centroid ---
        let cx = 0.0, cy = 0.0, cz = 0.0;
        for (let i = 0; i < params.shape.vertices.length; i++) {
            const v = params.shape.vertices[i];
            cx += v.x; cy += v.y; cz += v.z;
        }
        const center = vec3(cx / params.shape.vertices.length, cy / params.shape.vertices.length, cz / params.shape.vertices.length);
        const nfract_body = normalize(stressDir_body);

        // --- 3. Split mesh into two sets along plane ---
        const vA: Vec3[] = [], vB: Vec3[] = [];
        const fA: i32Vec3[] = [], fB: i32Vec3[] = [];

        for (let i = 0; i < params.shape.vertices.length; i++) {
            const v = params.shape.vertices[i];
            const side = dot(sub(v, center), nfract_body);
            if (side >= 0.0) vA.push(v); else vB.push(v);
        }

        for (let i = 0; i < params.shape.faces.length; i++) {
            const f = params.shape.faces[i];
            let countA = 0;
            const vx = params.shape.vertices[f.x];
            const vy = params.shape.vertices[f.y];
            const vz = params.shape.vertices[f.z];

            for (let i = 0; i < vA.length; i++) {
                if (vA[i].x === vx.x && vA[i].y === vx.y && vA[i].z === vx.z) countA++;
                if (vA[i].x === vy.x && vA[i].y === vy.y && vA[i].z === vy.z) countA++;
                if (vA[i].x === vz.x && vA[i].y === vz.y && vA[i].z === vz.z) countA++;
            }
            if (countA >= 2) fA.push(f); else fB.push(f);
        }

        const mesh1: Mesh3D = { vertices: StaticArray.fromArray(vA), faces: StaticArray.fromArray(fA) };
        const mesh2: Mesh3D = { vertices: StaticArray.fromArray(vB), faces: StaticArray.fromArray(fB) };

        // --- 4. Momentum-conserving separation (in world frame) ---
        const nfract_world = rotVec(R_bodyToWorld, nfract_body);
        const sepVel = scale(nfract_world, vmag * 0.02);

        const mass1 = params.mass_kg * 0.5;
        const mass2 = params.mass_kg * 0.5;

        const frag1: StaticArray<Asteroid> = compute({
            id: nextId(),
            rot_rad: 0.0,
            mass_kg: mass1,
            pos: params.pos,
            vel: add(params.vel, sepVel),
            shape: mesh1,
            spinAxis: params.spinAxis,
            rotr_radS: params.rotr_radS,
            strength_MPa: params.strength_MPa
        }, planet, dt);

        const frag2: StaticArray<Asteroid> = compute({
            id: nextId(),
            rot_rad: 0.0,
            mass_kg: mass2,
            pos: params.pos,
            vel: sub(params.vel, sepVel),
            shape: mesh2,
            spinAxis: params.spinAxis,
            rotr_radS: params.rotr_radS,
            strength_MPa: params.strength_MPa
        }, planet, dt);

        const result = new StaticArray<Asteroid>(frag1.length + frag2.length);
        for (let i = 0; i < frag1.length; i++) result[i] = frag1[i];
        for (let i = 0; i < frag2.length; i++) result[frag1.length + i] = frag2[i];
        return result;
    }

    // ---- 4. Erosion and rotation-aware scaling ----
    const vel_body = rotVec(R_worldToBody, normalize(params.vel)); // world->body
    for (let i = 0; i < params.shape.vertices.length; i++) {
        let v = params.shape.vertices[i];
        const bias = 1.0 - 0.05 * dot(normalize(v), vel_body); // leading side erodes faster
        const scaleFactor = Math.pow(params.mass_kg / massPrev, 1.0 / 3.0) * bias;
        params.shape.vertices[i] = scale(v, scaleFactor);
    }

    // ---- 5. Dynamics integration (RK4) ----
    function accel(planet: PlanetConsts, pos: Vec3, v: Vec3, m: f64): Vec3 {
        const geo = ecefToGeodetic({ pos, axis_m: planet.axis_m, ecc: planet.ecc, flat: planet.flat });
        const rho = airDensity(geo.alt);
        const vmag = norm(v);
        const dragMag = 0.5 * CD * rho * vmag * vmag * Math.pow(m / 2000.0, 2 / 3);
        const drag = scale(normalize(v), -dragMag / m);
        const gravity = scale(normalize(pos), -gravityAt({
            pos, eqtgrav_mS2: planet.eqtgrav_mS2,
            axis_m: planet.axis_m, ecc: planet.ecc, flat: planet.flat
        }));
        return add(gravity, drag);
    }

    const k1v = accel(planet, params.pos, params.vel, params.mass_kg);
    const k1r = params.vel;

    const k2v = accel(planet, add(params.pos, scale(k1r, dt * 0.5)), add(params.vel, scale(k1v, dt * 0.5)), params.mass_kg);
    const k2r = add(params.vel, scale(k1v, dt * 0.5));

    const k3v = accel(planet, add(params.pos, scale(k2r, dt * 0.5)), add(params.vel, scale(k2v, dt * 0.5)), params.mass_kg);
    const k3r = add(params.vel, scale(k2v, dt * 0.5));

    const k4v = accel(planet, add(params.pos, scale(k3r, dt)), add(params.vel, scale(k3v, dt)), params.mass_kg);
    const k4r = add(params.vel, scale(k3v, dt));

    params.vel = add(params.vel, scale(add(add(k1v, scale(add(k2v, k3v), 2.0)), k4v), dt / 6.0));
    params.pos = add(params.pos, scale(add(add(k1r, scale(add(k2r, k3r), 2.0)), k4r), dt / 6.0));

    // ---- 6. Rotation update ----
    params.rot_rad += params.rotr_radS * dt;
    if (params.rot_rad > 2.0 * Math.PI) params.rot_rad -= 2.0 * Math.PI;

    return [params];
}
class MeteorsimParams {
    planet!: PlanetConsts;
    orbital!: KeplerOrbit;
    asteroid!: AsteroidInitial;
    resol_ms!: i32;
}

export function meteorsim(params: MeteorsimParams): Array<SimulationFrame> {

    const dt = params.resol_ms as f64 / 1000.0; // timestep in seconds

    let movement = keplerToCartesian({ mu: params.planet.mu, orbital: params.orbital });


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


    // Integrate until hit ground or burn-up
    for (let step = 0; step < 200000; step++) {

        let keys = frags.keys();
        let state = new StaticArray<FragmentState>(keys.length);

        for (let i = 0; i < keys.length; i++) {
            let result = compute(frags.get(keys[i]), params.planet, dt)
            frags.delete(keys[i]);

            for (let j = 0; j < result.length; j++) {
                frags.set(result[j].id, result[j]);
            }
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
