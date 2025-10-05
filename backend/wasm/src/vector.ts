// Vector data structures and helpers

export class Vec3 { 
    x!: f64;
    y!: f64;
    z!: f64;
}

export class Mat3x3 {
    a!: Vec3;
    b!: Vec3;
    c!: Vec3;
}

// Vector Math Helpers
export function vec3(x: f64, y: f64, z: f64): Vec3 { return { x, y, z }; }

export function add(a: Vec3, b: Vec3): Vec3 { return { x: a.x + b.x, y: a.y + b.y, z: a.z + b.z }; }

export function sub(a: Vec3, b: Vec3): Vec3 { return { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z }; }

export function scale(a: Vec3, s: f64): Vec3 { return { x: a.x * s, y: a.y * s, z: a.z * s }; }

export function dot(a: Vec3, b: Vec3): f64 { return a.x*b.x + a.y*b.y + a.z*b.z; }

export function norm(a: Vec3): f64 { return Math.sqrt(dot(a,a)); }

export function normalize(a: Vec3): Vec3 { let n = norm(a); return scale(a, 1.0/n); }

// Matrix math helpers
export function rotationMat(axis: Vec3, ang_rad: f64): Mat3x3 {
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

export function rotVec(M: Mat3x3, v: Vec3): Vec3 {
    return vec3(
        M.a.x*v.x + M.a.y*v.y + M.a.z*v.z,
        M.b.x*v.x + M.b.y*v.y + M.b.z*v.z,
        M.c.x*v.x + M.c.y*v.y + M.c.z*v.z
    );
}
