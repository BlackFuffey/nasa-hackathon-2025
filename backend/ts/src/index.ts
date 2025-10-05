// backend/src/index.ts
import express from 'express';
import * as fs from 'fs';
import * as path from 'path';
import { instantiate } from '@assemblyscript/loader';


const app = express();
app.use(express.json());

let meteorsim: Function | null = null;


async function loadWasm(): Promise<Function> {
  if (meteorsim) return meteorsim;

  // Use .tmp/module.wasm in dev, build/module.wasm in prod
  const wasmPath = fs.existsSync(path.join(__dirname, '..', '.tmp', 'module.wasm'))
    ? path.join(__dirname, '..', '.tmp', 'module.wasm')
    : path.join(__dirname, '..', 'build', 'module.wasm');

  const wasmBytes = fs.readFileSync(wasmPath);
  const { exports } = await instantiate(wasmBytes);
  meteorsim = exports.meteorsim as Function;
  return meteorsim;
}


const EARTH = {
  mu: 3.986004418e14,
  radius: 6.371e6,
  axis_m: 6.378137e6,
  ecc: 0.08181919,
  rotr_radS: 7.292115e-5,
  eqtgrav_mS2: 9.7803267714,
  flat: 1 / 298.257223563,
};


function latLngToECEF(lat: number, lng: number): { x: number; y: number; z: number } {
  const φ = (lat * Math.PI) / 180;
  const λ = (lng * Math.PI) / 180;
  const a = 6378137;
  const e2 = 6.69437999014e-3;
  const N = a / Math.sqrt(1 - e2 * Math.sin(φ) * Math.sin(φ));
  return {
    x: (N) * Math.cos(φ) * Math.cos(λ),
    y: (N) * Math.cos(φ) * Math.sin(λ),
    z: (N * (1 - e2)) * Math.sin(φ),
  };
}


function generateSphere(radius: number, segments = 8): { vertices: { x: number; y: number; z: number }[]; faces: { 0: number; 1: number; 2: number }[] } {
  const vertices: { x: number; y: number; z: number }[] = [];
  const faces: { 0: number; 1: number; 2: number }[] = [];

  for (let i = 0; i <= segments; i++) {
    const phi = Math.PI * (i / segments);
    for (let j = 0; j <= segments; j++) {
      const theta = 2 * Math.PI * (j / segments);
      vertices.push({
        x: radius * Math.sin(phi) * Math.cos(theta),
        y: radius * Math.sin(phi) * Math.sin(theta),
        z: radius * Math.cos(phi),
      });
    }
  }

  for (let i = 0; i < segments; i++) {
    for (let j = 0; j < segments; j++) {
      const a = i * (segments + 1) + j;
      const b = a + 1;
      const c = a + (segments + 1);
      const d = c + 1;
      faces.push({ 0: a, 1: b, 2: c });
      faces.push({ 0: b, 1: d, 2: c });
    }
  }

  return { vertices, faces };
}


app.post('/simulate', async (req, res) => {
  const { lat, lng, diameter_m, velocity_kms, strength_MPa = 0.1 } = req.body;

  if (!lat || !lng || !diameter_m || !velocity_kms) {
    return res.status(400).json({
      error: 'Missing required parameters: lat, lng, diameter_m, velocity_kms',
    });
  }

  try {
    const pos = latLngToECEF(lat, lng);
    const speed = velocity_kms * 1000;
    const vel = { x: -pos.x, y: -pos.y, z: -pos.z };
    const norm = Math.hypot(vel.x, vel.y, vel.z);
    vel.x = (vel.x / norm) * speed;
    vel.y = (vel.y / norm) * speed;
    vel.z = (vel.z / norm) * speed;

    const mass = (4 / 3) * Math.PI * Math.pow(diameter_m / 2, 3) * 3000;
    const { vertices, faces } = generateSphere(diameter_m / 2, 6);

    const params = {
      planet: EARTH,
      orbital: {
        axis_m: 7e6,
        ecc: 0.9,
        incl_rad: 0.5,
        peri_rad: 0.3,
        lon_rad: 0.2,
        mean_rad: 0.1,
      },
      asteroid: {
        mass_kg: mass,
        rotr_radS: 0.1,
        spinAxis: { x: 0, y: 1, z: 0 },
        shape: { vertices, faces },
        strength_MPa,
      },
      resol_ms: 100,
    };

    const runSim = await loadWasm();
    const result = runSim(params);

    res.json(result);
  } catch (err: any) {
    console.error('❌ Simulation failed:', err.message);
    res.status(500).json({ error: 'Simulation failed', details: err.message });
  }
});

const PORT = process.env.PORT || 5000;
app.listen(PORT, () => {
  console.log(`🚀 Backend server running on http://localhost:${PORT}`);
});
