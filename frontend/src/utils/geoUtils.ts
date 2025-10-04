import { Location, Country, MeteorImpact } from './types';

// Calculate distance between two points on Earth (in kilometers)
export function calculateDistance(lat1: number, lng1: number, lat2: number, lng2: number): number {
  const R = 6371; // Earth's radius in kilometers
  const dLat = toRadians(lat2 - lat1);
  const dLng = toRadians(lng2 - lng1);
  
  const a = Math.sin(dLat / 2) * Math.sin(dLat / 2) +
    Math.cos(toRadians(lat1)) * Math.cos(toRadians(lat2)) *
    Math.sin(dLng / 2) * Math.sin(dLng / 2);
  
  const c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
  return R * c;
}

function toRadians(degrees: number): number {
  return degrees * (Math.PI / 180);
}

// Calculate if a location is within impact radius
export function isWithinImpactRadius(
  location: Location, 
  impact: MeteorImpact
): boolean {
  const distance = calculateDistance(
    location.lat, 
    location.lng, 
    impact.location.lat, 
    impact.location.lng
  );
  return distance <= impact.radius;
}

// Calculate affected countries for a meteor impact
export function calculateAffectedCountries(
  impact: MeteorImpact, 
  countries: Country[]
): string[] {
  return countries
    .filter(country => isWithinImpactRadius(country.coordinates, impact))
    .map(country => country.name);
}

// Generate hazard ring coordinates for visualization
export function generateHazardRing(
  center: Location, 
  radius: number, 
  segments: number = 64
): Location[] {
  const ring: Location[] = [];
  const latRad = toRadians(center.lat);
  const lngRad = toRadians(center.lng);
  const radiusRad = radius / 6371; // Convert km to radians
  
  for (let i = 0; i <= segments; i++) {
    const angle = (i * 2 * Math.PI) / segments;
    const lat = Math.asin(
      Math.sin(latRad) * Math.cos(radiusRad) +
      Math.cos(latRad) * Math.sin(radiusRad) * Math.cos(angle)
    );
    const lng = lngRad + Math.atan2(
      Math.sin(angle) * Math.sin(radiusRad) * Math.cos(latRad),
      Math.cos(radiusRad) - Math.sin(latRad) * Math.sin(lat)
    );
    
    ring.push({
      lat: toDegrees(lat),
      lng: toDegrees(lng)
    });
  }
  
  return ring;
}

function toDegrees(radians: number): number {
  return radians * (180 / Math.PI);
}

// Sample country data
export const SAMPLE_COUNTRIES: Country[] = [
  {
    name: 'United States',
    capital: 'Washington, D.C.',
    population: 331002651,
    area: 9833517,
    continent: 'North America',
    coordinates: { lat: 39.8283, lng: -98.5795 },
    code: 'US'
  },
  {
    name: 'China',
    capital: 'Beijing',
    population: 1439323776,
    area: 9596961,
    continent: 'Asia',
    coordinates: { lat: 35.8617, lng: 104.1954 },
    code: 'CN'
  },
  {
    name: 'India',
    capital: 'New Delhi',
    population: 1380004385,
    area: 3287263,
    continent: 'Asia',
    coordinates: { lat: 20.5937, lng: 78.9629 },
    code: 'IN'
  },
  {
    name: 'Brazil',
    capital: 'Brasília',
    population: 212559417,
    area: 8514877,
    continent: 'South America',
    coordinates: { lat: -14.2350, lng: -51.9253 },
    code: 'BR'
  },
  {
    name: 'Russia',
    capital: 'Moscow',
    population: 145934462,
    area: 17098242,
    continent: 'Europe/Asia',
    coordinates: { lat: 61.5240, lng: 105.3188 },
    code: 'RU'
  },
  {
    name: 'Canada',
    capital: 'Ottawa',
    population: 37742154,
    area: 9984670,
    continent: 'North America',
    coordinates: { lat: 56.1304, lng: -106.3468 },
    code: 'CA'
  },
  {
    name: 'Australia',
    capital: 'Canberra',
    population: 25499884,
    area: 7692024,
    continent: 'Oceania',
    coordinates: { lat: -25.2744, lng: 133.7751 },
    code: 'AU'
  },
  {
    name: 'Germany',
    capital: 'Berlin',
    population: 83783942,
    area: 357114,
    continent: 'Europe',
    coordinates: { lat: 51.1657, lng: 10.4515 },
    code: 'DE'
  },
  {
    name: 'France',
    capital: 'Paris',
    population: 65273511,
    area: 551695,
    continent: 'Europe',
    coordinates: { lat: 46.2276, lng: 2.2137 },
    code: 'FR'
  },
  {
    name: 'Japan',
    capital: 'Tokyo',
    population: 126476461,
    area: 377975,
    continent: 'Asia',
    coordinates: { lat: 36.2048, lng: 138.2529 },
    code: 'JP'
  },
  {
    name: 'United Kingdom',
    capital: 'London',
    population: 67886011,
    area: 242495,
    continent: 'Europe',
    coordinates: { lat: 55.3781, lng: -3.4360 },
    code: 'GB'
  },
  {
    name: 'Italy',
    capital: 'Rome',
    population: 60461826,
    area: 301340,
    continent: 'Europe',
    coordinates: { lat: 41.8719, lng: 12.5674 },
    code: 'IT'
  },
  {
    name: 'South Africa',
    capital: 'Cape Town',
    population: 59308690,
    area: 1221037,
    continent: 'Africa',
    coordinates: { lat: -30.5595, lng: 22.9375 },
    code: 'ZA'
  },
  {
    name: 'Mexico',
    capital: 'Mexico City',
    population: 128932753,
    area: 1964375,
    continent: 'North America',
    coordinates: { lat: 23.6345, lng: -102.5528 },
    code: 'MX'
  },
  {
    name: 'Argentina',
    capital: 'Buenos Aires',
    population: 45195774,
    area: 2780400,
    continent: 'South America',
    coordinates: { lat: -38.4161, lng: -63.6167 },
    code: 'AR'
  }
];
