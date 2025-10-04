import { OpenCageResponse, Location } from './types';

const OPENCAGE_API_KEY = process.env.REACT_APP_OPENAI_API_KEY;
const OPENCAGE_BASE_URL = 'https://api.opencagedata.com/geocode/v1/json';

// Cache for API responses to avoid rate limiting
const locationCache = new Map<string, Location>();

export async function reverseGeocode(lat: number, lng: number): Promise<Location | null> {
  // Check cache first
  const cacheKey = `${lat.toFixed(4)},${lng.toFixed(4)}`;
  if (locationCache.has(cacheKey)) {
    return locationCache.get(cacheKey)!;
  }

  if (!OPENCAGE_API_KEY) {
    console.warn('OpenCage API key not found. Please set REACT_APP_OPENAI_API_KEY in your .env file');
    return null;
  }

  try {
    const response = await fetch(
      `${OPENCAGE_BASE_URL}?q=${lat}+${lng}&key=${OPENCAGE_API_KEY}&no_annotations=1&limit=1`
    );

    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }

    const data: OpenCageResponse = await response.json();

    if (data.results && data.results.length > 0) {
      const result = data.results[0];
      const location: Location = {
        lat: result.geometry.lat,
        lng: result.geometry.lng,
        name: result.formatted,
        country: result.components.country,
        continent: result.components.continent
      };

      // Cache the result
      locationCache.set(cacheKey, location);
      return location;
    }

    return null;
  } catch (error) {
    console.error('Error fetching location data:', error);
    return null;
  }
}

// Get location name from coordinates
export async function getLocationName(lat: number, lng: number): Promise<string> {
  const location = await reverseGeocode(lat, lng);
  return location?.name || `Location: ${lat.toFixed(4)}, ${lng.toFixed(4)}`;
}

// Clear the location cache (useful for testing or memory management)
export function clearLocationCache(): void {
  locationCache.clear();
}
