// fetchNEOData.ts
import axios from "axios";

const apiKey = "qsce741pO9VyjytLVO3dSO9TudjRQfkPloSjAqSh";
const startDate = "2025-09-26";
const endDate = "2025-10-03";
const baseURL = `https://api.nasa.gov/neo/rest/v1/feed?start_date=${startDate}&end_date=${endDate}&api_key=${apiKey}`;

// Interfaces for NASA NEO API response
export interface Diameter {
  estimated_diameter_min: number;
  estimated_diameter_max: number;
}

export interface EstimatedDiameter {
  kilometers: Diameter;
  meters: Diameter;
  miles: Diameter;
  feet: Diameter;
}

export interface RelativeVelocity {
  kilometers_per_second: string;
  kilometers_per_hour: string;
  miles_per_hour: string;
}

export interface MissDistance {
  astronomical: string;
  lunar: string;
  kilometers: string;
  miles: string;
}

export interface CloseApproachData {
  close_approach_date: string;
  close_approach_date_full: string;
  epoch_date_close_approach: number;
  relative_velocity: RelativeVelocity;
  miss_distance: MissDistance;
  orbiting_body: string;
}

export interface NeoAsteroid {
  links: {
    self: string;
  };
  id: string;
  neo_reference_id: string;
  name: string;
  nasa_jpl_url: string;
  absolute_magnitude_h: number;
  estimated_diameter: EstimatedDiameter;
  is_potentially_hazardous_asteroid: boolean;
  close_approach_data: CloseApproachData[];
  is_sentry_object: boolean;
}

export interface NeoFeedResponse {
  links: {
    next: string;
    prev: string;
    self: string;
  };
  element_count: number;
  near_earth_objects: {
    [date: string]: NeoAsteroid[];
  };
}

export async function fetchNEOData(): Promise<NeoAsteroid[]> {
  try {
    const response = await axios.get<NeoFeedResponse>(baseURL);
    const neoMap = response.data.near_earth_objects;
    const neos = Object.values(neoMap).flat() as NeoAsteroid[];
    return neos;
  } catch (error) {
    console.error("Error fetching NEO data:", error);
    return [];
  }
}
