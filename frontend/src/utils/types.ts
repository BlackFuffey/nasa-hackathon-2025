// Geographic and location types
export interface Location {
  lat: number;
  lng: number;
  name?: string;
  country?: string;
  continent?: string;
}

export interface Country {
  name: string;
  capital: string;
  population: number;
  area: number;
  continent: string;
  coordinates: Location;
  code: string;
}

export interface MeteorImpact {
  id: string;
  location: Location;
  radius: number;
  intensity: number;
  timestamp: Date;
  affectedCountries: string[];
  estimatedDamage: number;
}

export interface GlobeClickEvent {
  lat: number;
  lng: number;
  point?: {
    lat: number;
    lng: number;
    size: number;
    color: string;
  };
}

// API response types
export interface OpenCageResponse {
  results: OpenCageResult[];
  status: {
    code: number;
    message: string;
  };
}

export interface OpenCageResult {
  components: {
    country?: string;
    country_code?: string;
    continent?: string;
    city?: string;
    town?: string;
    village?: string;
    state?: string;
    county?: string;
    postcode?: string;
    road?: string;
    house_number?: string;
  };
  formatted: string;
  geometry: {
    lat: number;
    lng: number;
  };
}

// Component prop types
export interface SimpleGlobeProps {
  width?: number;
  height?: number;
  onLocationClick?: (location: Location) => void;
  simulationMode?: boolean;
  impacts?: MeteorImpact[];
  centerOnLocation?: Location;
  simulationResult?: any[];
}

export interface GlobeSidePanelProps {
  selectedLocation?: Location;
  selectedCountry?: Country;
  simulationMode: boolean;
  onToggleSimulation: () => void;
  onClearImpacts: () => void;
  onResetView: () => void;
  onSearchCountry: (country: Country) => void;
  impacts: MeteorImpact[];
}

export interface CountryInfoModalProps {
  isOpen: boolean;
  onClose: () => void;
  location?: Location;
  country?: Country;
}
