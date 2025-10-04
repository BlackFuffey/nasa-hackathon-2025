import React, { useState, useCallback } from 'react';
import SimpleGlobe from './components/SimpleGlobe';
import GlobeSidePanel from './components/GlobeSidePanel';
import CountryInfoModal from './components/CountryInfoModal';
import { Location, Country, MeteorImpact } from './utils/types';
import { SAMPLE_COUNTRIES } from './utils/geoUtils';
import { reverseGeocode } from './utils/revgeo';
import './App.css';

function App() {
  const [selectedLocation, setSelectedLocation] = useState<Location | undefined>();
  const [selectedCountry, setSelectedCountry] = useState<Country | undefined>();
  const [simulationMode, setSimulationMode] = useState(false);
  const [impacts, setImpacts] = useState<MeteorImpact[]>([]);
  const [showModal, setShowModal] = useState(false);
  const [centerOnLocation, setCenterOnLocation] = useState<Location | undefined>();

  const handleLocationClick = useCallback(async (location: Location) => {
    console.log('🎯 ===== LOCATION CLICKED IN APP =====');
    console.log('🎯 Received location:', location);
    console.log('🎯 Location type:', typeof location);
    console.log('🎯 Location lat:', location.lat, 'type:', typeof location.lat);
    console.log('🎯 Location lng:', location.lng, 'type:', typeof location.lng);
    
    // Validate coordinates
    if (typeof location.lat !== 'number' || typeof location.lng !== 'number' || 
        isNaN(location.lat) || isNaN(location.lng)) {
      console.error('❌ Invalid coordinates received:', location);
      return;
    }
    
    console.log('✅ Coordinates are valid, proceeding...');
    setSelectedLocation(location);
    
    // Try to get more detailed location info
    console.log('🔍 Attempting reverse geocoding...');
    const detailedLocation = await reverseGeocode(location.lat, location.lng);
    if (detailedLocation) {
      console.log('✅ Detailed location from API:', detailedLocation);
      // Merge API data with original coordinates to preserve clicked location
      const mergedLocation = {
        ...detailedLocation,
        lat: location.lat,  // Keep original clicked coordinates
        lng: location.lng   // Keep original clicked coordinates
      };
      setSelectedLocation(mergedLocation);
    } else {
      console.log('⚠️ No detailed location from API');
      // Keep original location if API fails
      setSelectedLocation(location);
    }

    // Find matching country using multiple methods
    console.log('🔍 Looking for matching country...');
    let country: Country | null = null;
    let closestCountry: Country | null = null;
    let closestDistance = Infinity;
    
    // Method 1: Exact name match from API
    if (detailedLocation?.country) {
      country = SAMPLE_COUNTRIES.find(c => 
        c.name.toLowerCase() === detailedLocation.country!.toLowerCase()
      ) || null;
      if (country) {
        console.log('✅ Method 1 - Found by API country name:', country.name);
      }
    }
    
    // Method 2: Find closest country by coordinates
    if (!country) {
      console.log('🔍 Method 2 - Finding closest country by coordinates...');
      
      SAMPLE_COUNTRIES.forEach((c: Country) => {
        const distance = Math.sqrt(
          Math.pow(c.coordinates.lat - location.lat, 2) + 
          Math.pow(c.coordinates.lng - location.lng, 2)
        );
        if (distance < closestDistance) {
          closestDistance = distance;
          closestCountry = c;
        }
      });
      
      if (closestCountry && closestDistance < 50) { // Within 50 degrees
        country = closestCountry;
        console.log('✅ Method 2 - Found closest country:', (closestCountry as Country).name, 'Distance:', closestDistance);
      }
    }
    
    // Method 3: Check if coordinates are in known regions
    if (!country) {
      console.log('🔍 Method 3 - Checking known regions...');
      const { lat, lng } = location;
      
      // North America
      if (lat >= 15 && lat <= 85 && lng >= -180 && lng <= -50) {
        country = SAMPLE_COUNTRIES.find(c => c.name === 'United States') || null;
        console.log('✅ Method 3 - North America region, using US as default');
      }
      // Europe
      else if (lat >= 35 && lat <= 70 && lng >= -10 && lng <= 40) {
        country = SAMPLE_COUNTRIES.find(c => c.name === 'Germany') || null;
        console.log('✅ Method 3 - Europe region, using Germany as default');
      }
      // Asia
      else if (lat >= 5 && lat <= 55 && lng >= 70 && lng <= 180) {
        country = SAMPLE_COUNTRIES.find(c => c.name === 'China') || null;
        console.log('✅ Method 3 - Asia region, using China as default');
      }
      // South America
      else if (lat >= -60 && lat <= 15 && lng >= -85 && lng <= -30) {
        country = SAMPLE_COUNTRIES.find(c => c.name === 'Brazil') || null;
        console.log('✅ Method 3 - South America region, using Brazil as default');
      }
    }
    
    if (country) {
      console.log('✅ Found matching country:', country.name);
      setSelectedCountry(country);
    } else {
      console.log('⚠️ No matching country found for coordinates:', location);
    }

    // Show modal
    console.log('📱 Showing modal...');
    setShowModal(true);
    console.log('🎯 ===== LOCATION CLICK HANDLER COMPLETE =====');
  }, []);

  const handleToggleSimulation = useCallback(() => {
    setSimulationMode(!simulationMode);
  }, [simulationMode]);

  const handleClearImpacts = useCallback(() => {
    setImpacts([]);
  }, []);

  const handleResetView = useCallback(() => {
    // This would reset the globe view - implementation depends on globe instance
    console.log('Reset view clicked');
  }, []);

  const handleCloseModal = useCallback(() => {
    setShowModal(false);
  }, []);

  const handleSearchCountry = useCallback((country: Country) => {
    console.log('Searching for country:', country);
    setSelectedCountry(country);
    setSelectedLocation(country.coordinates);
    setCenterOnLocation(country.coordinates);
  }, []);

  return (
    <div className="App">
      <header className="App-header">
        <h1>NASA Space Apps 2025 - Meteor Impact Simulator</h1>
        <p>Interactive 3D Globe for Meteor Impact Analysis</p>
      </header>
      
      <div className="App-content">
        <main className="App-main">
          <SimpleGlobe 
            width={900} 
            height={700} 
            onLocationClick={handleLocationClick}
            simulationMode={simulationMode}
            impacts={impacts}
            centerOnLocation={centerOnLocation}
          />
        </main>
        
        <GlobeSidePanel
          selectedLocation={selectedLocation}
          selectedCountry={selectedCountry}
          simulationMode={simulationMode}
          onToggleSimulation={handleToggleSimulation}
          onClearImpacts={handleClearImpacts}
          onResetView={handleResetView}
          onSearchCountry={handleSearchCountry}
          impacts={impacts}
        />
      </div>

      <CountryInfoModal
        isOpen={showModal}
        onClose={handleCloseModal}
        location={selectedLocation}
        country={selectedCountry}
      />
    </div>
  );
}

export default App;
