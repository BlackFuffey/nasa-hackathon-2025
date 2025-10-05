import React, { useState, useCallback } from "react";
import SimpleGlobe from "./components/SimpleGlobe";
import GlobeSidePanel from "./components/GlobeSidePanel";
import CountryInfoModal from "./components/CountryInfoModal";
import { Location, Country, MeteorImpact } from "./utils/types";
import { SAMPLE_COUNTRIES } from "./utils/geoUtils";
import { reverseGeocode } from "./utils/revgeo";
import "./App.css";
import { AsteroidPanel } from "./components/AsteroidPanel";
import { NeoAsteroid } from "./data/fetch";

function App() {
  const [selectedLocation, setSelectedLocation] = useState<Location | undefined>();
  const [selectedCountry, setSelectedCountry] = useState<Country | undefined>();
  const [simulationMode, setSimulationMode] = useState(false);
  const [impacts, setImpacts] = useState<MeteorImpact[]>([]);
  const [showModal, setShowModal] = useState(false);
  const [centerOnLocation, setCenterOnLocation] = useState<Location | undefined>();
  const [selectedAsteroid, setSelectedAsteroid] = useState<NeoAsteroid | null>(null);
  const [asteroidSuccessMessage, setAsteroidSuccessMessage] = useState<string | null>(null);

  const handleLocationClick = useCallback(async (location: Location) => {
    setSelectedLocation(location);
    const detailedLocation = await reverseGeocode(location.lat, location.lng);
    if (detailedLocation) {
      setSelectedLocation({ ...detailedLocation, lat: location.lat, lng: location.lng });
    }
    let country: Country | null = null;
    if (detailedLocation?.country) {
      country = SAMPLE_COUNTRIES.find(
        (c) => c.name.toLowerCase() === detailedLocation.country!.toLowerCase()
      ) || null;
    }
    if (!country) {
      let closestCountry: Country | null = null;
      let closestDistance = Infinity;
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
      if (closestCountry && closestDistance < 50) {
        country = closestCountry;
      }
    }
    if (!country) {
      const { lat, lng } = location;
      if (lat >= 15 && lat <= 85 && lng >= -180 && lng <= -50) {
        country = SAMPLE_COUNTRIES.find(c => c.name === "United States") || null;
      } else if (lat >= 35 && lat <= 70 && lng >= -10 && lng <= 40) {
        country = SAMPLE_COUNTRIES.find(c => c.name === "Germany") || null;
      } else if (lat >= 5 && lat <= 55 && lng >= 70 && lng <= 180) {
        country = SAMPLE_COUNTRIES.find(c => c.name === "China") || null;
      } else if (lat >= -60 && lat <= 15 && lng >= -85 && lng <= -30) {
        country = SAMPLE_COUNTRIES.find(c => c.name === "Brazil") || null;
      }
    }
    if (country) setSelectedCountry(country);
    setShowModal(true);
  }, []);

  const handleToggleSimulation = useCallback(() => setSimulationMode((prev) => !prev), []);
  const handleClearImpacts = useCallback(() => setImpacts([]), []);
  const handleResetView = useCallback(() => setCenterOnLocation(undefined), []);
  const handleCloseModal = useCallback(() => setShowModal(false), []);
  const handleSearchCountry = useCallback((country: Country) => {
    setSelectedCountry(country);
    setSelectedLocation(country.coordinates);
    setCenterOnLocation(country.coordinates);
  }, []);

  const asteroidSelectorEnabled = !!selectedCountry || !!selectedLocation;
  const preparedData =
    selectedAsteroid && selectedCountry && selectedLocation
      ? {
        asteroid: selectedAsteroid,
        impactCoordinates: { lat: selectedLocation.lat, lng: selectedLocation.lng },
        country: selectedCountry,
      }
      : null;

  const handleAsteroidSelect = (asteroid: NeoAsteroid | null) => {
    setSelectedAsteroid(asteroid);
    if (asteroid) {
      setAsteroidSuccessMessage(`Asteroid "${asteroid.name}" selected!`);
      setTimeout(() => setAsteroidSuccessMessage(null), 2400); // hide after 2.4s
    }
  };

  return (
    <div className="App">
      <div className="App-layout">
        <div className="App-sidepanel">
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
        <div className="App-dashboard">
          <main className="App-main">
            <SimpleGlobe
              width={800}
              height={600}
              onLocationClick={handleLocationClick}
              simulationMode={simulationMode}
              impacts={impacts}
              centerOnLocation={centerOnLocation}
            />
          </main>
          <div className="App-sidebar">
            <AsteroidPanel
              onSelect={handleAsteroidSelect}
              selected={selectedAsteroid}
              disabled={!asteroidSelectorEnabled}
            />
          </div>
        </div>
      </div>
      <CountryInfoModal
        isOpen={showModal}
        onClose={handleCloseModal}
        location={selectedLocation}
        country={selectedCountry}
      />
      {preparedData && (
        <button
          className="submit-impact-btn"
          onClick={() => alert("Submit: " + JSON.stringify(preparedData))}
          style={{
            position: "fixed",
            bottom: "20px",
            right: "40px",
            padding: "12px 20px",
            background: "#27c3fc",
            color: "#001a2c",
            borderRadius: "10px",
            fontWeight: "bold",
            border: "none",
            cursor: "pointer",
            boxShadow: "0 2px 14px #37f3ff22",
          }}
        >
          Calculate Impact!
        </button>
      )}

      {asteroidSuccessMessage && (
        <div style={{
          position: "fixed",
          bottom: "35px",
          left: "50%",
          transform: "translateX(-50%)",
          background: "#18ff9e",
          color: "#003b27",
          padding: "16px 32px",
          fontWeight: "bold",
          borderRadius: "10px",
          fontSize: "1.1rem",
          boxShadow: "0 4px 24px #0fffc150",
          zIndex: 99999,
          transition: "opacity 0.4s"
        }}>
          {asteroidSuccessMessage}
        </div>
      )}
    </div>
  );
}

export default App;
