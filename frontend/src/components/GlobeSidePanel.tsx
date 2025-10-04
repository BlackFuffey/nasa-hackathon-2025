import React, { useState, useMemo } from 'react';
import { GlobeSidePanelProps, Country } from '../utils/types';
import { SAMPLE_COUNTRIES } from '../utils/geoUtils';
import './GlobeSidePanel.css';

const GlobeSidePanel: React.FC<GlobeSidePanelProps> = ({
  selectedLocation,
  selectedCountry,
  simulationMode,
  onToggleSimulation,
  onClearImpacts,
  onResetView,
  onSearchCountry,
  impacts
}) => {
  const [searchQuery, setSearchQuery] = useState('');
  const [showSearchResults, setShowSearchResults] = useState(false);

  // Filter countries based on search query
  const filteredCountries = useMemo(() => {
    if (!searchQuery.trim()) return [];
    
    return SAMPLE_COUNTRIES.filter(country =>
      country.name.toLowerCase().includes(searchQuery.toLowerCase()) ||
      country.capital.toLowerCase().includes(searchQuery.toLowerCase()) ||
      country.continent.toLowerCase().includes(searchQuery.toLowerCase())
    );
  }, [searchQuery]);

  const handleSearchChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const query = e.target.value;
    setSearchQuery(query);
    setShowSearchResults(query.length > 0);
  };

  const handleCountrySelect = (country: Country) => {
    onSearchCountry(country);
    setSearchQuery('');
    setShowSearchResults(false);
  };

  const handleSearchSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    if (filteredCountries.length > 0) {
      handleCountrySelect(filteredCountries[0]);
    }
  };
  return (
    <div className="side-panel">
      <div className="panel-header">
        <h2>Meteor Impact Simulator</h2>
        <p>Interactive 3D Globe Analysis</p>
      </div>

      <div className="panel-section">
        <h3>Search Countries</h3>
        <form onSubmit={handleSearchSubmit} className="search-form">
          <div className="search-container">
            <input
              type="text"
              placeholder="Search countries, capitals, or continents..."
              value={searchQuery}
              onChange={handleSearchChange}
              className="search-input"
            />
            <button type="submit" className="search-button">
              🔍
            </button>
          </div>
          
          {showSearchResults && (
            <div className="search-results">
              {filteredCountries.length > 0 ? (
                filteredCountries.slice(0, 5).map((country) => (
                  <div
                    key={country.code}
                    className="search-result-item"
                    onClick={() => handleCountrySelect(country)}
                  >
                    <div className="country-name">{country.name}</div>
                    <div className="country-details">
                      {country.capital} • {country.continent}
                    </div>
                  </div>
                ))
              ) : (
                <div className="no-results">No countries found</div>
              )}
            </div>
          )}
        </form>
      </div>

      <div className="panel-section">
        <h3>Simulation Controls</h3>
        <div className="control-buttons">
          <button 
            className={`simulation-toggle ${simulationMode ? 'active' : ''}`}
            onClick={onToggleSimulation}
          >
            {simulationMode ? 'Exit Simulation' : 'Start Simulation'}
          </button>
          
          {simulationMode && (
            <>
              <button 
                className="clear-button"
                onClick={onClearImpacts}
                disabled={impacts.length === 0}
              >
                Clear Impacts ({impacts.length})
              </button>
              
              <button 
                className="reset-button"
                onClick={onResetView}
              >
                Reset View
              </button>
            </>
          )}
        </div>
      </div>

      {selectedLocation && (
        <div className="panel-section">
          <h3>Selected Location</h3>
          <div className="location-info">
            <div className="info-item">
              <span className="label">Coordinates:</span>
              <span className="value">
                {selectedLocation.lat.toFixed(4)}, {selectedLocation.lng.toFixed(4)}
              </span>
            </div>
            {selectedLocation.name && (
              <div className="info-item">
                <span className="label">Name:</span>
                <span className="value">{selectedLocation.name}</span>
              </div>
            )}
            {selectedLocation.country && (
              <div className="info-item">
                <span className="label">Country:</span>
                <span className="value">{selectedLocation.country}</span>
              </div>
            )}
            {selectedLocation.continent && (
              <div className="info-item">
                <span className="label">Continent:</span>
                <span className="value">{selectedLocation.continent}</span>
              </div>
            )}
          </div>
        </div>
      )}

      {selectedCountry && (
        <div className="panel-section">
          <h3>Country Information</h3>
          <div className="country-info">
            <div className="country-header">
              <h4>{selectedCountry.name}</h4>
              <span className="country-code">{selectedCountry.code}</span>
            </div>
            <div className="info-grid">
              <div className="info-item">
                <span className="label">Capital:</span>
                <span className="value">{selectedCountry.capital}</span>
              </div>
              <div className="info-item">
                <span className="label">Population:</span>
                <span className="value">
                  {selectedCountry.population.toLocaleString()}
                </span>
              </div>
              <div className="info-item">
                <span className="label">Area:</span>
                <span className="value">
                  {selectedCountry.area.toLocaleString()} km²
                </span>
              </div>
              <div className="info-item">
                <span className="label">Continent:</span>
                <span className="value">{selectedCountry.continent}</span>
              </div>
            </div>
          </div>
        </div>
      )}

      {impacts.length > 0 && (
        <div className="panel-section">
          <h3>Meteor Impacts</h3>
          <div className="impacts-list">
            {impacts.map((impact, index) => (
              <div key={impact.id} className="impact-item">
                <div className="impact-header">
                  <span className="impact-number">#{index + 1}</span>
                  <span className="impact-location">
                    {impact.location.lat.toFixed(2)}, {impact.location.lng.toFixed(2)}
                  </span>
                </div>
                <div className="impact-details">
                  <div className="impact-stat">
                    <span>Radius:</span>
                    <span>{impact.radius} km</span>
                  </div>
                  <div className="impact-stat">
                    <span>Intensity:</span>
                    <span>{impact.intensity}/10</span>
                  </div>
                  <div className="impact-stat">
                    <span>Affected:</span>
                    <span>{impact.affectedCountries.length} countries</span>
                  </div>
                </div>
              </div>
            ))}
          </div>
        </div>
      )}

      <div className="panel-section">
        <h3>Instructions</h3>
        <div className="instructions">
          <p><strong>Navigation:</strong></p>
          <ul>
            <li>Drag to rotate the globe</li>
            <li>Scroll to zoom in/out</li>
            <li>Click anywhere for location info</li>
          </ul>
          
          {simulationMode && (
            <>
              <p><strong>Simulation Mode:</strong></p>
              <ul>
                <li>Click on ocean areas to place meteor impacts</li>
                <li>View hazard rings and affected areas</li>
                <li>Analyze potential damage zones</li>
              </ul>
            </>
          )}
        </div>
      </div>
    </div>
  );
};

export default GlobeSidePanel;
