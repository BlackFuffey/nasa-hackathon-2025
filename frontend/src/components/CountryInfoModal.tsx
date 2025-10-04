import React from 'react';
import { CountryInfoModalProps } from '../utils/types';
import './CountryInfoModal.css';

const CountryInfoModal: React.FC<CountryInfoModalProps> = ({
  isOpen,
  onClose,
  location,
  country
}) => {
  if (!isOpen) return null;

  return (
    <div className="modal-overlay" onClick={onClose}>
      <div className="modal-content" onClick={(e) => e.stopPropagation()}>
        <div className="modal-header">
          <h2>Location Details</h2>
          <button className="close-button" onClick={onClose}>
            ×
          </button>
        </div>

        <div className="modal-body">
          {location && (
            <div className="location-section">
              <h3>Coordinates</h3>
              <div className="coordinate-display">
                <div className="coord-item">
                  <span className="coord-label">Latitude:</span>
                  <span className="coord-value">
                    {typeof location.lat === 'number' && !isNaN(location.lat) 
                      ? location.lat.toFixed(6) + '°' 
                      : 'Invalid'}
                  </span>
                </div>
                <div className="coord-item">
                  <span className="coord-label">Longitude:</span>
                  <span className="coord-value">
                    {typeof location.lng === 'number' && !isNaN(location.lng) 
                      ? location.lng.toFixed(6) + '°' 
                      : 'Invalid'}
                  </span>
                </div>
              </div>
              
              {location.name && (
                <div className="location-name">
                  <h4>{location.name}</h4>
                </div>
              )}

              <div className="location-details">
                {location.country && (
                  <div className="detail-item">
                    <span className="detail-label">Country:</span>
                    <span className="detail-value">{location.country}</span>
                  </div>
                )}
                {location.continent && (
                  <div className="detail-item">
                    <span className="detail-label">Continent:</span>
                    <span className="detail-value">{location.continent}</span>
                  </div>
                )}
              </div>
            </div>
          )}

          {country && (
            <div className="country-section">
              <h3>Country Information</h3>
              <div className="country-header">
                <h4>{country.name}</h4>
                <span className="country-code">{country.code}</span>
              </div>

              <div className="country-stats">
                <div className="stat-grid">
                  <div className="stat-item">
                    <div className="stat-label">Capital</div>
                    <div className="stat-value">{country.capital}</div>
                  </div>
                  <div className="stat-item">
                    <div className="stat-label">Population</div>
                    <div className="stat-value">
                      {country.population.toLocaleString()}
                    </div>
                  </div>
                  <div className="stat-item">
                    <div className="stat-label">Area</div>
                    <div className="stat-value">
                      {country.area.toLocaleString()} km²
                    </div>
                  </div>
                  <div className="stat-item">
                    <div className="stat-label">Continent</div>
                    <div className="stat-value">{country.continent}</div>
                  </div>
                </div>

                <div className="population-density">
                  <div className="stat-label">Population Density</div>
                  <div className="stat-value">
                    {(country.population / country.area).toFixed(2)} people/km²
                  </div>
                </div>
              </div>
            </div>
          )}

          {!location && !country && (
            <div className="no-data">
              <p>No location data available</p>
            </div>
          )}
        </div>

        <div className="modal-footer">
          <button className="close-modal-button" onClick={onClose}>
            Close
          </button>
        </div>
      </div>
    </div>
  );
};

export default CountryInfoModal;
