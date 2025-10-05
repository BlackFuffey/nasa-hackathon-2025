import React, { useEffect, useState } from "react";
import { NeoAsteroid, fetchNEOData } from "../data/fetch"; // adjust import as needed
import "./AsteroidPanel.css";

interface AsteroidPanelProps {
  onSelect: (a: NeoAsteroid | null) => void;
  selected: NeoAsteroid | null;
  disabled?: boolean; // disables clicks if needed
}

export const AsteroidPanel: React.FC<AsteroidPanelProps> = ({
  onSelect,
  selected,
  disabled,
}) => {
  const [asteroids, setAsteroids] = useState<NeoAsteroid[]>([]);
  const [search, setSearch] = useState("");

  useEffect(() => {
    fetchNEOData().then((data) => setAsteroids(data));
  }, []);

  const filtered = asteroids
    .filter((a) =>
      a.name.toLowerCase().includes(search.toLowerCase()) ||
      a.id.includes(search)
    )
    .slice(0, 15);

  if (selected) {
    const sz = selected.estimated_diameter.meters;
    const vel = selected.close_approach_data[0]?.relative_velocity?.kilometers_per_second || "";
    return (
      <div className="asteroid-panel">
        <button className="back-btn" onClick={() => onSelect(null)} disabled={disabled}>
          ← Back to List
        </button>
        <div className="asteroid-card selected-card">
          <div className="asteroid-title">{selected.name}</div>
          <div><b>ID:</b> {selected.id}</div>
          <div><b>Diameter:</b> {sz.estimated_diameter_min.toFixed(0)}–{sz.estimated_diameter_max.toFixed(0)} m</div>
          <div><b>Velocity:</b> {Number(vel).toFixed(2)} km/s</div>
          <div>
            <b>Hazard:</b>{' '}
            <span className={selected.is_potentially_hazardous_asteroid ? "hazard" : "safe"}>
              {selected.is_potentially_hazardous_asteroid ? "Hazardous" : "Safe"}
            </span>
          </div>
          <button
            className="choose-asteroid-btn"
            onClick={() => onSelect(selected)}
            disabled={disabled}
          >
            Select This Asteroid
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="asteroid-panel">
      <h2 className="panel-title">Asteroid Database</h2>
      <input
        className="search-box"
        placeholder="Filter by name, size..."
        value={search}
        onChange={(e) => setSearch(e.target.value)}
        disabled={disabled}
      />
      <div className="asteroid-list">
        {filtered.map((a) => {
          const sz = a.estimated_diameter.meters;
          const vel = a.close_approach_data[0]?.relative_velocity?.kilometers_per_second || "";
          return (
            <div
              key={a.id}
              className="asteroid-card"
              onClick={() => !disabled && onSelect(a)}
              style={{ opacity: disabled ? 0.6 : 1, pointerEvents: disabled ? 'none' : 'auto' }}
            >
              <div className="asteroid-title">{a.name}</div>
              <div><b>ID:</b> {a.id}</div>
              <div><b>Size:</b> {sz.estimated_diameter_min.toFixed(0)}–{sz.estimated_diameter_max.toFixed(0)} m</div>
              <div><b>Vel:</b> {Number(vel).toFixed(2)} km/s</div>
              <span className={a.is_potentially_hazardous_asteroid ? "hazard" : "safe"}>
                {a.is_potentially_hazardous_asteroid ? "Hazardous" : "Safe"}
              </span>
            </div>
          );
        })}
      </div>
    </div>
  );
};
