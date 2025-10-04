import React, { useEffect, useRef } from 'react';
import Globe from 'globe.gl';
import { SimpleGlobeProps, Location } from '../utils/types';

const SimpleGlobe: React.FC<SimpleGlobeProps> = ({ 
  width = 800, 
  height = 600,
  onLocationClick,
  simulationMode = false,
  impacts = [],
  centerOnLocation
}) => {
  const globeRef = useRef<HTMLDivElement>(null);
  const globeInstance = useRef<any>(null);

  useEffect(() => {
    if (!globeRef.current) return;

    // Initialize the globe
    const globe = Globe()(globeRef.current)
      .globeImageUrl('//unpkg.com/three-globe/example/img/earth-blue-marble.jpg')
      .bumpImageUrl('//unpkg.com/three-globe/example/img/earth-topology.png')
      .backgroundImageUrl('//unpkg.com/three-globe/example/img/night-sky.png')
      .width(width)
      .height(height)
      .backgroundColor('rgba(0,0,0,0)')
      .showAtmosphere(true)
      .atmosphereColor('#00ccff')
      .atmosphereAltitude(0.15);

    // Store reference
    globeInstance.current = globe;

    // Note: Auto-rotation removed due to API compatibility issues
    // The globe can still be manually rotated by dragging

    // Add some sample data for visualization
    const sampleData = [
      { lat: 40.7128, lng: -74.0060, size: 0.5, color: '#ff6b6b' }, // New York
      { lat: 51.5074, lng: -0.1278, size: 0.5, color: '#4ecdc4' }, // London
      { lat: 35.6762, lng: 139.6503, size: 0.5, color: '#45b7d1' }, // Tokyo
      { lat: -33.8688, lng: 151.2093, size: 0.5, color: '#96ceb4' }, // Sydney
    ];

    // Add sample points
    globe.pointsData(sampleData)
      .pointColor('color')
      .pointAltitude(0.01)
      .pointRadius('size');

    // Add click handler for location lookup
    console.log('🔧 Setting up globe click handler...');
    
    // Try multiple event handling approaches
    globe.onGlobeClick((event: any) => {
      console.log('🌍 ===== GLOBE CLICKED =====');
      console.log('🌍 Raw event object:', event);
      console.log('🌍 Event type:', typeof event);
      console.log('🌍 Event keys:', event ? Object.keys(event) : 'No event');
      console.log('🌍 onLocationClick function:', onLocationClick);
      
      // Try different ways to extract coordinates
      let lat, lng;
      
      // Method 1: Direct properties
      if (event && typeof event.lat === 'number' && typeof event.lng === 'number' && !isNaN(event.lat) && !isNaN(event.lng)) {
        lat = event.lat;
        lng = event.lng;
        console.log('✅ Method 1 - Direct lat/lng:', { lat, lng });
      }
      // Method 2: Point object
      else if (event && event.point && typeof event.point.lat === 'number' && typeof event.point.lng === 'number' && !isNaN(event.point.lat) && !isNaN(event.point.lng)) {
        lat = event.point.lat;
        lng = event.point.lng;
        console.log('✅ Method 2 - Point coordinates:', { lat, lng });
      }
      // Method 3: Coordinates array
      else if (event && event.coordinates && Array.isArray(event.coordinates) && event.coordinates.length >= 2) {
        lng = event.coordinates[0];
        lat = event.coordinates[1];
        console.log('✅ Method 3 - Array coordinates:', { lat, lng });
      }
      // Method 4: Try to get coordinates from the globe itself
      else if (globeInstance.current) {
        console.log('🔄 Method 4 - Trying to get coordinates from globe instance...');
        try {
          // Try to get the current point of view or use a default location
          const pov = globeInstance.current.pointOfView();
          if (pov && typeof pov.lat === 'number' && typeof pov.lng === 'number') {
            lat = pov.lat;
            lng = pov.lng;
            console.log('✅ Method 4 - Got coordinates from POV:', { lat, lng });
          } else {
            // Use a default location (center of globe)
            lat = 0;
            lng = 0;
            console.log('✅ Method 4 - Using default coordinates:', { lat, lng });
          }
        } catch (error) {
          console.log('❌ Method 4 failed:', error);
          lat = 0;
          lng = 0;
          console.log('✅ Method 4 - Using fallback coordinates:', { lat, lng });
        }
      }
      else {
        console.warn('❌ All methods failed, using default coordinates');
        lat = 0;
        lng = 0;
      }
      
      console.log('🎯 Final coordinates:', { lat, lng });
      console.log('🎯 Coordinate types:', { latType: typeof lat, lngType: typeof lng });
      console.log('🎯 Are coordinates valid?', typeof lat === 'number' && typeof lng === 'number' && !isNaN(lat) && !isNaN(lng));
      
      if (onLocationClick && typeof lat === 'number' && typeof lng === 'number' && !isNaN(lat) && !isNaN(lng)) {
        console.log('📞 Calling onLocationClick with:', { lat, lng });
        const location: Location = { lat, lng };
        onLocationClick(location);
        console.log('✅ onLocationClick called successfully');
      } else {
        console.warn('❌ onLocationClick not available or invalid coordinates');
        console.log('❌ onLocationClick available:', !!onLocationClick);
        console.log('❌ lat valid:', typeof lat === 'number' && !isNaN(lat));
        console.log('❌ lng valid:', typeof lng === 'number' && !isNaN(lng));
      }
      
      console.log('🌍 ===== CLICK HANDLER END =====');
    });
    
    console.log('✅ Globe click handler set up complete');

    // Cleanup function
    return () => {
      if (globeInstance.current) {
        globeInstance.current._destructor?.();
      }
    };
  }, [width, height, onLocationClick]);

  // Update impacts when they change
  useEffect(() => {
    if (globeInstance.current && impacts.length > 0) {
      // Convert impacts to points data
      const impactPoints = impacts.map(impact => ({
        lat: impact.location.lat,
        lng: impact.location.lng,
        size: impact.radius / 100, // Scale down for visualization
        color: '#ff6b6b'
      }));

      globeInstance.current.pointsData(impactPoints)
        .pointColor('color')
        .pointAltitude(0.01)
        .pointRadius('size');
    }
  }, [impacts]);

  // Center globe on location when centerOnLocation changes
  useEffect(() => {
    if (globeInstance.current && centerOnLocation) {
      console.log('Centering globe on:', centerOnLocation);
      // Use the globe's pointOfView method to center on coordinates
      globeInstance.current.pointOfView({
        lat: centerOnLocation.lat,
        lng: centerOnLocation.lng,
        altitude: 2
      }, 1000); // 1 second animation
    }
  }, [centerOnLocation]);

  return (
    <div 
      ref={globeRef} 
      style={{ 
        width: `${width}px`, 
        height: `${height}px`,
        margin: '0 auto',
        borderRadius: '10px',
        overflow: 'hidden',
        boxShadow: '0 0 20px rgba(0, 204, 255, 0.3)'
      }} 
    />
  );
};

export default SimpleGlobe;
