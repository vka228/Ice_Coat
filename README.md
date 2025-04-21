# Quantifying the Dynamics of Bridge Formation in Crystal Growth Phenomena

<img src="/madia_gt/intro.gif" width="250" height="250"/> <img src="/madia_gt/demo_numergif.gif" width="250" height="250"/>


## Overview
This project simulates heat transfer and phase change (water to ice) in a 2D system containing cold pipes submerged in water. The numerical model solves the heat equation using finite difference methods, tracking temperature distribution and ice formation over time. The simulation visualizes the cooling process and records temperature evolution at specific monitoring points.

## Key Components

### Main Simulation (`numer.py`)
- Implements the core finite difference algorithm
- Manages temperature field calculations
- Handles phase change detection (water → ice)
- Tracks boundary conditions and material properties
- Generates visualization frames and data outputs

## Modeling Process

### Heat Transfer Physics
1. **Governing Equations**: Uses the heat equation with different thermal diffusivities for water and ice
   - ∂T/∂t = α∇²T where α is thermal diffusivity
   - Different α values for water (1.4×10⁻¹) and ice (1.6×10⁻¹)

2. **Numerical Method**:
   - Explicit finite difference scheme
   - Stability condition: Δt ≤ min(Δx²/(4α), Δy²/(4α))
   - Spatial discretization: 100×100 grid over 120mm domain

3. **Boundary Conditions**:
   - Outer boundaries maintained at constant temperature (275K)
   - Pipe surfaces fixed at 253K
   - Circular pipe geometry implemented with radius 20mm

### Phase Change Handling
1. **Detection**:
   - Tracks nodes where temperature crosses 273K (freezing point)
   - Maintains state map (`cond_map`) distinguishing water (0), ice (1), and pipe (2)

2. **Energy Considerations**:
   - Calculates latent heat of fusion (q_phase_change)
   - Uses material properties:
     - Density (ρ) = 10⁻⁶ (scaled)
     - Latent heat (λ) = 330,000 J/kg

3. **Interface Tracking**:
   - Identifies boundary nodes between ice and water
   - Marks these for visualization

### Data Collection
1. **Temperature Monitoring**:
   - Tracks three specific points (60,80), (60,95), (60,110)
   - Records temperature evolution over time
   - Outputs to text files (results_tp1.txt, etc.)

2. **Visualization**:
   - Color-mapped temperature fields
   - Time progression animation
   - Boundary markers for ice-water interface
   - Temperature-time plots for monitoring points

## Implementation Details

### Key Functions
- `setpipe()`: Initializes pipe geometry and temperature
- `update_cond_map()`: Updates material state (water/ice)
- `boundary()`: Detects ice-water interface nodes
- `gettemp()`: Records temperature at specified points
- `main()`: Core simulation loop

### Parameters
- Spatial domain: 120mm × 120mm
- Time duration: 6000 seconds
- Pipe temperature: 253K
- Ambient temperature: 275K
- Grid resolution: 100×100 nodes
- Thermal diffusivities: 
  - Water: 0.14 mm²/s
  - Ice: 0.16 mm²/s
