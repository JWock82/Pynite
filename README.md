# Pynite Flask Web Application

A comprehensive web-based finite element analysis application built with Flask and Pynite, featuring an intuitive interface for structural modeling and analysis.

## Features

- **Interactive 3D Visualization**: Built with Three.js for real-time model viewing
- **Complete FEA Workflow**: Create nodes, members, materials, sections, loads, and supports
- **Analysis Capabilities**: Linear analysis with support for P-Delta analysis
- **Results Visualization**: View displacements, reactions, and member force diagrams
- **Modern UI**: Clean, responsive interface with modal dialogs and real-time updates

## Installation

1. **Install Python Dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

2. **Install Node.js Dependencies** (optional, for development):
   ```bash
   npm install
   ```

## Running the Application

1. **Start the Flask Server**:
   ```bash
   python app.py
   ```

2. **Open your browser** and navigate to:
   ```
   http://localhost:5000
   ```

## Usage

### Creating a Model

1. Click "New Model" in the header
2. Enter a model name and click "Create"

### Building Your Structure

1. **Add Materials**: Define material properties (E, G, ν, ρ)
2. **Add Sections**: Define cross-sectional properties (A, Iy, Iz, J)
3. **Add Nodes**: Create nodes at specific coordinates
4. **Add Members**: Connect nodes with structural members
5. **Add Supports**: Define boundary conditions
6. **Add Loads**: Apply forces and moments

### Analysis

1. Click "Analyze" to run the finite element analysis
2. View results in the Results panel:
   - **Nodes Tab**: Displacements and reactions
   - **Members Tab**: Internal forces and moments
   - **Diagrams Tab**: Force and moment diagrams

### Visualization Options

- **Deformed Shape**: Show displaced structure
- **Show Loads**: Display applied loads
- **Show Reactions**: Display support reactions

## API Endpoints

The application provides a RESTful API for programmatic access:

- `POST /api/models` - Create new model
- `POST /api/models/{id}/nodes` - Add node
- `POST /api/models/{id}/materials` - Add material
- `POST /api/models/{id}/sections` - Add section
- `POST /api/models/{id}/members` - Add member
- `POST /api/models/{id}/supports` - Add support
- `POST /api/models/{id}/loads` - Add loads
- `POST /api/models/{id}/analyze` - Run analysis
- `GET /api/models/{id}/results/{combo}` - Get results

## Technology Stack

- **Backend**: Flask (Python)
- **Frontend**: Vanilla JavaScript, Three.js
- **FEA Engine**: Pynite
- **Visualization**: Three.js, Matplotlib
- **Styling**: Custom CSS with modern design principles

## Project Structure

```
├── app.py                 # Flask application
├── pynite_engine.py      # Pynite interface layer
├── requirements.txt      # Python dependencies
├── templates/
│   └── index.html       # Main application template
├── static/
│   ├── css/
│   │   └── styles.css   # Application styles
│   └── js/
│       ├── app.js       # Main application logic
│       ├── three-viewer.js  # 3D visualization
│       └── modals.js    # Modal management
└── Pynite/              # Pynite library files
```

## Contributing

This application demonstrates the integration of Pynite with modern web technologies. Feel free to extend it with additional features such as:

- Advanced load types
- Material nonlinearity
- Dynamic analysis
- Report generation
- Model import/export

## License

This project uses the same MIT license as Pynite.