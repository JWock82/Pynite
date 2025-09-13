"""
Pynite Engine - Backend interface for Pynite finite element analysis
"""

import uuid
import io
import base64
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, List, Any, Optional, Tuple
import sys
import os

# Add the current directory to Python path to import Pynite
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from Pynite import FEModel3D

class PyniteEngine:
    """
    Engine class that manages Pynite models and provides a clean interface
    for the Flask web application
    """
    
    def __init__(self):
        self.models: Dict[str, FEModel3D] = {}
        self.model_metadata: Dict[str, Dict[str, Any]] = {}
    
    def create_model(self, name: str) -> str:
        """Create a new finite element model"""
        model_id = str(uuid.uuid4())
        self.models[model_id] = FEModel3D()
        self.model_metadata[model_id] = {
            'name': name,
            'created': True,
            'analyzed': False,
            'nodes_count': 0,
            'members_count': 0,
            'materials_count': 0,
            'sections_count': 0
        }
        return model_id
    
    def add_node(self, model_id: str, name: str, x: float, y: float, z: float) -> None:
        """Add a node to the specified model"""
        if model_id not in self.models:
            raise ValueError(f"Model {model_id} not found")
        
        model = self.models[model_id]
        model.add_node(name, x, y, z)
        self.model_metadata[model_id]['nodes_count'] += 1
    
    def add_material(self, model_id: str, name: str, E: float, G: float, nu: float, rho: float) -> None:
        """Add a material to the specified model"""
        if model_id not in self.models:
            raise ValueError(f"Model {model_id} not found")
        
        model = self.models[model_id]
        model.add_material(name, E, G, nu, rho)
        self.model_metadata[model_id]['materials_count'] += 1
    
    def add_section(self, model_id: str, name: str, A: float, Iy: float, Iz: float, J: float) -> None:
        """Add a section to the specified model"""
        if model_id not in self.models:
            raise ValueError(f"Model {model_id} not found")
        
        model = self.models[model_id]
        model.add_section(name, A, Iy, Iz, J)
        self.model_metadata[model_id]['sections_count'] += 1
    
    def add_member(self, model_id: str, name: str, i_node: str, j_node: str, 
                   material: str, section: str) -> None:
        """Add a member to the specified model"""
        if model_id not in self.models:
            raise ValueError(f"Model {model_id} not found")
        
        model = self.models[model_id]
        model.add_member(name, i_node, j_node, material, section)
        self.model_metadata[model_id]['members_count'] += 1
    
    def add_support(self, model_id: str, node_name: str, support_dx: bool, support_dy: bool,
                   support_dz: bool, support_rx: bool, support_ry: bool, support_rz: bool) -> None:
        """Add support conditions to a node"""
        if model_id not in self.models:
            raise ValueError(f"Model {model_id} not found")
        
        model = self.models[model_id]
        model.def_support(node_name, support_dx, support_dy, support_dz, 
                         support_rx, support_ry, support_rz)
    
    def add_node_load(self, model_id: str, node_name: str, direction: str, 
                     magnitude: float, case: str = 'Case 1') -> None:
        """Add a nodal load to the model"""
        if model_id not in self.models:
            raise ValueError(f"Model {model_id} not found")
        
        model = self.models[model_id]
        model.add_node_load(node_name, direction, magnitude, case)
    
    def add_member_point_load(self, model_id: str, member_name: str, direction: str,
                             magnitude: float, location: float, case: str = 'Case 1') -> None:
        """Add a point load to a member"""
        if model_id not in self.models:
            raise ValueError(f"Model {model_id} not found")
        
        model = self.models[model_id]
        model.add_member_pt_load(member_name, direction, magnitude, location, case)
    
    def add_member_distributed_load(self, model_id: str, member_name: str, direction: str,
                                   w1: float, w2: float, x1: float, x2: float, case: str = 'Case 1') -> None:
        """Add a distributed load to a member"""
        if model_id not in self.models:
            raise ValueError(f"Model {model_id} not found")
        
        model = self.models[model_id]
        if x2 == -1:  # Full length
            model.add_member_dist_load(member_name, direction, w1, w2, case=case)
        else:
            model.add_member_dist_load(member_name, direction, w1, w2, x1, x2, case)
    
    def add_load_combo(self, model_id: str, name: str, factors: Dict[str, float]) -> None:
        """Add a load combination to the model"""
        if model_id not in self.models:
            raise ValueError(f"Model {model_id} not found")
        
        model = self.models[model_id]
        model.add_load_combo(name, factors)
    
    def analyze(self, model_id: str, analysis_type: str = 'linear') -> Dict[str, Any]:
        """Analyze the specified model"""
        if model_id not in self.models:
            raise ValueError(f"Model {model_id} not found")
        
        model = self.models[model_id]
        
        # Perform analysis based on type
        if analysis_type == 'linear':
            model.analyze_linear()
        elif analysis_type == 'pdelta':
            model.analyze_PDelta()
        else:
            model.analyze()
        
        # Mark as analyzed
        self.model_metadata[model_id]['analyzed'] = True
        
        # Return basic analysis summary
        return {
            'analysis_type': analysis_type,
            'load_combinations': list(model.load_combos.keys()),
            'nodes_count': len(model.nodes),
            'members_count': len(model.members),
            'converged': True
        }
    
    def get_results(self, model_id: str, combo_name: str) -> Dict[str, Any]:
        """Get analysis results for a specific load combination"""
        if model_id not in self.models:
            raise ValueError(f"Model {model_id} not found")
        
        model = self.models[model_id]
        
        # Node results
        node_results = {}
        for node_name, node in model.nodes.items():
            node_results[node_name] = {
                'displacements': {
                    'DX': float(node.DX.get(combo_name, 0)),
                    'DY': float(node.DY.get(combo_name, 0)),
                    'DZ': float(node.DZ.get(combo_name, 0)),
                    'RX': float(node.RX.get(combo_name, 0)),
                    'RY': float(node.RY.get(combo_name, 0)),
                    'RZ': float(node.RZ.get(combo_name, 0))
                },
                'reactions': {
                    'FX': float(node.RxnFX.get(combo_name, 0)),
                    'FY': float(node.RxnFY.get(combo_name, 0)),
                    'FZ': float(node.RxnFZ.get(combo_name, 0)),
                    'MX': float(node.RxnMX.get(combo_name, 0)),
                    'MY': float(node.RxnMY.get(combo_name, 0)),
                    'MZ': float(node.RxnMZ.get(combo_name, 0))
                }
            }
        
        # Member results
        member_results = {}
        for member_name, member in model.members.items():
            try:
                member_results[member_name] = {
                    'max_moment_Mz': float(member.max_moment('Mz', combo_name)),
                    'min_moment_Mz': float(member.min_moment('Mz', combo_name)),
                    'max_moment_My': float(member.max_moment('My', combo_name)),
                    'min_moment_My': float(member.min_moment('My', combo_name)),
                    'max_shear_Fy': float(member.max_shear('Fy', combo_name)),
                    'min_shear_Fy': float(member.min_shear('Fy', combo_name)),
                    'max_shear_Fz': float(member.max_shear('Fz', combo_name)),
                    'min_shear_Fz': float(member.min_shear('Fz', combo_name)),
                    'max_axial': float(member.max_axial(combo_name)),
                    'min_axial': float(member.min_axial(combo_name)),
                    'max_deflection_dy': float(member.max_deflection('dy', combo_name)),
                    'min_deflection_dy': float(member.min_deflection('dy', combo_name)),
                    'max_deflection_dz': float(member.max_deflection('dz', combo_name)),
                    'min_deflection_dz': float(member.min_deflection('dz', combo_name))
                }
            except Exception as e:
                member_results[member_name] = {'error': str(e)}
        
        return {
            'nodes': node_results,
            'members': member_results
        }
    
    def get_visualization_data(self, model_id: str, combo_name: str) -> Dict[str, Any]:
        """Get visualization data for the model"""
        if model_id not in self.models:
            raise ValueError(f"Model {model_id} not found")
        
        model = self.models[model_id]
        
        # Node data
        nodes = []
        for node_name, node in model.nodes.items():
            nodes.append({
                'name': node_name,
                'x': float(node.X),
                'y': float(node.Y),
                'z': float(node.Z),
                'displacement': {
                    'x': float(node.DX.get(combo_name, 0)),
                    'y': float(node.DY.get(combo_name, 0)),
                    'z': float(node.DZ.get(combo_name, 0))
                },
                'supports': {
                    'DX': bool(node.support_DX),
                    'DY': bool(node.support_DY),
                    'DZ': bool(node.support_DZ),
                    'RX': bool(node.support_RX),
                    'RY': bool(node.support_RY),
                    'RZ': bool(node.support_RZ)
                }
            })
        
        # Member data
        members = []
        for member_name, member in model.members.items():
            members.append({
                'name': member_name,
                'i_node': member.i_node.name,
                'j_node': member.j_node.name,
                'i_coords': [float(member.i_node.X), float(member.i_node.Y), float(member.i_node.Z)],
                'j_coords': [float(member.j_node.X), float(member.j_node.Y), float(member.j_node.Z)]
            })
        
        return {
            'nodes': nodes,
            'members': members
        }
    
    def generate_member_diagram(self, model_id: str, member_name: str, diagram_type: str,
                               direction: str, combo_name: str) -> str:
        """Generate a member diagram and return as base64 encoded image"""
        if model_id not in self.models:
            raise ValueError(f"Model {model_id} not found")
        
        model = self.models[model_id]
        member = model.members[member_name]
        
        # Create the plot
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Generate data points along the member
        n_points = 100
        x_array = np.linspace(0, member.L(), n_points)
        
        if diagram_type == 'moment':
            y_array = [member.moment(direction, x, combo_name) for x in x_array]
            ax.set_ylabel(f'Moment {direction}')
            ax.set_title(f'Moment Diagram - Member {member_name} - {combo_name}')
        elif diagram_type == 'shear':
            y_array = [member.shear(direction, x, combo_name) for x in x_array]
            ax.set_ylabel(f'Shear {direction}')
            ax.set_title(f'Shear Diagram - Member {member_name} - {combo_name}')
        elif diagram_type == 'deflection':
            y_array = [member.deflection(direction, x, combo_name) for x in x_array]
            ax.set_ylabel(f'Deflection {direction}')
            ax.set_title(f'Deflection Diagram - Member {member_name} - {combo_name}')
        elif diagram_type == 'axial':
            y_array = [member.axial(x, combo_name) for x in x_array]
            ax.set_ylabel('Axial Force')
            ax.set_title(f'Axial Force Diagram - Member {member_name} - {combo_name}')
        else:
            raise ValueError(f"Unknown diagram type: {diagram_type}")
        
        # Plot the diagram
        ax.plot(x_array, y_array, 'b-', linewidth=2)
        ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)
        ax.grid(True, alpha=0.3)
        ax.set_xlabel('Distance along member')
        
        # Convert plot to base64 string
        img_buffer = io.BytesIO()
        plt.savefig(img_buffer, format='png', dpi=150, bbox_inches='tight')
        img_buffer.seek(0)
        img_data = base64.b64encode(img_buffer.getvalue()).decode()
        plt.close(fig)
        
        return img_data
    
    def list_models(self) -> List[Dict[str, Any]]:
        """List all available models"""
        models_list = []
        for model_id, metadata in self.model_metadata.items():
            models_list.append({
                'id': model_id,
                'name': metadata['name'],
                'analyzed': metadata['analyzed'],
                'nodes_count': metadata['nodes_count'],
                'members_count': metadata['members_count']
            })
        return models_list
    
    def get_model_info(self, model_id: str) -> Dict[str, Any]:
        """Get detailed information about a model"""
        if model_id not in self.models:
            raise ValueError(f"Model {model_id} not found")
        
        model = self.models[model_id]
        metadata = self.model_metadata[model_id]
        
        # Get nodes info
        nodes_info = []
        for node_name, node in model.nodes.items():
            nodes_info.append({
                'name': node_name,
                'x': float(node.X),
                'y': float(node.Y),
                'z': float(node.Z),
                'supports': {
                    'DX': bool(node.support_DX),
                    'DY': bool(node.support_DY),
                    'DZ': bool(node.support_DZ),
                    'RX': bool(node.support_RX),
                    'RY': bool(node.support_RY),
                    'RZ': bool(node.support_RZ)
                }
            })
        
        # Get members info
        members_info = []
        for member_name, member in model.members.items():
            members_info.append({
                'name': member_name,
                'i_node': member.i_node.name,
                'j_node': member.j_node.name,
                'material': member.material.name,
                'section': member.section.name,
                'length': float(member.L())
            })
        
        # Get materials info
        materials_info = []
        for material_name, material in model.materials.items():
            materials_info.append({
                'name': material_name,
                'E': float(material.E),
                'G': float(material.G),
                'nu': float(material.nu),
                'rho': float(material.rho)
            })
        
        # Get sections info
        sections_info = []
        for section_name, section in model.sections.items():
            sections_info.append({
                'name': section_name,
                'A': float(section.A),
                'Iy': float(section.Iy),
                'Iz': float(section.Iz),
                'J': float(section.J)
            })
        
        # Get load combinations info
        load_combos_info = []
        for combo_name, combo in model.load_combos.items():
            load_combos_info.append({
                'name': combo_name,
                'factors': combo.factors
            })
        
        return {
            'id': model_id,
            'metadata': metadata,
            'nodes': nodes_info,
            'members': members_info,
            'materials': materials_info,
            'sections': sections_info,
            'load_combinations': load_combos_info
        }