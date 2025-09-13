from flask import Flask, request, jsonify, render_template, send_file
from flask_cors import CORS
import json
import io
import base64
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from typing import Dict, List, Any, Optional
import traceback
import os

# Import Pynite components
from pynite_engine import PyniteEngine

app = Flask(__name__)
CORS(app)

# Initialize the Pynite engine
engine = PyniteEngine()

@app.route('/')
def index():
    """Serve the main application page"""
    return render_template('index.html')

@app.route('/api/models', methods=['POST'])
def create_model():
    """Create a new finite element model"""
    try:
        data = request.get_json()
        model_name = data.get('name', 'Model1')
        
        model_id = engine.create_model(model_name)
        
        return jsonify({
            'success': True,
            'model_id': model_id,
            'message': f'Model "{model_name}" created successfully'
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e),
            'traceback': traceback.format_exc()
        }), 500

@app.route('/api/models/<model_id>/nodes', methods=['POST'])
def add_node(model_id):
    """Add a node to the model"""
    try:
        data = request.get_json()
        node_name = data['name']
        x = float(data['x'])
        y = float(data['y'])
        z = float(data['z'])
        
        engine.add_node(model_id, node_name, x, y, z)
        
        return jsonify({
            'success': True,
            'message': f'Node "{node_name}" added successfully'
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/models/<model_id>/materials', methods=['POST'])
def add_material(model_id):
    """Add a material to the model"""
    try:
        data = request.get_json()
        name = data['name']
        E = float(data['E'])
        G = float(data['G'])
        nu = float(data['nu'])
        rho = float(data['rho'])
        
        engine.add_material(model_id, name, E, G, nu, rho)
        
        return jsonify({
            'success': True,
            'message': f'Material "{name}" added successfully'
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/models/<model_id>/sections', methods=['POST'])
def add_section(model_id):
    """Add a section to the model"""
    try:
        data = request.get_json()
        name = data['name']
        A = float(data['A'])
        Iy = float(data['Iy'])
        Iz = float(data['Iz'])
        J = float(data['J'])
        
        engine.add_section(model_id, name, A, Iy, Iz, J)
        
        return jsonify({
            'success': True,
            'message': f'Section "{name}" added successfully'
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/models/<model_id>/members', methods=['POST'])
def add_member(model_id):
    """Add a member to the model"""
    try:
        data = request.get_json()
        name = data['name']
        i_node = data['i_node']
        j_node = data['j_node']
        material = data['material']
        section = data['section']
        
        engine.add_member(model_id, name, i_node, j_node, material, section)
        
        return jsonify({
            'success': True,
            'message': f'Member "{name}" added successfully'
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/models/<model_id>/supports', methods=['POST'])
def add_support(model_id):
    """Add support conditions to a node"""
    try:
        data = request.get_json()
        node_name = data['node']
        support_dx = data.get('support_DX', False)
        support_dy = data.get('support_DY', False)
        support_dz = data.get('support_DZ', False)
        support_rx = data.get('support_RX', False)
        support_ry = data.get('support_RY', False)
        support_rz = data.get('support_RZ', False)
        
        engine.add_support(model_id, node_name, support_dx, support_dy, support_dz, 
                          support_rx, support_ry, support_rz)
        
        return jsonify({
            'success': True,
            'message': f'Support added to node "{node_name}"'
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/models/<model_id>/loads', methods=['POST'])
def add_load(model_id):
    """Add loads to the model"""
    try:
        data = request.get_json()
        load_type = data['type']
        
        if load_type == 'node':
            node_name = data['node']
            direction = data['direction']
            magnitude = float(data['magnitude'])
            case = data.get('case', 'Case 1')
            
            engine.add_node_load(model_id, node_name, direction, magnitude, case)
            message = f'Node load added to "{node_name}"'
            
        elif load_type == 'member_point':
            member_name = data['member']
            direction = data['direction']
            magnitude = float(data['magnitude'])
            location = float(data['location'])
            case = data.get('case', 'Case 1')
            
            engine.add_member_point_load(model_id, member_name, direction, magnitude, location, case)
            message = f'Point load added to member "{member_name}"'
            
        elif load_type == 'member_distributed':
            member_name = data['member']
            direction = data['direction']
            w1 = float(data['w1'])
            w2 = float(data['w2'])
            x1 = float(data.get('x1', 0))
            x2 = float(data.get('x2', -1))  # -1 indicates full length
            case = data.get('case', 'Case 1')
            
            engine.add_member_distributed_load(model_id, member_name, direction, w1, w2, x1, x2, case)
            message = f'Distributed load added to member "{member_name}"'
        
        return jsonify({
            'success': True,
            'message': message
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/models/<model_id>/load_combos', methods=['POST'])
def add_load_combo(model_id):
    """Add a load combination to the model"""
    try:
        data = request.get_json()
        name = data['name']
        factors = data['factors']
        
        engine.add_load_combo(model_id, name, factors)
        
        return jsonify({
            'success': True,
            'message': f'Load combination "{name}" added successfully'
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/models/<model_id>/analyze', methods=['POST'])
def analyze_model(model_id):
    """Analyze the finite element model"""
    try:
        data = request.get_json()
        analysis_type = data.get('type', 'linear')
        
        results = engine.analyze(model_id, analysis_type)
        
        return jsonify({
            'success': True,
            'results': results,
            'message': 'Analysis completed successfully'
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e),
            'traceback': traceback.format_exc()
        }), 500

@app.route('/api/models/<model_id>/results/<combo_name>')
def get_results(model_id, combo_name):
    """Get analysis results for a specific load combination"""
    try:
        results = engine.get_results(model_id, combo_name)
        
        return jsonify({
            'success': True,
            'results': results
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/models/<model_id>/visualization')
def get_visualization(model_id):
    """Generate visualization data for the model"""
    try:
        combo_name = request.args.get('combo', 'Combo 1')
        viz_data = engine.get_visualization_data(model_id, combo_name)
        
        return jsonify({
            'success': True,
            'visualization': viz_data
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/models/<model_id>/diagram/<member_name>/<diagram_type>')
def get_member_diagram(model_id, member_name, diagram_type):
    """Generate member force/moment diagrams"""
    try:
        combo_name = request.args.get('combo', 'Combo 1')
        direction = request.args.get('direction', 'Mz')
        
        # Generate the diagram
        img_data = engine.generate_member_diagram(model_id, member_name, diagram_type, direction, combo_name)
        
        return jsonify({
            'success': True,
            'image': img_data,
            'content_type': 'image/png'
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/models')
def list_models():
    """List all available models"""
    try:
        models = engine.list_models()
        return jsonify({
            'success': True,
            'models': models
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/models/<model_id>')
def get_model_info(model_id):
    """Get detailed information about a model"""
    try:
        model_info = engine.get_model_info(model_id)
        return jsonify({
            'success': True,
            'model': model_info
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)