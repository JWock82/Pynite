/**
 * Main application JavaScript for Pynite FEA Web App
 */

class PyniteApp {
    constructor() {
        this.currentModelId = null;
        this.currentCombo = 'Combo 1';
        this.viewer = null;
        this.init();
    }

    init() {
        this.setupEventListeners();
        this.initializeViewer();
        this.updateStatus('Ready');
    }

    setupEventListeners() {
        // Header buttons
        document.getElementById('new-model-btn').addEventListener('click', () => this.showNewModelModal());
        document.getElementById('analyze-btn').addEventListener('click', () => this.analyzeModel());

        // Sidebar buttons
        document.getElementById('add-material-btn').addEventListener('click', () => this.showAddMaterialModal());
        document.getElementById('add-section-btn').addEventListener('click', () => this.showAddSectionModal());
        document.getElementById('add-node-btn').addEventListener('click', () => this.showAddNodeModal());
        document.getElementById('add-member-btn').addEventListener('click', () => this.showAddMemberModal());
        document.getElementById('add-node-load-btn').addEventListener('click', () => this.showAddNodeLoadModal());
        document.getElementById('add-member-load-btn').addEventListener('click', () => this.showAddMemberLoadModal());
        document.getElementById('add-support-btn').addEventListener('click', () => this.showAddSupportModal());
        document.getElementById('add-load-combo-btn').addEventListener('click', () => this.showAddLoadComboModal());

        // View controls
        document.getElementById('view-combo').addEventListener('change', (e) => {
            this.currentCombo = e.target.value;
            this.updateVisualization();
        });

        document.getElementById('show-deformed').addEventListener('change', () => this.updateVisualization());
        document.getElementById('show-loads').addEventListener('change', () => this.updateVisualization());
        document.getElementById('show-reactions').addEventListener('change', () => this.updateVisualization());

        // Results panel
        document.getElementById('close-results').addEventListener('click', () => this.hideResults());

        // Results tabs
        document.querySelectorAll('.tab-btn').forEach(btn => {
            btn.addEventListener('click', (e) => {
                document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
                e.target.classList.add('active');
                this.showResultsTab(e.target.dataset.tab);
            });
        });
    }

    initializeViewer() {
        this.viewer = new ThreeViewer('three-canvas');
    }

    async showNewModelModal() {
        const modal = new Modal('Create New Model', `
            <div class="form-group">
                <label for="model-name">Model Name</label>
                <input type="text" id="model-name" class="form-control" placeholder="Enter model name" value="New Model">
            </div>
        `, [
            { text: 'Cancel', class: 'btn-outline' },
            { text: 'Create', class: 'btn-primary', action: () => this.createModel() }
        ]);
        modal.show();
    }

    async createModel() {
        const modelName = document.getElementById('model-name').value;
        if (!modelName.trim()) {
            this.showError('Please enter a model name');
            return;
        }

        try {
            this.updateStatus('Creating model...');
            const response = await fetch('/api/models', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ name: modelName })
            });

            const result = await response.json();
            if (result.success) {
                this.currentModelId = result.model_id;
                this.updateModelInfo();
                this.enableAnalyzeButton();
                this.updateStatus(`Model "${modelName}" created`);
                Modal.close();
            } else {
                this.showError(result.error);
            }
        } catch (error) {
            this.showError('Failed to create model: ' + error.message);
        }
    }

    async showAddMaterialModal() {
        if (!this.currentModelId) {
            this.showError('Please create a model first');
            return;
        }

        const modal = new Modal('Add Material', `
            <div class="form-group">
                <label for="material-name">Material Name</label>
                <input type="text" id="material-name" class="form-control" placeholder="e.g., Steel, Concrete">
            </div>
            <div class="form-row">
                <div class="form-group">
                    <label for="material-E">Modulus of Elasticity (E)</label>
                    <input type="number" id="material-E" class="form-control" placeholder="29000" step="any">
                </div>
                <div class="form-group">
                    <label for="material-G">Shear Modulus (G)</label>
                    <input type="number" id="material-G" class="form-control" placeholder="11200" step="any">
                </div>
            </div>
            <div class="form-row">
                <div class="form-group">
                    <label for="material-nu">Poisson's Ratio (ν)</label>
                    <input type="number" id="material-nu" class="form-control" placeholder="0.3" step="any">
                </div>
                <div class="form-group">
                    <label for="material-rho">Density (ρ)</label>
                    <input type="number" id="material-rho" class="form-control" placeholder="0.000284" step="any">
                </div>
            </div>
        `, [
            { text: 'Cancel', class: 'btn-outline' },
            { text: 'Add Material', class: 'btn-primary', action: () => this.addMaterial() }
        ]);
        modal.show();
    }

    async addMaterial() {
        const data = {
            name: document.getElementById('material-name').value,
            E: parseFloat(document.getElementById('material-E').value),
            G: parseFloat(document.getElementById('material-G').value),
            nu: parseFloat(document.getElementById('material-nu').value),
            rho: parseFloat(document.getElementById('material-rho').value)
        };

        if (!data.name || isNaN(data.E) || isNaN(data.G) || isNaN(data.nu) || isNaN(data.rho)) {
            this.showError('Please fill in all material properties');
            return;
        }

        try {
            const response = await fetch(`/api/models/${this.currentModelId}/materials`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(data)
            });

            const result = await response.json();
            if (result.success) {
                this.updateMaterialsList();
                this.updateStatus(result.message);
                Modal.close();
            } else {
                this.showError(result.error);
            }
        } catch (error) {
            this.showError('Failed to add material: ' + error.message);
        }
    }

    async showAddSectionModal() {
        if (!this.currentModelId) {
            this.showError('Please create a model first');
            return;
        }

        const modal = new Modal('Add Section', `
            <div class="form-group">
                <label for="section-name">Section Name</label>
                <input type="text" id="section-name" class="form-control" placeholder="e.g., W12x26, Custom">
            </div>
            <div class="form-row">
                <div class="form-group">
                    <label for="section-A">Area (A)</label>
                    <input type="number" id="section-A" class="form-control" placeholder="7.65" step="any">
                </div>
                <div class="form-group">
                    <label for="section-Iy">Moment of Inertia Iy</label>
                    <input type="number" id="section-Iy" class="form-control" placeholder="17.3" step="any">
                </div>
            </div>
            <div class="form-row">
                <div class="form-group">
                    <label for="section-Iz">Moment of Inertia Iz</label>
                    <input type="number" id="section-Iz" class="form-control" placeholder="204" step="any">
                </div>
                <div class="form-group">
                    <label for="section-J">Torsional Constant (J)</label>
                    <input type="number" id="section-J" class="form-control" placeholder="0.3" step="any">
                </div>
            </div>
        `, [
            { text: 'Cancel', class: 'btn-outline' },
            { text: 'Add Section', class: 'btn-primary', action: () => this.addSection() }
        ]);
        modal.show();
    }

    async addSection() {
        const data = {
            name: document.getElementById('section-name').value,
            A: parseFloat(document.getElementById('section-A').value),
            Iy: parseFloat(document.getElementById('section-Iy').value),
            Iz: parseFloat(document.getElementById('section-Iz').value),
            J: parseFloat(document.getElementById('section-J').value)
        };

        if (!data.name || isNaN(data.A) || isNaN(data.Iy) || isNaN(data.Iz) || isNaN(data.J)) {
            this.showError('Please fill in all section properties');
            return;
        }

        try {
            const response = await fetch(`/api/models/${this.currentModelId}/sections`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(data)
            });

            const result = await response.json();
            if (result.success) {
                this.updateSectionsList();
                this.updateStatus(result.message);
                Modal.close();
            } else {
                this.showError(result.error);
            }
        } catch (error) {
            this.showError('Failed to add section: ' + error.message);
        }
    }

    async showAddNodeModal() {
        if (!this.currentModelId) {
            this.showError('Please create a model first');
            return;
        }

        const modal = new Modal('Add Node', `
            <div class="form-group">
                <label for="node-name">Node Name</label>
                <input type="text" id="node-name" class="form-control" placeholder="e.g., N1, Node1">
            </div>
            <div class="form-row-3">
                <div class="form-group">
                    <label for="node-x">X Coordinate</label>
                    <input type="number" id="node-x" class="form-control" placeholder="0" step="any">
                </div>
                <div class="form-group">
                    <label for="node-y">Y Coordinate</label>
                    <input type="number" id="node-y" class="form-control" placeholder="0" step="any">
                </div>
                <div class="form-group">
                    <label for="node-z">Z Coordinate</label>
                    <input type="number" id="node-z" class="form-control" placeholder="0" step="any">
                </div>
            </div>
        `, [
            { text: 'Cancel', class: 'btn-outline' },
            { text: 'Add Node', class: 'btn-primary', action: () => this.addNode() }
        ]);
        modal.show();
    }

    async addNode() {
        const data = {
            name: document.getElementById('node-name').value,
            x: parseFloat(document.getElementById('node-x').value) || 0,
            y: parseFloat(document.getElementById('node-y').value) || 0,
            z: parseFloat(document.getElementById('node-z').value) || 0
        };

        if (!data.name) {
            this.showError('Please enter a node name');
            return;
        }

        try {
            const response = await fetch(`/api/models/${this.currentModelId}/nodes`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(data)
            });

            const result = await response.json();
            if (result.success) {
                this.updateModelTree();
                this.updateVisualization();
                this.updateStatus(result.message);
                Modal.close();
            } else {
                this.showError(result.error);
            }
        } catch (error) {
            this.showError('Failed to add node: ' + error.message);
        }
    }

    async analyzeModel() {
        if (!this.currentModelId) {
            this.showError('No model to analyze');
            return;
        }

        try {
            this.updateStatus('Analyzing model...');
            document.getElementById('analyze-btn').disabled = true;

            const response = await fetch(`/api/models/${this.currentModelId}/analyze`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ type: 'linear' })
            });

            const result = await response.json();
            if (result.success) {
                this.updateStatus('Analysis completed');
                this.showResults();
                this.updateVisualization();
            } else {
                this.showError('Analysis failed: ' + result.error);
            }
        } catch (error) {
            this.showError('Analysis failed: ' + error.message);
        } finally {
            document.getElementById('analyze-btn').disabled = false;
        }
    }

    async showResults() {
        const resultsPanel = document.getElementById('results-panel');
        resultsPanel.classList.add('active');
        
        try {
            const response = await fetch(`/api/models/${this.currentModelId}/results/${this.currentCombo}`);
            const result = await response.json();
            
            if (result.success) {
                this.displayResults(result.results);
            } else {
                this.showError('Failed to get results: ' + result.error);
            }
        } catch (error) {
            this.showError('Failed to get results: ' + error.message);
        }
    }

    displayResults(results) {
        // Show node results by default
        this.showResultsTab('nodes', results);
    }

    showResultsTab(tab, results = null) {
        const content = document.getElementById('results-content');
        
        if (!results && this.currentModelId) {
            // Fetch results if not provided
            this.getResults().then(results => {
                if (results) this.showResultsTab(tab, results);
            });
            return;
        }

        switch (tab) {
            case 'nodes':
                content.innerHTML = this.generateNodeResultsTable(results.nodes);
                break;
            case 'members':
                content.innerHTML = this.generateMemberResultsTable(results.members);
                break;
            case 'diagrams':
                content.innerHTML = this.generateDiagramsView();
                break;
        }
    }

    generateNodeResultsTable(nodeResults) {
        let html = '<h4>Node Displacements & Reactions</h4>';
        html += '<table class="results-table">';
        html += '<thead><tr><th>Node</th><th>DX</th><th>DY</th><th>DZ</th><th>RX</th><th>RY</th><th>RZ</th></tr></thead>';
        html += '<tbody>';
        
        for (const [nodeName, data] of Object.entries(nodeResults)) {
            const disp = data.displacements;
            html += `<tr>
                <td>${nodeName}</td>
                <td>${disp.DX.toFixed(6)}</td>
                <td>${disp.DY.toFixed(6)}</td>
                <td>${disp.DZ.toFixed(6)}</td>
                <td>${disp.RX.toFixed(6)}</td>
                <td>${disp.RY.toFixed(6)}</td>
                <td>${disp.RZ.toFixed(6)}</td>
            </tr>`;
        }
        
        html += '</tbody></table>';
        
        html += '<h4 class="mt-2">Reactions</h4>';
        html += '<table class="results-table">';
        html += '<thead><tr><th>Node</th><th>FX</th><th>FY</th><th>FZ</th><th>MX</th><th>MY</th><th>MZ</th></tr></thead>';
        html += '<tbody>';
        
        for (const [nodeName, data] of Object.entries(nodeResults)) {
            const rxn = data.reactions;
            html += `<tr>
                <td>${nodeName}</td>
                <td>${rxn.FX.toFixed(3)}</td>
                <td>${rxn.FY.toFixed(3)}</td>
                <td>${rxn.FZ.toFixed(3)}</td>
                <td>${rxn.MX.toFixed(3)}</td>
                <td>${rxn.MY.toFixed(3)}</td>
                <td>${rxn.MZ.toFixed(3)}</td>
            </tr>`;
        }
        
        html += '</tbody></table>';
        return html;
    }

    generateMemberResultsTable(memberResults) {
        let html = '<h4>Member Forces & Moments</h4>';
        html += '<table class="results-table">';
        html += '<thead><tr><th>Member</th><th>Max Mz</th><th>Min Mz</th><th>Max Fy</th><th>Min Fy</th><th>Max Axial</th><th>Min Axial</th></tr></thead>';
        html += '<tbody>';
        
        for (const [memberName, data] of Object.entries(memberResults)) {
            if (data.error) {
                html += `<tr><td>${memberName}</td><td colspan="6">Error: ${data.error}</td></tr>`;
            } else {
                html += `<tr>
                    <td>${memberName}</td>
                    <td>${data.max_moment_Mz.toFixed(3)}</td>
                    <td>${data.min_moment_Mz.toFixed(3)}</td>
                    <td>${data.max_shear_Fy.toFixed(3)}</td>
                    <td>${data.min_shear_Fy.toFixed(3)}</td>
                    <td>${data.max_axial.toFixed(3)}</td>
                    <td>${data.min_axial.toFixed(3)}</td>
                </tr>`;
            }
        }
        
        html += '</tbody></table>';
        return html;
    }

    generateDiagramsView() {
        return `
            <h4>Member Diagrams</h4>
            <p>Select a member from the model tree to view force and moment diagrams.</p>
            <div id="diagram-content">
                <!-- Diagrams will be loaded here -->
            </div>
        `;
    }

    async getResults() {
        try {
            const response = await fetch(`/api/models/${this.currentModelId}/results/${this.currentCombo}`);
            const result = await response.json();
            return result.success ? result.results : null;
        } catch (error) {
            console.error('Failed to get results:', error);
            return null;
        }
    }

    hideResults() {
        document.getElementById('results-panel').classList.remove('active');
    }

    async updateVisualization() {
        if (!this.currentModelId) return;

        try {
            const response = await fetch(`/api/models/${this.currentModelId}/visualization?combo=${this.currentCombo}`);
            const result = await response.json();
            
            if (result.success && this.viewer) {
                this.viewer.updateModel(result.visualization);
            }
        } catch (error) {
            console.error('Failed to update visualization:', error);
        }
    }

    async updateModelTree() {
        if (!this.currentModelId) return;

        try {
            const response = await fetch(`/api/models/${this.currentModelId}`);
            const result = await response.json();
            
            if (result.success) {
                this.displayModelTree(result.model);
            }
        } catch (error) {
            console.error('Failed to update model tree:', error);
        }
    }

    displayModelTree(modelData) {
        const treeContainer = document.getElementById('model-tree');
        let html = `<div class="model-info"><strong>${modelData.metadata.name}</strong></div>`;
        
        // Nodes
        if (modelData.nodes.length > 0) {
            html += '<div class="tree-section"><h5>Nodes</h5>';
            modelData.nodes.forEach(node => {
                const supports = Object.values(node.supports).some(s => s) ? ' 🔒' : '';
                html += `<div class="node-item" data-node="${node.name}">
                    ${node.name} (${node.x.toFixed(1)}, ${node.y.toFixed(1)}, ${node.z.toFixed(1)})${supports}
                </div>`;
            });
            html += '</div>';
        }
        
        // Members
        if (modelData.members.length > 0) {
            html += '<div class="tree-section"><h5>Members</h5>';
            modelData.members.forEach(member => {
                html += `<div class="member-item" data-member="${member.name}">
                    ${member.name} (${member.i_node} → ${member.j_node})
                </div>`;
            });
            html += '</div>';
        }
        
        treeContainer.innerHTML = html;
    }

    updateMaterialsList() {
        // Implementation for updating materials list
    }

    updateSectionsList() {
        // Implementation for updating sections list
    }

    updateModelInfo() {
        const info = document.getElementById('model-info');
        if (this.currentModelId) {
            info.textContent = `Model: ${this.currentModelId.substring(0, 8)}...`;
        } else {
            info.textContent = '';
        }
    }

    enableAnalyzeButton() {
        document.getElementById('analyze-btn').disabled = false;
    }

    updateStatus(message) {
        document.getElementById('status-text').textContent = message;
    }

    showError(message) {
        console.error(message);
        this.updateStatus('Error: ' + message);
        // You could also show a toast notification here
    }

    // Placeholder methods for other modals
    showAddMemberModal() { this.showError('Add Member modal not implemented yet'); }
    showAddNodeLoadModal() { this.showError('Add Node Load modal not implemented yet'); }
    showAddMemberLoadModal() { this.showError('Add Member Load modal not implemented yet'); }
    showAddSupportModal() { this.showError('Add Support modal not implemented yet'); }
    showAddLoadComboModal() { this.showError('Add Load Combo modal not implemented yet'); }
}

// Initialize the application when the page loads
document.addEventListener('DOMContentLoaded', () => {
    window.app = new PyniteApp();
});