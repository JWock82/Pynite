/**
 * Modal management for the Pynite FEA Web App
 */

class Modal {
    constructor(title, content, actions = []) {
        this.title = title;
        this.content = content;
        this.actions = actions;
        this.overlay = document.getElementById('modal-overlay');
        this.contentContainer = document.getElementById('modal-content');
    }

    show() {
        this.render();
        this.overlay.classList.add('active');
        
        // Focus first input if available
        const firstInput = this.contentContainer.querySelector('input, select, textarea');
        if (firstInput) {
            setTimeout(() => firstInput.focus(), 100);
        }

        // Handle escape key
        this.escapeHandler = (e) => {
            if (e.key === 'Escape') {
                this.close();
            }
        };
        document.addEventListener('keydown', this.escapeHandler);

        // Handle overlay click
        this.overlayClickHandler = (e) => {
            if (e.target === this.overlay) {
                this.close();
            }
        };
        this.overlay.addEventListener('click', this.overlayClickHandler);
    }

    render() {
        let html = `
            <div class="modal-header">
                <h3>${this.title}</h3>
            </div>
            <div class="modal-body">
                ${this.content}
            </div>
        `;

        if (this.actions.length > 0) {
            html += '<div class="modal-actions">';
            this.actions.forEach((action, index) => {
                const btnClass = action.class || 'btn-outline';
                html += `<button class="btn ${btnClass}" data-action="${index}">${action.text}</button>`;
            });
            html += '</div>';
        }

        this.contentContainer.innerHTML = html;

        // Attach action handlers
        this.contentContainer.querySelectorAll('[data-action]').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const actionIndex = parseInt(e.target.dataset.action);
                const action = this.actions[actionIndex];
                if (action.action) {
                    action.action();
                } else {
                    this.close();
                }
            });
        });
    }

    close() {
        this.overlay.classList.remove('active');
        document.removeEventListener('keydown', this.escapeHandler);
        this.overlay.removeEventListener('click', this.overlayClickHandler);
    }

    static close() {
        const overlay = document.getElementById('modal-overlay');
        overlay.classList.remove('active');
    }
}

// Specific modal implementations
class AddNodeModal extends Modal {
    constructor(onSubmit) {
        const content = `
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
        `;

        const actions = [
            { text: 'Cancel', class: 'btn-outline' },
            { text: 'Add Node', class: 'btn-primary', action: onSubmit }
        ];

        super('Add Node', content, actions);
    }

    getData() {
        return {
            name: document.getElementById('node-name').value,
            x: parseFloat(document.getElementById('node-x').value) || 0,
            y: parseFloat(document.getElementById('node-y').value) || 0,
            z: parseFloat(document.getElementById('node-z').value) || 0
        };
    }
}

class AddMemberModal extends Modal {
    constructor(nodes, materials, sections, onSubmit) {
        const nodeOptions = nodes.map(node => `<option value="${node.name}">${node.name}</option>`).join('');
        const materialOptions = materials.map(mat => `<option value="${mat.name}">${mat.name}</option>`).join('');
        const sectionOptions = sections.map(sec => `<option value="${sec.name}">${sec.name}</option>`).join('');

        const content = `
            <div class="form-group">
                <label for="member-name">Member Name</label>
                <input type="text" id="member-name" class="form-control" placeholder="e.g., M1, Beam1">
            </div>
            <div class="form-row">
                <div class="form-group">
                    <label for="member-i-node">I-Node</label>
                    <select id="member-i-node" class="form-control">
                        <option value="">Select node...</option>
                        ${nodeOptions}
                    </select>
                </div>
                <div class="form-group">
                    <label for="member-j-node">J-Node</label>
                    <select id="member-j-node" class="form-control">
                        <option value="">Select node...</option>
                        ${nodeOptions}
                    </select>
                </div>
            </div>
            <div class="form-row">
                <div class="form-group">
                    <label for="member-material">Material</label>
                    <select id="member-material" class="form-control">
                        <option value="">Select material...</option>
                        ${materialOptions}
                    </select>
                </div>
                <div class="form-group">
                    <label for="member-section">Section</label>
                    <select id="member-section" class="form-control">
                        <option value="">Select section...</option>
                        ${sectionOptions}
                    </select>
                </div>
            </div>
        `;

        const actions = [
            { text: 'Cancel', class: 'btn-outline' },
            { text: 'Add Member', class: 'btn-primary', action: onSubmit }
        ];

        super('Add Member', content, actions);
    }

    getData() {
        return {
            name: document.getElementById('member-name').value,
            i_node: document.getElementById('member-i-node').value,
            j_node: document.getElementById('member-j-node').value,
            material: document.getElementById('member-material').value,
            section: document.getElementById('member-section').value
        };
    }
}

class AddSupportModal extends Modal {
    constructor(nodes, onSubmit) {
        const nodeOptions = nodes.map(node => `<option value="${node.name}">${node.name}</option>`).join('');

        const content = `
            <div class="form-group">
                <label for="support-node">Node</label>
                <select id="support-node" class="form-control">
                    <option value="">Select node...</option>
                    ${nodeOptions}
                </select>
            </div>
            <div class="form-group">
                <label>Translation Supports</label>
                <div class="form-row-3">
                    <label><input type="checkbox" id="support-dx"> DX</label>
                    <label><input type="checkbox" id="support-dy"> DY</label>
                    <label><input type="checkbox" id="support-dz"> DZ</label>
                </div>
            </div>
            <div class="form-group">
                <label>Rotation Supports</label>
                <div class="form-row-3">
                    <label><input type="checkbox" id="support-rx"> RX</label>
                    <label><input type="checkbox" id="support-ry"> RY</label>
                    <label><input type="checkbox" id="support-rz"> RZ</label>
                </div>
            </div>
            <div class="form-group">
                <label>Quick Support Types</label>
                <div class="form-row-3">
                    <button type="button" class="btn btn-sm btn-outline" onclick="this.setPinnedSupport()">Pinned</button>
                    <button type="button" class="btn btn-sm btn-outline" onclick="this.setFixedSupport()">Fixed</button>
                    <button type="button" class="btn btn-sm btn-outline" onclick="this.clearSupports()">Clear</button>
                </div>
            </div>
        `;

        const actions = [
            { text: 'Cancel', class: 'btn-outline' },
            { text: 'Add Support', class: 'btn-primary', action: onSubmit }
        ];

        super('Add Support', content, actions);
    }

    getData() {
        return {
            node: document.getElementById('support-node').value,
            support_DX: document.getElementById('support-dx').checked,
            support_DY: document.getElementById('support-dy').checked,
            support_DZ: document.getElementById('support-dz').checked,
            support_RX: document.getElementById('support-rx').checked,
            support_RY: document.getElementById('support-ry').checked,
            support_RZ: document.getElementById('support-rz').checked
        };
    }

    setPinnedSupport() {
        document.getElementById('support-dx').checked = true;
        document.getElementById('support-dy').checked = true;
        document.getElementById('support-dz').checked = true;
        document.getElementById('support-rx').checked = false;
        document.getElementById('support-ry').checked = false;
        document.getElementById('support-rz').checked = false;
    }

    setFixedSupport() {
        document.getElementById('support-dx').checked = true;
        document.getElementById('support-dy').checked = true;
        document.getElementById('support-dz').checked = true;
        document.getElementById('support-rx').checked = true;
        document.getElementById('support-ry').checked = true;
        document.getElementById('support-rz').checked = true;
    }

    clearSupports() {
        document.getElementById('support-dx').checked = false;
        document.getElementById('support-dy').checked = false;
        document.getElementById('support-dz').checked = false;
        document.getElementById('support-rx').checked = false;
        document.getElementById('support-ry').checked = false;
        document.getElementById('support-rz').checked = false;
    }
}

// Global functions for quick support setting (called from modal buttons)
window.setPinnedSupport = function() {
    document.getElementById('support-dx').checked = true;
    document.getElementById('support-dy').checked = true;
    document.getElementById('support-dz').checked = true;
    document.getElementById('support-rx').checked = false;
    document.getElementById('support-ry').checked = false;
    document.getElementById('support-rz').checked = false;
};

window.setFixedSupport = function() {
    document.getElementById('support-dx').checked = true;
    document.getElementById('support-dy').checked = true;
    document.getElementById('support-dz').checked = true;
    document.getElementById('support-rx').checked = true;
    document.getElementById('support-ry').checked = true;
    document.getElementById('support-rz').checked = true;
};

window.clearSupports = function() {
    document.getElementById('support-dx').checked = false;
    document.getElementById('support-dy').checked = false;
    document.getElementById('support-dz').checked = false;
    document.getElementById('support-rx').checked = false;
    document.getElementById('support-ry').checked = false;
    document.getElementById('support-rz').checked = false;
};