/**
 * Three.js viewer for 3D visualization of finite element models
 */

class ThreeViewer {
    constructor(canvasId) {
        this.canvas = document.getElementById(canvasId);
        this.scene = null;
        this.camera = null;
        this.renderer = null;
        this.controls = null;
        this.modelData = null;
        this.nodeObjects = [];
        this.memberObjects = [];
        
        this.init();
    }

    init() {
        this.setupScene();
        this.setupCamera();
        this.setupRenderer();
        this.setupControls();
        this.setupLighting();
        this.addGrid();
        this.animate();
        
        // Handle window resize
        window.addEventListener('resize', () => this.onWindowResize());
    }

    setupScene() {
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0xf8fafc);
    }

    setupCamera() {
        const container = this.canvas.parentElement;
        const aspect = container.clientWidth / container.clientHeight;
        this.camera = new THREE.PerspectiveCamera(75, aspect, 0.1, 1000);
        this.camera.position.set(10, 10, 10);
        this.camera.lookAt(0, 0, 0);
    }

    setupRenderer() {
        this.renderer = new THREE.WebGLRenderer({ 
            canvas: this.canvas,
            antialias: true 
        });
        this.renderer.setSize(this.canvas.parentElement.clientWidth, this.canvas.parentElement.clientHeight);
        this.renderer.shadowMap.enabled = true;
        this.renderer.shadowMap.type = THREE.PCFSoftShadowMap;
    }

    setupControls() {
        this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableDamping = true;
        this.controls.dampingFactor = 0.05;
        this.controls.screenSpacePanning = false;
        this.controls.minDistance = 1;
        this.controls.maxDistance = 100;
    }

    setupLighting() {
        // Ambient light
        const ambientLight = new THREE.AmbientLight(0x404040, 0.6);
        this.scene.add(ambientLight);

        // Directional light
        const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
        directionalLight.position.set(10, 10, 5);
        directionalLight.castShadow = true;
        directionalLight.shadow.mapSize.width = 2048;
        directionalLight.shadow.mapSize.height = 2048;
        this.scene.add(directionalLight);

        // Additional fill light
        const fillLight = new THREE.DirectionalLight(0xffffff, 0.3);
        fillLight.position.set(-10, -10, -5);
        this.scene.add(fillLight);
    }

    addGrid() {
        const gridHelper = new THREE.GridHelper(20, 20, 0x888888, 0xcccccc);
        gridHelper.rotateX(Math.PI / 2); // Rotate to XY plane
        this.scene.add(gridHelper);

        // Add axis helper
        const axesHelper = new THREE.AxesHelper(5);
        this.scene.add(axesHelper);
    }

    updateModel(modelData) {
        this.modelData = modelData;
        this.clearModel();
        this.renderNodes();
        this.renderMembers();
        this.fitCameraToModel();
    }

    clearModel() {
        // Remove existing model objects
        this.nodeObjects.forEach(obj => this.scene.remove(obj));
        this.memberObjects.forEach(obj => this.scene.remove(obj));
        this.nodeObjects = [];
        this.memberObjects = [];
    }

    renderNodes() {
        if (!this.modelData || !this.modelData.nodes) return;

        const nodeGeometry = new THREE.SphereGeometry(0.2, 16, 16);
        
        this.modelData.nodes.forEach(node => {
            // Choose color based on support conditions
            let color = 0x3b82f6; // Blue for free nodes
            if (Object.values(node.supports).some(s => s)) {
                color = 0xef4444; // Red for supported nodes
            }

            const nodeMaterial = new THREE.MeshLambertMaterial({ color });
            const nodeMesh = new THREE.Mesh(nodeGeometry, nodeMaterial);
            
            // Position the node
            nodeMesh.position.set(node.x, node.y, node.z);
            
            // Add displacement if available
            if (document.getElementById('show-deformed').checked && node.displacement) {
                const scale = 10; // Displacement scale factor
                nodeMesh.position.x += node.displacement.x * scale;
                nodeMesh.position.y += node.displacement.y * scale;
                nodeMesh.position.z += node.displacement.z * scale;
            }

            nodeMesh.userData = { type: 'node', name: node.name, data: node };
            this.scene.add(nodeMesh);
            this.nodeObjects.push(nodeMesh);

            // Add support symbols
            if (Object.values(node.supports).some(s => s)) {
                this.addSupportSymbol(node);
            }

            // Add node label
            this.addNodeLabel(node, nodeMesh.position);
        });
    }

    renderMembers() {
        if (!this.modelData || !this.modelData.members) return;

        this.modelData.members.forEach(member => {
            // Find the nodes
            const iNode = this.modelData.nodes.find(n => n.name === member.i_node);
            const jNode = this.modelData.nodes.find(n => n.name === member.j_node);
            
            if (!iNode || !jNode) return;

            // Create member geometry
            const points = [];
            let iPos = new THREE.Vector3(iNode.x, iNode.y, iNode.z);
            let jPos = new THREE.Vector3(jNode.x, jNode.y, jNode.z);

            // Apply displacements if showing deformed shape
            if (document.getElementById('show-deformed').checked) {
                const scale = 10;
                if (iNode.displacement) {
                    iPos.add(new THREE.Vector3(
                        iNode.displacement.x * scale,
                        iNode.displacement.y * scale,
                        iNode.displacement.z * scale
                    ));
                }
                if (jNode.displacement) {
                    jPos.add(new THREE.Vector3(
                        jNode.displacement.x * scale,
                        jNode.displacement.y * scale,
                        jNode.displacement.z * scale
                    ));
                }
            }

            points.push(iPos, jPos);
            
            const geometry = new THREE.BufferGeometry().setFromPoints(points);
            const material = new THREE.LineBasicMaterial({ 
                color: 0x64748b,
                linewidth: 3
            });
            
            const line = new THREE.Line(geometry, material);
            line.userData = { type: 'member', name: member.name, data: member };
            
            this.scene.add(line);
            this.memberObjects.push(line);

            // Add member label at midpoint
            const midpoint = new THREE.Vector3().addVectors(iPos, jPos).multiplyScalar(0.5);
            this.addMemberLabel(member, midpoint);
        });
    }

    addSupportSymbol(node) {
        // Create a simple support symbol (triangle for pinned, square for fixed)
        const hasTranslationSupport = node.supports.DX || node.supports.DY || node.supports.DZ;
        const hasRotationSupport = node.supports.RX || node.supports.RY || node.supports.RZ;

        if (hasTranslationSupport) {
            let geometry;
            if (hasRotationSupport) {
                // Fixed support - square
                geometry = new THREE.BoxGeometry(0.6, 0.6, 0.6);
            } else {
                // Pinned support - triangle
                geometry = new THREE.ConeGeometry(0.3, 0.6, 3);
            }

            const material = new THREE.MeshLambertMaterial({ color: 0x22c55e });
            const supportMesh = new THREE.Mesh(geometry, material);
            supportMesh.position.set(node.x, node.y - 0.5, node.z);
            
            this.scene.add(supportMesh);
            this.nodeObjects.push(supportMesh);
        }
    }

    addNodeLabel(node, position) {
        // Create text sprite for node label
        const canvas = document.createElement('canvas');
        const context = canvas.getContext('2d');
        canvas.width = 128;
        canvas.height = 64;
        
        context.fillStyle = 'rgba(255, 255, 255, 0.9)';
        context.fillRect(0, 0, canvas.width, canvas.height);
        context.fillStyle = 'black';
        context.font = '16px Arial';
        context.textAlign = 'center';
        context.fillText(node.name, canvas.width / 2, canvas.height / 2 + 6);
        
        const texture = new THREE.CanvasTexture(canvas);
        const spriteMaterial = new THREE.SpriteMaterial({ map: texture });
        const sprite = new THREE.Sprite(spriteMaterial);
        sprite.position.copy(position);
        sprite.position.y += 0.5;
        sprite.scale.set(1, 0.5, 1);
        
        this.scene.add(sprite);
        this.nodeObjects.push(sprite);
    }

    addMemberLabel(member, position) {
        // Create text sprite for member label
        const canvas = document.createElement('canvas');
        const context = canvas.getContext('2d');
        canvas.width = 128;
        canvas.height = 64;
        
        context.fillStyle = 'rgba(255, 255, 255, 0.8)';
        context.fillRect(0, 0, canvas.width, canvas.height);
        context.fillStyle = '#64748b';
        context.font = '14px Arial';
        context.textAlign = 'center';
        context.fillText(member.name, canvas.width / 2, canvas.height / 2 + 5);
        
        const texture = new THREE.CanvasTexture(canvas);
        const spriteMaterial = new THREE.SpriteMaterial({ map: texture });
        const sprite = new THREE.Sprite(spriteMaterial);
        sprite.position.copy(position);
        sprite.scale.set(0.8, 0.4, 1);
        
        this.scene.add(sprite);
        this.memberObjects.push(sprite);
    }

    fitCameraToModel() {
        if (!this.modelData || !this.modelData.nodes || this.modelData.nodes.length === 0) {
            return;
        }

        // Calculate bounding box
        const box = new THREE.Box3();
        this.modelData.nodes.forEach(node => {
            box.expandByPoint(new THREE.Vector3(node.x, node.y, node.z));
        });

        // Position camera to view the entire model
        const center = box.getCenter(new THREE.Vector3());
        const size = box.getSize(new THREE.Vector3());
        const maxDim = Math.max(size.x, size.y, size.z);
        
        const distance = maxDim * 2;
        this.camera.position.set(
            center.x + distance,
            center.y + distance,
            center.z + distance
        );
        this.camera.lookAt(center);
        this.controls.target.copy(center);
        this.controls.update();
    }

    animate() {
        requestAnimationFrame(() => this.animate());
        this.controls.update();
        this.renderer.render(this.scene, this.camera);
    }

    onWindowResize() {
        const container = this.canvas.parentElement;
        const width = container.clientWidth;
        const height = container.clientHeight;

        this.camera.aspect = width / height;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(width, height);
    }
}