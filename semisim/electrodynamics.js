// TODO
// Make graphics faster

class Electrodynamics {
    constructor() {

        // Initialize UI
        // Member declarations (unchanged, included for completeness)
        this.Ex = null;
        this.Ey = null;
        this.Hz = null;
        this.Hz_laplacian = null;
        this.rho_n = null;
        this.rho_p = null;
        this.rho_back = null;
        this.rho_abs = null;
        this.rho_free = null;
        this.mobility_factor = null;
        this.Jx_n = null;
        this.Jy_n = null;
        this.Jx_p = null;
        this.Jy_p = null;
        this.Jx_abs = null;
        this.Jy_abs = null;
        this.Jx_free = null;
        this.Jy_free = null;
        this.Dx = null;
        this.Dy = null;
        this.Bz = null;
        this.Sx = null;
        this.Sy = null;
        this.u = null;
        this.phi = null;
        this.G = null;
        this.F_n = null;
        this.F_p = null;
        this.grad_E0x_n = null;
        this.grad_E0y_n = null;
        this.grad_E0x_p = null;
        this.grad_E0y_p = null;
        this.grad_Fx_n = null;
        this.grad_Fy_n = null;
        this.grad_Fx_p = null;
        this.grad_Fy_p = null;
        this.Q = null;
        this.S = null;
        this.F = null;
        this.updateMiscFields = false;
        this.debug = null;
        this.materials = null;
        this.selection = null;
        this.clipboard = null;
        this.F0_n = null;
        this.F0_p = null;
        this.E0_n = null;
        this.E0_p = null;
        this.K = null;
        this.E_a = null;
        this.R = null;
        this.cmfx_n = null;
        this.cmfy_n = null;
        this.cmfx_p = null;
        this.cmfy_p = null;
        this.emfx = null;
        this.emfy = null;
        this.epsx = null;
        this.epsy = null;
        this.mu_z = null;
        this.conducting = null;
        this.conducting_x = null;
        this.conducting_y = null;
        this.absorptivity = null;
        this.absorptivity_x = null;
        this.absorptivity_y = null;
        this.MG_rho0 = null;
        this.MG_rho = null;
        this.MG_epsx = null;
        this.MG_epsy = null;
        this.MG_eps_avg = null;
        this.MG_phi1 = null;
        this.MG_phi2 = null;
        this.distance = null;
        this.visited = null;
        this.eps0 = 8.85e-12;
        this.mu0 = 1.257e-6;
        this.c = 1 / Math.sqrt(this.eps0 * this.mu0);
        this.c_squared = 1 / (this.eps0 * this.mu0);
        this.kT = 4.11e-21;
        this.beta = 1 / this.kT;
        this.e_charge = 1.6e-19;
        this.m_electron = 9e-31;
        this.mu_electron = 0.1400 * 1500;
        this.D_electron = this.mu_electron / (this.beta * this.e_charge);
        this.mu_hole = 0.5 * this.mu_electron;
        this.D_hole = this.mu_hole / (this.beta * this.e_charge);
        this.eVtoJ = this.e_charge;
        this.ni_semi = 1e16;
        this.W_semi = 4.7 * this.eVtoJ;
        this.E_b_semi = 1.12 * this.eVtoJ;
        this.ni_metal = 5e20;
        this.W_metal_default = this.W_semi;
        this.W_metal_high = this.W_semi + 0.3 * this.eVtoJ;
        this.W_metal_low = this.W_semi - 0.3 * this.eVtoJ;
        this.E_b_metal = this.E_b_semi;
        this.recombination_rate_semi = 1e7 / this.ni_semi;
        this.recombination_rate_metal = 1e8 / this.ni_semi;
        this.K_semi = this.ni_semi * this.ni_semi;
        this.recombination_cross_section = 1e-10;
        this.arrhenius_prefactor = this.recombination_cross_section * Math.sqrt(8 * this.kT / (Math.PI * this.m_electron / 2));
        this.E_a_semi = -Math.log(this.recombination_rate_semi / this.arrhenius_prefactor);
        this.E_a_metal = -Math.log(this.recombination_rate_metal / this.arrhenius_prefactor);
        this.q_n = -this.e_charge;
        this.q_p = this.e_charge;
        this.n_default_doping_concentration = 5e19;
        this.p_default_doping_concentration = 5e19;
        this.n_light_doping_concentration = 1e19;
        this.p_light_doping_concentration = 1e19;
        this.n_heavy_doping_concentration = 2.5e20;
        this.p_heavy_doping_concentration = 2.5e20;
        this.dielectric_eps_r = 5.0;
        this.ferromagnet_mu_r = 5.0;
        this.staticcharge_density = 10.0;
        this.E_sat = 5e5;
        this.junction_size = 3;
        this.voltageprobes = [];
        this.currentprobes = [];
        this.ground = null;
        this.default_width = 2.56e-5;
        this.width = 0;
        this.ds = 0;
        this.dt = 0;
        this.depth = 1e-3;
        this.default_resolution = 256; //Must be a power of 2
        this.WASM_array_padding = 1;
        this.nx = 0;
        this.ny = 0;
        this.absorber_width = 0;
        this.absorbing_coeff = 0;
        this.Hz_dissipation = 0;
        this.dt_maximum = 0;
        this.parity = -1;
        this.time = 0.0;
        this.stepnumber = 0;
        this.frame = 0;
        this.lastsimspeed = 0;
        this.error_detection_threshold = 1e-5;
        this.sign_violation = false;
        this.r = null;
        this.opts = null;
        this.help = null;
        this.screen = null;
        this.image_r = null;
        this.image_g = null;
        this.image_b = null;
        this.col_r = 0;
        this.col_g = 0;
        this.col_b = 0;
        this.alphaBG = 0;
        this.alphaFG = 0;
        this.rand = Math.random;
        this.bigfont = '15px sans-serif';
        this.regularfont = '12px sans-serif';
        this.scalefactor = 0;
        this.imgwidth = 0;
        this.imgheight = 0;
        this.targetframerate = 60;
        this.frameduration = 1000 / this.targetframerate;
        this.t5 = null;
        this.advanceframe = false;
        this.clear = false;
        this.reset = false;
        this.save = false;
        this.load = false;
        this.debugging = false;
        this.cut = false;
        this.copy = false;
        this.paste = false;
        this.delete = false;
        this.shift_down = false;
        this.ctrl_down = false;
        this.alt_down = false;
        this.pointerinfo = null;
        this.mouse_pressed = false;
        this.mouse_pressed_prev = false;
        this.modifier_pressed = false;
        this.moving_selection = false;
        this.dragging_selection = false;
        this.brush_changed = false;
        this.mousebutton = 0;
        this.mx = 0;
        this.my = 0;
        this.mx_start = 0;
        this.my_start = 0;
        this.mx_index = 0;
        this.my_index = 0;
        this.mx_start_index = 0;
        this.my_start_index = 0;
        this.mx_realspace = 0;
        this.my_realspace = 0;
        this.mxp_realspace = 0;
        this.myp_realspace = 0;
        this.mx_start_realspace = 0;
        this.my_start_realspace = 0;
        this.delta_mx_index = 0;
        this.delta_my_index = 0;
        this.EMF_selected = false;
        this.max_EMF = 5e5;
        this.prev_brush = Brush.DRAW;
        this.brushsize = 0;
        this.under_brush = null;
        this.selected = null;
        this.selected_EMF = null;
        this.text_x = 0;
        this.text_y = 0;
        this.texting = false;
        this.HAND_CURSOR = 'pointer';
        this.DEFAULT_CURSOR = 'default';
        this.infile = null;
        this.outfile = null;
        this.saveversion = 1;
        this.fileextension = '.semisim';

        // Initialize timers
        this.t4 = new Timer('Poisson constraint solver');
        this.t5 = new Timer('Graphics');
        this.t6 = new Timer('Iterate simulation');
        this.t7 = new Timer('Poisson potential solver');
        this.t8 = new Timer('Calc misc fields');
        this.t9 = new Timer('Debug');
    
        // Initialize variables
        this.prev_EMF_setting = 50;
        this.prev_boundary = BoundaryCondition.DISSIPATIVE;
        this.alphaFG = 1.0;
        this.alphaBG = 0.0;
        this.col_r = 1.0;
        this.col_g = 1.0;
        this.col_b = 1.0;
        this.rand = {
            seed: 4,
            nextFloat: function () {
                this.seed = (this.seed * 1664525 + 1013904223) & 0xFFFFFFFF;
                return this.seed / 0x100000000;
            },
            setSeed: function (seed) {
                this.seed = seed;
            }
        };
        this.selection = this.create2DArray(this.nx, this.ny, () => new Material());
        this.delta_mx_index = 0;
        this.delta_my_index = 0;
        this.moving_selection = false;
        this.texting = false;
        this.text_x = 0;
        this.text_y = 0;
        this.mouseX = 0;
        this.mouseY = 0;
        this.under_brush = this.create2DArray(this.nx, this.ny, () => false);
        this.debugging = false;
        this.texts = [];
        this.regularfont = '12px Arial';
        this.bigfont = '16px Arial';
        this.startingpath = '.';
        this.fileextension = '.json.gz';
        this.infile = null;
        this.outfile = null;
        this.saveversion = 1;
        this.mousebutton = 0;
        this.mx_start = 0;
        this.my_start = 0;
        this.advanceframe = false;
        this.shift_down = false;
        this.ctrl_down = false;
        this.cut = false;
        this.copy = false;
        this.paste = false;
        this.delete = false;
        this.prev_scalar_view = ScalarView.NONE;
        this.prev_vector_view = VectorView.NONE;
        this.alt_down = false;
        this.clear = false;
        this.reset = false;
        this.save = false;
        this.load = false;
        this.brush_changed = false;
        this.updateMiscFields = false;

        this.opts = document.getElementById('panel'); // Reference #panel from index.html

        // Initialize Canvas and RenderCanvas
        this.canvas = document.getElementById('canvas');
        this.ctx = this.canvas.getContext('2d');
        this.r = new RenderCanvas(this, this.canvas);
        this.canvas.tabIndex = 0; // Make canvas focusable for key events

        // UI controls
        this.gui_simspeed = document.getElementById('gui_simspeed');
        this.gui_simspeed_2 = document.getElementById('gui_simspeed_2');
        this.gui_paused = document.getElementById('gui_paused');
        this.gui_reset = document.getElementById('gui_reset');
        this.gui_resetall = document.getElementById('gui_resetall');
        this.gui_save = document.getElementById('gui_save');
        this.gui_open = document.getElementById('gui_open');
        this.gui_help = document.getElementById('gui_help');
        this.gui_editdesc = document.getElementById('gui_editdesc');
        this.gui_material = document.getElementById('gui_material');
        this.gui_material.innerHTML = MaterialType_values
            .filter(m => m !== MaterialType.ABSORBER)
            .map(m => `<option value="${m.name}">${m.name}</option>`).join('');
        this.gui_material.selectedIndex = 0;
        this.gui_brush = document.getElementById('gui_brush');
        this.gui_brush.innerHTML = Brush_values
            .map(b => `<option value="${b.name}">${b.name}</option>`).join('');
        this.gui_brush.selectedIndex = 0;
        this.gui_brush_1 = document.getElementById('gui_brush_1');
        this.gui_brush_1.innerHTML = Object.values(BrushShape)
            .map(bs => `<option value="${bs.name}">${bs.name}</option>`).join('');
        this.gui_brush_1.selectedIndex = 1;
        this.gui_brushsize = document.getElementById('gui_brushsize');
        this.gui_brush_highlight = document.getElementById('gui_brush_highlight');
        this.gui_brush_highlight_label = document.getElementById('gui_brush_highlight_label');
        this.gui_view = document.getElementById('gui_view');
        this.gui_view.innerHTML = Object.values(ScalarView)
            .filter(v => v !== ScalarView.DEBUG)
            .map(v => `<option value="${v.name}">${v.name}</option>`).join('');
        this.gui_view.selectedIndex = 3;
        this.gui_view_vec = document.getElementById('gui_view_vec');
        this.gui_view_vec.innerHTML = Object.values(VectorView)
            .map(v => `<option value="${v.name}">${v.name}</option>`).join('');
        this.gui_view_vec.selectedIndex = 1;
        this.gui_bc = document.getElementById('gui_bc');
        this.gui_bc.innerHTML = Object.values(BoundaryCondition)
            .map(bc => `<option value="${bc.name}">${bc.name}</option>`).join('');
        this.gui_bc.selectedIndex = 0;
        this.gui_parameter2 = document.getElementById('gui_parameter2');
        this.gui_parameter2_text = document.getElementById('direction-label');
        this.gui_parameter3 = document.getElementById('gui_parameter3');
        this.gui_parameter3_text = document.getElementById('gui_parameter3_text');
        this.gui_brightness = document.getElementById('gui_brightness');
        this.gui_elem_colors = document.getElementById('gui_elem_colors');
        this.gui_tooltip = document.getElementById('gui_tooltip');
        this.gui_brightness_vec = document.getElementById('gui_brightness_vec');
        this.gui_view_vec_mode = document.getElementById('gui_view_vec_mode');
        this.gui_view_vec_mode.innerHTML = Object.values(VectorMode)
            .map(vm => `<option value="${vm.name}">${vm.name}</option>`).join('');
        this.gui_view_vec_mode.selectedIndex = 0;
        this.gui_text_bg = document.getElementById('gui_text_bg');
        this.gui_stepsizelbl = document.getElementById('step-size-label');
        this.gui_stepslbl = document.getElementById('steps-frame-label');
        this.gui_lblBrushSize = document.getElementById('brush-size-label');
        this.textPane = document.getElementById('textPane');
        this.help = document.getElementById('help');
        this.gui_close_help = document.getElementById('gui_close_help');

        // Initialize labels
        this.gui_stepsizelbl.textContent = `Step size: ${this.gui_simspeed.value}`;
        this.gui_lblBrushSize.textContent = `Brush size: ${this.gui_brushsize.value}`;
        this.gui_parameter3_text.textContent = '';

        // Event listeners
        this.canvas.addEventListener('mousedown', this.mousePressed.bind(this));
        this.canvas.addEventListener('mousemove', this.mouseMoved.bind(this));
        this.canvas.addEventListener('mouseup', this.mouseReleased.bind(this));
        this.canvas.addEventListener('mousemove', this.mouseDragged.bind(this));
        this.canvas.addEventListener('wheel', this.mouseWheelMoved.bind(this));
        this.canvas.addEventListener('keydown', this.keyPressed.bind(this));
        this.canvas.addEventListener('keyup', this.keyReleased.bind(this));
        this.canvas.addEventListener('keypress', this.keyTyped.bind(this));
        this.canvas.addEventListener('contextmenu', event => event.preventDefault());
        this.gui_reset.addEventListener('click', this.actionPerformed.bind(this));
        this.gui_resetall.addEventListener('click', this.actionPerformed.bind(this));
        this.gui_save.addEventListener('click', this.actionPerformed.bind(this));
        this.gui_open.addEventListener('click', this.actionPerformed.bind(this));
        this.gui_help.addEventListener('click', this.actionPerformed.bind(this));
        this.gui_close_help.addEventListener('click', this.actionPerformed.bind(this));
        this.gui_editdesc.addEventListener('click', this.actionPerformed.bind(this));
        this.gui_material.addEventListener('change', this.actionPerformed.bind(this));
        this.gui_brush.addEventListener('change', this.actionPerformed.bind(this));
        this.gui_brush_1.addEventListener('change', this.actionPerformed.bind(this));
        this.gui_brushsize.addEventListener('input', this.actionPerformed.bind(this));
        this.gui_brush_highlight.addEventListener('change', this.actionPerformed.bind(this));
        this.gui_view.addEventListener('change', this.actionPerformed.bind(this));
        this.gui_view_vec.addEventListener('change', this.actionPerformed.bind(this));
        this.gui_bc.addEventListener('change', this.actionPerformed.bind(this));
        this.gui_parameter2.addEventListener('input', this.actionPerformed.bind(this));
        this.gui_parameter3.addEventListener('input', this.actionPerformed.bind(this));
        this.gui_brightness.addEventListener('input', this.actionPerformed.bind(this));
        this.gui_elem_colors.addEventListener('change', this.actionPerformed.bind(this));
        this.gui_tooltip.addEventListener('change', this.actionPerformed.bind(this));
        this.gui_brightness_vec.addEventListener('input', this.actionPerformed.bind(this));
        this.gui_view_vec_mode.addEventListener('change', this.actionPerformed.bind(this));
        this.gui_text_bg.addEventListener('change', this.actionPerformed.bind(this));

        // Initialize grid
        this.detect64Bit();
        this.initializeGrid(this.default_width, this.default_resolution);

        this.canvas.tabIndex = 0;

        // Initialize Font7x5
        this.Font7x5 = new Font7x5();

        var initial_file = location.search.split('file=')[1];
        this.readFile(initial_file);
    }

    detect64Bit() {
        console.warn('64-bit check skipped; JavaScript assumes 64-bit environment');
    }

    resetFields(resetall) {
        for (let i = 0; i < this.nx; i++) {
            for (let j = 0; j < this.ny; j++) {
                if (resetall) {
                    if (!this.materials[i][j]) this.materials[i][j] = new Material();
                    else this.materials[i][j].erase();

                    if (!this.selection[i][j]) this.selection[i][j] = new Material();
                    else this.selection[i][j].erase();

                    if (!this.clipboard[i][j]) this.clipboard[i][j] = new Material();
                    else this.clipboard[i][j].erase();

                    this.F0_n[i][j] = 0;
                    this.F0_p[i][j] = 0;
                    this.F_n[i][j] = 0;
                    this.F_p[i][j] = 0;
                    this.E0_n[i][j] = 0;
                    this.E0_p[i][j] = 0;
                    this.K[i][j] = 0;
                    this.E_a[i][j] = 0;
                    this.R[i][j] = 0;

                    this.cmfx_n[i][j] = 0;
                    this.cmfy_n[i][j] = 0;
                    this.cmfx_p[i][j] = 0;
                    this.cmfy_p[i][j] = 0;

                    this.emfx[i][j] = 0.0;
                    this.emfy[i][j] = 0.0;

                    this.epsx[i][j] = this.eps0;
                    this.epsy[i][j] = this.eps0;
                    this.mu_z[i][j] = this.mu0;

                    this.conducting[i][j] = 0;
                    this.conducting_x[i][j] = 0;
                    this.conducting_y[i][j] = 0;

                    this.absorptivity[i][j] = 0;
                    this.absorptivity_x[i][j] = 0;
                    this.absorptivity_y[i][j] = 0;
                }

                this.Ex[i][j] = 0.0;
                this.Ey[i][j] = 0.0;
                this.Bz[i][j] = 0.0;
                this.Hz_laplacian[i][j] = 0.0;

                this.rho_abs[i][j] = 0.0;
                this.rho_n[i][j] = 0.0;
                this.rho_p[i][j] = 0.0;
                this.rho_back[i][j] = 0.0;
                this.rho_free[i][j] = 0.0;
                this.mobility_factor[i][j] = 0.0;

                this.Jx_abs[i][j] = 0.0;
                this.Jy_abs[i][j] = 0.0;
                this.Jx_n[i][j] = 0.0;
                this.Jy_n[i][j] = 0.0;
                this.Jx_p[i][j] = 0.0;
                this.Jy_p[i][j] = 0.0;
                this.Jx_free[i][j] = 0.0;
                this.Jy_free[i][j] = 0.0;

                this.MG_phi1[i][j] = 0.0;
                this.MG_phi2[i][j] = 0.0;

                this.distance[i][j] = Number.MAX_SAFE_INTEGER;
                this.visited[i][j] = false;

                this.Dx[i][j] = 0.0;
                this.Dy[i][j] = 0.0;
                this.Hz[i][j] = 0.0;
                this.Sx[i][j] = 0.0;
                this.Sy[i][j] = 0.0;
                this.u[i][j] = 0.0;
                this.phi[i][j] = 0.0;

                this.G[i][j] = 0.0;
                this.F_n[i][j] = 0.0;
                this.F_p[i][j] = 0.0;
                this.grad_E0x_n[i][j] = 0.0;
                this.grad_E0y_n[i][j] = 0.0;
                this.grad_E0x_p[i][j] = 0.0;
                this.grad_E0y_p[i][j] = 0.0;
                this.grad_Fx_n[i][j] = 0.0;
                this.grad_Fy_n[i][j] = 0.0;
                this.grad_Fx_p[i][j] = 0.0;
                this.grad_Fy_p[i][j] = 0.0;
                this.Q[i][j] = 0.0;
                this.S[i][j] = 0.0;
                this.F[i][j] = 0.0;
                this.debug[i][j] = 0.0;

                this.selected[i][j] = false;
                this.selected_EMF[i][j] = false;
            }
        }

        if (resetall) {
            this.voltageprobes = [];
            this.currentprobes = [];
            this.ground = null;
            // opts.setTitle replaced with no-op; HTML title is static
        }

        this.constructBoundary();
        this.initializeAllMaterials();
        this.updateAllMaterials();
        this.checkCFL();
    }

    constructBoundary() {
        this.absorber_width = Math.floor(0.05 * this.nx);
        this.absorbing_coeff = 50 * this.c / this.width;
        const max_stretch = 10;
        const coeff = Math.log(max_stretch);

        // Placeholder: gui_bc.value should match BoundaryCondition name
        const selectedBC = Object.values(BoundaryCondition).find(bc => bc.name === this.gui_bc.value) || BoundaryCondition.DISSIPATIVE;

        if (selectedBC === BoundaryCondition.DISSIPATIVE) {
            for (let i = 0; i < this.nx; i++) {
                for (let j = 0; j < this.ny; j++) {
                    const dx = Math.max(0, this.absorber_width - Math.min(i, this.nx - 1 - i)) / this.absorber_width;
                    const dy = Math.max(0, this.absorber_width - Math.min(j, this.ny - 1 - j)) / this.absorber_width;
                    const depth = Math.sqrt(dx * dx + dy * dy);
                    if (depth > 0) {
                        const stretchfactor = Math.exp(coeff * depth);

                        this.materials[i][j].erase();
                        this.materials[i][j].type = MaterialType.ABSORBER;
                        this.materials[i][j].eps_r = stretchfactor;
                        this.materials[i][j].mu_r = stretchfactor;
                        this.materials[i][j].absorptivity = 1;
                    }
                }
            }
        } else {
            for (let i = 0; i < this.nx; i++) {
                for (let j = 0; j < this.ny; j++) {
                    if (this.materials[i][j].type === MaterialType.ABSORBER) {
                        this.materials[i][j].erase();
                    }
                }
            }
        }
    }



    iterateSimulation() {
        this.t6.start();

        this.stepnumber++;

        Module._iterateSimulation(
            this.nx, this.ny, this.ny + 1, this.dt, this.ds,
            this.Hz_dissipation,
            this.absorbing_coeff,
            this.mu_electron,
            this.mu_hole,
            this.D_electron,
            this.D_hole,
            this.q_n,
            this.q_p,
            this.E_sat,

            this.ptr_Hz,
            this.ptr_Hz_laplacian,
            this.ptr_Ex,
            this.ptr_Ey,
            this.ptr_mu_z,
            this.ptr_absorptivity,
            this.ptr_epsx,
            this.ptr_epsy,
            this.ptr_absorptivity_x,
            this.ptr_absorptivity_y,
            this.ptr_emfx,
            this.ptr_emfy,
            this.ptr_cmfx_n,
            this.ptr_cmfx_p,
            this.ptr_cmfy_n,
            this.ptr_cmfy_p,
            this.ptr_mobility_factor,

            this.ptr_rho_n,
            this.ptr_rho_p,
            this.ptr_rho_abs,
            this.ptr_rho_back,
            this.ptr_rho_free,
            this.ptr_R,
            this.ptr_K,
            this.ptr_conducting,

            this.ptr_Jx_n,
            this.ptr_Jx_p,
            this.ptr_Jx_abs,
            this.ptr_Jy_n,
            this.ptr_Jy_p,
            this.ptr_Jy_abs,
            this.ptr_conducting_x,
            this.ptr_conducting_y);

        this.t6.stop();

        if (this.stepnumber % 100 === 0) {
            this.t4.start();
            this.multigridSolve(true, false);
            this.t4.stop();
        }

        this.time += this.dt;
        this.advanceframe = false;
    }

    multigridSolve(computeE, computePhi) {
        Module._multigridSolve(this.nx, this.ny, this.ny+this.WASM_array_padding, this.log2_resolution, computeE, computePhi,
            this.ptr_MG_phi1,
            this.ptr_MG_phi2,
            this.ptr_MG_rho,
            this.ptr_MG_rho0,
            this.ptr_MG_epsx,
            this.ptr_MG_epsy,
            this.ptr_MG_eps_avg,
            this.ptr_epsx,
            this.ptr_epsy,
            this.ptr_Ex,
            this.ptr_Ey,
            this.ptr_phi,
            this.ptr_rho_free,
            this.ds,
            this.width);

        if (computePhi)
            this.copyWasmHeapTo2DArray(this.ptr_phi, this.nx, this.ny, this.phi);
        if (computeE) {
            this.copyWasmHeapTo2DArray(this.ptr_Ex, this.nx, this.ny, this.Ex);
            this.copyWasmHeapTo2DArray(this.ptr_Ey, this.nx, this.ny, this.Ey);
        }
    }

    // Existing methods (unchanged, included for completeness)
    detect64Bit() {
        console.warn('64-bit check skipped; JavaScript assumes 64-bit environment');
    }

    run() {
        try {
            const iterationmultiplier = parseInt(this.gui_simspeed_2.value);

            if (this.clear) {
                this.resetFields(false);
                this.multigridSolve(true, false);
                this.time = 0.0;
                this.clear = false;
            }

            if (this.reset) {
                if (confirm('Do you wish to reset the entire simulation?')) {
                    this.resetFields(true);
                    this.time = 0.0;
                }
                this.reset = false;
            }

            if (parseInt(this.gui_simspeed.value) !== this.lastsimspeed) {
                this.lastsimspeed = parseInt(this.gui_simspeed.value);
                this.dt = this.dt_maximum * (this.lastsimspeed / 20.0);
            }

            this.handleMouseInput();

            if (!this.gui_paused.checked || this.advanceframe) {
                for (let i = 0; i < iterationmultiplier; i++) {
                    this.iterateSimulation();
                }
                this.downloadDataFromWasmHeap();

                if (this.frame % 2 === 0) {
                    this.calcMiscFields(true);
                } else {
                    this.calcMiscFields(false);
                }

                this.frame++;
            } else if (this.updateMiscFields) {
                this.calcMiscFields(true);
            }

            if (this.save) {
                this.writeFile();
                this.save = false;
            }

            if (this.load) {
                this.readFile();
                this.load = false;
            }

            this.r.render();
        } catch (e) {
            alert(`Error: ${e.message}`);
            throw e;
        }

        requestAnimationFrame(this.run.bind(this));
    }


    initalizePointers() {
        this.ptr_Hz = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_Hz_laplacian = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_Ex = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_Ey = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_mu_z = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_absorptivity = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_epsx = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_epsy = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_absorptivity_x = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_absorptivity_y = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_emfx = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_emfy = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_cmfx_n = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_cmfx_p = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_cmfy_n = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_cmfy_p = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_mobility_factor = this.allocateHeapArray(this.nx, this.ny);

        this.ptr_rho_n = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_rho_p = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_rho_abs = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_rho_back = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_rho_free = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_R = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_K = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_conducting = this.allocateHeapArray(this.nx, this.ny);

        this.ptr_Jx_n = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_Jx_p = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_Jx_abs = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_Jy_n = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_Jy_p = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_Jy_abs = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_conducting_x = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_conducting_y = this.allocateHeapArray(this.nx, this.ny);

        this.ptr_MG_phi1 = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_MG_phi2 = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_MG_rho = this.allocate3DHeapArray(this.nx, this.ny);
        this.ptr_MG_rho0 = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_MG_epsx = this.allocate3DHeapArray(this.nx, this.ny);
        this.ptr_MG_epsy = this.allocate3DHeapArray(this.nx, this.ny);
        this.ptr_MG_eps_avg = this.allocateHeapArray(this.nx, this.ny);
        this.ptr_phi = this.allocateHeapArray(this.nx, this.ny);
    }

    uploadDataToWasmHeap() {
        this.copy2DArrayToWasmHeap(this.Hz, this.nx, this.ny, this.ptr_Hz);
        this.copy2DArrayToWasmHeap(this.Hz_laplacian, this.nx, this.ny, this.ptr_Hz_laplacian);
        this.copy2DArrayToWasmHeap(this.Ex, this.nx, this.ny, this.ptr_Ex);
        this.copy2DArrayToWasmHeap(this.Ey, this.nx, this.ny, this.ptr_Ey);
        this.copy2DArrayToWasmHeap(this.mu_z, this.nx, this.ny, this.ptr_mu_z);
        this.copy2DArrayToWasmHeap(this.absorptivity, this.nx, this.ny, this.ptr_absorptivity);
        this.copy2DArrayToWasmHeap(this.epsx, this.nx, this.ny, this.ptr_epsx);
        this.copy2DArrayToWasmHeap(this.epsy, this.nx, this.ny, this.ptr_epsy);
        this.copy2DArrayToWasmHeap(this.absorptivity_x, this.nx, this.ny, this.ptr_absorptivity_x);
        this.copy2DArrayToWasmHeap(this.absorptivity_y, this.nx, this.ny, this.ptr_absorptivity_y);
        this.copy2DArrayToWasmHeap(this.emfx, this.nx, this.ny, this.ptr_emfx);
        this.copy2DArrayToWasmHeap(this.emfy, this.nx, this.ny, this.ptr_emfy);
        this.copy2DArrayToWasmHeap(this.cmfx_n, this.nx, this.ny, this.ptr_cmfx_n);
        this.copy2DArrayToWasmHeap(this.cmfx_p, this.nx, this.ny, this.ptr_cmfx_p);
        this.copy2DArrayToWasmHeap(this.cmfy_n, this.nx, this.ny, this.ptr_cmfy_n);
        this.copy2DArrayToWasmHeap(this.cmfy_p, this.nx, this.ny, this.ptr_cmfy_p);
        this.copy2DArrayToWasmHeap(this.mobility_factor, this.nx, this.ny, this.ptr_mobility_factor);
        this.copy2DArrayToWasmHeap(this.rho_n, this.nx, this.ny, this.ptr_rho_n);
        this.copy2DArrayToWasmHeap(this.rho_p, this.nx, this.ny, this.ptr_rho_p);
        this.copy2DArrayToWasmHeap(this.rho_abs, this.nx, this.ny, this.ptr_rho_abs);
        this.copy2DArrayToWasmHeap(this.rho_back, this.nx, this.ny, this.ptr_rho_back);
        this.copy2DArrayToWasmHeap(this.rho_free, this.nx, this.ny, this.ptr_rho_free);
        this.copy2DArrayToWasmHeap(this.R, this.nx, this.ny, this.ptr_R);
        this.copy2DArrayToWasmHeap(this.K, this.nx, this.ny, this.ptr_K);
        this.copy2DArrayToWasmHeap(this.conducting, this.nx, this.ny, this.ptr_conducting);
        this.copy2DArrayToWasmHeap(this.Jx_n, this.nx, this.ny, this.ptr_Jx_n);
        this.copy2DArrayToWasmHeap(this.Jx_p, this.nx, this.ny, this.ptr_Jx_p);
        this.copy2DArrayToWasmHeap(this.Jx_abs, this.nx, this.ny, this.ptr_Jx_abs);
        this.copy2DArrayToWasmHeap(this.Jy_n, this.nx, this.ny, this.ptr_Jy_n);
        this.copy2DArrayToWasmHeap(this.Jy_p, this.nx, this.ny, this.ptr_Jy_p);
        this.copy2DArrayToWasmHeap(this.Jy_abs, this.nx, this.ny, this.ptr_Jy_abs);
        this.copy2DArrayToWasmHeap(this.conducting_x, this.nx, this.ny, this.ptr_conducting_x);
        this.copy2DArrayToWasmHeap(this.conducting_y, this.nx, this.ny, this.ptr_conducting_y);

        this.copy2DArrayToWasmHeap(this.MG_phi1, this.nx, this.ny, this.ptr_MG_phi1);
        this.copy2DArrayToWasmHeap(this.MG_phi2, this.nx, this.ny, this.ptr_MG_phi2);
        this.copy3DArrayToWasmHeap(this.MG_rho, this.nx, this.ny, this.ptr_MG_rho);
        this.copy2DArrayToWasmHeap(this.MG_rho0, this.nx, this.ny, this.ptr_MG_rho0);
        this.copy3DArrayToWasmHeap(this.MG_epsx, this.nx, this.ny, this.ptr_MG_epsx);
        this.copy3DArrayToWasmHeap(this.MG_epsy, this.nx, this.ny, this.ptr_MG_epsy);
        this.copy2DArrayToWasmHeap(this.MG_eps_avg, this.nx, this.ny, this.ptr_MG_eps_avg);
        this.copy2DArrayToWasmHeap(this.phi, this.nx, this.ny, this.ptr_phi);
    }

    downloadDataFromWasmHeap() {
        this.copyWasmHeapTo2DArray(this.ptr_Hz, this.nx, this.ny, this.Hz);
        this.copyWasmHeapTo2DArray(this.ptr_Ex, this.nx, this.ny, this.Ex);
        this.copyWasmHeapTo2DArray(this.ptr_Ey, this.nx, this.ny, this.Ey);
        this.copyWasmHeapTo2DArray(this.ptr_rho_n, this.nx, this.ny, this.rho_n);
        this.copyWasmHeapTo2DArray(this.ptr_rho_p, this.nx, this.ny, this.rho_p);
        this.copyWasmHeapTo2DArray(this.ptr_rho_abs, this.nx, this.ny, this.rho_abs);
        this.copyWasmHeapTo2DArray(this.ptr_rho_back, this.nx, this.ny, this.rho_back);
        this.copyWasmHeapTo2DArray(this.ptr_rho_free, this.nx, this.ny, this.rho_free);
        this.copyWasmHeapTo2DArray(this.ptr_R, this.nx, this.ny, this.R);
        this.copyWasmHeapTo2DArray(this.ptr_Jx_n, this.nx, this.ny, this.Jx_n);
        this.copyWasmHeapTo2DArray(this.ptr_Jx_p, this.nx, this.ny, this.Jx_p);
        this.copyWasmHeapTo2DArray(this.ptr_Jx_abs, this.nx, this.ny, this.Jx_abs);
        this.copyWasmHeapTo2DArray(this.ptr_Jy_n, this.nx, this.ny, this.Jy_n);
        this.copyWasmHeapTo2DArray(this.ptr_Jy_p, this.nx, this.ny, this.Jy_p);
        this.copyWasmHeapTo2DArray(this.ptr_Jy_abs, this.nx, this.ny, this.Jy_abs);
    }

    allocateHeapArray(nx, ny, type = Float64Array) {
        const bytesPerElement = type.BYTES_PER_ELEMENT;
        const ptr = Module._malloc((nx+this.WASM_array_padding) * (ny+this.WASM_array_padding) * bytesPerElement);
        //const ptr = Module._aligned_alloc(64, nx * ny * bytesPerElement)
        return ptr;
    }

    allocate3DHeapArray(nx, ny, type = Float64Array) {
        const bytesPerElement = type.BYTES_PER_ELEMENT;
        const ptr = Module._malloc((this.log2_resolution+1) * (nx+this.WASM_array_padding) * (ny+this.WASM_array_padding) * bytesPerElement);
        //const ptr = Module._aligned_alloc(64, (this.log2_resolution+1) * nx * ny * bytesPerElement);
        return ptr;
    }

    copy2DArrayToWasmHeap(js2DArray, nx, ny, ptr, type = Float64Array) {
        const flatArray = new type((nx+this.WASM_array_padding) * (ny+this.WASM_array_padding));
        for (let i = 0; i < nx; i++) {
            for (let j = 0; j < ny; j++) {
                flatArray[i * (ny+this.WASM_array_padding) + j] = js2DArray[i][j];
            }
        }

        const bytesPerElement = type.BYTES_PER_ELEMENT;

        // Select correct heap view
        const heap = (type === Float32Array) ? Module.HEAPF32 :
            (type === Float64Array) ? Module.HEAPF64 :
                (type === Int32Array) ? Module.HEAP32 :
                    (type === Uint8Array) ? Module.HEAPU8 :
                        (() => { throw new Error("Unsupported type"); })();

        heap.set(flatArray, ptr / bytesPerElement);
        return ptr;
    }

    copyWasmHeapTo2DArray(ptr, nx, ny, result, type = Float64Array) {
        const bytesPerElement = type.BYTES_PER_ELEMENT;

        const heap = (type === Float32Array) ? Module.HEAPF32 :
            (type === Float64Array) ? Module.HEAPF64 :
                (type === Int32Array) ? Module.HEAP32 :
                    (type === Uint8Array) ? Module.HEAPU8 :
                        (() => { throw new Error("Unsupported type"); })();

        const flat = new type(heap.buffer, ptr, (nx+this.WASM_array_padding) * (ny+this.WASM_array_padding));

        for (let i = 0; i < nx; i++) {
            for (let j = 0; j < ny; j++) {
                result[i][j] = flat[i * (ny+this.WASM_array_padding) + j];
            }
        }
    }

    copy3DArrayToWasmHeap(js3DArray, nx, ny, ptr, type = Float64Array) {
        const flatArray = new type((this.log2_resolution+1) * (nx+this.WASM_array_padding) * (ny+this.WASM_array_padding));
        for (let k = 0; k < (this.log2_resolution+1); k++) {
            for (let i = 0; i < nx; i++) {
                for (let j = 0; j < ny; j++) {
                    flatArray[k * (nx+this.WASM_array_padding) * (ny+this.WASM_array_padding) + i * (ny+this.WASM_array_padding) + j] = js3DArray[k][i][j];
                }
            }
        }

        const bytesPerElement = type.BYTES_PER_ELEMENT;

        // Select correct heap view
        const heap = (type === Float32Array) ? Module.HEAPF32 :
            (type === Float64Array) ? Module.HEAPF64 :
                (type === Int32Array) ? Module.HEAP32 :
                    (type === Uint8Array) ? Module.HEAPU8 :
                        (() => { throw new Error("Unsupported type"); })();

        heap.set(flatArray, ptr / bytesPerElement);
        return ptr;
    }

    copyWasmHeapTo3DArray(ptr, nx, ny, result, type = Float64Array) {
        const bytesPerElement = type.BYTES_PER_ELEMENT;

        const heap = (type === Float32Array) ? Module.HEAPF32 :
            (type === Float64Array) ? Module.HEAPF64 :
                (type === Int32Array) ? Module.HEAP32 :
                    (type === Uint8Array) ? Module.HEAPU8 :
                        (() => { throw new Error("Unsupported type"); })();

        const flat = new type(heap.buffer, ptr, (this.log2_resolution+1) * (nx+this.WASM_array_padding) * (ny+this.WASM_array_padding));

        for (let k = 0; k < (this.log2_resolution+1); k++) {
            for (let i = 0; i < nx; i++) {
                for (let j = 0; j < ny; j++) {
                    result[i][j] = flat[k * (nx+this.WASM_array_padding) * (ny+this.WASM_array_padding) + i * (ny+this.WASM_array_padding) + j];
                }
            }
        }
    }

    create2DArray(nx, ny, type = Float64Array, initializer = () => 0) {
        const arr = new Array(nx);
        for (let i = 0; i < nx; i++) {
            arr[i] = new type(ny);
            if (initializer() !== 0) {
                for (let j = 0; j < ny; j++) {
                    arr[i][j] = initializer();
                }
            }
        }
        return arr;
    }

    create3DArray(nz, nx, ny, type = Float64Array) {
        const arr = new Array(nz);
        for (let k = 0; k < nz; k++) {
            arr[k] = this.create2DArray(nx, ny, type);
        }
        return arr;
    }

    initializeGrid(width, resolution) {
        this.ds = width / resolution;
        this.nx = resolution;
        this.ny = resolution;
        this.width = this.nx * this.ds;
        this.log2_resolution = Math.round(Math.log(resolution) / Math.log(2));
        console.assert((1 << this.log2_resolution) === resolution, "resolution must be a power of 2");
        this.dt_maximum = 0.9 * this.ds / (Math.sqrt(2) * this.c);
        this.Hz_dissipation = 0.01 * this.ds * this.ds / this.dt_maximum;

        this.Ex = this.create2DArray(this.nx, this.ny);
        this.Ey = this.create2DArray(this.nx, this.ny);
        this.Bz = this.create2DArray(this.nx, this.ny);
        this.Hz_laplacian = this.create2DArray(this.nx, this.ny);
        this.rho_abs = this.create2DArray(this.nx, this.ny);
        this.rho_n = this.create2DArray(this.nx, this.ny);
        this.rho_p = this.create2DArray(this.nx, this.ny);
        this.rho_back = this.create2DArray(this.nx, this.ny);
        this.rho_free = this.create2DArray(this.nx, this.ny);
        this.mobility_factor = this.create2DArray(this.nx, this.ny);
        this.Jx_abs = this.create2DArray(this.nx, this.ny);
        this.Jy_abs = this.create2DArray(this.nx, this.ny);
        this.Jx_n = this.create2DArray(this.nx, this.ny);
        this.Jy_n = this.create2DArray(this.nx, this.ny);
        this.Jx_p = this.create2DArray(this.nx, this.ny);
        this.Jy_p = this.create2DArray(this.nx, this.ny);
        this.Jx_free = this.create2DArray(this.nx, this.ny);
        this.Jy_free = this.create2DArray(this.nx, this.ny);
        this.materials = this.create2DArray(this.nx, this.ny, Array, () => new Material());
        this.selection = this.create2DArray(this.nx, this.ny, Array, () => new Material());
        this.clipboard = this.create2DArray(this.nx, this.ny, Array, () => new Material());
        this.K = this.create2DArray(this.nx, this.ny);
        this.F0_n = this.create2DArray(this.nx, this.ny);
        this.F0_p = this.create2DArray(this.nx, this.ny);
        this.F_n = this.create2DArray(this.nx, this.ny);
        this.F_p = this.create2DArray(this.nx, this.ny);
        this.E0_n = this.create2DArray(this.nx, this.ny);
        this.E0_p = this.create2DArray(this.nx, this.ny);
        this.E_a = this.create2DArray(this.nx, this.ny);
        this.R = this.create2DArray(this.nx, this.ny);
        this.cmfx_n = this.create2DArray(this.nx, this.ny);
        this.cmfy_n = this.create2DArray(this.nx, this.ny);
        this.cmfx_p = this.create2DArray(this.nx, this.ny);
        this.cmfy_p = this.create2DArray(this.nx, this.ny);
        this.conducting = this.create2DArray(this.nx, this.ny, Int32Array);
        this.conducting_x = this.create2DArray(this.nx, this.ny, Int32Array);
        this.conducting_y = this.create2DArray(this.nx, this.ny, Int32Array);
        this.absorptivity = this.create2DArray(this.nx, this.ny);
        this.absorptivity_x = this.create2DArray(this.nx, this.ny);
        this.absorptivity_y = this.create2DArray(this.nx, this.ny);
        this.emfx = this.create2DArray(this.nx, this.ny);
        this.emfy = this.create2DArray(this.nx, this.ny);
        this.epsx = this.create2DArray(this.nx, this.ny);
        this.epsy = this.create2DArray(this.nx, this.ny);
        this.mu_z = this.create2DArray(this.nx, this.ny);
        this.Dx = this.create2DArray(this.nx, this.ny);
        this.Dy = this.create2DArray(this.nx, this.ny);
        this.Hz = this.create2DArray(this.nx, this.ny);
        this.Sx = this.create2DArray(this.nx, this.ny);
        this.Sy = this.create2DArray(this.nx, this.ny);
        this.u = this.create2DArray(this.nx, this.ny);
        this.phi = this.create2DArray(this.nx, this.ny);
        this.G = this.create2DArray(this.nx, this.ny);
        this.grad_E0x_n = this.create2DArray(this.nx, this.ny);
        this.grad_E0y_n = this.create2DArray(this.nx, this.ny);
        this.grad_E0x_p = this.create2DArray(this.nx, this.ny);
        this.grad_E0y_p = this.create2DArray(this.nx, this.ny);
        this.grad_Fx_n = this.create2DArray(this.nx, this.ny);
        this.grad_Fy_n = this.create2DArray(this.nx, this.ny);
        this.grad_Fx_p = this.create2DArray(this.nx, this.ny);
        this.grad_Fy_p = this.create2DArray(this.nx, this.ny);
        this.Q = this.create2DArray(this.nx, this.ny);
        this.S = this.create2DArray(this.nx, this.ny);
        this.F = this.create2DArray(this.nx, this.ny);
        this.debug = this.create2DArray(this.nx, this.ny);
        this.MG_rho0 = this.create2DArray(this.nx, this.ny);
        this.MG_rho = this.create3DArray(this.log2_resolution+1, this.nx, this.ny);
        this.MG_epsx = this.create3DArray(this.log2_resolution+1, this.nx, this.ny);
        this.MG_epsy = this.create3DArray(this.log2_resolution+1, this.nx, this.ny);
        this.MG_eps_avg = this.create2DArray(this.nx, this.ny);
        this.MG_phi1 = this.create2DArray(this.nx, this.ny);
        this.MG_phi2 = this.create2DArray(this.nx, this.ny);
        this.distance = this.create2DArray(this.nx, this.ny, Int32Array);
        this.visited = this.create2DArray(this.nx, this.ny, Array, () => false);
        this.under_brush = this.create2DArray(this.nx, this.ny, Array, () => false);
        this.selected = this.create2DArray(this.nx, this.ny, Array, () => false);
        this.selected_EMF = this.create2DArray(this.nx, this.ny, Array, () => false);

        this.initalizePointers();

        this.voltageprobes = [];
        this.currentprobes = [];

        this.resetFields(true);
        this.multigridSolve(true, false);
        this.calcMiscFields(false);

        this.text_x = 0;
        this.text_y = 0;
        this.texting = false;

        this.image_r = this.create2DArray(this.nx, this.ny, Float32Array);
        this.image_g = this.create2DArray(this.nx, this.ny, Float32Array);
        this.image_b = this.create2DArray(this.nx, this.ny, Float32Array);
        this.scalefactor = Math.floor(768 / this.ny);
        this.imgwidth = Math.floor(this.scalefactor * this.nx);
        this.imgheight = Math.floor(this.scalefactor * this.ny);

        this.canvas.width = this.imgwidth;
        this.canvas.height = this.imgheight;
    }

    updateAllMaterials() {
        for (let i = 1; i < this.nx - 1; i++) {
            for (let j = 1; j < this.ny - 1; j++) {
                this.mu_z[i][j] = 0.25 * this.mu0 * (
                    this.materials[i][j].mu_r +
                    this.materials[i + 1][j].mu_r +
                    this.materials[i][j + 1].mu_r +
                    this.materials[i + 1][j + 1].mu_r
                );
                this.absorptivity[i][j] = 0.25 * (
                    this.materials[i][j].absorptivity +
                    this.materials[i + 1][j].absorptivity +
                    this.materials[i][j + 1].absorptivity +
                    this.materials[i + 1][j + 1].absorptivity
                );
            }
        }

        for (let i = 0; i < this.nx; i++) {
            for (let j = 0; j < this.ny; j++) {
                this.rho_back[i][j] = this.materials[i][j].rho_back;
                this.conducting[i][j] = this.materials[i][j].conducting*this.materials[i][j].activated;
            }
        }

        for (let i = 0; i < this.nx - 1; i++) {
            for (let j = 1; j < this.ny - 1; j++) {
                this.conducting_x[i][j] = Math.min(this.materials[i + 1][j].conducting*this.materials[i + 1][j].activated, this.materials[i][j].conducting*this.materials[i][j].activated);
                this.emfx[i][j] = 0.5 * (
                    this.materials[i + 1][j].emf * Math.cos(this.materials[i + 1][j].emf_direction) * this.materials[i + 1][j].activated +
                    this.materials[i][j].emf * Math.cos(this.materials[i][j].emf_direction) * this.materials[i][j].activated
                );
                this.epsx[i][j] = this.eps0 * 0.5 * (this.materials[i + 1][j].eps_r + this.materials[i][j].eps_r);
                this.absorptivity_x[i][j] = Math.min(this.materials[i + 1][j].absorptivity, this.materials[i][j].absorptivity);
            }
        }

        for (let i = 1; i < this.nx - 1; i++) {
            for (let j = 0; j < this.ny - 1; j++) {
                this.conducting_y[i][j] = Math.min(this.materials[i][j + 1].conducting*this.materials[i][j + 1].activated, this.materials[i][j].conducting*this.materials[i][j].activated);
                this.emfy[i][j] = 0.5 * (
                    this.materials[i][j + 1].emf * Math.sin(this.materials[i][j + 1].emf_direction) * this.materials[i][j + 1].activated +
                    this.materials[i][j].emf * Math.sin(this.materials[i][j].emf_direction) * this.materials[i][j].activated
                );
                this.epsy[i][j] = this.eps0 * 0.5 * (this.materials[i][j + 1].eps_r + this.materials[i][j].eps_r);
                this.absorptivity_y[i][j] = Math.min(this.materials[i][j + 1].absorptivity, this.materials[i][j].absorptivity);
            }
        }

        this.computeChemicalForces();

        for (let i = 0; i < this.nx; i++) {
            for (let j = 0; j < this.ny; j++) {
                if (this.K[i][j] > 0 || this.materials[i][j].type == MaterialType.SWITCH) {
                    if (this.rho_n[i][j] === 0 && this.rho_p[i][j] === 0) {
                        this.rho_n[i][j] = this.calcEquilibriumElectronCharge(this.rho_back[i][j], this.K[i][j]);
                        this.rho_p[i][j] = this.calcEquilibriumHoleCharge(this.rho_back[i][j], this.K[i][j]);
                    }
                } else {
                    this.rho_n[i][j] = 0;
                    this.rho_p[i][j] = 0;
                }

                if (this.materials[i][j].absorptivity === 0) {
                    this.rho_abs[i][j] = 0;
                }

                this.rho_free[i][j] = this.rho_abs[i][j] + this.rho_n[i][j] + this.rho_p[i][j] + this.rho_back[i][j];
            }
        }

        for (let i = 0; i < this.nx - 1; i++) {
            for (let j = 0; j < this.ny; j++) {
                this.Jx_n[i][j] *= this.conducting_x[i][j];
                this.Jx_p[i][j] *= this.conducting_x[i][j];
                this.Jx_abs[i][j] *= (this.materials[i][j].absorptivity > 0) ? 1 : 0;
            }
        }

        for (let i = 0; i < this.nx; i++) {
            for (let j = 0; j < this.ny - 1; j++) {
                this.Jy_n[i][j] *= this.conducting_y[i][j];
                this.Jy_p[i][j] *= this.conducting_y[i][j];
                this.Jy_abs[i][j] *= (this.absorptivity_y[i][j] > 0) ? 1 : 0;
            }
        }

        this.prescaleDielectric();
        this.uploadDataToWasmHeap();
    }

    updateJustEMFs() {
        for (let i = 0; i < this.nx - 1; i++) {
            for (let j = 1; j < this.ny - 1; j++) {
                this.emfx[i][j] = 0.5 * (
                    this.materials[i + 1][j].emf * Math.cos(this.materials[i + 1][j].emf_direction) * this.materials[i + 1][j].activated +
                    this.materials[i][j].emf * Math.cos(this.materials[i][j].emf_direction) * this.materials[i][j].activated
                );
            }
        }

        for (let i = 1; i < this.nx - 1; i++) {
            for (let j = 0; j < this.ny - 1; j++) {
                this.emfy[i][j] = 0.5 * (
                    this.materials[i][j + 1].emf * Math.sin(this.materials[i][j + 1].emf_direction) * this.materials[i][j + 1].activated +
                    this.materials[i][j].emf * Math.sin(this.materials[i][j].emf_direction) * this.materials[i][j].activated
                );
            }
        }
        this.uploadDataToWasmHeap();
    }

    computeFreeEnergy(i, j, ni, W, B) {
        const K = ni * ni;

        this.E0_n[i][j] = -W + 0.5 * B;
        this.E0_p[i][j] = W + 0.5 * B;

        const TS_n = (Math.log(K) / this.beta + this.E0_n[i][j] + this.E0_p[i][j]) / 2.0;
        const TS_p = TS_n;

        this.F0_n[i][j] = this.E0_n[i][j] - TS_n;
        this.F0_p[i][j] = this.E0_p[i][j] - TS_p;
    }

    calcEquilibriumElectronCharge(rho_back, K) {
        const B = rho_back / this.e_charge;
        return -this.e_charge * 0.5 * (B + Math.sqrt(B * B + 4 * K));
    }

    calcEquilibriumHoleCharge(rho_back, K) {
        const B = -rho_back / this.e_charge;
        return this.e_charge * 0.5 * (B + Math.sqrt(B * B + 4 * K));
    }

    markNeighborhood(i1, j1, max_distance) {
        for (let i = i1 - max_distance; i <= i1 + max_distance; i++) {
            for (let j = j1 - max_distance; j <= j1 + max_distance; j++) {
                if (i >= 0 && j >= 0 && i < this.nx && j < this.ny) {
                    this.distance[i][j] = max_distance+1;
                    this.visited[i][j] = false;
                }
            }
        }
        this.distance[i1][j1] = 0;

        while (true) {
            let smallest_length = max_distance+1;
            let smallest_i = 0;
            let smallest_j = 0;
            for (let i = i1 - max_distance; i <= i1 + max_distance; i++) {
                for (let j = j1 - max_distance; j <= j1 + max_distance; j++) {
                    if (i >= 0 && j >= 0 && i < this.nx && j < this.ny) {
                        if (this.distance[i][j] < smallest_length && !this.visited[i][j] && this.conducting[i][j] === 1) {
                            smallest_length = this.distance[i][j];
                            smallest_i = i;
                            smallest_j = j;
                        }
                    }
                }
            }

            if (smallest_length === max_distance+1) break;

            for (let di = -1; di <= 1; di++) {
                for (let dj = -1; dj <= 1; dj++) {
                    if (Math.abs(di) + Math.abs(dj) > 0 &&
                        smallest_i + di >= 0 && smallest_j + dj >= 0 &&
                        smallest_i + di < this.nx && smallest_j + dj < this.ny &&
                        this.distance[smallest_i][smallest_j] < max_distance &&
                        this.conducting[smallest_i + di][smallest_j + dj] === 1) {
                        this.distance[smallest_i + di][smallest_j + dj] = Math.min(
                            this.distance[smallest_i + di][smallest_j + dj],
                            this.distance[smallest_i][smallest_j] + 1
                        );
                    }
                }
            }

            this.visited[smallest_i][smallest_j] = true;
        }
    }

    computeChemicalForces() {
        for (let i = 0; i < this.nx; i++) {
            for (let j = 0; j < this.ny; j++) {
                if (this.conducting[i][j] === 1) {
                    this.computeFreeEnergy(i, j, this.materials[i][j].ni, this.materials[i][j].W, this.materials[i][j].Eb);
                } else {
                    this.F0_n[i][j] = 0;
                    this.F0_p[i][j] = 0;
                    this.E0_n[i][j] = 0;
                    this.E0_p[i][j] = 0;
                }
                this.E_a[i][j] = this.materials[i][j].Ea;
            }
        }

        // Smooth free energy for stability
        const F_n_tmp = this.copyArray(this.F0_n);
        const F_p_tmp = this.copyArray(this.F0_p);
        const E_n_tmp = this.copyArray(this.E0_n);
        const E_p_tmp = this.copyArray(this.E0_p);
        const Ea_tmp = this.copyArray(this.E_a);

        for (let i = 0; i < this.nx; i++) {
            for (let j = 0; j < this.ny; j++) {
                if (this.conducting[i][j] === 1) {
                    let F_n_sum = 0;
                    let F_p_sum = 0;
                    let E_n_sum = 0;
                    let E_p_sum = 0;
                    let Ea_sum = 0;
                    let neighbors = 0;

                    this.markNeighborhood(i, j, this.junction_size);

                    for (let di = -this.junction_size; di <= this.junction_size; di++) {
                        for (let dj = -this.junction_size; dj <= this.junction_size; dj++) {
                            if (i + di >= 0 && j + dj >= 0 && i + di < this.nx && j + dj < this.ny &&
                                this.conducting[i + di][j + dj] === 1 && this.visited[i + di][j + dj]) {
                                F_p_sum += F_p_tmp[i + di][j + dj];
                                F_n_sum += F_n_tmp[i + di][j + dj];
                                E_p_sum += E_p_tmp[i + di][j + dj];
                                E_n_sum += E_n_tmp[i + di][j + dj];
                                Ea_sum += Ea_tmp[i + di][j + dj];
                                neighbors += 1;
                                this.distance[i + di][j + dj] = Number.MAX_SAFE_INTEGER;
                                this.visited[i + di][j + dj] = false;
                            }
                        }
                    }
                    this.F0_n[i][j] = F_n_sum / neighbors;
                    this.F0_p[i][j] = F_p_sum / neighbors;
                    this.E0_n[i][j] = E_n_sum / neighbors;
                    this.E0_p[i][j] = E_p_sum / neighbors;
                    this.E_a[i][j] = Ea_sum / neighbors;
                    this.K[i][j] = Math.exp(-this.beta * (this.F0_n[i][j] + this.F0_p[i][j]));
                    this.R[i][j] = this.arrhenius_prefactor * Math.exp(-this.E_a[i][j]);
                } else {
                    this.E_a[i][j] = 0;
                    this.K[i][j] = 0;
                    this.R[i][j] = 0;
                }
            }
        }

        for (let i = 0; i < this.nx - 1; i++) {
            for (let j = 1; j < this.ny - 1; j++) {
                this.cmfx_n[i][j] = -this.conducting_x[i][j] * (this.F0_n[i + 1][j] - this.F0_n[i][j]) / this.ds;
                this.cmfx_p[i][j] = -this.conducting_x[i][j] * (this.F0_p[i + 1][j] - this.F0_p[i][j]) / this.ds;
            }
        }

        for (let i = 1; i < this.nx - 1; i++) {
            for (let j = 0; j < this.ny - 1; j++) {
                this.cmfy_n[i][j] = -this.conducting_y[i][j] * (this.F0_n[i][j + 1] - this.F0_n[i][j]) / this.ds;
                this.cmfy_p[i][j] = -this.conducting_y[i][j] * (this.F0_p[i][j + 1] - this.F0_p[i][j]) / this.ds;
            }
        }
    }

    copyArray(array) {
        const newarray = new Array(array.length);
        for (let i = 0; i < array.length; i++) {
            newarray[i] = array[i].slice();
        }
        return newarray;
    }

    initializeAllMaterials() {
        for (let i = 0; i < this.nx; i++) {
            for (let j = 0; j < this.ny; j++) {
                this.initializeMaterial(i, j);
            }
        }
    }

    initializeMaterial(i, j, material = this.materials[i][j].type) {
        if (i < 0 || j < 0 || i >= this.nx || j >= this.ny || this.materials[i][j].modified) {
            return;
        }

        this.materials[i][j].type = material;

        if (material === MaterialType.DIELECTRIC) {
            this.materials[i][j].eps_r = this.dielectric_eps_r;
        } else if (material === MaterialType.FERROMAGNET) {
            this.materials[i][j].mu_r = this.ferromagnet_mu_r;
        } else if (material === MaterialType.POS_CHARGE) {
            this.materials[i][j].rho_back = this.staticcharge_density;
        } else if (material === MaterialType.NEG_CHARGE) {
            this.materials[i][j].rho_back = -this.staticcharge_density;
        }

        if (MaterialType.isConducting(material)) {
            this.materials[i][j].conducting = 1;
            this.materials[i][j].ni = this.ni_metal;
            this.materials[i][j].W = this.W_metal_default;
            this.materials[i][j].Eb = this.E_b_metal;
            this.materials[i][j].Ea = this.E_a_metal;

            if (material === MaterialType.METAL_HIGH_W) {
                this.materials[i][j].W = this.W_metal_high;
            } else if (material === MaterialType.METAL_LOW_W) {
                this.materials[i][j].W = this.W_metal_low;
            } else if (material === MaterialType.METAL_HIGH_C) {
                this.materials[i][j].ni = 2.5 * this.ni_metal;
            } else if (material === MaterialType.METAL_LOW_C) {
                this.materials[i][j].ni = 0.25 * this.ni_metal;
            }
        }

        if (MaterialType.isSemiconducting(material)) {
            this.materials[i][j].conducting = 1;
            this.materials[i][j].semiconducting = 1;
            this.materials[i][j].ni = this.ni_semi;
            this.materials[i][j].W = this.W_semi;
            this.materials[i][j].Eb = this.E_b_semi;
            this.materials[i][j].Ea = this.E_a_semi;

            if (material === MaterialType.SEMI_P_TYPE) {
                this.materials[i][j].rho_back = -this.p_default_doping_concentration * this.e_charge;
            } else if (material === MaterialType.SEMI_N_TYPE) {
                this.materials[i][j].rho_back = this.n_default_doping_concentration * this.e_charge;
            } else if (material === MaterialType.SEMI_HEAVY_P_TYPE) {
                this.materials[i][j].rho_back = -this.p_heavy_doping_concentration * this.e_charge;
            } else if (material === MaterialType.SEMI_HEAVY_N_TYPE) {
                this.materials[i][j].rho_back = this.n_heavy_doping_concentration * this.e_charge;
            } else if (material === MaterialType.SEMI_LIGHT_P_TYPE) {
                this.materials[i][j].rho_back = -this.p_light_doping_concentration * this.e_charge;
            } else if (material === MaterialType.SEMI_LIGHT_N_TYPE) {
                this.materials[i][j].rho_back = this.n_light_doping_concentration * this.e_charge;
            }
        }
    }

    handleMouseInput() {
        const brush = Object.values(Brush).find(b => b.name === this.gui_brush.value) || Brush.DRAW;
        const brushshape = Object.values(BrushShape).find(bs => bs.name === this.gui_brush_1.value) || BrushShape.CIRCLE;

        let pressing = false;
        let releasing = false;

        if (this.mouse_pressed) {
            if (!this.mouse_pressed_prev) {
                pressing = true;
                this.canvas.focus();
                if (Brush.isMaterialModifyingBrush(brush) || brush === Brush.SELECT) {
                    this.gui_paused.checked = true;
                }
            }
        } else {
            if (this.mouse_pressed_prev) {
                releasing = true;
            }
        }
        this.mouse_pressed_prev = this.mouse_pressed;

        if (!this.mouse_pressed && !releasing) {
            if (this.ctrl_down || this.shift_down) {
                if (!this.modifier_pressed && Brush.isMaterialModifyingBrush(brush)) {
                    this.modifier_pressed = true;
                    this.prev_brush = brush;

                    if (this.ctrl_down) {
                        this.gui_brush.value = Brush.FILL.name;
                    } else if (this.shift_down) {
                        this.gui_brush.value = Brush.LINE.name;
                    }
                }
            } else {
                if (this.modifier_pressed) {
                    this.modifier_pressed = false;
                    this.gui_brush.value = Brush.DRAW.name;
                }
            }
        }

        this.mx_realspace = Math.round((this.mx - 1) / this.scalefactor - 0.5) * this.ds;
        this.my_realspace = Math.round((this.my - 1) / this.scalefactor - 0.5) * this.ds;

        this.mx_start_realspace = Math.round((this.mx_start - 1) / this.scalefactor - 0.5) * this.ds;
        this.my_start_realspace = Math.round((this.my_start - 1) / this.scalefactor - 0.5) * this.ds;

        this.mx_index = Math.round((this.mx - 1) / this.scalefactor - 0.5);
        this.my_index = Math.round((this.my - 1) / this.scalefactor - 0.5);

        this.mx_start_index = Math.round((this.mx_start - 1) / this.scalefactor - 0.5);
        this.my_start_index = Math.round((this.my_start - 1) / this.scalefactor - 0.5);

        if (this.mx_index < 0) this.mx_index = 0;
        if (this.my_index < 0) this.my_index = 0;
        if (this.mx_index >= this.nx) this.mx_index = this.nx - 1;
        if (this.my_index >= this.ny) this.my_index = this.ny - 1;

        if (this.mx_start_index < 0) this.mx_start_index = 0;
        if (this.my_start_index < 0) this.my_start_index = 0;
        if (this.mx_start_index >= this.nx) this.mx_start_index = this.nx - 1;
        if (this.my_start_index >= this.ny) this.my_start_index = this.ny - 1;

        this.gui_stepsizelbl.textContent = `Step size: ${this.getSI(this.dt, 's')}`;
        this.gui_stepslbl.textContent = `Steps/frame: ${this.gui_simspeed_2.value}`;

        this.brushsize = (this.width / 500) * (Math.pow(10.0, this.gui_brushsize.value / 500.0) + this.gui_brushsize.value / 100.0);
        this.gui_lblBrushSize.textContent = `Brush size: ${Math.ceil(this.brushsize / this.ds)}`;

        if (!Brush.isMaterialModifyingBrush(brush)) {
            this.gui_parameter2.style.display = 'none';
            this.gui_parameter2_text.style.display = 'none';
            this.gui_parameter2_text.textContent = '';
        }

        if (Brush.isMaterialModifyingBrush(brush) && brush !== Brush.FILL) {
			this.gui_brush_1.style.display = 'inline';
			this.gui_brush_highlight.style.display = 'inline';
			this.gui_brush_highlight_label.style.display = 'inline';
			this.gui_brushsize.style.display = 'inline';
			this.gui_lblBrushSize.style.display = 'inline';
		} else {
			this.gui_brush_1.style.display = 'none';
			this.gui_brush_highlight.style.display = 'none';
			this.gui_brush_highlight_label.style.display = 'none';
			this.gui_brushsize.style.display = 'none';
			this.gui_lblBrushSize.style.display = 'none';
		}
		

		if (Brush.isMaterialModifyingBrush(brush) && brush !== Brush.ERASE) {
			this.gui_material.style.display = 'inline';
		} else {
			this.gui_material.style.display = 'none';
		}

        if (this.brush_changed) {
            if (!(brush === Brush.SELECT || brush === Brush.FLOODSELECT)) {
                for (let i = 0; i < this.nx; i++) {
                    for (let j = 0; j < this.ny; j++) {
                        this.selected[i][j] = false;
                    }
                }
            }

            if (brush !== Brush.INTERACT) {
                for (let i = 0; i < this.nx; i++) {
                    for (let j = 0; j < this.ny; j++) {
                        this.selected_EMF[i][j] = false;
                    }
                }
                this.EMF_selected = false;
            }

            this.brush_changed = false;
        }

        let update = false;

        if (this.cut || this.copy) {
            let i_min = this.nx - 1;
            let j_min = this.ny - 1;
            for (let i = 0; i < this.nx; i++) {
                for (let j = 0; j < this.ny; j++) {
                    this.clipboard[i][j].erase();
                    if (this.selected[i][j] && this.materials[i][j].type !== MaterialType.VACUUM) {
                        if (i < i_min) i_min = i;
                        if (j < j_min) j_min = j;
                    }
                }
            }
            for (let i = 0; i < this.nx; i++) {
                for (let j = 0; j < this.ny; j++) {
                    if (this.selected[i][j] && this.materials[i][j].type !== MaterialType.VACUUM) {
                        this.clipboard[i - i_min][j - j_min] = this.materials[i][j].clone();
                        if (this.cut) {
                            this.materials[i][j].erase();
                        }
                    }
                    if (this.cut) {
                        this.selected[i][j] = false;
                    }
                }
            }

            if (this.cut) update = true;

            this.cut = false;
            this.copy = false;
        }

        if (this.paste) {
            for (let i = 0; i < this.nx; i++) {
                for (let j = 0; j < this.ny; j++) {
                    this.selected[i][j] = false;
                }
            }

            for (let i = 0; i < this.nx; i++) {
                for (let j = 0; j < this.ny; j++) {
                    this.selection[i][j] = this.clipboard[i][j].clone();
                }
            }
            this.moving_selection = true;
            this.dragging_selection = false;
            this.paste = false;
            update = true;
        }

        if (this.delete) {
            for (let i = 0; i < this.nx; i++) {
                for (let j = 0; j < this.ny; j++) {
                    if (this.selected[i][j] && this.materials[i][j].type !== MaterialType.VACUUM) {
                        this.materials[i][j].erase();
                    }
                }
            }
            this.delete = false;
            update = true;
        }

        switch (brush) {
            case Brush.DRAW:
            case Brush.LINE:
            case Brush.REPLACE:
            case Brush.ERASE:
            case Brush.FILL:
                if ((this.mousebutton === 1 || this.alt_down) && pressing) {
                    this.gui_material.value = this.materials[this.mx_index][this.my_index].type.name;
                }

                let angle = 0;
                let mat = Object.values(MaterialType).find(m => m.name === this.gui_material.value) || MaterialType.VACUUM;

                if (mat === MaterialType.EMF) {
                    this.gui_parameter2.style.display = 'inline';
                    this.gui_parameter2_text.style.display = 'inline';

                    const directionval = Math.floor(this.gui_parameter2.value / 6);
                    if (directionval === 0) {
                        this.gui_parameter2_text.textContent = 'EMF direction: Up';
                        angle = -Math.PI / 2;
                    } else if (directionval === 1) {
                        this.gui_parameter2_text.textContent = 'EMF direction: Right';
                        angle = 0;
                    } else if (directionval === 2) {
                        this.gui_parameter2_text.textContent = 'EMF direction: Down';
                        angle = Math.PI / 2;
                    } else if (directionval === 3) {
                        this.gui_parameter2_text.textContent = 'EMF direction: Left';
                        angle = Math.PI;
                    }
                } else {
                    this.gui_parameter2.style.display = 'none';
                    this.gui_parameter2_text.style.display = 'none';
                    this.gui_parameter2_text.textContent = '';
                }

                if (this.mousebutton === 2 || brush === Brush.ERASE) {
                    mat = MaterialType.VACUUM;
                }

                if (!(this.mousebutton === 1 || this.alt_down)) {
                    if (brush === Brush.LINE) {
                        if (releasing) {
                            this.drawMaterialLine(this.mx_start_realspace, this.my_start_realspace, this.mx_realspace, this.my_realspace, brush, brushshape, mat, this.brushsize, angle);
                        }
                    } else if (brush === Brush.FILL) {
                        if (pressing) {
                            this.floodFillSet(this.mx_index, this.my_index, this.materials[this.mx_index][this.my_index].type, mat, angle);
                        }
                    } else if (this.mouse_pressed) {
                        this.drawMaterialLine(this.mxp_realspace, this.myp_realspace, this.mx_realspace, this.my_realspace, brush, brushshape, mat, this.brushsize, angle);
                    }
                }

                if (Brush.isBrushShapeImportant(brush)) {
                    for (let i = 0; i < this.nx; i++) {
                        for (let j = 0; j < this.ny; j++) {
                            if (this.gui_brush_highlight.checked) {
                                const cx = i * this.ds;
                                const cy = j * this.ds;

                                const px = cx - this.mx_realspace;
                                const py = cy - this.my_realspace;
                                let r = 0;

                                if (brushshape === BrushShape.CIRCLE) {
                                    r = Math.sqrt(px * px + py * py);
                                } else if (brushshape === BrushShape.SQUARE) {
                                    r = Math.max(Math.abs(px), Math.abs(py));
                                }
                                this.under_brush[i][j] = (r <= this.brushsize);
                            } else {
                                this.under_brush[i][j] = false;
                            }
                        }
                    }
                }
                break;

            case Brush.INTERACT:
                if (this.materials[this.mx_index][this.my_index].type === MaterialType.EMF || this.materials[this.mx_index][this.my_index].type === MaterialType.SWITCH) {
                    this.canvas.style.cursor = this.HAND_CURSOR;
                } else {
                    this.canvas.style.cursor = this.DEFAULT_CURSOR;
                }

                if (pressing) {
                    const turn_on_EMF = !this.selected_EMF[this.mx_index][this.my_index];

                    for (let i = 0; i < this.nx; i++) {
                        for (let j = 0; j < this.ny; j++) {
                            this.selected_EMF[i][j] = false;
                        }
                    }
                    this.EMF_selected = false;

                    if (this.materials[this.mx_index][this.my_index].type === MaterialType.EMF && turn_on_EMF) {
                        this.floodFillSelectEMF(this.mx_index, this.my_index, true);
                        this.EMF_selected = true;
                        const setting = Math.round(50 * this.materials[this.mx_index][this.my_index].emf / this.max_EMF);
                        this.gui_parameter3.value = setting;
                        this.prev_EMF_setting = setting;
                    }

                    if (this.materials[this.mx_index][this.my_index].type == MaterialType.SWITCH) {
                        this.floodFillToggleSwitch(this.mx_index, this.my_index, 1 - this.materials[this.mx_index][this.my_index].activated);
                        update = true;
                    }
                }

            case Brush.FLOODSELECT:
            case Brush.SELECT:
                if (pressing) {
                    if (brush === Brush.FLOODSELECT && !this.moving_selection) {
                        this.floodFillSelect(this.mx_index, this.my_index, this.materials[this.mx_index][this.my_index].type, !this.selected[this.mx_index][this.my_index]);
                    } else if (this.moving_selection && !this.dragging_selection) {
                        for (let i = 0; i < this.nx; i++) {
                            for (let j = 0; j < this.ny; j++) {
                                const si = i - this.delta_mx_index;
                                const sj = j - this.delta_my_index;
                                if (si >= 0 && sj >= 0 && si < this.nx && sj < this.ny && this.selection[si][sj].type !== MaterialType.VACUUM) {
                                    this.materials[i][j].erase();
                                    this.materials[i][j] = this.selection[si][sj].clone();
                                    this.selected[i][j] = true;
                                }
                            }
                        }
                        this.moving_selection = false;
                        this.dragging_selection = false;
                    } else if (!this.moving_selection && this.selected[this.mx_index][this.my_index]) {
                        for (let i = 0; i < this.nx; i++) {
                            for (let j = 0; j < this.ny; j++) {
                                this.selection[i][j].erase();
                                if (this.selected[i][j]) {
                                    this.selection[i][j] = this.materials[i][j].clone();
                                    this.selected[i][j] = false;
                                    this.materials[i][j].erase();
                                }
                            }
                        }
                        this.moving_selection = true;
                        this.dragging_selection = true;
                        this.delta_mx_index = 0;
                        this.delta_my_index = 0;
                    }
                } else if (this.mouse_pressed) {
                    if (brush !== Brush.FLOODSELECT) {
                        if (this.dragging_selection) {
                            this.delta_mx_index = this.mx_index - this.mx_start_index;
                            this.delta_my_index = this.my_index - this.my_start_index;
                        } else {
                            const mx0 = Math.min(this.mx_start_index, this.mx_index);
                            const my0 = Math.min(this.my_start_index, this.my_index);
                            const mx1 = Math.max(this.mx_start_index, this.mx_index);
                            const my1 = Math.max(this.my_start_index, this.my_index);

                            for (let i = 0; i < this.nx; i++) {
                                for (let j = 0; j < this.ny; j++) {
                                    this.selected[i][j] = (i >= mx0 && i <= mx1 && j >= my0 && j <= my1);
                                }
                            }
                        }
                    }
                } else if (releasing) {
                    if (brush !== Brush.FLOODSELECT) {
                        this.delta_mx_index = this.mx_index - this.mx_start_index;
                        this.delta_my_index = this.my_index - this.my_start_index;
                        if (this.dragging_selection) {
                            for (let i = 0; i < this.nx; i++) {
                                for (let j = 0; j < this.ny; j++) {
                                    const si = i - this.delta_mx_index;
                                    const sj = j - this.delta_my_index;
                                    if (si >= 0 && sj >= 0 && si < this.nx && sj < this.ny && this.selection[si][sj].type !== MaterialType.VACUUM) {
                                        this.materials[i][j].erase();
                                        this.materials[i][j] = this.selection[si][sj].clone();
                                        this.selected[i][j] = true;
                                    }
                                }
                            }
                            this.moving_selection = false;
                            this.dragging_selection = false;
                        } else {
                            if (this.delta_mx_index === 0 && this.delta_my_index === 0) {
                                for (let i = 0; i < this.nx; i++) {
                                    for (let j = 0; j < this.ny; j++) {
                                        this.selected[i][j] = false;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    this.delta_mx_index = this.mx_index;
                    this.delta_my_index = this.my_index;
                }
                break;

            case Brush.CURRENT:
                if (pressing) {
                    const p = new CurrentProbe();
                    p.x1 = this.mx_start_index;
                    p.y1 = this.my_start_index;
                    p.x2 = this.mx_index;
                    p.y2 = this.my_index;
                    this.currentprobes.push(p);
                } else if (this.mouse_pressed) {
                    const p = this.currentprobes[this.currentprobes.length - 1];
                    p.x2 = this.mx_index;
                    p.y2 = this.my_index;
                }
                break;

            case Brush.VOLTAGE:
                if (pressing) {
                    const p = new VoltageProbe();
                    p.x = this.mx_start_index;
                    p.y = this.my_start_index;
                    this.voltageprobes.push(p);
                } else if (this.mouse_pressed) {
                    const p = this.voltageprobes[this.voltageprobes.length - 1];
                    p.x = this.mx_index;
                    p.y = this.my_index;
                }
                break;

            case Brush.DELETEPROBE:
                this.canvas.style.cursor = this.DEFAULT_CURSOR;
                let i = 0;
                while (i < this.voltageprobes.length) {
                    const p = this.voltageprobes[i];
                    if (this.length(p.x - this.mx_index, p.y - this.my_index) < 3) {
                        this.canvas.style.cursor = this.HAND_CURSOR;
                        if (pressing) {
                            this.voltageprobes.splice(i, 1);
                            i--;
                        }
                    }
                    i++;
                }

                i = 0;
                while (i < this.currentprobes.length) {
                    const p = this.currentprobes[i];
                    if (this.length(p.x1 - this.mx_index, p.y1 - this.my_index) < 3 || this.length(p.x2 - this.mx_index, p.y2 - this.my_index) < 3) {
                        this.canvas.style.cursor = this.HAND_CURSOR;
                        if (pressing) {
                            this.currentprobes.splice(i, 1);
                            i--;
                        }
                    }
                    i++;
                }

                if (this.ground && this.length(this.ground.x - this.mx_index, this.ground.y - this.my_index) < 3) {
                    this.canvas.style.cursor = this.HAND_CURSOR;
                    if (pressing) {
                        this.ground = null;
                    }
                }
                break;

            case Brush.TEXT:
                this.canvas.style.cursor = this.HAND_CURSOR;
                if (this.mouse_pressed) {
                    this.texting = true;
                    this.text_x = this.mx_index;
                    this.text_y = this.my_index;
                }
                break;

            case Brush.GROUND:
                if (pressing) {
                    if (!this.ground) {
                        this.ground = new VoltageProbe();
                    }
                    this.ground.x = this.mx_start_index;
                    this.ground.y = this.my_start_index;
                } else if (this.mouse_pressed) {
                    this.ground.x = this.mx_index;
                    this.ground.y = this.my_index;
                }
                break;
        }

        if (brush !== Brush.TEXT) {
            this.texting = false;
        }

        this.setEMFs();

        const selectedBC = Object.values(BoundaryCondition).find(bc => bc.name === this.gui_bc.value) || BoundaryCondition.DISSIPATIVE;
        if (releasing || selectedBC !== this.prev_boundary || update) {
            this.constructBoundary();
            this.updateAllMaterials();
            this.multigridSolve(true, false);
        }

        this.prev_boundary = selectedBC;

        this.mxp_realspace = this.mx_realspace;
        this.myp_realspace = this.my_realspace;
    }

    drawMaterialLine(x1, y1, x2, y2, brush, brushshape, mat, brushsize, EMF_angle) {
        const a = new Vector(0, 0);
        const b = new Vector(0, 0);
        const p = new Vector(0, 0);
        const ab = new Vector(0, 0);
        for (let i = 0; i < this.nx; i++) {
            for (let j = 0; j < this.ny; j++) {
                const cx = i * this.ds;
                const cy = j * this.ds;

                a.initialize(x1, y1);
                b.initialize(x2, y2);
                p.initialize(cx, cy);
                p.addmult(a, -1);
                ab.copyFrom(b);
                ab.addmult(a, -1);

                let l2 = ab.dot(ab);
                if (l2 === 0) l2 = 1;
                const t = Math.max(0, Math.min(1, p.dot(ab) / l2)); // clamp
                ab.scalarmult(t);
                p.addmult(ab, -1);
                let r = 0;
                if (brushshape === BrushShape.CIRCLE) {
                    r = Math.sqrt(p.dot(p));
                } else if (brushshape === BrushShape.SQUARE) {
                    r = Math.max(Math.abs(p.x), Math.abs(p.y));
                }
                if (r <= brushsize) {
                    if (mat === MaterialType.VACUUM) {
                        this.materials[i][j].erase();
                    } else if (this.materials[i][j].type === MaterialType.VACUUM || brush === Brush.REPLACE) {
                        this.materials[i][j].erase();
                        this.initializeMaterial(i, j, mat);
                        if (mat === MaterialType.EMF) this.materials[i][j].emf_direction = EMF_angle;
                    }
                }
            }
        }
    }

    setEMFs() {
        const EMF_setting = parseInt(this.gui_parameter3.value);
        const new_EMF = this.max_EMF * EMF_setting / 50.0;

        const brush = Object.values(Brush).find(b => b.name === this.gui_brush.value) || Brush.DRAW;
        if (brush === Brush.INTERACT && this.EMF_selected) {
            this.gui_parameter3.style.display = 'inline';
            this.gui_parameter3_text.style.display = 'inline';
            this.gui_parameter3_text.textContent = `EMF: ${this.getSI(new_EMF, 'V/m')}`;
        } else {
            this.gui_parameter3.style.display = 'none';
            this.gui_parameter3_text.style.display = 'none';
            this.gui_parameter3_text.textContent = '';
        }

        if (EMF_setting !== this.prev_EMF_setting && this.EMF_selected) {
            for (let i = 0; i < this.nx; i++) {
                for (let j = 0; j < this.ny; j++) {
                    if (this.materials[i][j].type === MaterialType.EMF && this.selected_EMF[i][j]) {
                        this.materials[i][j].emf = new_EMF;
                    }
                }
            }
            this.updateJustEMFs();
        }

        this.prev_EMF_setting = EMF_setting;
    }

    checkCFL() {
        console.log(`Wave equation CFL ratio = ${(Math.sqrt(2) * this.c) / (this.ds / this.dt_maximum)}`);
        console.log(`Electron diffusion CFL ratio = ${this.dt_maximum / (this.ds * this.ds / (4 * this.D_electron))}`);
        console.log(`Hole diffusion CFL ratio = ${this.dt_maximum / (this.ds * this.ds / (4 * this.D_hole))}`);
    }

    calcMiscFields(updatePhi) {
        const view_scalar = Object.values(ScalarView).find(v => v.name === this.gui_view.value) || ScalarView.POTENTIAL;

        if (updatePhi) {
            this.t7.start();
            this.multigridSolve(false, true);
            this.t7.stop();
        }

        this.t8.start();

        this.sign_violation = false;
        for (let i = 0; i < this.nx; i++) {
            for (let j = 0; j < this.ny; j++) {
                if (this.conducting[i][j] === 0) {
                    this.F_n[i][j] = NaN;
                    this.F_p[i][j] = NaN;
                    this.F[i][j] = NaN;
                } else {
                    this.F_n[i][j] = this.F0_n[i][j] + this.kT * Math.log(this.rho_n[i][j] / this.q_n);
                    this.F_p[i][j] = this.F0_p[i][j] + this.kT * Math.log(this.rho_p[i][j] / this.q_p);

                    const probe_n = -this.e_charge * this.ni_metal;
                    const probe_p = this.e_charge * this.ni_metal;
                    const sigma_n = this.mu_electron * (-this.rho_n[i][j]);
                    const sigma_p = this.mu_hole * this.rho_p[i][j];

                    this.F[i][j] = (sigma_n * (this.F_n[i][j] / this.q_n + this.phi[i][j]) + sigma_p * (this.F_p[i][j] / this.q_p + this.phi[i][j]))
                        / (sigma_n + sigma_p) - this.W_semi / this.eVtoJ;
                }

                this.G[i][j] = this.conducting[i][j] * this.R[i][j] * (this.K[i][j] - this.rho_n[i][j] * this.rho_p[i][j] / (this.q_n * this.q_p));

                if (this.rho_n[i][j] > this.error_detection_threshold || this.rho_p[i][j] < -this.error_detection_threshold) {
                    this.sign_violation = true;
                }
            }
        }

        for (let i = 0; i < this.nx - 1; i++) {
            for (let j = 1; j < this.ny - 1; j++) {
                this.Jx_free[i][j] = this.Jx_abs[i][j] + this.Jx_n[i][j] + this.Jx_p[i][j];
                this.Dx[i][j] = this.epsx[i][j] * this.Ex[i][j];
                this.Sy[i][j] = -this.Ex[i][j] * 0.5 * (this.Hz[i][j] + this.Hz[i][j - 1]);
            }
        }

        for (let i = 1; i < this.nx - 1; i++) {
            for (let j = 0; j < this.ny - 1; j++) {
                this.Jy_free[i][j] = this.Jy_abs[i][j] + this.Jy_n[i][j] + this.Jy_p[i][j];
                this.Dy[i][j] = this.epsy[i][j] * this.Ey[i][j];
                this.Sx[i][j] = this.Ey[i][j] * 0.5 * (this.Hz[i][j] + this.Hz[i - 1][j]);
            }
        }

        for (let i = 1; i < this.nx - 2; i++) {
            for (let j = 1; j < this.ny - 2; j++) {
                this.Bz[i][j] = this.Hz[i][j] * this.mu_z[i][j];
                this.u[i][j] = 0.25 * (this.Ex[i][j] * this.Dx[i][j] + this.Ex[i][j + 1] * this.Dx[i][j + 1])
                    + 0.25 * (this.Ey[i][j] * this.Dy[i][j] + this.Ey[i + 1][j] * this.Dy[i + 1][j])
                    + 0.5 * this.Hz[i][j] * this.Bz[i][j];
            }
        }

        if (view_scalar === ScalarView.HEAT) {
            for (let i = 0; i < this.nx - 1; i++) {
                for (let j = 1; j < this.ny - 1; j++) {
                    this.grad_E0x_n[i][j] = -this.conducting_x[i][j] * (this.E0_n[i + 1][j] - this.E0_n[i][j]) / this.ds;
                    this.grad_E0x_p[i][j] = -this.conducting_x[i][j] * (this.E0_p[i + 1][j] - this.E0_p[i][j]) / this.ds;
                }
            }

            for (let i = 1; i < this.nx - 1; i++) {
                for (let j = 0; j < this.ny - 1; j++) {
                    this.grad_E0y_n[i][j] = -this.conducting_y[i][j] * (this.E0_n[i][j + 1] - this.E0_n[i][j]) / this.ds;
                    this.grad_E0y_p[i][j] = -this.conducting_y[i][j] * (this.E0_p[i][j + 1] - this.E0_p[i][j]) / this.ds;
                }
            }

            for (let i = 1; i < this.nx - 1; i++) {
                for (let j = 1; j < this.ny - 1; j++) {
                    const n_contrib = -this.G[i][j] * this.E0_n[i][j];
                    const jn_contrib = 0.5 * (
                        this.grad_E0x_n[i - 1][j] * this.Jx_n[i - 1][j] +
                        this.grad_E0x_n[i][j] * this.Jx_n[i][j] +
                        this.grad_E0y_n[i][j - 1] * this.Jy_n[i][j - 1] +
                        this.grad_E0y_n[i][j] * this.Jy_n[i][j]
                    ) / this.q_n;

                    const p_contrib = -this.G[i][j] * this.E0_p[i][j];
                    const jp_contrib = 0.5 * (
                        this.grad_E0x_p[i - 1][j] * this.Jx_p[i - 1][j] +
                        this.grad_E0x_p[i][j] * this.Jx_p[i][j] +
                        this.grad_E0y_p[i][j - 1] * this.Jy_p[i][j - 1] +
                        this.grad_E0y_p[i][j] * this.Jy_p[i][j]
                    ) / this.q_p;

                    const ohm_contrib = 0.5 * (
                        this.Ex[i - 1][j] * this.Jx_free[i - 1][j] +
                        this.Ex[i][j] * this.Jx_free[i][j] +
                        this.Ey[i][j - 1] * this.Jy_free[i][j - 1] +
                        this.Ey[i][j] * this.Jy_free[i][j]
                    );

                    this.Q[i][j] = n_contrib + jn_contrib + p_contrib + jp_contrib + ohm_contrib;
                }
            }
        }

        if (view_scalar === ScalarView.ENTROPY) {
            for (let i = 0; i < this.nx - 1; i++) {
                for (let j = 1; j < this.ny - 1; j++) {
                    this.grad_Fx_n[i][j] = this.conducting_x[i][j] === 1 ? -(this.F_n[i + 1][j] - this.F_n[i][j]) / this.ds : 0;
                    this.grad_Fx_p[i][j] = this.conducting_x[i][j] === 1 ? -(this.F_p[i + 1][j] - this.F_p[i][j]) / this.ds : 0;
                }
            }

            for (let i = 1; i < this.nx - 1; i++) {
                for (let j = 0; j < this.ny - 1; j++) {
                    this.grad_Fy_n[i][j] = this.conducting_y[i][j] === 1 ? -(this.F_n[i][j + 1] - this.F_n[i][j]) / this.ds : 0;
                    this.grad_Fy_p[i][j] = this.conducting_y[i][j] === 1 ? -(this.F_p[i][j + 1] - this.F_p[i][j]) / this.ds : 0;
                }
            }

            for (let i = 1; i < this.nx - 1; i++) {
                for (let j = 1; j < this.ny - 1; j++) {
                    const n_contrib = -this.G[i][j] * this.F_n[i][j];
                    const jn_contrib = 0.5 * (
                        this.grad_Fx_n[i - 1][j] * this.Jx_n[i - 1][j] +
                        this.grad_Fx_n[i][j] * this.Jx_n[i][j] +
                        this.grad_Fy_n[i][j - 1] * this.Jy_n[i][j - 1] +
                        this.grad_Fy_n[i][j] * this.Jy_n[i][j]
                    ) / this.q_n;

                    const p_contrib = -this.G[i][j] * this.F_p[i][j];
                    const jp_contrib = 0.5 * (
                        this.grad_Fx_p[i - 1][j] * this.Jx_p[i - 1][j] +
                        this.grad_Fx_p[i][j] * this.Jx_p[i][j] +
                        this.grad_Fy_p[i][j - 1] * this.Jy_p[i][j - 1] +
                        this.grad_Fy_p[i][j] * this.Jy_p[i][j]
                    ) / this.q_p;

                    const ohm_contrib = 0.5 * (
                        this.Ex[i - 1][j] * this.Jx_free[i - 1][j] +
                        this.Ex[i][j] * this.Jx_free[i][j] +
                        this.Ey[i][j - 1] * this.Jy_free[i][j - 1] +
                        this.Ey[i][j] * this.Jy_free[i][j]
                    );

                    this.S[i][j] = n_contrib + jn_contrib + p_contrib + jp_contrib + ohm_contrib;
                }
            }
        }

        for (const p of this.voltageprobes) {
            p.potential = this.F[p.x][p.y];
        }

        if (this.ground) {
            this.ground.potential = this.F[this.ground.x][this.ground.y];
        }

        for (const p of this.currentprobes) {
            p.current = this.calcCurrent(p.x1, p.y1, p.x2, p.y2);
        }

        this.updateMiscFields = false;

        this.t8.stop();
    }

    calcCurrent(x0, y0, x1, y1) {
        let dy = y1 - y0;
        let dx = x1 - x0;
        let t = 0.5;
        let J = 0;

        if (Math.abs(dx) > Math.abs(dy)) {
            const m = dy / dx;
            t += y0;
            dx = dx < 0 ? -1 : 1;
            const m_dx = m * dx;
            while (x0 !== x1) {
                const x0_prev = x0;
                const t_prev = t;

                x0 += dx;
                t += m_dx;

                J += this.accumCurrent(x0_prev, Math.floor(t_prev), x0, Math.floor(t), this.Jx_free, this.Jy_free);
            }
        } else {
            const m = dx / dy;
            t += x0;
            dy = dy < 0 ? -1 : 1;
            const m_dy = m * dy;
            while (y0 !== y1) {
                const y0_prev = y0;
                const t_prev = t;

                y0 += dy;
                t += m_dy;

                J += this.accumCurrent(Math.floor(t_prev), y0_prev, Math.floor(t), y0, this.Jx_free, this.Jy_free);
            }
        }

        return J;
    }

    accumCurrent(x0, y0, x1, y1, Jx, Jy) {
        const dx = x1 - x0;
        const dy = y1 - y0;
        if (dx === 1 && dy === 0) {
            return Math.abs(Jy[x0 + 1][y0]) * this.ds;
        } else if (dx === -1 && dy === 0) {
            return Math.abs(Jy[x0][y0]) * this.ds;
        } else if (dx === 0 && dy === 1) {
            return Math.abs(Jx[x0][y0 + 1]) * this.ds;
        } else if (dx === 0 && dy === -1) {
            return Math.abs(Jx[x0][y0]) * this.ds;
        } else if (dx === 1 && dy === 1) {
            return (Math.abs(Jx[x0][y0 + 1]) + Math.abs(Jy[x0 + 1][y0 + 1])) * this.ds;
        } else if (dx === 1 && dy === -1) {
            return (Math.abs(Jx[x0][y0]) + Math.abs(Jy[x0 + 1][y0 - 1])) * this.ds;
        } else if (dx === -1 && dy === 1) {
            return (Math.abs(Jx[x0][y0 + 1]) + Math.abs(Jy[x0][y0 + 1])) * this.ds;
        } else if (dx === -1 && dy === -1) {
            return (Math.abs(Jx[x0][y0]) + Math.abs(Jy[x0][y0 - 1])) * this.ds;
        }
        return 0;
    }

    prescaleDielectric() {
        this.downscale_x_vector(this.epsx, this.MG_epsx, this.log2_resolution, this.nx, this.ny);
        this.downscale_y_vector(this.epsy, this.MG_epsy, this.log2_resolution, this.nx, this.ny);
    }

    bilinearinterp(array, x, y) {
        let xfloor = Math.floor(x);
        let yfloor = Math.floor(y);
        let fx = x - xfloor;
        let fy = y - yfloor;
        if (Math.abs(x - Math.round(x)) < 1e-2 || Math.abs(y - Math.round(y)) < 1e-2) {
            let i = Math.round(x);
            let j = Math.round(y);
            if (i < 0) i = 0;
            if (j < 0) j = 0;
            if (i >= this.nx) i = this.nx - 1;
            if (j >= this.ny) j = this.ny - 1;
            return array[i][j];
        }

        if (xfloor < 0) {
            xfloor = 0;
            fx = 0.0;
        } else if (xfloor >= this.nx - 1) {
            xfloor = this.nx - 2;
            fx = 1.0;
        }
        if (yfloor < 0) {
            yfloor = 0;
            fy = 0.0;
        } else if (yfloor >= this.ny - 1) {
            yfloor = this.ny - 2;
            fy = 1.0;
        }
        const va = array[xfloor][yfloor] * (1.0 - fx) + array[xfloor + 1][yfloor] * fx;
        const vb = array[xfloor][yfloor + 1] * (1.0 - fx) + array[xfloor + 1][yfloor + 1] * fx;

        return va * (1.0 - fy) + vb * fy;
    }

    downscale_x_vector(source, dest, steps, nx, ny) {
        // Copy initial data
        for (let i = 0; i < nx - 1; i++) {
            for (let j = 0; j < ny; j++) {
                dest[steps][i][j] = source[i][j];
            }
        }
    
        let nx_d = nx;
        let ny_d = ny;
    
        for (let k = steps - 1; k >= 0; k--) {
            nx_d = Math.floor(nx_d / 2);
            ny_d = Math.floor(ny_d / 2);
    
            for (let i = 0; i < nx_d - 1; i++) {
                for (let j = 0; j < ny_d; j++) {
                    dest[k][i][j] = 0.125 * (
                        dest[k + 1][2 * i][2 * j] + dest[k + 1][2 * i][2 * j + 1] +
                        2 * dest[k + 1][2 * i + 1][2 * j] + 2 * dest[k + 1][2 * i + 1][2 * j + 1] +
                        dest[k + 1][2 * i + 2][2 * j] + dest[k + 1][2 * i + 2][2 * j + 1]
                    );
                }
            }
        }
    }
    
    downscale_y_vector(source, dest, steps, nx, ny) {
        // Copy initial data
        for (let i = 0; i < nx; i++) {
            for (let j = 0; j < ny - 1; j++) {
                dest[steps][i][j] = source[i][j];
            }
        }
    
        let nx_d = nx;
        let ny_d = ny;
    
        for (let k = steps - 1; k >= 0; k--) {
            nx_d = Math.floor(nx_d / 2);
            ny_d = Math.floor(ny_d / 2);
    
            for (let i = 0; i < nx_d; i++) {
                for (let j = 0; j < ny_d - 1; j++) {
                    dest[k][i][j] = 0.125 * (
                        dest[k + 1][2 * i][2 * j] + dest[k + 1][2 * i + 1][2 * j] +
                        2 * dest[k + 1][2 * i][2 * j + 1] + 2 * dest[k + 1][2 * i + 1][2 * j + 1] +
                        dest[k + 1][2 * i][2 * j + 2] + dest[k + 1][2 * i + 1][2 * j + 2]
                    );
                }
            }
        }
    }

    floodFillSet(i, j, old_mat, new_mat, EMF_angle) {
        if (old_mat === new_mat) return;

        const queue = [];
        queue.push(new FloodFillCoordinate(i, j));

        while (queue.length > 0) {
            const coord = queue.shift();
            if (coord.i >= 0 && coord.i < this.nx && coord.j >= 0 && coord.j < this.ny &&
                this.materials[coord.i][coord.j].type === old_mat &&
                this.materials[coord.i][coord.j].type !== new_mat) {
                this.materials[coord.i][coord.j].erase();
                this.initializeMaterial(coord.i, coord.j, new_mat);
                if (new_mat === MaterialType.EMF) this.materials[coord.i][coord.j].emf_direction = EMF_angle;
                queue.push(new FloodFillCoordinate(coord.i - 1, coord.j));
                queue.push(new FloodFillCoordinate(coord.i + 1, coord.j));
                queue.push(new FloodFillCoordinate(coord.i, coord.j - 1));
                queue.push(new FloodFillCoordinate(coord.i, coord.j + 1));
            }
        }
    }

    floodFillSelect(i, j, mat, select) {
        const queue = [];
        queue.push(new FloodFillCoordinate(i, j));

        while (queue.length > 0) {
            const coord = queue.shift();
            if (coord.i >= 0 && coord.i < this.nx && coord.j >= 0 && coord.j < this.ny &&
                this.materials[coord.i][coord.j].type === mat &&
                this.selected[coord.i][coord.j] !== select) {
                this.selected[coord.i][coord.j] = select;
                queue.push(new FloodFillCoordinate(coord.i - 1, coord.j));
                queue.push(new FloodFillCoordinate(coord.i + 1, coord.j));
                queue.push(new FloodFillCoordinate(coord.i, coord.j - 1));
                queue.push(new FloodFillCoordinate(coord.i, coord.j + 1));
            }
        }
    }

    floodFillSelectEMF(i, j, select) {
        const queue = [];
        queue.push(new FloodFillCoordinate(i, j));

        while (queue.length > 0) {
            const coord = queue.shift();
            if (coord.i >= 0 && coord.i < this.nx && coord.j >= 0 && coord.j < this.ny &&
                this.materials[coord.i][coord.j].type === MaterialType.EMF &&
                this.selected_EMF[coord.i][coord.j] !== select) {
                this.selected_EMF[coord.i][coord.j] = select;
                queue.push(new FloodFillCoordinate(coord.i - 1, coord.j));
                queue.push(new FloodFillCoordinate(coord.i + 1, coord.j));
                queue.push(new FloodFillCoordinate(coord.i, coord.j - 1));
                queue.push(new FloodFillCoordinate(coord.i, coord.j + 1));
            }
        }
    }

    floodFillToggleSwitch(i, j, active) {
        const queue = [];
        queue.push(new FloodFillCoordinate(i, j));

        while (queue.length > 0) {
            const coord = queue.shift();
            if (coord.i >= 0 && coord.i < this.nx && coord.j >= 0 && coord.j < this.ny &&
                this.materials[coord.i][coord.j].type === MaterialType.SWITCH &&
                this.materials[coord.i][coord.j].activated !== active) {
                this.materials[coord.i][coord.j].activated = active;
                queue.push(new FloodFillCoordinate(coord.i - 1, coord.j));
                queue.push(new FloodFillCoordinate(coord.i + 1, coord.j));
                queue.push(new FloodFillCoordinate(coord.i, coord.j - 1));
                queue.push(new FloodFillCoordinate(coord.i, coord.j + 1));
            }
        }
    }

    drawPixelRectangle(x, y, w, h) {
        for (let i = x; i < x + w; i++) {
            for (let j = y; j < y + h; j++) {
                this.setPixel(i, j);
            }
        }
    }

    setPixel(i, j) {
        if (i < 0 || j < 0 || i >= this.nx || j >= this.ny) return;

        this.image_r[i][j] = this.image_r[i][j] * this.alphaBG + this.col_r * this.alphaFG;
        this.image_g[i][j] = this.image_g[i][j] * this.alphaBG + this.col_g * this.alphaFG;
        this.image_b[i][j] = this.image_b[i][j] * this.alphaBG + this.col_b * this.alphaFG;
    }

    stampPixelData() {
        this.t9.start();
        const scansize = this.nx * this.scalefactor;
        const imageData = this.ctx.createImageData(scansize, this.ny * this.scalefactor);
        const data = imageData.data;

        for (let j = 0; j < this.ny; j++) {
            for (let i = 0; i < this.nx; i++) {
                const scale = 1.0 / Math.max(this.image_r[i][j], this.image_g[i][j], this.image_b[i][j], 1.0);
                const r = this.clamp(Math.floor(256 * this.image_r[i][j] * scale), 0, 255);
                const g = this.clamp(Math.floor(256 * this.image_g[i][j] * scale), 0, 255);
                const b = this.clamp(Math.floor(256 * this.image_b[i][j] * scale), 0, 255);
    
                for (let dy = 0; dy < this.scalefactor; dy++) {
                    for (let dx = 0; dx < this.scalefactor; dx++) {
                        const x = i * this.scalefactor + dx;
                        const y = j * this.scalefactor + dy;
                        const idx = (x + y * scansize) * 4;
                        data[idx] = r;
                        data[idx + 1] = g;
                        data[idx + 2] = b;
                        data[idx + 3] = 255;
                    }
                }
            }
        }

        this.ctx.putImageData(imageData, 0, 0);
        this.t9.stop();
    }

    drawLine(x0, y0, x1, y1, draw_starting_point) {
        const ctx = this.ctx;
        const scansize = this.scalefactor * this.nx;
        const scanheight = this.scalefactor * this.ny;

        // Verify coordinates are within bounds
        if (x0 < 0 || x0 >= scansize || y0 < 0 || y0 >= scanheight ||
            x1 < 0 || x1 >= scansize || y1 < 0 || y1 >= scanheight) {
            return; // Skip invalid lines
        }

        // Map scaled coordinates to canvas pixels (assuming canvas.width = scansize)
        const px0 = x0;
        const py0 = y0;
        const px1 = x1;
        const py1 = y1;

        // Draw line with alpha blending
        ctx.strokeStyle = `rgb(${Math.floor(this.col_r * 255)}, ${Math.floor(this.col_g * 255)}, ${Math.floor(this.col_b * 255)})`;
        ctx.globalAlpha = this.alphaFG;
        ctx.globalCompositeOperation = 'source-over'; // Mimics alphaBG blending
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(px0, py0);
        ctx.lineTo(px1, py1);
        ctx.stroke();

        // Draw starting point
        if (draw_starting_point) {
            ctx.fillStyle = `rgb(${Math.floor(this.col_r * 255)}, ${Math.floor(this.col_g * 255)}, ${Math.floor(this.col_b * 255)})`;
            ctx.beginPath();
            ctx.arc(px0, py0, 2, 0, 2 * Math.PI);
            ctx.fill();
        }

        ctx.globalAlpha = 1.0; // Reset alpha
        ctx.globalCompositeOperation = 'source-over'; // Reset composite mode
    }

    drawPixelLine(x0, y0, x1, y1) {
        let dy = y1 - y0;
        let dx = x1 - x0;
        let t = 0.5;

        this.setPixel(x0, y0);

        if (Math.abs(dx) > Math.abs(dy)) {
            const m = dy / dx;
            t += y0;
            dx = dx < 0 ? -1 : 1;
            const m_dx = m * dx;
            while (x0 !== x1) {
                x0 += dx;
                t += m_dx;
                this.setPixel(x0, Math.floor(t));
            }
        } else {
            const m = dx / dy;
            t += x0;
            dy = dy < 0 ? -1 : 1;
            const m_dy = m * dy;
            while (y0 !== y1) {
                y0 += dy;
                t += m_dy;
                this.setPixel(Math.floor(t), y0);
            }
        }
    }

    setColor(r, g, b) {
        this.col_r = r / 255;
        this.col_g = g / 255;
        this.col_b = b / 255;
    }

    setColorFloat(r, g, b) {
        if (!Number.isFinite(r + g + b)) {
            this.col_r = 0;
            this.col_g = 0;
            this.col_b = 0;
            return;
        }
        const scale = 1.0 / Math.max(r, g, b, 1.0);
        this.col_r = Math.max(r * scale, 0);
        this.col_g = Math.max(g * scale, 0);
        this.col_b = Math.max(b * scale, 0);
    }

    max(x, y, z, w) {
        if (w === undefined) {
            return Math.max(Math.max(x, y), z);
        }
        return Math.max(Math.max(Math.max(x, y), z), w);
    }

    min(x, y, z, w) {
        if (w === undefined) {
            return Math.min(Math.min(x, y), z);
        }
        return Math.min(Math.min(Math.min(x, y), z), w);
    }

    clamp(val, min, max) {
        if (typeof val === 'number' && isNaN(val)) return min;
        return Math.max(min, Math.min(max, val));
    }

    setalphaBG(alpha) {
        this.alphaBG = alpha;
    }

    setalphaFG(alpha) {
        this.alphaFG = alpha;
    }

    render() {
        this.t5.start();

        //this.ctx.imageSmoothingEnabled = false;

        const scalingconstant = 10.0 * Math.pow(10.0, parseInt(this.gui_brightness.value) / 10.0);

        // Clear image buffers
        for (let i = 0; i < this.nx; i++) {
            for (let j = 0; j < this.ny; j++) {
                this.image_r[i][j] = 0;
                this.image_g[i][j] = 0;
                this.image_b[i][j] = 0;
            }
        }

        this.clearStrings();

        // Draw pixels
        this.setalphaBG(0);
        this.setalphaFG(1);

        const brush = Object.values(Brush).find(b => b.name === this.gui_brush.value) || Brush.DRAW;

        if (this.gui_elem_colors.checked) {
            for (let i = 0; i < this.nx; i++) {
                for (let j = 0; j < this.ny; j++) {
                    this.setColor(
                        this.materials[i][j].type.color_r,
                        this.materials[i][j].type.color_g,
                        this.materials[i][j].type.color_b
                    );
                    this.setPixel(i, j);
                }
            }
        } else {
            for (let i = 0; i < this.nx; i++) {
                for (let j = 0; j < this.ny; j++) {
                    const gray = this.materials[i][j].type.color_grayscale;
                    this.setColor(gray, gray, gray);
                    this.setPixel(i, j);
                }
            }
        }

        if (this.moving_selection) {
            if (this.gui_elem_colors.checked) {
                for (let i = 0; i < this.nx; i++) {
                    for (let j = 0; j < this.ny; j++) {
                        const si = i - this.delta_mx_index;
                        const sj = j - this.delta_my_index;
                        if (si >= 0 && sj >= 0 && si < this.nx && sj < this.ny &&
                            this.selection[si][sj].type !== MaterialType.VACUUM) {
                            this.setColor(
                                this.selection[si][sj].type.color_r,
                                this.selection[si][sj].type.color_g,
                                this.selection[si][sj].type.color_b
                            );
                            this.setPixel(i, j);
                        }
                    }
                }
            } else {
                for (let i = 0; i < this.nx; i++) {
                    for (let j = 0; j < this.ny; j++) {
                        const si = i - this.delta_mx_index;
                        const sj = j - this.delta_my_index;
                        if (si >= 0 && sj >= 0 && si < this.nx && sj < this.ny &&
                            this.selection[si][sj].type !== MaterialType.VACUUM) {
                            const gray952 = this.selection[si][sj].type.color_grayscale;
                            this.setColor(gray952, gray952, gray952);
                            this.setPixel(i, j);
                        }
                    }
                }
            }
        }

        const scalarView = Object.values(ScalarView).find(v => v.name === this.gui_view.value) || ScalarView.NONE;
        if (scalarView !== ScalarView.NONE) {
            this.setalphaBG(1.0);
            this.setalphaFG(1.0);
            for (let i = 1; i < this.nx - 1; i++) {
                for (let j = 1; j < this.ny - 1; j++) {
                    let v, vx, vy, rc, bc, gc, it;
                    switch (scalarView) {
                        case ScalarView.NONE:
                            this.setColorFloat(0, 0, 0);
                            break;
                        case ScalarView.B_FIELD:
                            v = this.parity * 0.25 * (
                                this.Bz[i][j] + this.Bz[i - 1][j] +
                                this.Bz[i][j - 1] + this.Bz[i - 1][j - 1]
                            ) * 1e5 * scalingconstant;
                            this.setColorFloat(v, Math.abs(v), -v);
                            break;
                        case ScalarView.E_FIELD:
                            vx = 0.5 * (this.Ex[i][j] + this.Ex[i][j + 1]) * scalingconstant / 1e5;
                            vy = 0.5 * (this.Ey[i][j] + this.Ey[i + 1][j]) * scalingconstant / 1e5;
                            this.setColorFloat(0, this.length(vx, vy), 0);
                            break;
                        case ScalarView.H_FIELD:
                            v = this.parity * this.Hz[i][j] * this.mu0 * 1e5 * scalingconstant;
                            this.setColorFloat(v, Math.abs(v), -v);
                            break;
                        case ScalarView.CURRENT:
                            vx = 0.5 * (this.Jx_free[i][j] + this.Jx_free[i - 1][j]) * scalingconstant;
                            vy = 0.5 * (this.Jy_free[i][j] + this.Jy_free[i][j - 1]) * scalingconstant;
                            this.setColorFloat(0, this.length(vx, vy) / 1e7, 0);
                            break;
                        case ScalarView.POTENTIAL:
                            v = this.phi[i][j] * scalingconstant;
                            this.setColorFloat(v, 0, -v);
                            break;
                        case ScalarView.CHARGE:
                            v = this.rho_free[i][j] * scalingconstant;
                            rc = Math.min(Math.max(v, 0), 1);
                            bc = Math.min(Math.max(-v, 0), 1);
                            gc = Math.min(rc, bc);
                            this.setalphaBG(1 - 0.5 * gc);
                            this.setColorFloat(rc, gc, bc);
                            break;
                        case ScalarView.BACKGROUND_CHARGE:
                            v = this.rho_back[i][j] * scalingconstant;
                            this.setColorFloat(v, 0, -v);
                            break;
                        case ScalarView.ELECTRON_CHARGE:
                            v = this.rho_n[i][j] * scalingconstant;
                            this.setColorFloat(v, 0, -v);
                            break;
                        case ScalarView.HOLE_CHARGE:
                            v = this.rho_p[i][j] * scalingconstant;
                            this.setColorFloat(v, 0, -v);
                            break;
                        case ScalarView.COMBINED_CHARGE:
                            rc = this.clamp(0.2 * Math.log(this.rho_p[i][j] * scalingconstant), 0, 1);
                            bc = this.clamp(0.2 * Math.log(-this.rho_n[i][j] * scalingconstant), 0, 1);
                            gc = Math.min(rc, bc);
                            it = Math.max(rc, bc);
                            this.setalphaBG(1 - 0.5 * gc);
                            this.setalphaFG(0.6 * it);
                            this.setColorFloat(rc, gc, bc);
                            break;
                        case ScalarView.ENERGY:
                            v = this.u[i][j] * scalingconstant;
                            this.setColorFloat(0, v, 0);
                            break;
                        case ScalarView.ENTROPY:
                            v = this.S[i][j] * scalingconstant / 1e12;
                            this.setColorFloat(v, Math.abs(v), -v);
                            break;
                        case ScalarView.HEAT:
                            v = this.Q[i][j] * scalingconstant / 1e12;
                            this.setColorFloat(v, Math.abs(v), -v);
                            break;
                        case ScalarView.ELECTRON_POTENTIAL:
                            v = this.conducting[i][j] * (
                                this.F_n[i][j] / this.q_n + this.phi[i][j] - this.W_semi / this.eVtoJ
                            ) * scalingconstant;
                            this.setColorFloat(v, Math.abs(v), -v);
                            break;
                        case ScalarView.HOLE_POTENTIAL:
                            v = this.conducting[i][j] * (
                                this.F_p[i][j] / this.q_p + this.phi[i][j] - this.W_semi / this.eVtoJ
                            ) * scalingconstant;
                            this.setColorFloat(v, Math.abs(v), -v);
                            break;
                        case ScalarView.DEBUG:
                            v = (this.visited[i][j] ? 1 : 0) * scalingconstant;
                            this.setColorFloat(v, 0, -v);
                            break;
                        case ScalarView.RECOMBINATION:
                            v = -this.G[i][j] / 1e30 * scalingconstant;
                            this.setColorFloat(v, Math.abs(v), -v);
                            break;
                        case ScalarView.AVERAGE_POTENTIAL:
                            v = this.F[i][j] * scalingconstant;
                            this.setColorFloat(v, Math.abs(v), -v);
                            break;
                        case ScalarView.LIGHT:
                            v = -this.materials[i][j].semiconducting * this.G[i][j] / 1e30 * scalingconstant;
                            this.setColorFloat(v, v, v);
                            break;
                    }
                    this.setPixel(i, j);
                }
            }
        }

        const highlight = this.gui_brush_highlight.checked && Brush.isBrushShapeImportant(brush);
        for (let i = 0; i < this.nx; i++) {
            for (let j = 0; j < this.ny; j++) {
                if (this.materials[i][j].type === MaterialType.EMF && i > 0 && j > 0 && i < this.nx - 1 && j < this.ny - 1) {
                    this.setalphaBG(0.25);
                    this.setalphaFG(0.75);

                    let offset = 0;
                    if (this.selected_EMF[i][j]) {
                        offset = 60 * (2 * ((i + j) % 2) - 1);
                    }
                    if (this.materials[i + 1][j].type !== MaterialType.EMF ||
                        this.materials[i - 1][j].type !== MaterialType.EMF ||
                        this.materials[i][j + 1].type !== MaterialType.EMF ||
                        this.materials[i][j - 1].type !== MaterialType.EMF) {
                        offset = -30;
                    }

                    const delta_r = MaterialType.EMF.color_r + offset;
                    const delta_g = MaterialType.EMF.color_g + offset;
                    const delta_b = MaterialType.EMF.color_b + offset;
                    this.setColor(delta_r, delta_g, delta_b);
                    this.setPixel(i, j);
                } else if (!this.materials[i][j].activated) {
                    this.setalphaBG(0.25);
                    this.setalphaFG(0.75);
                    this.setColor(0, 0, 0);
                    this.setPixel(i, j);
                }

                if (this.selected[i][j] || (highlight && this.under_brush[i][j])) {
                    this.setalphaBG(0.75);
                    this.setalphaFG(0.25);

                    const s = this.selected[i][j] ? 1 : 0;
                    const h = highlight && this.under_brush[i][j] ? 1 : 0;
                    const delta_r = 256 * s + 256 * h;
                    const delta_g = 100 * s + 256 * h;
                    const delta_b = 256 * s + 256 * h;

                    this.setColor(delta_r, delta_g, delta_b);
                    this.setPixel(i, j);
                }
            }
        }

        this.setalphaBG(1);
        this.setalphaFG(1);
        this.setColorFloat(0.7, 0.7, 0.7);

        if (brush === Brush.LINE && this.mouse_pressed) {
            this.drawPixelLine(this.mx_start_index, this.my_start_index, this.mx_index, this.my_index);
        }

        this.setalphaFG(1.0);
        this.setColorFloat(0.5, 1.0, 1.0);

        for (const p of this.currentprobes) {
            this.drawPixelRectangle(p.x1 - 1, p.y1 - 1, 3, 3);
            this.drawPixelRectangle(p.x2 - 1, p.y2 - 1, 3, 3);
        }

        for (const p of this.voltageprobes) {
            this.drawPixelRectangle(p.x - 1, p.y - 1, 3, 3);
        }

        if (this.ground) {
            this.drawPixelRectangle(this.ground.x - 1, this.ground.y - 1, 3, 3);
        }

        this.setalphaFG(0.3);
        this.setColorFloat(0.5, 1.0, 1.0);

        for (const p of this.currentprobes) {
            this.drawPixelLine(p.x1, p.y1, p.x2, p.y2);
        }

        this.setalphaFG(0.8);
        this.setColorFloat(1.0, 1.0, 1.0);

        if (this.texting) {
            this.drawPixelLine(this.text_x, this.text_y, this.text_x, this.text_y + 7);
        }

        this.stampPixelData();

        // Draw vectors
        let vectorView = Object.values(VectorView).find(v => v.name === this.gui_view_vec.value) || VectorView.NONE;
        if (vectorView !== VectorView.NONE) {
            this.setalphaBG(1.0);

            const arrowlength = 10.0 / this.scalefactor;
            const vectorscalingconstant = 0.01 * Math.pow(10.0, parseInt(this.gui_brightness_vec.value) / 5.0);
            const vector_display_mode = Object.values(VectorMode).find(vm => vm.name === this.gui_view_vec_mode.value) || VectorMode.ARROWS;

            this.rand.setSeed(4);

            const density = 75;
            let randomness = 0;
            if (vector_display_mode === VectorMode.ARROWS) {
                randomness = 0.5;
            } else if (vector_display_mode === VectorMode.LINES) {
                randomness = 0.75;
            }

            let vf_x = null;
            let vf_y = null;
            switch (vectorView) {
                case VectorView.D_FIELD:
                    vf_x = this.Dx;
                    vf_y = this.Dy;
                    break;
                case VectorView.E_FIELD:
                    vf_x = this.Ex;
                    vf_y = this.Ey;
                    break;
                case VectorView.ELECTRON_CURRENT:
                    vf_x = this.Jx_n;
                    vf_y = this.Jy_n;
                    break;
                case VectorView.HOLE_CURRENT:
                    vf_x = this.Jx_p;
                    vf_y = this.Jy_p;
                    break;
                case VectorView.TOTAL_CURRENT:
                    vf_x = this.Jx_free;
                    vf_y = this.Jy_free;
                    break;
                case VectorView.POYNTING:
                    vf_x = this.Sx;
                    vf_y = this.Sy;
                    break;
                case VectorView.EMF:
                    vf_x = this.emfx;
                    vf_y = this.emfy;
                    break;
            }

            const ctr = new Vector(0, 0);
            const arrow = new Vector(0, 0);
            const tip1 = new Vector(0, 0);
            const tip2 = new Vector(0, 0);
            const body1 = new Vector(0, 0);
            const body2 = new Vector(0, 0);

            for (let i = 0; i < density; i++) {
                for (let j = 0; j < density; j++) {
                    const x = this.nx * (i + randomness * (this.rand.nextFloat() - 0.5)) / density;
                    const y = this.ny * (j + randomness * (this.rand.nextFloat() - 0.5)) / density;

                    if (vector_display_mode === VectorMode.LINES && vf_x !== null) {
                        for (let sign = -1; sign <= 1; sign += 2) {
                            let prevx = x;
                            let prevy = y;
                            let dx, dy;

                            for (let k = 0; k < 10; k++) {
                                dx = this.bilinearinterp(vf_x, prevx - 0.5, prevy);
                                dy = this.bilinearinterp(vf_y, prevx, prevy - 0.5);

                                const fieldmagnitude = Math.sqrt(dx * dx + dy * dy);
                                this.setalphaFG(0.1 * Math.sqrt(1 / (0.1 * k * k + 1.0) * vectorscalingconstant * fieldmagnitude));

                                if (fieldmagnitude !== 0) {
                                    dx /= fieldmagnitude;
                                    dy /= fieldmagnitude;
                                }

                                this.setColorFloat(1.0, 1.0, 1.0);

                                const nextx = prevx + dx * arrowlength * 0.25 * sign;
                                const nexty = prevy + dy * arrowlength * 0.25 * sign;

                                this.drawLine(
                                    Math.floor((prevx + 0.5) * this.scalefactor),
                                    Math.floor((prevy + 0.5) * this.scalefactor),
                                    Math.floor((nextx + 0.5) * this.scalefactor),
                                    Math.floor((nexty + 0.5) * this.scalefactor),
                                    false
                                );

                                prevx = nextx;
                                prevy = nexty;
                            }
                        }
                    } else if (vector_display_mode === VectorMode.ARROWS && vf_x !== null) {
                        ctr.x = x + 0.5;
                        ctr.y = y + 0.5;

                        arrow.x = this.bilinearinterp(vf_x, x - 0.5, y);
                        arrow.y = this.bilinearinterp(vf_y, x, y - 0.5);

                        const fieldmagnitude = Math.max(0.1, vectorscalingconstant * Math.sqrt(arrow.dot(arrow)));
                        arrow.normalize();
                        tip1.copyFrom(arrow);
                        tip2.copyFrom(arrow);
                        tip1.rotate(Math.PI * 5.0 / 6.0);
                        tip2.rotate(Math.PI * 7.0 / 6.0);

                        body1.copyFrom(ctr);
                        body1.addmult(arrow, -0.5 * arrowlength);
                        body2.copyFrom(ctr);
                        body2.addmult(arrow, 0.5 * arrowlength);
                        tip1.scalarmult(0.4 * arrowlength);
                        tip1.add(body2);
                        tip2.scalarmult(0.4 * arrowlength);
                        tip2.add(body2);

                        this.setColorFloat(fieldmagnitude, fieldmagnitude, fieldmagnitude);
                        this.setalphaFG(0.1 * Math.sqrt(fieldmagnitude));
                        this.drawLine(
                            Math.floor(body1.x * this.scalefactor),
                            Math.floor(body1.y * this.scalefactor),
                            Math.floor(body2.x * this.scalefactor),
                            Math.floor(body2.y * this.scalefactor),
                            false
                        );
                        this.drawLine(
                            Math.floor(body2.x * this.scalefactor),
                            Math.floor(body2.y * this.scalefactor),
                            Math.floor(tip1.x * this.scalefactor),
                            Math.floor(tip1.y * this.scalefactor),
                            false
                        );
                        this.drawLine(
                            Math.floor(body2.x * this.scalefactor),
                            Math.floor(body2.y * this.scalefactor),
                            Math.floor(tip2.x * this.scalefactor),
                            Math.floor(tip2.y * this.scalefactor),
                            false
                        );
                    }
                }
            }
        }

        // Draw text
        this.ctx.font = '12px Arial';
        this.ctx.textAlign = 'left';
        this.ctx.textBaseline = 'alphabetic';

        for (const p of this.voltageprobes) {
            const potential = this.ground ? p.potential - this.ground.potential : p.potential;
            this.drawStringWithBackgroundAndBorder(
                `V = ${this.getSI(potential, 'V')}`,
                p.x * this.scalefactor - 5,
                p.y * this.scalefactor - 12,
                this.ctx
            );
        }

        if (this.ground) {
            this.drawStringWithBackgroundAndBorder(
                `Ground = ${this.getSI(0, 'V')}`,
                this.ground.x * this.scalefactor - 5,
                this.ground.y * this.scalefactor - 12,
                this.ctx
            );
        }

        for (const p of this.currentprobes) {
            const xa = 0.5 * (p.x1 + p.x2) * this.scalefactor;
            const ya = 0.5 * (p.y1 + p.y2) * this.scalefactor;
            let dx = p.x2 - p.x1;
            let dy = p.y2 - p.y1;
            const len = this.length(dx, dy);
            dx = dx / len;
            dy = dy / len;
            if (Math.abs(dx) > Math.abs(dy)) {
                dx = -Math.abs(dx);
            } else {
                dy = -2 * Math.abs(dy);
            }
            this.drawStringWithBackgroundAndBorder(
                `I = ${this.getSI(p.current * this.depth, 'A')}`,
                Math.floor(xa - 8 * dy) - 5,
                Math.floor(ya + 12 * dx) + 5,
                this.ctx
            );
        }

        this.drawStringBackgrounds(this.ctx);
        this.drawStrings(this.ctx);
        this.startNewStringLayer();

        {
            // Tooltip
            let mx_t = this.mx_index;
            let my_t = this.my_index;
            let mi = this.mx_index;
            let mj = this.my_index;

            if (mi < 0) mi = 0;
            if (mi >= this.nx - 1) mi = this.nx - 2;
            if (mj < 0) mj = 0;
            if (mj >= this.ny - 1) mj = this.ny - 2;

            const mat = this.materials[mi][mj];

            let vspacing = 12;
            let voffset = 1 + this.my;
            let hoffset = 5 + this.mx + 15;

            if (this.gui_tooltip.checked) {
                if (voffset + 22 * vspacing > this.ny * this.scalefactor) {
                    voffset = voffset - ((voffset + 22 * vspacing) - this.ny * this.scalefactor);
                }
                if (hoffset + 120 > this.ny * this.scalefactor) {
                    hoffset = hoffset - ((hoffset + 120) - this.ny * this.scalefactor);
                }
            }

            const name = `Material: ${mat.type.name}${mat.modified ? ' (Modified)' : ''}`;
            this.drawBigStringWithBackground(name, hoffset, voffset + 1 * vspacing, this.ctx);

            if (this.gui_tooltip.checked) {
                voffset += 3;
                let line = 2;
                this.drawTwoColumnString('E', this.getSI(this.length(
                    this.bilinearinterp(this.Ex, mx_t - 0.5, my_t),
                    this.bilinearinterp(this.Ey, mx_t, my_t - 0.5)
                ), 'V/m'), hoffset, voffset + line * vspacing, this.ctx); line++;
                this.drawTwoColumnString('B', this.getSI(this.parity * this.bilinearinterp(this.Bz, mx_t - 0.5, my_t - 0.5), 'T'), hoffset, voffset + line * vspacing, this.ctx); line++;
                this.drawTwoColumnString('φ', this.getSI(this.bilinearinterp(this.phi, mx_t, my_t), 'V'), hoffset, voffset + line * vspacing, this.ctx); line++;
                this.drawTwoColumnString('ℰ', this.getSI(mat.emf, 'V/m'), hoffset, voffset + line * vspacing, this.ctx); line++;
                this.drawTwoColumnString('ε/ε₀', this.getSI(mat.eps_r, ''), hoffset, voffset + line * vspacing, this.ctx); line++;
                this.drawTwoColumnString('μ/μ₀', this.getSI(mat.mu_r, ''), hoffset, voffset + line * vspacing, this.ctx); line++;
                this.drawTwoColumnString('ρₙ', this.getSI(this.bilinearinterp(this.rho_n, mx_t, my_t), 'C/m³'), hoffset, voffset + line * vspacing, this.ctx); line++;
                this.drawTwoColumnString('ρₚ', this.getSI(this.bilinearinterp(this.rho_p, mx_t, my_t), 'C/m³'), hoffset, voffset + line * vspacing, this.ctx); line++;
                this.drawTwoColumnString('ρ₀', this.getSI(this.bilinearinterp(this.rho_back, mx_t, my_t), 'C/m³'), hoffset, voffset + line * vspacing, this.ctx); line++;
                this.drawTwoColumnString('ρ', this.getSI(this.bilinearinterp(this.rho_free, mx_t, my_t), 'C/m³'), hoffset, voffset + line * vspacing, this.ctx); line++;
                this.drawTwoColumnString('Jₙ', this.getSI(this.length(
                    this.bilinearinterp(this.Jx_n, mx_t - 0.5, my_t),
                    this.bilinearinterp(this.Jy_n, mx_t, my_t - 0.5)
                ), 'A/m²'), hoffset, voffset + line * vspacing, this.ctx); line++;
                this.drawTwoColumnString('Jₚ', this.getSI(this.length(
                    this.bilinearinterp(this.Jx_p, mx_t - 0.5, my_t),
                    this.bilinearinterp(this.Jy_p, mx_t, my_t - 0.5)
                ), 'A/m²'), hoffset, voffset + line * vspacing, this.ctx); line++;
                this.drawTwoColumnString('J', this.getSI(this.length(
                    this.bilinearinterp(this.Jx_free, mx_t - 0.5, my_t),
                    this.bilinearinterp(this.Jy_free, mx_t, my_t - 0.5)
                ), 'A/m²'), hoffset, voffset + line * vspacing, this.ctx); line++;
                this.drawTwoColumnString('Fₙ', this.getSI(
                    this.bilinearinterp(this.F_n, mx_t, my_t) / this.q_n +
                    this.bilinearinterp(this.phi, mx_t, my_t) - this.W_semi / this.eVtoJ, 'V'
                ), hoffset, voffset + line * vspacing, this.ctx); line++;
                this.drawTwoColumnString('Fₚ', this.getSI(
                    this.bilinearinterp(this.F_p, mx_t, my_t) / this.q_p +
                    this.bilinearinterp(this.phi, mx_t, my_t) - this.W_semi / this.eVtoJ, 'V'
                ), hoffset, voffset + line * vspacing, this.ctx); line++;
                this.drawTwoColumnString('F', this.getSI(this.bilinearinterp(this.F, mx_t, my_t), 'V'), hoffset, voffset + line * vspacing, this.ctx); line++;
                this.drawTwoColumnString('CMFₙ', this.getSI(this.length(
                    this.bilinearinterp(this.cmfx_n, mx_t - 0.5, my_t),
                    this.bilinearinterp(this.cmfy_n, mx_t, my_t - 0.5)
                ) / this.q_n, 'V/m'), hoffset, voffset + line * vspacing, this.ctx); line++;
                this.drawTwoColumnString('CMFₚ', this.getSI(this.length(
                    this.bilinearinterp(this.cmfx_p, mx_t - 0.5, my_t),
                    this.bilinearinterp(this.cmfy_n, mx_t, my_t - 0.5)
                ) / this.q_p, 'V/m'), hoffset, voffset + line * vspacing, this.ctx); line++;
                this.drawTwoColumnString('x', this.getSI(mx_t * this.ds, 'm'), hoffset, voffset + line * vspacing, this.ctx); line++;
                this.drawTwoColumnString('y', this.getSI(this.width - (my_t + 1) * this.ds, 'm'), hoffset, voffset + line * vspacing, this.ctx); line++;
            }
        }

        this.drawStringBackgrounds(this.ctx);
        this.drawStrings(this.ctx);
        this.startNewStringLayer();

        // Status text
        let vspacing = 13;
        let voffset = 3;
        let hoffset = 5;
        let line = 1;
        this.drawStringWithBackground(
            `Time: ${this.getSI(this.time, 's')}`,
            hoffset,
            voffset + line * vspacing,
            this.ctx
        ); line++;
        if (this.gui_paused.checked) {
            this.drawStringWithBackground(
                'Paused',
                hoffset,
                voffset + line * vspacing,
                this.ctx
            ); line++;
        }
        if (this.sign_violation) {
            this.drawStringWithBackground(
                'Warning: Numerical instability detected. Please decrease timestep.',
                hoffset,
                voffset + line * vspacing,
                this.ctx
            ); line++;
        }
        if (this.debugging) {
            const total = window.performance.memory.totalJSHeapSize;
            const used = window.performance.memory.usedJSHeapSize;
            this.drawStringWithBackground(
                `Used memory ${this.getSI(used, 'B')}`,
                hoffset,
                voffset + line * vspacing,
                this.ctx
            ); line++;
            this.drawStringWithBackground(
                `Total memory ${this.getSI(total, 'B')}`,
                hoffset,
                voffset + line * vspacing,
                this.ctx
            ); line++;
            this.drawStringWithBackground(
                `${this.t4.name} ${this.getSI(this.t4.time, 's')}`,
                hoffset,
                voffset + line * vspacing,
                this.ctx
            ); line++;
            this.drawStringWithBackground(
                `${this.t5.name} ${this.getSI(this.t5.avgtime, 's')}`,
                hoffset,
                voffset + line * vspacing,
                this.ctx
            ); line++;
            this.drawStringWithBackground(
                `${this.t6.name} ${this.getSI(this.t6.avgtime * parseInt(this.gui_simspeed_2.value), 's')}`,
                hoffset,
                voffset + line * vspacing,
                this.ctx
            ); line++;
            this.drawStringWithBackground(
                `${this.t7.name} ${this.getSI(this.t7.avgtime, 's')}`,
                hoffset,
                voffset + line * vspacing,
                this.ctx
            ); line++;
            this.drawStringWithBackground(
                `${this.t8.name} ${this.getSI(this.t8.avgtime, 's')}`,
                hoffset,
                voffset + line * vspacing,
                this.ctx
            ); line++;
            this.drawStringWithBackground(
                `${this.t9.name} ${this.getSI(this.t9.avgtime, 's')}`,
                hoffset,
                voffset + line * vspacing,
                this.ctx
            ); line++;
        }

        this.drawStringBackgrounds(this.ctx);
        this.drawStrings(this.ctx);

        // Brush outline
        if (!this.gui_brush_highlight.checked) {
            const r = Math.floor(this.scalefactor * this.brushsize / this.ds);
            const brushshape = this.gui_brush_1.selectedIndex;
            if (Brush.isBrushShapeImportant(brush)) {
                this.ctx.strokeStyle = 'rgb(50, 50, 50)';
                this.ctx.beginPath();
                if (brushshape === 0) {
                    this.ctx.arc(this.mx, this.my, r, 0, 2 * Math.PI);
                } else {
                    this.ctx.rect(this.mx - r, this.my - r, 2 * r - 2, 2 * r - 2);
                }
                this.ctx.stroke();

                this.ctx.strokeStyle = 'rgb(200, 200, 200)';
                this.ctx.beginPath();
                if (brushshape === 0) {
                    this.ctx.arc(this.mx, this.my, r + 1, 0, 2 * Math.PI);
                } else {
                    this.ctx.rect(this.mx - r - 1, this.my - r - 1, 2 * r, 2 * r);
                }
                this.ctx.stroke();
            }
        }

        this.t5.stop();
    }

    clearStrings() {
        this.texts = [];
    }

    drawStringBackgrounds(ctx) {
        if (!this.gui_text_bg.checked) return;

        // First pass: Draw gray border
        for (const text of this.texts) {
            if (text.hasBackground) {
                ctx.font = text.big ? this.bigfont : this.regularfont;
                const metrics = ctx.measureText(text.text);
                const width = Math.max(text.minwidth || 0, metrics.width + 8);
                const height = (text.big ? 16 : 12) + 4; // Approximate font height
                const x = text.x - 3;
                const y = text.y - height + 6;

                ctx.fillStyle = 'gray';
                ctx.fillRect(x - 2, y - 4, width + 3, height + 4);
            }
        }

        // Second pass: Draw black background
        for (const text of this.texts) {
            if (text.hasBackground) {
                ctx.font = text.big ? this.bigfont : this.regularfont;
                const metrics = ctx.measureText(text.text);
                const width = Math.max(text.minwidth || 0, metrics.width + 8);
                const height = (text.big ? 16 : 12) + 4;
                const x = text.x - 3;
                const y = text.y - height + 6;

                ctx.fillStyle = 'black';
                ctx.fillRect(x, y - 2, width - 1, height);
            }
        }
    }

    drawStrings(ctx) {
        for (const text of this.texts) {
            ctx.font = text.big ? this.bigfont : this.regularfont;
            ctx.textAlign = 'left';
            ctx.textBaseline = 'alphabetic';

            ctx.fillStyle = 'rgba(64, 64, 64, 1)';
            ctx.fillText(text.text, text.x + 1, text.y + 1);
            ctx.fillStyle = 'white';
            ctx.fillText(text.text, text.x, text.y);
        }
    }

    startNewStringLayer() {
        this.texts = [];
    }

    drawString(str1, x, y, ctx) {
        this.texts.push(new Text(str1, x, y, false, false, false));
    }

    drawStringWithBackground(str1, x, y, ctx) {
        this.texts.push(new Text(str1, x, y, false, true, false));
    }

    drawStringWithBackgroundAndBorder(str1, x, y, ctx) {
        this.texts.push(new Text(str1, x, y, false, true, true));
    }

    drawBigString(str1, x, y, ctx) {
        this.texts.push(new Text(str1, x, y, true, false, false));
    }

    drawBigStringWithBackground(str1, x, y, ctx) {
        this.texts.push(new Text(str1, x, y, true, true, false));
    }

    drawTwoColumnString(str1, str2, x, y, ctx) {
        const text1 = new Text(str1.padEnd(10), x, y, false, true, true);
        const text2 = new Text(str2, x + 40, y, false, true, true);
        text2.minwidth = 80;
        this.texts.push(text1, text2);
    }

    length(x, y) {
        return Math.sqrt(x * x + y * y);
    }

    getSI(quantity, unit) {
        if (!Number.isFinite(quantity)) {
            return `${quantity} ${unit}`;
        }

        const mag = Math.abs(quantity);
        const precision = (x) => x.toFixed(2);
        if (mag < 1e-18) {
            return `0 ${unit}`;
        } else if (mag < 1e-15) {
            return `${precision(quantity * 1e15)} f${unit}`;
        } else if (mag < 1e-12) {
            return `${precision(quantity * 1e15)} f${unit}`;
        } else if (mag < 1e-9) {
            return `${precision(quantity * 1e12)} p${unit}`;
        } else if (mag < 1e-6) {
            return `${precision(quantity * 1e9)} n${unit}`;
        } else if (mag < 1e-3) {
            return `${precision(quantity * 1e6)} μ${unit}`;
        } else if (mag < 1) {
            return `${precision(quantity * 1e3)} m${unit}`;
        } else if (mag < 1e3) {
            return `${precision(quantity)} ${unit}`;
        } else if (mag < 1e6) {
            return `${precision(quantity * 1e-3)} k${unit}`;
        } else if (mag < 1e9) {
            return `${precision(quantity * 1e-6)} M${unit}`;
        } else if (mag < 1e12) {
            return `${precision(quantity * 1e-9)} G${unit}`;
        } else {
            return `${precision(quantity * 1e-12)} T${unit}`;
        }
    }

    async readFile(url = null) {
        const loadingDialog = document.createElement('div');
        loadingDialog.style.position = 'fixed';
        loadingDialog.style.top = '50%';
        loadingDialog.style.left = '50%';
        loadingDialog.style.transform = 'translate(-50%, -50%)';
        loadingDialog.style.background = 'white';
        loadingDialog.style.padding = '20px';
        loadingDialog.style.border = '1px solid black';
        loadingDialog.textContent = 'Loading file, please wait...';

        try {
            // Create file input element
            let arrayBuffer = null;

            if (url == null) {
                const input = document.createElement('input');
                input.type = 'file';
                input.accept = this.fileextension;
    
                // Simulate file chooser dialog
                const filePromise = new Promise((resolve) => {
                    input.addEventListener('change', () => {
                        resolve(input.files[0] || null);
                    });
                    input.click();
                });
    
                this.infile = await filePromise;
                if (!this.infile) {
                    return false;
                }
            }

            // Show loading dialog
            document.body.appendChild(loadingDialog);

            if (url == null)
            {
                arrayBuffer = await this.infile.arrayBuffer();
            } else {
                const response = await fetch(url);
                arrayBuffer = await response.arrayBuffer();
            }

            // Read and parse file
            const decompressed = pako.ungzip(new Uint8Array(arrayBuffer), { to: 'string' }); // Decompress GZIP using pako
            const data = JSON.parse(decompressed);

            if (!data.version || data.version !== 1) {
                throw new Error('Invalid file version');
            }

            // Initialize grid
            this.initializeGrid(this.width, this.default_resolution);

            // Load data
            for (const [key, value] of Object.entries(data)) {
                switch (key) {
                    case 'time':
                        this.time = value;
                        break;
                    case 'gui_paused':
                        this.gui_paused.checked = value;
                        break;
                    case 'gui_tooltip':
                        this.gui_tooltip.checked = value;
                        break;
                    case 'gui_text_bg':
                        this.gui_text_bg.checked = value;
                        break;
                    case 'gui_view':
                        this.gui_view.selectedIndex = value;
                        break;
                    case 'gui_view_vec':
                        this.gui_view_vec.selectedIndex = value;
                        break;
                    case 'gui_view_vec_mode':
                        this.gui_view_vec_mode.selectedIndex = value;
                        break;
                    case 'gui_simspeed':
                        this.gui_simspeed.value = value;
                        break;
                    case 'gui_simspeed_2':
                        this.gui_simspeed_2.value = value;
                        break;
                    case 'gui_brightness':
                        this.gui_brightness.value = value;
                        break;
                    case 'gui_brightness_vec':
                        this.gui_brightness_vec.value = value;
                        break;
                    case 'gui_elem_colors':
                        this.gui_elem_colors.checked = value;
                        break;
                    case 'gui_bc':
                        this.gui_bc.selectedIndex = value;
                        break;
                    case 'description':
                        this.textPane.value = value;
                        break;
                    case 'ex':
                        this.Ex = this.validateArraySize(value);
                        break;
                    case 'ey':
                        this.Ey = this.validateArraySize(value);
                        break;
                    case 'hz':
                        this.Hz = this.validateArraySize(value);
                        break;
                    case 'rho_c':
                        this.rho_abs = this.validateArraySize(value);
                        break;
                    case 'rho_n':
                        this.rho_n = this.validateArraySize(value);
                        break;
                    case 'rho_p':
                        this.rho_p = this.validateArraySize(value);
                        break;
                    case 'rho_back':
                        this.rho_back = this.validateArraySize(value);
                        break;
                    case 'rho_free':
                        this.rho_free = this.validateArraySize(value);
                        break;
                    case 'jx_c':
                        this.Jx_abs = this.validateArraySize(value);
                        break;
                    case 'jy_c':
                        this.Jy_abs = this.validateArraySize(value);
                        break;
                    case 'jx_n':
                        this.Jx_n = this.validateArraySize(value);
                        break;
                    case 'jy_n':
                        this.Jy_n = this.validateArraySize(value);
                        break;
                    case 'jx_p':
                        this.Jx_p = this.validateArraySize(value);
                        break;
                    case 'jy_p':
                        this.Jy_p = this.validateArraySize(value);
                        break;
                    case 'materials':
                        this.materials = this.create2DArray(this.nx, this.ny, Material, () => new Material());
                        for (let i = 0; i < Math.min(this.nx, data.materials.length); i++) {
                            for (let j = 0; j < Math.min(this.ny, data.materials[i].length); j++) {
                                const savedMaterial = data.materials[i][j];
                                if (savedMaterial && savedMaterial.type) {
                                    const typeKey = savedMaterial.type;
                                    if (MaterialType[typeKey]) {
                                        const material = new Material();
                                        material.type = MaterialType[typeKey]; // Map to MaterialType object
                                        material.modified = savedMaterial.modified || false;
                                        material.activated = savedMaterial.activated || 1;
                                        material.conducting = savedMaterial.conducting || 0;
                                        material.semiconducting = savedMaterial.semiconducting || 0;
                                        material.emf = savedMaterial.emf || 0.0;
                                        material.emf_direction = savedMaterial.emf_direction || 0.0;
                                        material.eps_r = savedMaterial.eps_r || 1.0;
                                        material.mu_r = savedMaterial.mu_r || 1.0;
                                        material.rho_back = savedMaterial.rho_back || 0.0;
                                        material.ni = savedMaterial.ni || 0;
                                        material.W = savedMaterial.W || 0;
                                        material.Eb = savedMaterial.Eb || 0;
                                        material.Ea = savedMaterial.Ea || 0;
                                        material.absorptivity = savedMaterial.absorptivity || 0.0;
                                        this.materials[i][j] = material;
                                    } else {
                                        console.warn(`Unknown material type: ${typeKey} at [${i}][${j}]`);
                                        this.materials[i][j] = new Material(); // Default to VACUUM
                                    }
                                }
                            }
                        }
                        break;
                    case 'voltageprobes':
                        this.voltageprobes = value.map(item => Object.assign(new VoltageProbe(), item));
                        break;
                    case 'currentprobes':
                        this.currentprobes = value.map(item => Object.assign(new CurrentProbe(), item));
                        break;
                    case 'ground':
                        this.ground = value ? Object.assign(new VoltageProbe(), value) : null;
                        break;
                    default:
                        break;
                }
            }

            this.textPane.readOnly = true;
            this.textPane.scrollTop = 0;
            this.constructBoundary();
            this.initializeAllMaterials();
            this.updateAllMaterials();
            this.calcMiscFields(true);

            document.body.removeChild(loadingDialog);

            if (url != null) {
                document.title = `Brandon's semiconductor simulator - ${url}`;
            } else {
                document.title = `Brandon's semiconductor simulator - ${this.infile.name}`;
            }

            return true;
        } catch (error) {
            console.error('Error loading file:', error);
            alert('Error: Unable to load file.');
            document.body.removeChild(loadingDialog);

            return false;
        }
    }

    validateArraySize(array) {
        if (array === null || array === undefined) {
            return Array(this.nx).fill().map(() => Array(this.ny).fill(0));
        }

        if (array.length !== this.nx) {
            throw new Error('Array size mismatch.');
        }

        for (let i = 0; i < array.length; i++) {
            if (array[i].length !== this.ny) {
                throw new Error('Array size mismatch.');
            }
        }

        return array;
    }

    async writeFile() {
        // Show saving dialog
        const savingDialog = document.createElement('div');
        savingDialog.style.position = 'fixed';
        savingDialog.style.top = '50%';
        savingDialog.style.left = '50%';
        savingDialog.style.transform = 'translate(-50%, -50%)';
        savingDialog.style.background = 'white';
        savingDialog.style.padding = '20px';
        savingDialog.style.border = '1px solid black';
        savingDialog.textContent = 'Saving file, please wait...';
        
        try {
            // Prompt for filename
            let filename = prompt('Enter filename to save:', `simulation${this.fileextension}`);
            if (!filename) {
                return false;
            }
            if (!filename.endsWith(this.fileextension)) {
                filename += this.fileextension;
            }

            // Check for overwrite (simulated)
            if (this.outfile && this.outfile.name === filename) {
                if (!confirm('A file with that name already exists. Do you wish to overwrite it?')) {
                    return false;
                }
            }
            this.outfile = { name: filename };

            document.body.appendChild(savingDialog);

            // Prepare data
            const data = {
                version: this.saveversion,
                time: this.time,
                gui_paused: this.gui_paused.checked,
                gui_tooltip: this.gui_tooltip.checked,
                gui_text_bg: this.gui_text_bg.checked,
                gui_view: this.gui_view.selectedIndex,
                gui_view_vec: this.gui_view_vec.selectedIndex,
                gui_view_vec_mode: this.gui_view_vec_mode.selectedIndex,
                gui_simspeed: parseInt(this.gui_simspeed.value),
                gui_simspeed_2: parseInt(this.gui_simspeed_2.value),
                gui_brightness: parseInt(this.gui_brightness.value),
                gui_brightness_vec: parseInt(this.gui_brightness_vec.value),
                gui_elem_colors: this.gui_elem_colors.checked,
                gui_bc: this.gui_bc.selectedIndex,
                description: this.textPane.value,
                ex: this.Ex,
                ey: this.Ey,
                hz: this.Hz,
                rho_c: this.rho_abs,
                rho_n: this.rho_n,
                rho_p: this.rho_p,
                rho_back: this.rho_back,
                rho_free: this.rho_free,
                jx_c: this.Jx_abs,
                jy_c: this.Jy_abs,
                jx_n: this.Jx_n,
                jy_n: this.Jy_n,
                jx_p: this.Jx_p,
                jy_p: this.Jy_p,
                materials: this.materials,
                voltageprobes: this.voltageprobes,
                currentprobes: this.currentprobes,
                ground: this.ground
            };

            // Serialize and compress
            const json = JSON.stringify(data);
            let blob;
            if (this.fileextension === '.json.gz') {
                const compressed = pako.gzip(json);
                blob = new Blob([compressed], { type: 'application/gzip' });
            } else {
                blob = new Blob([json], { type: 'application/json' });
            }

            // Trigger download
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = filename;
            a.click();
            URL.revokeObjectURL(url);

            document.body.removeChild(savingDialog);
            document.title = `Brandon's semiconductor simulator - ${filename}`;

            return true;
        } catch (error) {
            console.error('Error saving file:', error);
            alert('Error: Unable to save file.');
            document.body.removeChild(savingDialog);
            return false;
        }
    }

    mouseClicked(event) {
        // No-op
    }

    mouseEntered(event) {
        // No-op
    }

    mouseExited(event) {
        // No-op
    }

    mousePressed(event) {
        this.mouse_pressed = true;
        this.mousebutton = event.button; //  + 1; // JavaScript: 0=left, 1=middle, 2=right; Java: 1=left, 2=middle, 3=right
        this.mx = event.offsetX;
        this.my = event.offsetY;
        this.mx_start = event.offsetX;
        this.my_start = event.offsetY;
    }

    mouseReleased(event) {
        this.mouse_pressed = false;
    }

    mouseDragged(event) {
        this.mx = event.offsetX;
        this.my = event.offsetY;
    }

    mouseMoved(event) {
        this.mx = event.offsetX;
        this.my = event.offsetY;
    }

    key_pause() {
        this.gui_paused.checked = !this.gui_paused.checked;
    }

    key_frame() {
        this.advanceframe = true;
    }

    key_dbg() {
        this.debugging = !this.debugging;
        Timer.allEnabled = this.debugging;
    }

    key_changebrush() {
        this.gui_brush_1.selectedIndex = (this.gui_brush_1.selectedIndex + 1) % 2;
    }

    key_shift() {
        this.shift_down = true;
    }

    key_shift_up() {
        this.shift_down = false;
    }

    key_ctrl() {
        this.ctrl_down = true;
    }

    key_ctrl_up() {
        this.ctrl_down = false;
    }

    key_cut() {
        this.cut = true;
    }

    key_copy() {
        this.copy = true;
    }

    key_paste() {
        this.paste = true;
    }

    key_delete() {
        this.delete = true;
    }

    key_color() {
        this.gui_elem_colors.checked = !this.gui_elem_colors.checked;
    }

    key_scalar_view() {
        const currentView = Object.values(ScalarView)[this.gui_view.selectedIndex];
        if (currentView === ScalarView.NONE) {
            this.gui_view.value = this.prev_scalar_view.name;
        } else {
            this.prev_scalar_view = currentView;
            this.gui_view.value = ScalarView.NONE.name;
        }
    }

    key_vector_view() {
        const currentView = Object.values(VectorView)[this.gui_view_vec.selectedIndex];
        if (currentView === VectorView.NONE) {
            this.gui_view_vec.value = this.prev_vector_view.name;
        } else {
            this.prev_vector_view = currentView;
            this.gui_view_vec.value = VectorView.NONE.name;
        }
    }

    key_tooltip() {
        this.gui_tooltip.checked = !this.gui_tooltip.checked;
    }

    key_textbg() {
        this.gui_text_bg.checked = !this.gui_text_bg.checked;
    }

    key_alt() {
        this.alt_down = true;
    }

    key_alt_up() {
        this.alt_down = false;
    }

    mouseWheelMoved(event) {
        event.preventDefault();
        const rotation = event.deltaY > 0 ? 1 : -1; // Normalize to +/-1
        this.gui_brushsize.value = parseInt(this.gui_brushsize.value) - 10 * rotation;
    }

    actionPerformed(event) {
        const source = event.target;
        if (source === this.gui_reset) {
            this.clear = true;
        } else if (source === this.gui_resetall) {
            this.reset = true;
        } else if (source === this.gui_save) {
            this.save = true;
        } else if (source === this.gui_open) {
            this.load = true;
        } else if (source === this.gui_help) {
            this.help.style.display = 'block';
        } else if (source === this.gui_close_help) {
            this.help.style.display = '';
        } else if (source === this.gui_editdesc) {
            this.textPane.readOnly = !this.textPane.readOnly;
        } else if (source === this.gui_view) {
            this.updateMiscFields = true;
        } else if (source === this.gui_view_vec) {
            this.updateMiscFields = true;
        } else if (source === this.gui_brush) {
            this.brush_changed = true;
        }
    }

    keyTyped(event) {
        // Do nothing
    }

    keyPressed(e) {
        e.preventDefault(); // Prevent browser actions (e.g., scrolling)

        const key = e.key.toUpperCase();
        const ctrlOrMeta = e.ctrlKey || e.metaKey;
        const shift = e.shiftKey;
        const alt = e.altKey;

        if (!this.texting) {
            if (key === 'P' || key === ' ') {
                this.key_pause();
            } else if (key === 'F') {
                this.key_frame();
            } else if (key === 'D') {
                this.key_dbg();
            } else if (key === 'Q') {
                this.key_changebrush();
            } else if (key === 'X' && ctrlOrMeta) {
                this.key_cut();
            } else if (key === 'C' && ctrlOrMeta) {
                this.key_copy();
            } else if (key === 'V' && ctrlOrMeta) {
                this.key_paste();
            } else if (key === 'BACKSPACE' || key === 'DELETE') {
                this.key_delete();
            } else if (key === 'C') {
                this.key_color();
            } else if (key === 'S') {
                this.key_scalar_view();
            } else if (key === 'V') {
                this.key_vector_view();
            } else if (key === 'T') {
                this.key_tooltip();
            } else if (key === 'G') {
                this.key_textbg();
            }
        }

        if (key === 'SHIFT' && shift) {
            this.key_shift();
        } else if (key === 'CONTROL' && ctrlOrMeta) {
            this.key_ctrl();
        } else if (key === 'META' && ctrlOrMeta) {
            this.key_ctrl();
        } else if (key === 'ALT' && alt) {
            this.key_alt();
        }

        if (this.texting) {
            const char = e.key;
            if (this.Font7x5.characterExists(char)) {
                for (let i = 0; i < 5; i++) {
                    for (let j = 0; j < 7; j++) {
                        if (
                            this.text_x + i + 1 >= 0 &&
                            this.text_x + i + 1 < this.nx &&
                            this.text_y + j >= 0 &&
                            this.text_y + j < this.ny &&
                            this.Font7x5.getPixel(char, 4 - i, j) === 1 &&
                            this.materials[this.text_x + i + 1][this.text_y + j].type === MaterialType.VACUUM
                        ) {
                            this.initializeMaterial(this.text_x + i + 1, this.text_y + j, MaterialType.DECO);
                        }
                    }
                }
                this.text_x += 6;
                this.updateAllMaterials();
            }
            if (e.key === 'Backspace' || e.key === 'Delete') {
                this.text_x -= 6;
                for (let i = 0; i < 5; i++) {
                    for (let j = 0; j < 7; j++) {
                        if (
                            this.text_x + i + 1 >= 0 &&
                            this.text_x + i + 1 < this.nx &&
                            this.text_y + j >= 0 &&
                            this.text_y + j < this.ny &&
                            this.materials[this.text_x + i + 1][this.text_y + j].type === MaterialType.DECO
                        ) {
                            this.materials[this.text_x + i + 1][this.text_y + j].erase();
                        }
                    }
                }
                this.updateAllMaterials();
            }
        }
    }

    keyReleased(e) {
        const key = e.key.toUpperCase();
        if (key === 'SHIFT') {
            this.key_shift_up();
        } else if (key === 'CONTROL' || key === 'META') {
            this.key_ctrl_up();
        } else if (key === 'ALT') {
            this.key_alt_up();
        }
    }
}
