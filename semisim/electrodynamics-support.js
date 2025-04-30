const MaterialType = Object.freeze({
    EMF: { name: "Voltage source (Adjustable)", color_r: 230, color_g: 216, color_b: 46, color_grayscale: 230, toString: function() { return `Material: ${this.name}`; } },
    SWITCH: { name: "Switch", color_r: 194, color_g: 194, color_b: 194, color_grayscale: 120, toString: function() { return `Material: ${this.name}`; } },
    METAL: { name: "Metal", color_r: 153, color_g: 153, color_b: 153, color_grayscale: 120, toString: function() { return `Material: ${this.name}`; } },
    METAL_HIGH_C: { name: "Conductive metal", color_r: 191, color_g: 191, color_b: 191, color_grayscale: 120, toString: function() { return `Material: ${this.name}`; } },
    METAL_LOW_C: { name: "Resistive metal", color_r: 94, color_g: 94, color_b: 94, color_grayscale: 120, toString: function() { return `Material: ${this.name}`; } },
    METAL_HIGH_W: { name: "High workfunction metal", color_r: 163, color_g: 116, color_b: 116, color_grayscale: 120, toString: function() { return `Material: ${this.name}`; } },
    METAL_LOW_W: { name: "Low workfunction metal", color_r: 116, color_g: 121, color_b: 163, color_grayscale: 120, toString: function() { return `Material: ${this.name}`; } },
    SEMI: { name: "Intrinsic semiconductor", color_r: 207, color_g: 161, color_b: 212, color_grayscale: 120, toString: function() { return `Material: ${this.name}`; } },
    SEMI_P_TYPE: { name: "P-type semiconductor", color_r: 191, color_g: 74, color_b: 34, color_grayscale: 120, toString: function() { return `Material: ${this.name}`; } },
    SEMI_N_TYPE: { name: "N-type semiconductor", color_r: 84, color_g: 123, color_b: 191, color_grayscale: 120, toString: function() { return `Material: ${this.name}`; } },
    SEMI_HEAVY_P_TYPE: { name: "Heavily doped P-type semiconductor", color_r: 204, color_g: 41, color_b: 41, color_grayscale: 120, toString: function() { return `Material: ${this.name}`; } },
    SEMI_HEAVY_N_TYPE: { name: "Heavily doped N-type semiconductor", color_r: 39, color_g: 52, color_b: 194, color_grayscale: 120, toString: function() { return `Material: ${this.name}`; } },
    SEMI_LIGHT_P_TYPE: { name: "Lightly doped P-type semiconductor", color_r: 201, color_g: 131, color_b: 73, color_grayscale: 120, toString: function() { return `Material: ${this.name}`; } },
    SEMI_LIGHT_N_TYPE: { name: "Lightly doped N-type semiconductor", color_r: 137, color_g: 188, color_b: 204, color_grayscale: 120, toString: function() { return `Material: ${this.name}`; } },
    DIELECTRIC: { name: "Dielectric", color_r: 81, color_g: 171, color_b: 51, color_grayscale: 120, toString: function() { return `Material: ${this.name}`; } },
    FERROMAGNET: { name: "Ferromagnet", color_r: 116, color_g: 50, color_b: 117, color_grayscale: 120, toString: function() { return `Material: ${this.name}`; } },
    DECO: { name: "Decoration", color_r: 255, color_g: 255, color_b: 255, color_grayscale: 255, toString: function() { return `Material: ${this.name}`; } },
    ABSORBER: { name: "Absorber", color_r: 50, color_g: 50, color_b: 50, color_grayscale: Math.round(0.7 * Math.max(50, 50)), toString: function() { return `Material: ${this.name}`; } },
    VACUUM: { name: "Vacuum", color_r: 20, color_g: 20, color_b: 20, color_grayscale: Math.round(0.7 * Math.max(20, 20)), toString: function() { return `Material: ${this.name}`; } },

    isConducting(material) {
        return [MaterialType.EMF, MaterialType.SWITCH, MaterialType.METAL, MaterialType.METAL_HIGH_W, MaterialType.METAL_LOW_W, MaterialType.METAL_HIGH_C, MaterialType.METAL_LOW_C].includes(material);
    },

    isSemiconducting(material) {
        return [MaterialType.SEMI_P_TYPE, MaterialType.SEMI_N_TYPE, MaterialType.SEMI, MaterialType.SEMI_HEAVY_P_TYPE, MaterialType.SEMI_HEAVY_N_TYPE, MaterialType.SEMI_LIGHT_P_TYPE, MaterialType.SEMI_LIGHT_N_TYPE].includes(material);
    }
});


MaterialType_values = [
    MaterialType.EMF,
    MaterialType.SWITCH,
    MaterialType.METAL,
    MaterialType.METAL_HIGH_C,
    MaterialType.METAL_LOW_C,
    MaterialType.METAL_HIGH_W,
    MaterialType.METAL_LOW_W,
    MaterialType.SEMI,
    MaterialType.SEMI_P_TYPE,
    MaterialType.SEMI_N_TYPE,
    MaterialType.SEMI_HEAVY_P_TYPE,
    MaterialType.SEMI_HEAVY_N_TYPE,
    MaterialType.SEMI_LIGHT_P_TYPE,
    MaterialType.SEMI_LIGHT_N_TYPE,
    MaterialType.DIELECTRIC,
    MaterialType.FERROMAGNET,
    MaterialType.DECO,
    MaterialType.ABSORBER,
    MaterialType.VACUUM
];

class Material {
    constructor() {
        this.type = MaterialType.VACUUM;
        this.modified = false;
        this.activated = 1;
        this.conducting = 0;
        this.semiconducting = 0;
        this.emf = 0.0;
        this.emf_direction = 0.0;
        this.eps_r = 1.0;
        this.mu_r = 1.0;
        this.rho_back = 0.0;
        this.ni = 0;
        this.W = 0;
        this.Eb = 0;
        this.Ea = 0;
        this.absorptivity = 0.0;
    }

    erase() {
        this.type = MaterialType.VACUUM;
        this.modified = false;
        this.activated = 1;
        this.conducting = 0;
        this.semiconducting = 0;
        this.emf = 0.0;
        this.emf_direction = 0.0;
        this.eps_r = 1.0;
        this.mu_r = 1.0;
        this.rho_back = 0.0;
        this.ni = 0;
        this.W = 0;
        this.Eb = 0;
        this.Ea = 0;
        this.absorptivity = 0;
    }

    clone() {
        const cloned = new Material();
        cloned.type = this.type;
        cloned.modified = this.modified;
        cloned.activated = this.activated;
        cloned.conducting = this.conducting;
        cloned.semiconducting = this.semiconducting;
        cloned.emf = this.emf;
        cloned.emf_direction = this.emf_direction;
        cloned.eps_r = this.eps_r;
        cloned.mu_r = this.mu_r;
        cloned.rho_back = this.rho_back;
        cloned.ni = this.ni;
        cloned.W = this.W;
        cloned.Eb = this.Eb;
        cloned.Ea = this.Ea;
        cloned.absorptivity = this.absorptivity;
        return cloned;
    }
}

class VoltageProbe {
    constructor(x = 0, y = 0, potential = 0) {
        this.x = x;
        this.y = y;
        this.potential = potential;
    }
}

class CurrentProbe {
    constructor(x1 = 0, y1 = 0, x2 = 0, y2 = 0, current = 0) {
        this.x1 = x1;
        this.y1 = y1;
        this.x2 = x2;
        this.y2 = y2;
        this.current = current;
    }
}

class RenderCanvas {
    constructor(electrodynamics, canvasElement) {
        this.parent = electrodynamics;
        this.canvas = canvasElement;
        this.ctx = canvasElement.getContext('2d');
    }

    render() {
        this.parent.render();
        if (this.parent.screen) {
            this.ctx.putImageData(this.parent.screen, 0, 0);
        }
    }
}

class Timer {
    static allEnabled = false;

    constructor(name, enabled = true) {
        this.name = name;
        this.enabled = enabled;
        this.tstart = 0;
        this.outputavg = false;
        this.avgtime = 0;
        this.time = 0;
    }

    start() {
        if (this.enabled) {
            this.tstart = performance.now();
        }
    }

    disableOutput() {
        this.enabled = false;
    }

    enableOutput() {
        this.enabled = true;
    }

    stop(msg = '') {
        if (Timer.allEnabled && this.enabled) {
            const tend = performance.now();
            const diff = tend - this.tstart;
            this.time = diff / 1000;
            this.avgtime = this.avgtime * 0.95 + this.time * 0.05;
            if (msg) {
                console.log(`${this.name}: ${this.time}s`);
            }
        }
    }
}

class Vector {
    constructor(x = 0, y = 0) {
        this.x = x;
        this.y = y;
    }

    copy() {
        return new Vector(this.x, this.y);
    }

    copyFrom(b) {
        this.x = b.x;
        this.y = b.y;
    }

    initialize(x, y) {
        this.x = x;
        this.y = y;
    }

    add(b) {
        this.x += b.x;
        this.y += b.y;
    }

    scalarmult(c) {
        this.x *= c;
        this.y *= c;
    }

    addmult(b, c) {
        this.x += b.x * c;
        this.y += b.y * c;
    }

    rotate(theta) {
        const xf = this.x * Math.cos(theta) + this.y * Math.sin(theta);
        const yf = -this.x * Math.sin(theta) + this.y * Math.cos(theta);
        this.x = xf;
        this.y = yf;
    }

    normalize() {
        const magnitude = Math.sqrt(this.x * this.x + this.y * this.y);
        if (magnitude !== 0) {
            this.x /= magnitude;
            this.y /= magnitude;
        }
    }

    dot(b) {
        return this.x * b.x + this.y * b.y;
    }
}

const Brush = Object.freeze({
    INTERACT: { name: "Interact", toString: function() { return this.name; } },
    DRAW: { name: "Draw", toString: function() { return this.name; } },
    VOLTAGE: { name: "Add voltage probe", toString: function() { return this.name; } },
    CURRENT: { name: "Add current probe", toString: function() { return this.name; } },
    GROUND: { name: "Add ground", toString: function() { return this.name; } },
    DELETEPROBE: { name: "Delete probe", toString: function() { return this.name; } },
    REPLACE: { name: "Replace", toString: function() { return this.name; } },
    LINE: { name: "Line", toString: function() { return this.name; } },
    FILL: { name: "Fill", toString: function() { return this.name; } },
    ERASE: { name: "Eraser", toString: function() { return this.name; } },
    SELECT: { name: "Select & Move", toString: function() { return this.name; } },
    FLOODSELECT: { name: "Select region", toString: function() { return this.name; } },
    TEXT: { name: "Add text", toString: function() { return this.name; } },

    isMaterialModifyingBrush(brush) {
        return [Brush.DRAW, Brush.LINE, Brush.REPLACE, Brush.ERASE, Brush.FILL].includes(brush);
    },

    isBrushShapeImportant(brush) {
        return [Brush.DRAW, Brush.LINE, Brush.REPLACE, Brush.ERASE].includes(brush);
    }
});

Brush_values = [
    Brush.INTERACT,
    Brush.DRAW,
    Brush.VOLTAGE,
    Brush.CURRENT,
    Brush.GROUND,
    Brush.DELETEPROBE,
    Brush.REPLACE,
    Brush.LINE,
    Brush.FILL,
    Brush.ERASE,
    Brush.SELECT,
    Brush.FLOODSELECT,
    Brush.TEXT
];

const BrushShape = Object.freeze({
    CIRCLE: { name: "Circle brush", toString: function() { return this.name; } },
    SQUARE: { name: "Square brush", toString: function() { return this.name; } }
});

const ScalarView = Object.freeze({
    NONE: { name: "No scalar overlay", toString: function() { return this.name; } },
    E_FIELD: { name: "View E field magnitude", toString: function() { return this.name; } },
    B_FIELD: { name: "View B field", toString: function() { return this.name; } },
    CHARGE: { name: "View ρ: Net charge density", toString: function() { return this.name; } },
    CURRENT: { name: "View J: Total current magnitude", toString: function() { return this.name; } },
    H_FIELD: { name: "View H field", toString: function() { return this.name; } },
    POTENTIAL: { name: "View φ: Electric scalar potential", toString: function() { return this.name; } },
    ENERGY: { name: "View u: Electromagnetic energy density", toString: function() { return this.name; } },
    ELECTRON_CHARGE: { name: "View ρₙ: Electron charge density", toString: function() { return this.name; } },
    HOLE_CHARGE: { name: "View ρₚ: Hole charge density", toString: function() { return this.name; } },
    COMBINED_CHARGE: { name: "View: Combined electron+hole charge density", toString: function() { return this.name; } },
    BACKGROUND_CHARGE: { name: "View ρ₀: Background charge density", toString: function() { return this.name; } },
    HEAT: { name: "View Q: Heat dissipation", toString: function() { return this.name; } },
    ENTROPY: { name: "View s: Entropy generation (Free energy dissipation)", toString: function() { return this.name; } },
    ELECTRON_POTENTIAL: { name: "View Fₙ: Electron chemical potential", toString: function() { return this.name; } },
    HOLE_POTENTIAL: { name: "View Fₚ: Hole chemical potential", toString: function() { return this.name; } },
    AVERAGE_POTENTIAL: { name: "View F: Average electrochemical potential", toString: function() { return this.name; } },
    RECOMBINATION: { name: "View R: Recombination rate", toString: function() { return this.name; } },
    LIGHT: { name: "View: Emitted light", toString: function() { return this.name; } },
    DEBUG: { name: "Debug", toString: function() { return this.name; } }
});

const VectorView = Object.freeze({
    NONE: { name: "No vector overlay", toString: function() { return this.name; } },
    E_FIELD: { name: "View E field", toString: function() { return this.name; } },
    D_FIELD: { name: "View D field", toString: function() { return this.name; } },
    ELECTRON_CURRENT: { name: "View Jₙ: Electron current", toString: function() { return this.name; } },
    HOLE_CURRENT: { name: "View Jₚ: Hole current", toString: function() { return this.name; } },
    TOTAL_CURRENT: { name: "View J: Total current", toString: function() { return this.name; } },
    EMF: { name: "View ℰ: External electromotive force", toString: function() { return this.name; } },
    POYNTING: { name: "View S: Poynting vector", toString: function() { return this.name; } }
});

const VectorMode = Object.freeze({
    ARROWS: { name: "Show vectors", toString: function() { return this.name; } },
    LINES: { name: "Show lines", toString: function() { return this.name; } }
});

const BoundaryCondition = Object.freeze({
    DISSIPATIVE: { name: "Absorbing boundary", toString: function() { return this.name; } },
    CONDUCTING: { name: "Conducting boundary", toString: function() { return this.name; } }
});

class FloodFillCoordinate {
    constructor(i, j) {
        this.i = i;
        this.j = j;
    }
}

class Text {
    constructor(text, x, y, big, hasBackground, hasBorder) {
        this.text = text;
        this.x = x;
        this.y = y;
        this.big = big;
        this.hasBackground = hasBackground;
        this.hasBorder = hasBorder;
        this.minwidth = 0;
    }
}

class Font7x5 {
    constructor() {
        this.FONT_MAP = {
        'A': new Uint8Array([0b01110, 0b10001, 0b10001, 0b11111, 0b10001, 0b10001, 0b10001]),
        'B': new Uint8Array([0b11110, 0b10001, 0b10001, 0b11110, 0b10001, 0b10001, 0b11110]),
        'C': new Uint8Array([0b01110, 0b10001, 0b10000, 0b10000, 0b10000, 0b10001, 0b01110]),
        'D': new Uint8Array([0b11100, 0b10010, 0b10001, 0b10001, 0b10001, 0b10010, 0b11100]),
        'E': new Uint8Array([0b11111, 0b10000, 0b10000, 0b11110, 0b10000, 0b10000, 0b11111]),
        'F': new Uint8Array([0b11111, 0b10000, 0b10000, 0b11110, 0b10000, 0b10000, 0b10000]),
        'G': new Uint8Array([0b01110, 0b10001, 0b10000, 0b10111, 0b10001, 0b10001, 0b01110]),
        'H': new Uint8Array([0b10001, 0b10001, 0b10001, 0b11111, 0b10001, 0b10001, 0b10001]),
        'I': new Uint8Array([0b01110, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b01110]),
        'J': new Uint8Array([0b00001, 0b00001, 0b00001, 0b00001, 0b10001, 0b10001, 0b01110]),
        'K': new Uint8Array([0b10001, 0b10010, 0b10100, 0b11000, 0b10100, 0b10010, 0b10001]),
        'L': new Uint8Array([0b10000, 0b10000, 0b10000, 0b10000, 0b10000, 0b10000, 0b11111]),
        'M': new Uint8Array([0b10001, 0b11011, 0b10101, 0b10101, 0b10001, 0b10001, 0b10001]),
        'N': new Uint8Array([0b10001, 0b10001, 0b11001, 0b10101, 0b10011, 0b10001, 0b10001]),
        'O': new Uint8Array([0b01110, 0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01110]),
        'P': new Uint8Array([0b11110, 0b10001, 0b10001, 0b11110, 0b10000, 0b10000, 0b10000]),
        'Q': new Uint8Array([0b01110, 0b10001, 0b10001, 0b10001, 0b10101, 0b10010, 0b01101]),
        'R': new Uint8Array([0b11110, 0b10001, 0b10001, 0b11110, 0b10100, 0b10010, 0b10001]),
        'S': new Uint8Array([0b01111, 0b10000, 0b10000, 0b01110, 0b00001, 0b00001, 0b11110]),
        'T': new Uint8Array([0b11111, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100]),
        'U': new Uint8Array([0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01110]),
        'V': new Uint8Array([0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01010, 0b00100]),
        'W': new Uint8Array([0b10001, 0b10001, 0b10001, 0b10101, 0b10101, 0b10101, 0b01010]),
        'X': new Uint8Array([0b10001, 0b10001, 0b01010, 0b00100, 0b01010, 0b10001, 0b10001]),
        'Y': new Uint8Array([0b10001, 0b10001, 0b01010, 0b00100, 0b00100, 0b00100, 0b00100]),
        'Z': new Uint8Array([0b11111, 0b00001, 0b00010, 0b00100, 0b01000, 0b10000, 0b11111]),

        'a': new Uint8Array([0b00000, 0b00000, 0b01110, 0b00001, 0b01111, 0b10001, 0b01111]),
        'b': new Uint8Array([0b10000, 0b10000, 0b10110, 0b11001, 0b10001, 0b10001, 0b11110]),
        'c': new Uint8Array([0b00000, 0b00000, 0b01110, 0b10000, 0b10000, 0b10001, 0b01110]),
        'd': new Uint8Array([0b00001, 0b00001, 0b01101, 0b10011, 0b10001, 0b10001, 0b01111]),
        'e': new Uint8Array([0b00000, 0b00000, 0b01110, 0b10001, 0b11111, 0b10000, 0b01110]),
        'f': new Uint8Array([0b00110, 0b01001, 0b01000, 0b11100, 0b01000, 0b01000, 0b01000]),
        'g': new Uint8Array([0b00000, 0b01111, 0b10001, 0b10001, 0b01111, 0b00001, 0b01110]),
        'h': new Uint8Array([0b10000, 0b10000, 0b10110, 0b11001, 0b10001, 0b10001, 0b10001]),
        'i': new Uint8Array([0b00100, 0b00000, 0b01100, 0b00100, 0b00100, 0b00100, 0b01110]),
        'j': new Uint8Array([0b00010, 0b00000, 0b00110, 0b00010, 0b00010, 0b10010, 0b01100]),
        'k': new Uint8Array([0b10000, 0b10000, 0b10010, 0b10100, 0b11000, 0b10100, 0b10010]),
        'l': new Uint8Array([0b11000, 0b01000, 0b01000, 0b01000, 0b01000, 0b01001, 0b00110]),
        'm': new Uint8Array([0b00000, 0b00000, 0b11010, 0b10101, 0b10101, 0b10101, 0b10101]),
        'n': new Uint8Array([0b00000, 0b00000, 0b10110, 0b11001, 0b10001, 0b10001, 0b10001]),
        'o': new Uint8Array([0b00000, 0b00000, 0b01110, 0b10001, 0b10001, 0b10001, 0b01110]),
        'p': new Uint8Array([0b00000, 0b00000, 0b11110, 0b10001, 0b10001, 0b11110, 0b10000]),
        'q': new Uint8Array([0b00000, 0b00000, 0b01111, 0b10001, 0b10001, 0b01111, 0b00001]),
        'r': new Uint8Array([0b00000, 0b00000, 0b10110, 0b11001, 0b10000, 0b10000, 0b10000]),
        's': new Uint8Array([0b00000, 0b00000, 0b01111, 0b10000, 0b01110, 0b00001, 0b11110]),
        't': new Uint8Array([0b01000, 0b01000, 0b11100, 0b01000, 0b01000, 0b01001, 0b00110]),
        'u': new Uint8Array([0b00000, 0b00000, 0b10001, 0b10001, 0b10001, 0b10011, 0b01101]),
        'v': new Uint8Array([0b00000, 0b00000, 0b10001, 0b10001, 0b10001, 0b01010, 0b00100]),
        'w': new Uint8Array([0b00000, 0b00000, 0b10001, 0b10001, 0b10101, 0b10101, 0b01010]),
        'x': new Uint8Array([0b00000, 0b00000, 0b10001, 0b01010, 0b00100, 0b01010, 0b10001]),
        'y': new Uint8Array([0b00000, 0b00000, 0b10001, 0b10001, 0b01111, 0b00001, 0b01110]),
        'z': new Uint8Array([0b00000, 0b00000, 0b11111, 0b00010, 0b00100, 0b01000, 0b11111]),
        
        '0': new Uint8Array([0b01110, 0b10001, 0b10011, 0b10101, 0b11001, 0b10001, 0b01110]),
        '1': new Uint8Array([0b00100, 0b01100, 0b00100, 0b00100, 0b00100, 0b00100, 0b01110]),
        '2': new Uint8Array([0b01110, 0b10001, 0b00001, 0b00010, 0b00100, 0b01000, 0b11111]),
        '3': new Uint8Array([0b11111, 0b00010, 0b00100, 0b00010, 0b00001, 0b10001, 0b01110]),
        '4': new Uint8Array([0b00010, 0b00110, 0b01010, 0b10010, 0b11111, 0b00010, 0b00010]),
        '5': new Uint8Array([0b11111, 0b10000, 0b11110, 0b00001, 0b00001, 0b10001, 0b01110]),
        '6': new Uint8Array([0b00110, 0b01000, 0b10000, 0b11110, 0b10001, 0b10001, 0b01110]),
        '7': new Uint8Array([0b11111, 0b00001, 0b00010, 0b00100, 0b01000, 0b10000, 0b10000]),
        '8': new Uint8Array([0b01110, 0b10001, 0b10001, 0b01110, 0b10001, 0b10001, 0b01110]),
        '9': new Uint8Array([0b01110, 0b10001, 0b10001, 0b01111, 0b00001, 0b00010, 0b01100]),

        ' ': new Uint8Array([0b00000, 0b00000, 0b00000, 0b00000, 0b00000, 0b00000, 0b00000]),
        '!': new Uint8Array([0b00100, 0b00100, 0b00100, 0b00100, 0b00000, 0b00000, 0b00100]),
        '"': new Uint8Array([0b01010, 0b01010, 0b00000, 0b00000, 0b00000, 0b00000, 0b00000]),
        '#': new Uint8Array([0b01010, 0b11111, 0b01010, 0b01010, 0b11111, 0b01010, 0b01010]),
        '$': new Uint8Array([0b00100, 0b01111, 0b10100, 0b01110, 0b00101, 0b11110, 0b00100]),
        '%': new Uint8Array([0b11000, 0b11001, 0b00010, 0b00100, 0b01000, 0b10011, 0b00011]),
        '&': new Uint8Array([0b01000, 0b10100, 0b10100, 0b01000, 0b10101, 0b10010, 0b01101]),
        '\'': new Uint8Array([0b00100, 0b00100, 0b00000, 0b00000, 0b00000, 0b00000, 0b00000]),
        '(': new Uint8Array([0b00010, 0b00100, 0b01000, 0b01000, 0b01000, 0b00100, 0b00010]),
        ')': new Uint8Array([0b01000, 0b00100, 0b00010, 0b00010, 0b00010, 0b00100, 0b01000]),
        '*': new Uint8Array([0b00000, 0b00100, 0b10101, 0b01110, 0b10101, 0b00100, 0b00000]),
        '+': new Uint8Array([0b00000, 0b00100, 0b00100, 0b11111, 0b00100, 0b00100, 0b00000]),
        ',': new Uint8Array([0b00000, 0b00000, 0b00000, 0b00000, 0b00100, 0b00100, 0b01000]),
        '-': new Uint8Array([0b00000, 0b00000, 0b00000, 0b11111, 0b00000, 0b00000, 0b00000]),
        '.': new Uint8Array([0b00000, 0b00000, 0b00000, 0b00000, 0b00000, 0b00100, 0b00100]),
        '/': new Uint8Array([0b00001, 0b00010, 0b00100, 0b01000, 0b10000, 0b00000, 0b00000]),
        ':': new Uint8Array([0b00000, 0b00100, 0b00100, 0b00000, 0b00100, 0b00100, 0b00000]),
        ',': new Uint8Array([0b00000, 0b00100, 0b00100, 0b00000, 0b00100, 0b00100, 0b01000]),
        '<': new Uint8Array([0b00010, 0b00100, 0b01000, 0b10000, 0b01000, 0b00100, 0b00010]),
        '=': new Uint8Array([0b00000, 0b11111, 0b00000, 0b11111, 0b00000, 0b00000, 0b00000]),
        '>': new Uint8Array([0b10000, 0b01000, 0b00100, 0b00010, 0b00100, 0b01000, 0b10000]),
        '?': new Uint8Array([0b01110, 0b10001, 0b00001, 0b00010, 0b00100, 0b00000, 0b00100]),
        '@': new Uint8Array([0b01110, 0b10001, 0b00001, 0b01101, 0b10101, 0b10101, 0b01110]),
        '[': new Uint8Array([0b01110, 0b01000, 0b01000, 0b01000, 0b01000, 0b01000, 0b01110]),
        '\\': new Uint8Array([0b10000, 0b01000, 0b00100, 0b00010, 0b00001, 0b00000, 0b00000]),
        ']': new Uint8Array([0b01110, 0b00010, 0b00010, 0b00010, 0b00010, 0b00010, 0b01110]),
        '^': new Uint8Array([0b00100, 0b01010, 0b10001, 0b00000, 0b00000, 0b00000, 0b00000]),
        '_': new Uint8Array([0b00000, 0b00000, 0b00000, 0b00000, 0b00000, 0b00000, 0b11111]),
        '`': new Uint8Array([0b01000, 0b00100, 0b00010, 0b00000, 0b00000, 0b00000, 0b00000]),
        '{': new Uint8Array([0b00010, 0b00100, 0b00100, 0b01000, 0b00100, 0b00100, 0b00010]),
        '|': new Uint8Array([0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100]),
        '}': new Uint8Array([0b01000, 0b00100, 0b00100, 0b00010, 0b00100, 0b00100, 0b01000]),
        '~': new Uint8Array([0b00000, 0b00000, 0b00000, 0b01001, 0b10110, 0b00000, 0b00000])
        };
    }

    characterExists(c) {
        const bitmap = this.FONT_MAP[c];
        if (!bitmap)
            return false;
        return true;
    }

    getPixel(c, i, j) {
        const bitmap = this.FONT_MAP[c];
        if (!bitmap || j < 0 || j >= bitmap.length || i < 0 || i >= 5) return 0;
        const row = bitmap[j];
        return (row >> i) & 1;
    }
}

