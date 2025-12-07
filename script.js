// Memristor Models
class BiolekMemristor {
    constructor(mu_v = 1e-9, D = 1e-8, R_ON = 100, R_OFF = 16e3, w_init = null, p = 1) {
        this.mu_v = mu_v;
        this.D = D;
        this.R_ON = R_ON;
        this.R_OFF = R_OFF;
        this.w = w_init !== null ? Math.max(0, Math.min(w_init, this.D)) : D / 2;
        this.p = p;
        this._updateMemristance();
    }

    _updateMemristance() {
        this.M = this.R_ON * (this.w / this.D) + this.R_OFF * (1 - (this.w / this.D));
        this.W = 1.0 / this.M;
    }

    stepFunction(i) {
        return i >= 0 ? 1.0 : 0.0;
    }

    windowFunction(x, i) {
        return 1.0 - Math.pow(x - this.stepFunction(-i), 2 * this.p);
    }

    updateState(V, dt) {
        const i = this.getCurrent(V);
        const f_w = this.windowFunction(this.w / this.D, i);
        const dw_dt = this.mu_v * (this.R_ON / this.D) * i * f_w;
        this.w += dw_dt * dt;
        this.w = Math.max(0, Math.min(this.w, this.D));
        this._updateMemristance();
    }

    getCurrent(V) {
        return V / this.M;
    }
}

class JoglekarMemristor {
    constructor(mu_v = 1e-9, D = 1e-8, R_ON = 100, R_OFF = 16e3, w_init = null, p = 1) {
        this.mu_v = mu_v;
        this.D = D;
        this.R_ON = R_ON;
        this.R_OFF = R_OFF;
        this.w = w_init !== null ? Math.max(0, Math.min(w_init, this.D)) : D / 2;
        this.p = p;
        this._updateMemristance();
    }

    _updateMemristance() {
        this.M = this.R_ON * (this.w / this.D) + this.R_OFF * (1 - (this.w / this.D));
        this.W = 1.0 / this.M;
    }

    windowFunction(x) {
        return 1 - Math.pow(2 * x - 1, 2 * this.p);
    }

    updateState(V, dt) {
        const i = this.getCurrent(V);
        const f_w = this.windowFunction(this.w / this.D);
        const dw_dt = this.mu_v * (this.R_ON / this.D) * i * f_w;
        this.w += dw_dt * dt;
        this.w = Math.max(0, Math.min(this.w, this.D));
        this._updateMemristance();
    }

    getCurrent(V) {
        return V / this.M;
    }
}

class LinearIonDriftMemristor {
    constructor(mu_v = 1e-9, D = 1e-8, R_ON = 100, R_OFF = 16e3, w_init = null) {
        this.mu_v = mu_v;
        this.D = D;
        this.R_ON = R_ON;
        this.R_OFF = R_OFF;
        this.w = w_init !== null ? Math.max(0, Math.min(w_init, this.D)) : D / 2;
        this._updateMemristance();
    }

    _updateMemristance() {
        this.M = this.R_ON * (this.w / this.D) + this.R_OFF * (1 - (this.w / this.D));
        this.W = 1.0 / this.M;
    }

    updateState(V, dt) {
        const i = this.getCurrent(V);
        const dw_dt = (this.mu_v * this.R_ON / this.D) * i;
        this.w += dw_dt * dt;
        this.w = Math.max(0, Math.min(this.w, this.D));
        this._updateMemristance();
    }

    getCurrent(V) {
        return V / this.M;
    }
}

class VTEAMMemristor {
    constructor(k_off = 5e-4, k_on = -10, alpha_off = 3, alpha_on = 1, w_off = 3e-9, w_on = 0, w_init = 0, a_off = 0.8, a_on = 0.2, w_c = 0.12, u_off = 0.5, u_on = -0.5, R_on = 100, R_off = 2.5e3) {
        this.k_off = k_off; this.k_on = k_on;
        this.alpha_off = alpha_off; this.alpha_on = alpha_on;
        this.w_off = w_off; this.w_on = w_on;
        this.w = Math.min(this.w_off, Math.max(this.w_on, w_init));
        this.a_off = a_off; this.a_on = a_on; this.w_c = w_c;
        this.u_off = u_off; this.u_on = u_on;
        this.R_on = R_on; this.R_off = R_off;
        this._updateMemristance();
    }

    f_off(w) {
        return Math.exp(-Math.exp((w - this.a_off) / this.w_c));
    }
    f_on(w) {
        return Math.exp(-Math.exp(-(w - this.a_on) / this.w_c));
    }

    _updateMemristance() {
        const λ = Math.log(this.R_off / this.R_on);
        this.M = this.R_on * Math.exp((λ / (this.w_off - this.w_on)) * (this.w - this.w_on));
        this.W = 1 / this.M;
    }

    updateState(u, dt) {
        let dw_dt = 0;
        if (0 < this.u_off && u > this.u_off) {
            dw_dt = this.k_off * Math.pow(u / this.u_off - 1, this.alpha_off) * this.f_off(this.w);
        } else if (u < this.u_on && this.u_on < 0) {
            dw_dt = this.k_on * Math.pow(u / this.u_on - 1, this.alpha_on) * this.f_on(this.w);
        }
        this.w = Math.min(this.w_off, Math.max(this.w_on, this.w + dw_dt * dt));
        this._updateMemristance();
    }

    getCurrent(u) {
        const λ = Math.log(this.R_off / this.R_on);
        const expTerm = Math.exp((-λ / (this.w_off - this.w_on)) * (this.w - this.w_on));
        return (u / this.R_on) * expTerm;
    }
}

class MMSMemristor {
    constructor(R_on = 500, R_off = 1500, U_on = 0.27, U_off = 0.27, tau = 1e-4, T = 298.5, x_init = 0) {
        this.R_on = R_on;
        this.R_off = R_off;
        this.U_on = U_on;
        this.U_off = U_off;
        this.tau = tau;
        this.T = T;
        this.q = 1.602176634e-19;
        this.k = 1.380649e-23;
        this.x = Math.max(0, Math.min(x_init, 1));
        this._updateMemristance();
    }

    _updateMemristance() {
        this.W = this.x / this.R_on + (1.0 - this.x) / this.R_off;
        this.M = 1.0 / this.W;
    }

    updateState(V, dt) {
        const alpha = dt / this.tau;
        const beta = this.q / (this.k * this.T);
        const P_on = alpha / (1 + Math.exp(-beta * (V - this.U_on)));
        const P_off = alpha * (1 - 1 / (1 + Math.exp(-beta * (V + this.U_off))));
        const N_on = P_on * (1 - this.x);
        const N_off = P_off * this.x;
        this.x += N_on - N_off;
        this.x = Math.max(0, Math.min(this.x, 1));
        this._updateMemristance();
    }

    getCurrent(V) {
        return V * this.W;
    }
}

class YakopcicMemristor {
    constructor(A_p = 4000, A_n = 4000, U_p = 0.5, U_n = 0.5, alpha_p = 1, alpha_n = 5, x_p = 0.3, x_n = 0.3, a1 = 0.17, a2 = 0.17, b = 0.05, x_init = 0, x_on = 0) {
        this.A_p = A_p; this.A_n = A_n;
        this.U_p = U_p; this.U_n = U_n;
        this.alpha_p = alpha_p; this.alpha_n = alpha_n;
        this.x_p = x_p; this.x_n = x_n;
        this.a1 = a1; this.a2 = a2;
        this.b = b;
        this.x_on = x_on;
        this.x = Math.min(1, Math.max(this.x_on, x_init));
        this._updateMemristance();
    }

    g(u) {
        if (u > this.U_p) return this.A_p * (Math.exp(u) - Math.exp(this.U_p));
        if (u < -this.U_n) return -this.A_n * (Math.exp(-u) - Math.exp(this.U_n));
        return 0;
    }

    f_p(x) {
        if (x >= this.x_p) {
            const wp = (this.x_p - x) / (1 - this.x_p) + 1;
            return Math.exp(-this.alpha_p * (x - this.x_p)) * wp;
        }
        return 1;
    }
    f_n(x) {
        if (x <= (1 - this.x_n)) {
            const wn = x / (1 - this.x_n);
            return Math.exp(this.alpha_n * (x + this.x_n - 1)) * wn;
        }
        return 1;
    }

    f(x, u) {
        return u >= 0 ? this.f_p(x) : this.f_n(x);
    }

    _updateMemristance() {
        this.W = this.a1 * this.x * Math.sinh(this.b);
        this.M = 1 / this.W;
    }

    updateState(u, dt) {
        const dx = this.g(u) * this.f(this.x, u);
        this.x = Math.min(1, Math.max(this.x_on, this.x + dx * dt));
        this._updateMemristance();
    }

    getCurrent(u) {
        if (u >= 0) return this.a1 * this.x * Math.sinh(this.b * u);
        else return this.a2 * this.x * Math.sinh(this.b * u);
    }
}

// Waveform generators
function generateWaveform(type, t, frequency, amplitude) {
    const omega = 2 * Math.PI * frequency;

    switch (type) {
        case 'sine':
            return t.map(time => amplitude * Math.sin(omega * time));
        case 'square':
            return t.map(time => amplitude * Math.sign(Math.sin(omega * time)));
        case 'triangle':
            return t.map(time => {
                const phase = (omega * time) % (2 * Math.PI);
                if (phase < Math.PI) {
                    return amplitude * (2 * phase / Math.PI - 1);
                } else {
                    return amplitude * (3 - 2 * phase / Math.PI);
                }
            });
        case 'sawtooth':
            return t.map(time => {
                const phase = (omega * time) % (2 * Math.PI);
                return amplitude * (2 * phase / (2 * Math.PI) - 1);
            });
        default:
            return t.map(time => amplitude * Math.sin(omega * time));
    }
}

// Simulation function
function simulateMemristor(memristor, V_seq, dt) {
    const I_mem = [];

    for (let i = 0; i < V_seq.length; i++) {
        const current = memristor.getCurrent(V_seq[i]);
        I_mem.push(current);
        memristor.updateState(V_seq[i], dt);
    }

    return I_mem;
}

// Theme Management
function getSystemTheme() {
    return window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
}

function initializeTheme() {
    const savedTheme = localStorage.getItem('theme');
    const theme = savedTheme || getSystemTheme();
    document.documentElement.setAttribute('data-theme', theme);
    updateThemeButton(theme);

    // Listen for system theme changes only if user hasn't set a preference
    if (!savedTheme) {
        window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', (e) => {
            if (!localStorage.getItem('theme')) {
                const newTheme = e.matches ? 'dark' : 'light';
                document.documentElement.setAttribute('data-theme', newTheme);
                updateThemeButton(newTheme);
                if (window.simulatorUI) {
                    window.simulatorUI.runSimulation();
                }
            }
        });
    }
}

function toggleTheme() {
    const currentTheme = document.documentElement.getAttribute('data-theme');
    const newTheme = currentTheme === 'dark' ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', newTheme);
    localStorage.setItem('theme', newTheme);
    updateThemeButton(newTheme);

    // Re-render the chart with new theme colors
    if (window.simulatorUI) {
        window.simulatorUI.runSimulation();
    }
}

function updateThemeButton(theme) {
    const button = document.getElementById('themeToggle');
    if (button) {
        const icon = button.querySelector('.icon');
        const text = button.querySelector('.theme-text');
        if (theme === 'dark') {
            icon.className = 'fas fa-moon icon';
            text.textContent = 'Dark';
        } else {
            icon.className = 'fas fa-sun icon';
            text.textContent = 'Light';
        }
    }
}

// UI Management
class SimulatorUI {
    constructor() {
        this.initializeElements();
        this.setupEventListeners();
        this.updateAllDisplays();
        this.initializeChart();
        this.updateModelControls();
        this.runSimulation();
    }

    initializeElements() {
        // Get all UI elements
        this.elements = {
            modelSelect: document.getElementById('modelSelect'),
            waveformSelect: document.getElementById('waveformSelect'),
            frequency: document.getElementById('frequency'),
            amplitude: document.getElementById('amplitude'),
            muV: document.getElementById('muV'),
            D: document.getElementById('D'),
            rOn: document.getElementById('rOn'),
            rOff: document.getElementById('rOff'),
            uOn: document.getElementById('uOn'),
            uOff: document.getElementById('uOff'),
            tau: document.getElementById('tau'),
            temperature: document.getElementById('temperature'),
            p: document.getElementById('p'),
            resetBtn: document.getElementById('resetBtn'),
            statusText: document.getElementById('statusText'),
            pControl: document.getElementById('pControl'),
            vtKOff: document.getElementById('vtKOff'),
            vtKOn: document.getElementById('vtKOn'),
            vtAlphaOff: document.getElementById('vtAlphaOff'),
            vtAlphaOn: document.getElementById('vtAlphaOn'),
            vtWOff: document.getElementById('vtWOff'),
            vtWOn: document.getElementById('vtWOn'),
            vtAOff: document.getElementById('vtAOff'),
            vtAOn: document.getElementById('vtAOn'),
            vtWC: document.getElementById('vtWC'),
            vtUOff: document.getElementById('vtUOff'),
            vtUOn: document.getElementById('vtUOn'),
            ycAp: document.getElementById('ycAp'),
            ycAn: document.getElementById('ycAn'),
            ycUp: document.getElementById('ycUp'),
            ycUn: document.getElementById('ycUn'),
            ycAlphaP: document.getElementById('ycAlphaP'),
            ycAlphaN: document.getElementById('ycAlphaN'),
            ycXp: document.getElementById('ycXp'),
            ycXn: document.getElementById('ycXn'),
            ycA1: document.getElementById('ycA1'),
            ycA2: document.getElementById('ycA2'),
            ycB: document.getElementById('ycB'),
            ycXOn: document.getElementById('ycXOn'),
            duration: document.getElementById('duration'),
            timeStep: document.getElementById('timeStep')
        };

        // Value display elements
        this.valueDisplays = {
            freqValue: document.getElementById('freqValue'),
            ampValue: document.getElementById('ampValue'),
            muVValue: document.getElementById('muVValue'),
            DValue: document.getElementById('DValue'),
            rOnValue: document.getElementById('rOnValue'),
            rOffValue: document.getElementById('rOffValue'),
            uOnValue: document.getElementById('uOnValue'),
            uOffValue: document.getElementById('uOffValue'),
            tauValue: document.getElementById('tauValue'),
            temperatureValue: document.getElementById('temperatureValue'),
            pValue: document.getElementById('pValue'),
            vtKOffValue: document.getElementById('vtKOffValue'),
            vtKOnValue: document.getElementById('vtKOnValue'),
            vtAlphaOffValue: document.getElementById('vtAlphaOffValue'),
            vtAlphaOnValue: document.getElementById('vtAlphaOnValue'),
            vtWOffValue: document.getElementById('vtWOffValue'),
            vtWOnValue: document.getElementById('vtWOnValue'),
            vtAOffValue: document.getElementById('vtAOffValue'),
            vtAOnValue: document.getElementById('vtAOnValue'),
            vtWCValue: document.getElementById('vtWCValue'),
            vtUOffValue: document.getElementById('vtUOffValue'),
            vtUOnValue: document.getElementById('vtUOnValue'),
            ycApValue: document.getElementById('ycApValue'),
            ycAnValue: document.getElementById('ycAnValue'),
            ycUpValue: document.getElementById('ycUpValue'),
            ycUnValue: document.getElementById('ycUnValue'),
            ycAlphaPValue: document.getElementById('ycAlphaPValue'),
            ycAlphaNValue: document.getElementById('ycAlphaNValue'),
            ycXpValue: document.getElementById('ycXpValue'),
            ycXnValue: document.getElementById('ycXnValue'),
            ycA1Value: document.getElementById('ycA1Value'),
            ycA2Value: document.getElementById('ycA2Value'),
            ycBValue: document.getElementById('ycBValue'),
            ycXOnValue: document.getElementById('ycXOnValue'),
            durationValue: document.getElementById('durationValue'),
            timeStepValue: document.getElementById('timeStepValue')
        };
    }

    setupEventListeners() {
        // Theme toggle
        const themeToggle = document.getElementById('themeToggle');
        if (themeToggle) {
            themeToggle.addEventListener('click', toggleTheme);
        }

        // Model selection
        this.elements.modelSelect.addEventListener('change', () => {
            this.updateModelControls();
            this.runSimulation();
        });

        // Real-time parameter updates
        Object.keys(this.elements).forEach(key => {
            if (this.elements[key].type === 'range' || this.elements[key].tagName === 'SELECT') {
                this.elements[key].addEventListener('input', () => {
                    this.updateAllDisplays();
                    this.runSimulation();
                });
            }
        });

        // Buttons
        this.elements.resetBtn.addEventListener('click', () => this.resetSimulation());
    }

    updateModelControls() {
        const m = this.elements.modelSelect.value;

        // Hide all model-specific parameter rows
        document.querySelectorAll(
            '.param-linear, .param-biolek, .param-joglekar, .param-vteam, .param-mms, .param-yakopcic'
        ).forEach(el => el.style.display = 'none');

        // Show only those belonging to the selected model
        document.querySelectorAll('.param-' + m)
            .forEach(el => el.style.display = 'flex');
    }

    updateAllDisplays() {
        // Update frequency display
        const freq = parseFloat(this.elements.frequency.value);
        this.valueDisplays.freqValue.textContent = freq >= 1000 ? `${(freq / 1000).toFixed(0)}k` : freq.toString();

        // Update amplitude display
        this.valueDisplays.ampValue.textContent = parseFloat(this.elements.amplitude.value).toFixed(1);

        // Update scientific notation displays
        this.valueDisplays.muVValue.textContent = parseFloat(this.elements.muV.value).toExponential(0);
        this.valueDisplays.DValue.textContent = parseFloat(this.elements.D.value).toExponential(0);
        this.valueDisplays.timeStepValue.textContent = parseFloat(this.elements.timeStep.value).toExponential(0);

        // Update resistance displays
        const rOff = parseFloat(this.elements.rOff.value);
        this.valueDisplays.rOnValue.textContent = this.elements.rOn.value;
        this.valueDisplays.rOffValue.textContent = rOff >= 1000 ? `${(rOff / 1000).toFixed(0)}k` : rOff.toString();

        // Update MMS model displays
        this.valueDisplays.uOnValue.textContent = parseFloat(this.elements.uOn.value).toFixed(2);
        this.valueDisplays.uOffValue.textContent = parseFloat(this.elements.uOff.value).toFixed(2);
        this.valueDisplays.tauValue.textContent = parseFloat(this.elements.tau.value).toExponential(0);
        this.valueDisplays.temperatureValue.textContent = parseFloat(this.elements.temperature.value).toFixed(1);

        // Update VTEAM model displays
        this.valueDisplays.vtKOffValue.textContent = parseFloat(this.elements.vtKOff.value).toExponential(1);
        this.valueDisplays.vtKOnValue.textContent = parseFloat(this.elements.vtKOn.value).toFixed(1);
        this.valueDisplays.vtAlphaOffValue.textContent = parseFloat(this.elements.vtAlphaOff.value).toFixed(1);
        this.valueDisplays.vtAlphaOnValue.textContent = parseFloat(this.elements.vtAlphaOn.value).toFixed(1);
        this.valueDisplays.vtWOffValue.textContent = parseFloat(this.elements.vtWOff.value).toExponential(1);
        this.valueDisplays.vtWOnValue.textContent = parseFloat(this.elements.vtWOn.value).toExponential(1);
        this.valueDisplays.vtAOffValue.textContent = parseFloat(this.elements.vtAOff.value).toFixed(2);
        this.valueDisplays.vtAOnValue.textContent = parseFloat(this.elements.vtAOn.value).toFixed(2);
        this.valueDisplays.vtWCValue.textContent = parseFloat(this.elements.vtWC.value).toFixed(2);
        this.valueDisplays.vtUOffValue.textContent = parseFloat(this.elements.vtUOff.value).toFixed(2);
        this.valueDisplays.vtUOnValue.textContent = parseFloat(this.elements.vtUOn.value).toFixed(2);

        // Update Yakopcic model displays
        this.valueDisplays.ycApValue.textContent = parseFloat(this.elements.ycAp.value).toFixed(0);
        this.valueDisplays.ycAnValue.textContent = parseFloat(this.elements.ycAn.value).toFixed(0);
        this.valueDisplays.ycUpValue.textContent = parseFloat(this.elements.ycUp.value).toFixed(2);
        this.valueDisplays.ycUnValue.textContent = parseFloat(this.elements.ycUn.value).toFixed(2);
        this.valueDisplays.ycAlphaPValue.textContent = parseFloat(this.elements.ycAlphaP.value).toFixed(1);
        this.valueDisplays.ycAlphaNValue.textContent = parseFloat(this.elements.ycAlphaN.value).toFixed(1);
        this.valueDisplays.ycXpValue.textContent = parseFloat(this.elements.ycXp.value).toFixed(2);
        this.valueDisplays.ycXnValue.textContent = parseFloat(this.elements.ycXn.value).toFixed(2);
        this.valueDisplays.ycA1Value.textContent = parseFloat(this.elements.ycA1.value).toFixed(2);
        this.valueDisplays.ycA2Value.textContent = parseFloat(this.elements.ycA2.value).toFixed(2);
        this.valueDisplays.ycBValue.textContent = parseFloat(this.elements.ycB.value).toFixed(3);
        this.valueDisplays.ycXOnValue.textContent = parseFloat(this.elements.ycXOn.value).toFixed(2);

        // Update other displays
        this.valueDisplays.pValue.textContent = this.elements.p.value;
        this.valueDisplays.durationValue.textContent = parseFloat(this.elements.duration.value).toFixed(1);
    }

    getChartColors() {
        const isDark = document.documentElement.getAttribute('data-theme') === 'dark';
        return {
            gridColor: isDark ? '#334155' : '#e2e8f0',
            zerolineColor: isDark ? '#475569' : '#cbd5e0',
            backgroundColor: isDark ? '#1e293b' : '#ffffff',
            textColor: isDark ? '#f1f5f9' : '#2d3748',
            lineColor: isDark ? '#818cf8' : '#4f46e5'
        };
    }

    initializeChart() {
        const colors = this.getChartColors();
        const layout = {
            title: false,
            xaxis: {
                title: 'Voltage (V)',
                gridcolor: colors.gridColor,
                zeroline: true,
                zerolinecolor: colors.zerolineColor,
                color: colors.textColor
            },
            yaxis: {
                title: 'Current (A)',
                gridcolor: colors.gridColor,
                zeroline: true,
                zerolinecolor: colors.zerolineColor,
                color: colors.textColor
            },
            plot_bgcolor: colors.backgroundColor,
            paper_bgcolor: colors.backgroundColor,
            font: {
                family: 'Inter, sans-serif',
                size: 12,
                color: colors.textColor
            },
            margin: { t: 20, r: 20, b: 60, l: 80 }
        };

        const config = {
            responsive: true,
            displayModeBar: false
        };

        Plotly.newPlot('chart', [], layout, config);
    }

    runSimulation() {
        this.elements.statusText.textContent = 'Running simulation...';

        setTimeout(() => {
            try {
                // Get parameters
                const modelType = this.elements.modelSelect.value;
                const waveform = this.elements.waveformSelect.value;
                const frequency = parseFloat(this.elements.frequency.value);
                const amplitude = parseFloat(this.elements.amplitude.value);
                const mu_v = parseFloat(this.elements.muV.value);
                const D = parseFloat(this.elements.D.value);
                const R_ON = parseFloat(this.elements.rOn.value);
                const R_OFF = parseFloat(this.elements.rOff.value);
                const U_ON = parseFloat(this.elements.uOn.value);
                const U_OFF = parseFloat(this.elements.uOff.value);
                const tau = parseFloat(this.elements.tau.value);
                const temperature = parseFloat(this.elements.temperature.value);
                const k_off = parseFloat(this.elements.vtKOff.value);
                const k_on = parseFloat(this.elements.vtKOn.value);
                const alpha_off = parseFloat(this.elements.vtAlphaOff.value);
                const alpha_on = parseFloat(this.elements.vtAlphaOn.value);
                const w_off = parseFloat(this.elements.vtWOff.value);
                const w_on = parseFloat(this.elements.vtWOn.value);
                const a_off = parseFloat(this.elements.vtAOff.value);
                const a_on = parseFloat(this.elements.vtAOn.value);
                const w_c = parseFloat(this.elements.vtWC.value);
                const u_off = parseFloat(this.elements.vtUOff.value);
                const u_on = parseFloat(this.elements.vtUOn.value);
                const p = parseFloat(this.elements.p.value);
                const A_p = parseFloat(this.elements.ycAp.value);
                const A_n = parseFloat(this.elements.ycAn.value);
                const U_p = parseFloat(this.elements.ycUp.value);
                const U_n = parseFloat(this.elements.ycUn.value);
                const alpha_p = parseFloat(this.elements.ycAlphaP.value);
                const alpha_n = parseFloat(this.elements.ycAlphaN.value);
                const x_p = parseFloat(this.elements.ycXp.value);
                const x_n = parseFloat(this.elements.ycXn.value);
                const a1 = parseFloat(this.elements.ycA1.value);
                const a2 = parseFloat(this.elements.ycA2.value);
                const b = parseFloat(this.elements.ycB.value);
                const x_on = parseFloat(this.elements.ycXOn.value);
                const duration = parseFloat(this.elements.duration.value) * 1e-6; // Convert to seconds
                const dt = parseFloat(this.elements.timeStep.value) * 1e-6; // Convert to seconds

                // Create memristor
                let memristor;
                switch (modelType) {
                    case 'biolek':
                        memristor = new BiolekMemristor(mu_v, D, R_ON, R_OFF, null, p);
                        break;
                    case 'joglekar':
                        memristor = new JoglekarMemristor(mu_v, D, R_ON, R_OFF, null, p);
                        break;
                    case 'linear':
                        memristor = new LinearIonDriftMemristor(mu_v, D, R_ON, R_OFF);
                        break;
                    case 'vteam':
                        memristor = new VTEAMMemristor(k_off, k_on, alpha_off, alpha_on, w_off, w_on, null, a_off, a_on, w_c, u_off, u_on, R_ON, R_OFF);
                        break;
                    case 'mms':
                        memristor = new MMSMemristor(R_ON, R_OFF, U_ON, U_OFF, tau, temperature);
                        break;
                    case 'yakopcic':
                        memristor = new YakopcicMemristor(A_p, A_n, U_p, U_n, alpha_p, alpha_n, x_p, x_n, a1, a2, b, 0.11, x_on);
                        break;
                }

                // Generate time vector
                const numPoints = Math.floor(duration / dt);
                const t = Array.from({ length: numPoints }, (_, i) => i * dt);

                // Generate voltage sequence
                const V_seq = generateWaveform(waveform, t, frequency, amplitude);

                // Run simulation
                const I_mem = simulateMemristor(memristor, V_seq, dt);

                // Plot results
                this.plotHysteresis(V_seq, I_mem);

                this.elements.statusText.textContent = `Simulation complete (${numPoints} points)`;
            } catch (error) {
                console.error('Simulation error:', error);
                this.elements.statusText.textContent = 'Simulation error';
            }
        }, 10);
    }

    plotHysteresis(V_seq, I_mem) {
        const colors = this.getChartColors();
        const trace = {
            x: V_seq,
            y: I_mem,
            type: 'scatter',
            mode: 'lines',
            line: {
                color: colors.lineColor,
                width: 2
            },
            name: 'I-V Curve'
        };

        const layout = {
            title: false,
            xaxis: {
                title: 'Voltage (V)',
                gridcolor: colors.gridColor,
                zeroline: true,
                zerolinecolor: colors.zerolineColor,
                color: colors.textColor
            },
            yaxis: {
                title: 'Current (A)',
                gridcolor: colors.gridColor,
                zeroline: true,
                zerolinecolor: colors.zerolineColor,
                color: colors.textColor
            },
            plot_bgcolor: colors.backgroundColor,
            paper_bgcolor: colors.backgroundColor,
            font: {
                family: 'Inter, sans-serif',
                size: 12,
                color: colors.textColor
            },
            margin: { t: 20, r: 20, b: 60, l: 80 },
            showlegend: false
        };

        Plotly.newPlot('chart', [trace], layout, { responsive: true, displayModeBar: false });
    }

    resetSimulation() {
        // Reset all controls to default values
        Object.values(this.elements).forEach(el => {
            if (el.tagName === 'INPUT' && el.type === 'range') {
                el.value = el.defaultValue;
            }
            else if (el.tagName === 'SELECT') {
                // restore the originally selected option(s)
                Array.from(el.options).forEach(opt => {
                    opt.selected = opt.defaultSelected;
                });
            }
        });

        this.updateAllDisplays();
        this.updateModelControls();
        this.initializeChart();
        this.elements.statusText.textContent = 'Reset complete - Ready to simulate';
        this.runSimulation();
    }
}

// Initialize the application
document.addEventListener('DOMContentLoaded', () => {
    initializeTheme();
    window.simulatorUI = new SimulatorUI();
});
