/* Latch top-10 viewer: the ten best latch designs of the PC search, one at a time, in 3D.
 *
 * Mounted onto every [data-latch-top10] element (README: proteins/latch_top10/). Built like the biosensor viewer
 * (latch-viewer.js: same card, chips, tools and Mol* 5.13 vendored in ../molstar/, fetched once the figure comes near
 * the viewport; load gv-viewer.css and latch-viewer.css first, then latch-top10.css).
 *
 * Data (built by proteins/latch_top10/tools/prep_top10.py):
 *   designs.json              per design: card facts and one camera frame over both states
 *   dNN_closed.pdb, _open.pdb the design with the latch shut / swung open; the B-factor column holds the part
 *                             (0 binder, 1 latch, 2 target-binding domain, 3 PD-L1), the theme colours by it
 * They load from this script's folder, or from data-base="…/" on the mount element; data-binder="…" renames
 * the binder in the legend (default "GvpA binder").
 *
 * One click on a design button shows that design in the state chosen by the chips (closed, open, both: the open
 * latch is then drawn faint over the closed one). Camera, spin and the chosen state carry over to the next design.
 */
(function () {
  'use strict';

  var HERE = (document.currentScript && document.currentScript.src || '').replace(/[^/]*$/, '');
  var MOLSTAR = HERE + '../molstar/';

  var ICON = {
    reset: '<path d="M4 12a8 8 0 1 0 2.4-5.7" /><path d="M4 4v4.5h4.5" />',
    spin: '<ellipse cx="12" cy="12" rx="9" ry="4" /><path d="M12 3v18" stroke-dasharray="2 2.5" />' +
          '<path d="m16.5 14.6 2.6 1-.8 2.7" />',
    expand: '<path d="M4 9V4h5M15 4h5v5M20 15v5h-5M9 20H4v-5" />',
    shrink: '<path d="M9 4v5H4M20 9h-5V4M15 20v-5h5M4 15h5v5" />'
  };
  function icon(name) {
    return '<svg class="gvv-icon" viewBox="0 0 24 24" aria-hidden="true" focusable="false" fill="none" ' +
           'stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round">' +
           ICON[name] + '</svg>';
  }

  /* what the chips offer, and what each shows */
  var VIEWS = [
    { key: 'closed', label: 'Closed', show: ['closed'] },
    { key: 'open', label: 'Open', show: ['open'] },
    { key: 'both', label: 'Both', show: ['closed', 'ghost'] }
  ];
  var LEGEND = [['--lv-binder', 'binder'], ['--lv-latch', 'latch'], ['--lv-tbd', 'target-binding domain'],
                ['--lv-target', 'PD-L1']];

  function el(tag, cls, html) {
    var e = document.createElement(tag);
    if (cls) e.className = cls;
    if (html != null) e.innerHTML = html;
    return e;
  }
  function stillMotion() {
    var m = document.documentElement.getAttribute('data-motion') || 'auto';
    return m === 'reduced' || (m !== 'full' && matchMedia('(prefers-reduced-motion: reduce)').matches);
  }
  var colourCtx = document.createElement('canvas').getContext('2d');
  function colour(host, name) {
    var v = getComputedStyle(host).getPropertyValue(name).trim();
    colourCtx.fillStyle = '#000';
    colourCtx.fillStyle = v;
    var s = colourCtx.fillStyle;
    if (s.charAt(0) === '#') return parseInt(s.slice(1, 7), 16);
    var m = s.match(/[\d.]+/g) || [0, 0, 0];
    return (Math.round(m[0]) << 16) | (Math.round(m[1]) << 8) | Math.round(m[2]);
  }
  function loadScript(src) {
    return new Promise(function (ok, fail) {
      if (window.molstar) return ok();
      var s = document.createElement('script');
      s.src = src; s.async = true;
      s.onload = function () { ok(); };
      s.onerror = function () { fail(new Error('could not load ' + src)); };
      document.head.appendChild(s);
    });
  }
  function loadStyle(href) {
    if (document.querySelector('link[data-molstar]')) return;
    var l = document.createElement('link');
    l.rel = 'stylesheet'; l.href = href; l.setAttribute('data-molstar', '');
    document.head.appendChild(l);
  }
  function webglOk() {
    try {
      var c = document.createElement('canvas');
      return !!(c.getContext('webgl2') || c.getContext('webgl'));
    } catch (e) { return false; }
  }

  /** "phi72(+61)+phi76(-90)" -> "φ72 +61° · φ76 −90°" */
  function hingeText(moves) {
    return moves.split('+').reduce(function (acc, part) {   // a "+61" also contains a plus: re-join the pieces
      if (/^(phi|psi)/.test(part)) acc.push(part); else acc[acc.length - 1] += '+' + part;
      return acc;
    }, []).map(function (m) {
      var x = m.match(/^(phi|psi)(\d+)\(([+-]\d+)\)$/);
      return x ? (x[1] === 'phi' ? 'φ' : 'ψ') + x[2] + ' ' + x[3].replace('-', '−') + '°' : m;
    }).join(' · ');
  }

  function mount(host) {
    var base = host.getAttribute('data-base') || HERE;
    var state = { design: 0, view: 0, spin: false };

    /* ---- markup ---- */
    host.innerHTML = '';
    host.classList.add('gvv', 'lv', 'lt');
    var picker = el('ol', 'lt-picker');
    picker.setAttribute('aria-label', 'Top 10 latch designs');
    var stage = el('div', 'gvv-stage');
    var canvas = el('div', 'gvv-canvas');
    var card = el('div', 'gvv-card');
    card.setAttribute('aria-live', 'polite');
    var tools = el('div', 'gvv-tools');
    var chips = el('div', 'gvv-chips');
    var legend = el('ul', 'gvv-legend');
    var hint = el('div', 'gvv-hint', 'Use <kbd>Ctrl</kbd> + scroll to zoom');
    var veil = el('div', 'gvv-veil',
      '<span class="gvv-pulse" aria-hidden="true"></span><span class="gvv-veil-text">Preparing the 3D model…</span>');
    stage.append(canvas, tools, chips, legend, hint, veil);
    host.append(picker, stage, card);

    function toolButton(name, label, fn) {
      var b = el('button', 'gvv-tool', icon(name));
      b.type = 'button';
      b.setAttribute('aria-label', label);
      b.setAttribute('data-tool', name);
      b.title = label;
      b.addEventListener('click', fn);
      tools.appendChild(b);
      return b;
    }
    toolButton('reset', 'Reset view', function () { frame(designs[state.design], 700); });
    var bSpin = toolButton('spin', 'Rotate', function () { setSpin(!state.spin); });
    bSpin.setAttribute('aria-pressed', 'false');
    var bFull = null;
    if (host.requestFullscreen) {
      bFull = toolButton('expand', 'Full screen', function () {
        if (document.fullscreenElement === host) document.exitFullscreen();
        else host.requestFullscreen().catch(function () {});
      });
    }

    VIEWS.forEach(function (v, i) {
      var c = el('button', 'gvv-chip', '<span>' + v.label + '</span>');
      c.type = 'button';
      c.addEventListener('click', function () { state.view = i; paintView(); apply(); frame(designs[state.design], 500); });
      chips.appendChild(c);
    });

    var binderName = host.getAttribute('data-binder') || 'GvpA binder';
    LEGEND.forEach(function (k) {
      legend.appendChild(el('li', 'gvv-key',
        '<span class="gvv-swatch" style="--gvv-swatch: var(' + k[0] + ')"></span>' +
        '<span class="gvv-key-label">' + (k[0] === '--lv-binder' ? binderName : k[1]) + '</span>'));
    });

    host.addEventListener('keydown', function (e) {
      if (e.target.closest && e.target.closest('input, textarea')) return;
      if (e.key === 'ArrowRight' && state.design < designs.length - 1) { go(state.design + 1); e.preventDefault(); }
      if (e.key === 'ArrowLeft' && state.design > 0) { go(state.design - 1); e.preventDefault(); }
    });

    var hintTimer;
    stage.addEventListener('wheel', function (e) {
      if (e.ctrlKey || e.metaKey || document.fullscreenElement === host) return;
      e.stopPropagation();
      hint.classList.add('is-on');
      clearTimeout(hintTimer);
      hintTimer = setTimeout(function () { hint.classList.remove('is-on'); }, 1400);
    }, { capture: true, passive: true });
    if (/Mac|iPhone|iPad/.test(navigator.platform)) hint.innerHTML = 'Use <kbd>⌘</kbd> + scroll to zoom';

    var designs = [], plugin, lib, V3, COLORS, built = {};

    /* the design buttons and card come from designs.json, which also arrives before Mol* */
    var jsonP = fetch(base + 'designs.json').then(function (r) {
      if (!r.ok) throw new Error('designs.json: HTTP ' + r.status);
      return r.json();
    });
    jsonP.then(function (j) {
      designs = j.designs;
      designs.forEach(function (d, i) {
        var li = el('li', 'lt-item');
        var b = el('button', 'lt-btn',
          '<span class="lt-ring">' + d.rank + '</span><span class="lt-id">' + d.id + '</span>');
        b.type = 'button';
        b.title = 'Design ' + d.rank + ' (' + d.id + ')';
        b.addEventListener('click', function () { go(i); });
        li.appendChild(b);
        picker.appendChild(li);
      });
      paintDesign(0);
      paintView();
    }).catch(function (err) { console.error('[latch-top10]', err); fail('The design list could not be loaded.'); });

    if (!webglOk()) return fail('This 3D figure needs WebGL, which this browser has switched off.');

    var started = false;
    function start() {
      if (started) return;
      started = true;
      boot().catch(function (err) {
        console.error('[latch-top10]', err);
        fail('The 3D model could not be loaded. Reload the page to try again.');
      });
    }
    if ('IntersectionObserver' in window) {
      var io = new IntersectionObserver(function (es) {
        if (es.some(function (x) { return x.isIntersecting; })) { io.disconnect(); start(); }
      }, { rootMargin: '600px 0px' });
      io.observe(host);
    } else start();

    function fail(msg) {
      veil.classList.add('is-error');
      veil.querySelector('.gvv-veil-text').textContent = msg;
      veil.hidden = false;
    }

    async function boot() {
      loadStyle(MOLSTAR + 'molstar.css');
      await Promise.all([loadScript(MOLSTAR + 'molstar.js'), jsonP]);
      lib = window.molstar.lib;
      V3 = lib.math.LinearAlgebra.Vec3;
      var viewer = await window.molstar.Viewer.create(canvas, {
        extensions: [],
        layoutIsExpanded: false,
        layoutShowControls: false,
        layoutShowRemoteState: false,
        layoutShowSequence: false,
        layoutShowLog: false,
        layoutShowLeftPanel: false,
        collapseLeftPanel: true,
        viewportShowExpand: false,
        viewportShowControls: false,
        viewportShowSettings: false,
        viewportShowSelectionMode: false,
        viewportShowAnimation: false,
        viewportShowTrajectoryControls: false,
        viewportShowReset: false,
        viewportShowScreenshotControls: false,
        viewportShowToggleFullscreen: false,
        viewportFocusBehavior: 'disabled',
        volumeStreamingDisabled: true
      });
      plugin = viewer.plugin;
      COLORS = readColours();
      styleCanvas();

      await go(state.design, true);
      veil.hidden = true;
      host.classList.add('is-ready');
      host.__pv = { plugin: plugin, lib: lib, built: built };       // what static/dev/dev-mode.js works on
      host.dispatchEvent(new CustomEvent('pv:ready'));

      if (window.ResizeObserver) new ResizeObserver(function () { viewer.handleResize(); }).observe(stage);
      document.addEventListener('fullscreenchange', function () {
        var full = document.fullscreenElement === host;
        if (bFull) {
          bFull.innerHTML = icon(full ? 'shrink' : 'expand');
          bFull.setAttribute('aria-label', full ? 'Leave full screen' : 'Full screen');
          bFull.title = bFull.getAttribute('aria-label');
        }
        setTimeout(function () { viewer.handleResize(); frame(designs[state.design], 0); }, 60);
      });

      /* build the neighbours while the reader looks at this design, so the next click is instant */
      var idle = window.requestIdleCallback || function (f) { setTimeout(f, 300); };
      idle(async function () {
        for (var i = 0; i < designs.length; i++) {
          try { await ensure(i); } catch (e) { console.warn('[latch-top10] prebuild', i, e); }
        }
      });
    }

    function readColours() {
      var c = {};
      ['--lv-binder', '--lv-latch', '--lv-tbd', '--lv-target', '--gv-outline', '--gv-highlight', '--gv-bg']
        .forEach(function (n) { c[n] = colour(host, n); });
      return c;
    }

    function styleCanvas() {
      var p = plugin.canvas3d.props;
      plugin.canvas3d.setProps({
        transparentBackground: true,
        checkeredTransparentBackground: false,
        camera: Object.assign({}, p.camera, {
          manualReset: true,
          helper: Object.assign({}, p.camera.helper, { axes: { name: 'off', params: {} } }) }),
        renderer: Object.assign({}, p.renderer, {
          backgroundColor: COLORS['--gv-bg'],
          highlightColor: COLORS['--gv-highlight'],
          selectColor: COLORS['--gv-highlight'],
          highlightStrength: 0.45,
          celSteps: 4,
          ambientIntensity: 0.55,
          exposure: 1.05
        }),
        postprocessing: Object.assign({}, p.postprocessing, {
          occlusion: { name: 'on', params: {
            samples: 32, multiScale: { name: 'off', params: {} }, radius: 5, bias: 0.8, blurKernelSize: 15,
            blurDepthBias: 0.5, resolutionScale: 1, color: COLORS['--gv-outline'], transparentThreshold: 0.4 } },
          outline: { name: 'on', params: {
            scale: 1, threshold: 0.33, color: COLORS['--gv-outline'], includeTransparent: true } }
        }),
        cameraFog: { name: 'on', params: { intensity: 12 } }
      });
    }

    /* ---- one design = closed + open structure, three representations (closed, open, open as a faint ghost) ---- */
    var LOOK = { celShaded: true, material: { metalness: 0, roughness: 1, bumpiness: 0 } };
    function cartoon(alpha) {
      // the files carry the part in the B-factor column; the theme runs high to low along its list, so PD-L1 (3) comes first
      return { type: 'cartoon', typeParams: Object.assign({ sizeFactor: 0.34, alpha: alpha }, LOOK),
        color: 'uncertainty', colorParams: { domain: [0, 3],
          list: { kind: 'interpolate', colors: [COLORS['--lv-target'], COLORS['--lv-tbd'], COLORS['--lv-latch'], COLORS['--lv-binder']] } } };
    }

    async function buildState(tag, which, withGhost) {
      var data = await plugin.builders.data.download(
        { url: base + tag + '_' + which + '.pdb', isBinary: false, label: tag + ' ' + which }, { state: { isGhost: true } });
      var traj = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      var struct = await plugin.builders.structure.createStructure(await plugin.builders.structure.createModel(traj));
      var add = function (props) {
        return plugin.builders.structure.representation.addRepresentation(struct, props, { initialState: { isHidden: true } })
          .then(function (s) { return s.ref; });
      };
      var out = {}; out[which] = await add(cartoon(1));
      if (withGhost) out.ghost = await add(cartoon(0.3));
      return out;
    }

    var pending = {};
    function ensure(i) {
      var tag = designs[i].tag;
      if (!pending[tag]) {
        pending[tag] = Promise.all([buildState(tag, 'closed', false), buildState(tag, 'open', true)]).then(function (r) {
          built[tag] = Object.assign({}, r[0], r[1]);
          if (state.design === i) show(i);
          return built[tag];
        });
      }
      return pending[tag];
    }

    function setVisible(ref, on) {
      var cell = ref && plugin.state.data.cells.get(ref);
      if (cell && !!cell.state.isHidden === on) plugin.state.data.updateCellState(ref, { isHidden: !on });
    }

    /** show design i in the chosen view, hide every other design and the rest of this one */
    function show(i) {
      var want = VIEWS[state.view].show;
      Object.keys(built).forEach(function (tag) {
        var mine = tag === designs[i].tag;
        ['closed', 'open', 'ghost'].forEach(function (part) {
          setVisible(built[tag][part], mine && want.indexOf(part) >= 0);
        });
      });
      plugin.canvas3d.requestDraw();
    }

    function apply() {
      if (!plugin) return;
      if (built[designs[state.design].tag]) show(state.design);
    }

    /* ---- camera ----
       The card sits under the stage, so the whole stage is for the model: its half-extent across the screen
       must fit the stage's width, its half-extent up the screen the height minus room for the chips and legend. */
    function frame(d, ms) {
      if (!plugin || !d) return;
      var f = d.frame, e = f[VIEWS[state.view].key === 'both' ? 'both' : VIEWS[state.view].key];
      var aspect = stage.clientWidth / Math.max(stage.clientHeight, 1);
      var r = Math.max(e.h / (0.94 * aspect), e.v / 0.8, e.z / 0.9);
      var cam = plugin.canvas3d.camera;
      var snap = cam.getInvariantFocus(V3.create.apply(null, e.target), r, V3.create.apply(null, f.up),
                                       V3.create.apply(null, f.dir));
      snap.radiusMax = 700;
      cam.setState(snap, stillMotion() ? 0 : ms);
      plugin.canvas3d.requestDraw();
    }

    function setSpin(on) {
      state.spin = on;
      bSpin.setAttribute('aria-pressed', on ? 'true' : 'false');
      if (!plugin) return;
      var tb = plugin.canvas3d.props.trackball;
      plugin.canvas3d.setProps({ trackball: Object.assign({}, tb, {
        animate: on ? { name: 'spin', params: { speed: 0.25, axis: [0, -1, 0] } } : { name: 'off', params: {} } }) });
    }

    /* ---- the buttons and card (they work before Mol* has arrived) ---- */
    function paintView() {
      Array.prototype.forEach.call(chips.children, function (c, i) {
        c.setAttribute('aria-pressed', i === state.view ? 'true' : 'false');
      });
    }

    function paintDesign(i) {
      var d = designs[i];
      Array.prototype.forEach.call(picker.querySelectorAll('.lt-btn'), function (b, j) {
        if (j === i) b.setAttribute('aria-current', 'true'); else b.removeAttribute('aria-current');
      });
      var end = d.terminus === 'c' ? 'C' : 'N';
      card.innerHTML = '<span class="gvv-eyebrow">Design ' + d.rank + ' of ' + designs.length + ' · ' + d.id + '</span>' +
        '<strong class="gvv-title">Hinge ' + hingeText(d.moves) + '</strong>' +
        '<span class="gvv-text">A ' + d.n_latch + '-residue latch on the ' + end + ' terminus of the ' + d.n_binder +
        '-residue binder, ' + (d.k ? 'joined to the target-binding domain by a ' + d.k + '-residue GS linker (not drawn)'
                                   : 'fused directly to the target-binding domain') +
        '. Closed, PD-L1 on the domain runs into the binder; open, the swung block keeps ' + d.clearance.toFixed(1) +
        ' Å clearance (search score ' + d.quality.toFixed(2) + '). Drag to turn it.</span>';
    }

    async function go(i, first) {
      if (!first && i === state.design) return;
      state.design = i;
      paintDesign(i);
      if (!plugin) return;
      host.classList.add('is-busy');
      try {
        await ensure(i);
        if (state.design !== i) return;        // clicked on past it while it built
        show(i);
        frame(designs[i], first ? 0 : 700);
      } finally {
        if (state.design === i) host.classList.remove('is-busy');
      }
    }
  }

  function init() {
    Array.prototype.forEach.call(document.querySelectorAll('[data-latch-top10]'), mount);
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
