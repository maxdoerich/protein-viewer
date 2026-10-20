/* GvpC-repeat viewer: two steps. 1: four GvpA and one GvpC repeat, coloured by sequence. 2: the camera flies onto
 * GvpC, GvpA turns see-through, and the GvpC residues take the colour of their job (anchor, other interface, rest).
 *
 * Mounted onto every [data-gvpc-viewer] element. Mol* 5.13 is vendored in ../molstar/; the layout, text box, legend and
 * buttons are the gas-vesicle viewer's (gv-viewer.css, loaded first), the colours are --gc-* in gvpc-viewer.css.
 *
 * Data (built by proteins/gvpc_animation/tools/prep_viewer.py), from this script's folder or data-base="…/":
 *   gvpc-assembly.cif  assemblies "gvpa" (chains A-D) and "gvpc" (chain F); B-factor of chain F = the residue's job
 *   gvpc-viewer.json   the two camera frames
 */
(function () {
  'use strict';

  var HERE = (document.currentScript && document.currentScript.src || '').replace(/[^/]*$/, '');
  var MOLSTAR = HERE + '../molstar/';

  var ICON = {
    // four GvpA side by side, GvpC lying across them
    context: '<rect x="3" y="12" width="3.6" height="8" rx="1.8" /><rect x="8.2" y="12" width="3.6" height="8" rx="1.8" />' +
             '<rect x="13.4" y="12" width="3.6" height="8" rx="1.8" /><rect x="18.6" y="12" width="2.4" height="8" rx="1.2" />' +
             '<path d="M2.5 7.5h19" stroke-width="2.4" />',
    gvpc: '<path d="M2.5 12c1.5-5 3.5-5 4.5 0s3 5 4.5 0 3-5 4.5 0 3 5 4.5 0" stroke-width="2.4" />',
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

  /* `show` names components: the three cartoons made in boot(). */
  var GVPA_SEE_THROUGH = 0.3;
  var STEPS = [
    { key: 'context', label: 'Context', icon: 'context', frame: 'all',
      title: 'GvpC on four GvpA',
      text: 'One repeat of GvpC (33 residues) lies on four neighbouring GvpA of the gas vesicle wall. Every sequence has ' +
            'its own colour: GvpA (four copies) and GvpC.',
      show: ['gvpa', 'seq'], alpha: { gvpa: 1 },
      legend: [['--gc-gvpa', 'GvpA'], ['--gc-gvpc', 'GvpC']] },
    { key: 'function', label: 'GvpC residues', icon: 'gvpc', frame: 'gvpc',
      title: 'What each GvpC residue does',
      text: 'The camera moves onto GvpC and GvpA turns see-through. The anchor residues (4, 11, 12, 15, 19, 22 of the ' +
            'repeat) are dark blue, the other interface residues (1, 5, 16, 33) light blue, all remaining residues orange.',
      show: ['gvpa', 'func'], alpha: { gvpa: GVPA_SEE_THROUGH },
      legend: [['--gc-anchor', 'anchor residues'], ['--gc-site', 'other interface residues'], ['--gc-rest', 'other residues']] }
  ];

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
  /** a custom property (resolved to its token) as 0xRRGGBB */
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
  function ease(k) { return k < .5 ? 4 * k * k * k : 1 - Math.pow(-2 * k + 2, 3) / 2; }   // in-out cubic
  function nextFrame() { return new Promise(function (r) { requestAnimationFrame(function () { r(); }); }); }

  function mount(host) {
    var base = host.getAttribute('data-base') || HERE;
    var state = { step: -1, token: 0, busy: Promise.resolve() };

    /* ---- markup ---- */
    host.innerHTML = '';
    host.classList.add('gvv', 'gc');
    var stage = el('div', 'gvv-stage');
    var canvas = el('div', 'gvv-canvas');
    var card = el('div', 'gvv-card');
    card.setAttribute('aria-live', 'polite');
    var tools = el('div', 'gvv-tools');
    var legend = el('ul', 'gvv-legend');
    var hint = el('div', 'gvv-hint', 'Use <kbd>Ctrl</kbd> + scroll to zoom');
    var veil = el('div', 'gvv-veil',
      '<span class="gvv-pulse" aria-hidden="true"></span><span class="gvv-veil-text">Preparing the 3D model…</span>');
    stage.append(canvas, tools, legend, hint, veil);

    var steps = el('ol', 'gvv-steps');
    steps.setAttribute('aria-label', 'Steps');
    STEPS.forEach(function (s, i) {
      var li = el('li', 'gvv-step');
      var b = el('button', 'gvv-step-btn',
        '<span class="gvv-ring">' + icon(s.icon) + '</span>' +
        '<span class="gvv-step-label"><span class="gvv-step-n">' + (i + 1) + '</span>' + s.label + '</span>');
      b.type = 'button';
      b.title = s.title;
      b.addEventListener('click', function () { go(i); });
      li.appendChild(b);
      steps.appendChild(li);
    });
    host.append(stage, card, steps);

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
    toolButton('reset', 'Reset view', function () { frame(STEPS[Math.max(state.step, 0)], 700); });
    var bSpin = toolButton('spin', 'Rotate', function () { setSpin(!state.spin); });
    bSpin.setAttribute('aria-pressed', 'false');
    var bFull = null;
    if (host.requestFullscreen) {
      bFull = toolButton('expand', 'Full screen', function () {
        if (document.fullscreenElement === host) document.exitFullscreen();
        else host.requestFullscreen().catch(function () {});
      });
    }

    host.addEventListener('keydown', function (e) {          // arrows move between steps while focus is in the figure
      if (e.target.closest && e.target.closest('input, textarea')) return;
      if (e.key === 'ArrowRight' && state.step < STEPS.length - 1) { go(state.step + 1); e.preventDefault(); }
      if (e.key === 'ArrowLeft' && state.step > 0) { go(state.step - 1); e.preventDefault(); }
    });

    var hintTimer;                                           // the page keeps scrolling over the figure; Ctrl/⌘ + wheel zooms
    stage.addEventListener('wheel', function (e) {
      if (e.ctrlKey || e.metaKey || document.fullscreenElement === host) return;
      e.stopPropagation();
      hint.classList.add('is-on');
      clearTimeout(hintTimer);
      hintTimer = setTimeout(function () { hint.classList.remove('is-on'); }, 1400);
    }, { capture: true, passive: true });
    if (/Mac|iPhone|iPad/.test(navigator.platform)) hint.innerHTML = 'Use <kbd>⌘</kbd> + scroll to zoom';

    paintStep(0);                                            // the text reads correctly before Mol* arrives

    function fail(msg) {
      veil.classList.add('is-error');
      veil.querySelector('.gvv-veil-text').textContent = msg;
      veil.hidden = false;
    }
    if (!webglOk()) return fail('This 3D figure needs WebGL, which this browser has switched off.');

    var started = false;                                     // Mol* is fetched once the figure nears the viewport
    function start() {
      if (started) return;
      started = true;
      boot().catch(function (err) {
        console.error('[gvpc-viewer]', err);
        fail('The 3D model could not be loaded. Reload the page to try again.');
      });
    }
    if ('IntersectionObserver' in window) {
      var io = new IntersectionObserver(function (es) {
        if (es.some(function (x) { return x.isIntersecting; })) { io.disconnect(); start(); }
      }, { rootMargin: '600px 0px' });
      io.observe(host);
    } else start();

    var viewer, plugin, lib, manifest, COLORS, built = {};

    async function boot() {
      loadStyle(MOLSTAR + 'molstar.css');
      var manifestP = fetch(base + 'gvpc-viewer.json').then(function (r) {
        if (!r.ok) throw new Error('gvpc-viewer.json: HTTP ' + r.status);
        return r.json();
      });
      await loadScript(MOLSTAR + 'molstar.js');
      lib = window.molstar.lib;
      viewer = await window.molstar.Viewer.create(canvas, {
        extensions: [], layoutIsExpanded: false, layoutShowControls: false, layoutShowRemoteState: false,
        layoutShowSequence: false, layoutShowLog: false, layoutShowLeftPanel: false, collapseLeftPanel: true,
        viewportShowExpand: false, viewportShowControls: false, viewportShowSettings: false,
        viewportShowSelectionMode: false, viewportShowAnimation: false, viewportShowTrajectoryControls: false,
        viewportShowReset: false, viewportShowScreenshotControls: false, viewportShowToggleFullscreen: false,
        viewportFocusBehavior: 'disabled', volumeStreamingDisabled: true
      });
      plugin = viewer.plugin;
      manifest = await manifestP;
      COLORS = {};
      ['--gc-gvpa', '--gc-gvpc', '--gc-anchor', '--gc-site', '--gc-rest', '--gv-outline', '--gv-highlight', '--gv-bg']
        .forEach(function (n) { COLORS[n] = colour(host, n); });
      styleCanvas();

      var data = await plugin.builders.data.download(
        { url: base + 'gvpc-assembly.cif', isBinary: false, label: 'GvpC repeat' }, { state: { isGhost: true } });
      var traj = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
      var model = await plugin.builders.structure.createModel(traj);
      var gvpa = await plugin.builders.structure.createStructure(model, { name: 'assembly', params: { id: 'gvpa' } });
      var gvpc = await plugin.builders.structure.createStructure(model, { name: 'assembly', params: { id: 'gvpc' } });
      built.gvpa = await rep(gvpa, uniform(COLORS['--gc-gvpa']));
      built.seq = await rep(gvpc, uniform(COLORS['--gc-gvpc']));
      // chain F's B-factor is the residue's job: 2 anchor, 1 other interface, 0 the rest. Mol*'s uncertainty scale runs
      // backwards (low value = last colour), hence anchor first.
      built.func = await rep(gvpc, { color: 'uncertainty', colorParams: { domain: [0, 2],
        list: { kind: 'interpolate', colors: [COLORS['--gc-anchor'], COLORS['--gc-site'], COLORS['--gc-rest']] } } });

      var first = Math.max(state.step, 0);
      state.step = -1;
      await go(first, true);
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
        setTimeout(function () { viewer.handleResize(); frame(STEPS[state.step], 0); }, 60);
      });
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
            scale: 1, threshold: 0.33, color: COLORS['--gv-outline'], includeTransparent: true } }   // the outline stays while GvpA fades
        }),
        cameraFog: { name: 'on', params: { intensity: 12 } }
      });
    }

    /* ---- components: three cartoons, built once, then only shown, hidden or faded ---- */
    var LOOK = { celShaded: true, material: { metalness: 0, roughness: 1, bumpiness: 0 } };
    function uniform(c) { return { color: 'uniform', colorParams: { value: c } }; }
    async function rep(structure, colourProps) {
      var sel = await plugin.builders.structure.representation.addRepresentation(structure,
        Object.assign({ type: 'cartoon', typeParams: Object.assign({ sizeFactor: 0.32 }, LOOK) }, colourProps),
        { initialState: { isHidden: true } });
      return sel.ref;
    }
    function setVisible(name, on) {
      var cell = plugin.state.data.cells.get(built[name]);
      if (cell && !!cell.state.isHidden === on) plugin.state.data.updateCellState(built[name], { isHidden: !on });
    }
    /* the opacity is the representation's alphaFactor: Mol* draws with alpha x alphaFactor, no state rebuild per frame.
       Never exactly 0: Mol* leaves a fully transparent part out of the visible bounds, which can clamp the camera. */
    var cur = {};
    function setAlpha(name, a) {
      a = Math.max(a, 0.01);
      var cell = plugin.state.data.cells.get(built[name]);
      if (!cell || !cell.obj || (cur[name] == null ? 1 : cur[name]) === a) return;
      cur[name] = a;
      cell.obj.data.repr.setState({ alphaFactor: a });
      plugin.canvas3d.update(cell.obj.data.repr);
    }
    function alphaOf(name) { return cur[name] == null ? 1 : cur[name]; }

    /* ---- camera ---- */
    function frame(step, ms) {
      if (!plugin || !manifest) return;
      var Vec3 = lib.math.LinearAlgebra.Vec3, cam = plugin.canvas3d.camera;
      var f = manifest.frames[step.frame];
      var dir = Vec3.create.apply(null, f.dir), up = Vec3.create.apply(null, f.up), target = Vec3.create.apply(null, f.target);
      var right = Vec3.cross(Vec3(), dir, up);
      var aspect = cam.viewport.width / Math.max(cam.viewport.height, 1);
      var radius = Math.max(f.halfHeight, f.halfWidth / aspect);   // whichever of width and height limits
      if (getComputedStyle(card).position === 'absolute') {   // from 768 up the text box sits over the scene's top-left
        radius *= 1.15;
        Vec3.scaleAndAdd(target, target, right, -0.2 * radius);
        Vec3.scaleAndAdd(target, target, up, 0.1 * radius);
      }
      var snap = cam.getInvariantFocus(target, radius, up, dir);
      snap.radiusMax = manifest.frames.all.radius * 4;
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

    /* ---- the step UI (works before and after the 3D part loads) ---- */
    function paintStep(i) {
      var s = STEPS[i];
      Array.prototype.forEach.call(steps.querySelectorAll('.gvv-step-btn'), function (b, j) {
        if (j === i) b.setAttribute('aria-current', 'step'); else b.removeAttribute('aria-current');
        b.parentNode.classList.toggle('is-done', j < i);
      });
      steps.style.setProperty('--gvv-progress', String(i / (STEPS.length - 1)));
      card.innerHTML = '<span class="gvv-eyebrow">Step ' + (i + 1) + ' of ' + STEPS.length + '</span>' +
        '<strong class="gvv-title">' + s.title + '</strong><span class="gvv-text">' + s.text + '</span>';
      legend.innerHTML = '';
      s.legend.forEach(function (k) {
        legend.appendChild(el('li', 'gvv-key',
          '<span class="gvv-swatch" style="--gvv-swatch: var(' + k[0] + ')"></span><span class="gvv-key-label">' + k[1] + '</span>'));
      });
    }

    /* ---- step 1 -> 2: the camera flies onto GvpC (1.2 s) while GvpA fades to see-through, then the job colours fade in
       over the sequence colour. Going back just switches. ---- */
    var FLY = 1200, FADE_A = 900, FADE_C = 800;
    async function go(i, first) {
      if (i === state.step) return state.busy;
      var from = state.step, s = STEPS[i], token = ++state.token;
      state.step = i;
      paintStep(i);
      if (!plugin) return state.busy;
      host.classList.add('is-busy');
      state.busy = (async function () {
        var move = !first && i > from && !stillMotion();
        Object.keys(built).forEach(function (n) {
          var on = s.show.indexOf(n) >= 0;
          if (on && !(move && n === 'func')) setVisible(n, true);
          else if (!on && !(move && n === 'seq')) setVisible(n, false);    // seq stays under the fading colours
          if (!move) setAlpha(n, on && s.alpha[n] != null ? s.alpha[n] : 1);
        });
        frame(s, first ? 0 : 1100);
        if (!move) return;
        setAlpha('func', 0);
        setVisible('func', true);
        var a0 = alphaOf('gvpa'), t0 = performance.now(), k;
        do {
          if (token !== state.token) return;
          var t = performance.now() - t0;
          setAlpha('gvpa', a0 + (s.alpha.gvpa - a0) * ease(Math.min(t / FADE_A, 1)));
          setAlpha('func', ease(Math.min(Math.max(t - FADE_A, 0) / FADE_C, 1)));
          plugin.canvas3d.requestDraw();
          k = t < FADE_A + FADE_C;
          if (k) await nextFrame();
        } while (k);
        setVisible('seq', false);
      })().finally(function () {
        if (state.step === i) host.classList.remove('is-busy');
      });
      return state.busy;
    }
  }

  function init() {
    Array.prototype.forEach.call(document.querySelectorAll('[data-gvpc-viewer]'), function (h) {
      if (!h.__gcv) { h.__gcv = true; mount(h); }
    });
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
