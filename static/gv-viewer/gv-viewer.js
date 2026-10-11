/* Gas-vesicle assembly viewer: the interactive 3D figure on the Model page.
 *
 * Mounted by the ::custom block "gv-assembly" in content/Project/Model.md onto every
 * [data-gv-viewer] element. Renders with Mol* 5.13 (MIT, vendored in ../molstar/, see
 * NOTICE.md), which is only fetched once the figure scrolls near the viewport.
 *
 * Data (built by proteins/gv_animation/tools/prep_viewer.py):
 *   gv-assembly.cif  GvpA + one GvpC repeat, every copy an mmCIF assembly operator
 *   emd-14238.mrc    the EMD-14238 map, cropped to the wall
 *   gv-viewer.json   isovalue + camera frames
 * They load from this script's folder, or from data-base="…/" on the mount element
 * (the iGEM CDN folder, once uploaded).
 *
 * Colours are not chosen here: each one is a --gv-* property in gv-viewer.css that
 * points at a site token, read back at start-up.
 */
(function () {
  'use strict';

  var HERE = (document.currentScript && document.currentScript.src || '').replace(/[^/]*$/, '');
  var MOLSTAR = HERE + '../molstar/';

  /* ---------------------------------------------------------------- icons -- */
  /* 24 x 24, drawn with currentColor strokes so they take the button's colour. */
  var ICON = {
    // one GvpA: a helix, then the beta-hairpin
    monomer: '<path d="M4.5 4c3 .6 3 2.4 0 3s-3 2.4 0 3 3 2.4 0 3 -3 2.4 0 3" />' +
             '<path d="M11 20V8.5a3 3 0 0 1 6 0V16" /><path d="M14.5 13.5 17 16.5l2.5-3" />',
    // a patch of wall: staggered rows of subunits = ribs
    rib: '<rect x="2.5" y="4" width="5.5" height="4" rx="2" /><rect x="9.25" y="4" width="5.5" height="4" rx="2" />' +
         '<rect x="16" y="4" width="5.5" height="4" rx="2" /><rect x="5.9" y="10" width="5.5" height="4" rx="2" />' +
         '<rect x="12.6" y="10" width="5.5" height="4" rx="2" /><rect x="2.5" y="16" width="5.5" height="4" rx="2" />' +
         '<rect x="9.25" y="16" width="5.5" height="4" rx="2" /><rect x="16" y="16" width="5.5" height="4" rx="2" />',
    // the whole vesicle: cylinder with conical tips, ribs across it
    vesicle: '<path d="M2 12 6.5 6.5h11L22 12l-4.5 5.5h-11Z" /><path d="M9 6.5v11M12 6.5v11M15 6.5v11" />',
    // vesicle wrapped by GvpC
    gvpc: '<path d="M2 12 6.5 6.5h11L22 12l-4.5 5.5h-11Z" />' +
          '<path d="M5 14.5 9.5 5M9.5 19l5-14M14.5 19 19 9.5" stroke-width="2.4" />',
    reset: '<path d="M4 12a8 8 0 1 0 2.4-5.7" /><path d="M4 4v4.5h4.5" />',
    spin: '<ellipse cx="12" cy="12" rx="9" ry="4" /><path d="M12 3v18" stroke-dasharray="2 2.5" />' +
          '<path d="m16.5 14.6 2.6 1-.8 2.7" />',
    expand: '<path d="M4 9V4h5M15 4h5v5M20 15v5h-5M9 20H4v-5" />',
    shrink: '<path d="M9 4v5H4M20 9h-5V4M15 20v-5h5M4 15h5v5" />',
    density: '<path d="M4 7.5 12 3l8 4.5v9L12 21l-8-4.5Z" /><path d="M4 7.5 12 12l8-4.5M12 12v9" stroke-dasharray="2 2" />',
    hydro: '<path d="M12 3.5c3.5 4.2 6 7.4 6 10.5a6 6 0 0 1-12 0c0-3.1 2.5-6.3 6-10.5Z" /><path d="M6 14h12" />' +
           '<path d="M6.6 16.5A6 6 0 0 0 17.4 16.5Z" fill="currentColor" stroke="none" />'
  };

  function icon(name) {
    return '<svg class="gvv-icon" viewBox="0 0 24 24" aria-hidden="true" focusable="false" fill="none" ' +
           'stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round">' +
           ICON[name] + '</svg>';
  }

  /* ---------------------------------------------------------------- steps -- */
  /* What each step shows. Text follows the team's captions in
     proteins/gv_animation/captions.txt. `show` names components (see COMPONENTS). */
  var ALL_STEPS = [
    { key: 'monomer', label: 'GvpA', icon: 'monomer', frame: 'monomer',
      title: 'One GvpA',
      text: 'The wall of a gas vesicle is one protein layer thick, built from thousands of copies of ' +
            'this one small protein. GvpA folds into an α–β–β–α motif: two helices and a β-hairpin.',
      show: ['heroSS'],
      toggle: { key: 'hydro', label: 'Hydropathy', icon: 'hydro', on: false, show: ['heroHydro'], hide: ['heroSS'],
                text: 'Coloured by hydropathy (Kyte–Doolittle): the side facing the gas is strongly hydrophobic, ' +
                      'the outer surface hydrophilic. Gas passes through, water is kept out.' },
      legend: [['--gv-helix', 'α-helix'], ['--gv-strand', 'β-strand'], ['--gv-loop', 'loop']],
      legendAlt: 'hydro' },
    { key: 'rib', label: 'Ribs', icon: 'rib', frame: 'rib',
      title: 'GvpA packs into ribs',
      text: 'Copies of GvpA pack side by side into ribs, which gives the wall its corrugated shape. ' +
            'The pale skin is the cryo-EM density of the wall (EMD-14238) around the fitted model (PDB 7R1C).',
      show: ['heroFlat', '@patch'],
      toggle: { key: 'density', label: 'Density', icon: 'density', on: true, show: ['density'], hide: [] },
      legend: [['--gv-hero', 'the GvpA from step 1'], ['--gv-rib-a', 'its rib'], ['--gv-rib-b', 'the next ribs'],
               ['--gv-density', 'EMD-14238 density', 'skin']] },
    { key: 'vesicle', label: 'Vesicle', icon: 'vesicle', frame: 'vesicle',
      title: '1730 GvpA build one vesicle',
      text: 'A vesicle starts as a small bicone and grows at its central seam, where new ribs are added and push ' +
            'the two halves apart. Each half is one continuous, low-pitch helix of GvpA; the halves have opposite polarity.',
      show: ['heroBlob', '@shell'],
      legend: [['--gv-half-a', 'helical half 1'], ['--gv-half-b', 'helical half 2'], ['--gv-hero', 'the GvpA from step 1']] },
    { key: 'gvpc', label: '+ GvpC', icon: 'gvpc', frame: 'vesicle',
      title: 'GvpC reinforces the shell',
      text: 'GvpC binds to the outside along the ribs and forms a helical cage that reinforces the shell. ' +
            'Shown is a predicted full-length GvpC (AlphaFold Server), placed where GvpC sits in the 8GBS model, ' +
            'on the cylinder only: a straight helix does not follow the cones. Its coverage is illustrative.',
      show: ['heroPlain', '@shell', 'gvpc'],
      legend: [['--gv-gvpc', 'GvpC'], ['--gv-half-a', 'GvpA, half 1'], ['--gv-half-b', 'GvpA, half 2']] }
  ];

  /* The 4-GvpA preview ([data-gv-viewer][data-preview="hydro4"]): four neighbours of one rib in
     hydropathy colours, seen edge-on with the vesicle's outside up and its gas-filled inside down. */
  var SIDE_OFFSET = 32;      // Å from the wall's centre to where the Outside / Inside labels are anchored
  var HYDRO4 = { key: 'hydro4', label: 'GvpA', icon: 'monomer', frame: 'hydro4',
    title: 'Four GvpA, outside and inside',
    text: 'Four neighbouring GvpA of one rib, coloured by hydropathy (Kyte–Doolittle). The wall is one layer thick: ' +
          'its outer face is hydrophilic and meets the water, its inner face is strongly hydrophobic and meets the gas. ' +
          'Gas passes through, water is kept out.',
    show: ['@hydro4'], legend: [], legendAlt: 'hydro' };

  /* -------------------------------------------------------------- helpers -- */
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
  /** a --gv-* custom property (resolved to its token) as 0xRRGGBB */
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

  /* ================================================================ viewer == */
  function mount(host) {
    var base = host.getAttribute('data-base') || HERE;
    var preview = host.getAttribute('data-preview') === 'hydro4';
    var STEPS = preview ? [HYDRO4] : ALL_STEPS;
    var state = { step: -1, toggles: {}, busy: Promise.resolve(), token: 0 };
    STEPS.forEach(function (s) { if (s.toggle) state.toggles[s.toggle.key] = s.toggle.on; });

    /* ---- markup ---- */
    host.innerHTML = '';
    host.classList.add('gvv');
    var stage = el('div', 'gvv-stage');
    var canvas = el('div', 'gvv-canvas');
    var card = el('div', 'gvv-card'), cardEl = card;
    card.setAttribute('aria-live', 'polite');
    var tools = el('div', 'gvv-tools');
    var chips = el('div', 'gvv-chips');
    var legend = el('ul', 'gvv-legend');
    var hint = el('div', 'gvv-hint', 'Use <kbd>Ctrl</kbd> + scroll to zoom');
    var veil = el('div', 'gvv-veil',
      '<span class="gvv-pulse" aria-hidden="true"></span><span class="gvv-veil-text">Preparing the 3D model…</span>');
    stage.append(canvas, tools, chips, legend, hint, veil);
    var sides = [];
    if (preview) {
      host.classList.add('gvv--preview');      // which way is out: placed in 3D by the camera listener below
      sides = [el('div', 'gvv-side gvv-side--out', '<strong>Outside</strong><span>water · hydrophilic</span>'),
               el('div', 'gvv-side gvv-side--in', '<strong>Inside</strong><span>gas · hydrophobic</span>')];
      stage.append.apply(stage, sides);
    }

    var steps = el('ol', 'gvv-steps');
    steps.setAttribute('aria-label', 'Assembly steps');
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
    host.append(stage, steps);
    if (!preview) stage.after(card);       // the preview has no text box
    steps.hidden = STEPS.length < 2;

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
    var bReset = toolButton('reset', 'Reset view', function () { frame(STEPS[Math.max(state.step, 0)], 700); });
    var bSpin = toolButton('spin', 'Rotate', function () { setSpin(!state.spin); });
    bSpin.setAttribute('aria-pressed', 'false');
    var bFull = null;
    if (host.requestFullscreen) {
      bFull = toolButton('expand', 'Full screen', function () {
        if (document.fullscreenElement === host) document.exitFullscreen();
        else host.requestFullscreen().catch(function () {});
      });
    }

    /* keyboard: arrows move between steps while focus is in the figure */
    host.addEventListener('keydown', function (e) {
      if (e.target.closest && e.target.closest('input, textarea')) return;
      if (e.key === 'ArrowRight' && state.step < STEPS.length - 1) { go(state.step + 1); e.preventDefault(); }
      if (e.key === 'ArrowLeft' && state.step > 0) { go(state.step - 1); e.preventDefault(); }
    });

    /* the page keeps scrolling over the figure; Ctrl/⌘ + wheel (and trackpad pinch) zooms */
    var hintTimer;
    stage.addEventListener('wheel', function (e) {
      if (e.ctrlKey || e.metaKey || document.fullscreenElement === host) return;
      e.stopPropagation();
      hint.classList.add('is-on');
      clearTimeout(hintTimer);
      hintTimer = setTimeout(function () { hint.classList.remove('is-on'); }, 1400);
    }, { capture: true, passive: true });
    if (/Mac|iPhone|iPad/.test(navigator.platform)) hint.innerHTML = 'Use <kbd>⌘</kbd> + scroll to zoom';

    /* step UI without the 3D part: the page reads correctly before Mol* arrives */
    paintStep(0);

    if (!webglOk()) return fail('This 3D figure needs WebGL, which this browser has switched off.');

    /* ---- load Mol* when the figure comes near ---- */
    var started = false;
    function start() {
      if (started) return;
      started = true;
      boot().catch(function (err) {
        console.error('[gv-viewer]', err);
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

    var viewer, plugin, lib, manifest, shell, model, built = {}, COLORS;

    async function boot() {
      loadStyle(MOLSTAR + 'molstar.css');
      function json(name) {
        return fetch(base + name).then(function (r) {
          if (!r.ok) throw new Error(name + ': HTTP ' + r.status);
          return r.json();
        });
      }
      var manifestP = json('gv-viewer.json'), shellP = json('gv-shell.json');
      await loadScript(MOLSTAR + 'molstar.js');
      lib = window.molstar.lib;
      viewer = await window.molstar.Viewer.create(canvas, {
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
      manifest = await manifestP;
      shell = await shellP;
      COLORS = readColours();
      styleCanvas();

      var data = await plugin.builders.data.download(
        { url: base + 'gv-assembly.cif', isBinary: false, label: 'Gas vesicle' }, { state: { isGhost: true } });
      var traj = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
      model = await plugin.builders.structure.createModel(traj);

      var first = Math.max(state.step, 0);
      state.step = -1;
      await go(first, true);
      veil.hidden = true;
      host.classList.add('is-ready');
      host.__pv = { plugin: plugin, lib: lib, built: built };       // what static/dev/dev-mode.js works on
      host.dispatchEvent(new CustomEvent('pv:ready'));

      if (preview) {         // each label rides a 3D point just outside / inside the wall, so it follows every turn
        var cam = plugin.canvas3d.camera, h4 = manifest.hydro4, Vec3 = lib.math.LinearAlgebra.Vec3, pt = Vec3();
        var place = function () {
          var w = stage.clientWidth, h = stage.clientHeight;       // CSS pixels: cam.viewport is in device pixels
          sides.forEach(function (e, i) {
            var k = i ? -SIDE_OFFSET : SIDE_OFFSET;
            Vec3.set(pt, h4.target[0] + h4.up[0] * k, h4.target[1] + h4.up[1] * k, h4.target[2] + h4.up[2] * k);
            Vec3.transformMat4(pt, pt, cam.projectionView);
            var x = (pt[0] * .5 + .5) * w, y = (1 - (pt[1] * .5 + .5)) * h, m = 24;
            e.style.left = Math.min(Math.max(x, m), w - m) + 'px';
            e.style.top = Math.min(Math.max(y, m), h - m) + 'px';
          });
        };
        plugin.canvas3d.didDraw.subscribe(place);    // after every draw, when projectionView is current
        place();
      }
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

      /* build what the later steps need while the reader looks at step 1, so a click is instant */
      var idle = window.requestIdleCallback || function (f) { setTimeout(f, 300); };
      if (!preview) idle(async function () {
        var later = ['heroHydro', 'heroFlat'].concat(group('@patch'), ['density', 'heroBlob'], group('@shell'), group('@pool'),
                                                  ['heroPlain', 'gvpc']);
        for (var name of later) {
          try { await ensure(name); } catch (e) { console.warn('[gv-viewer] prebuild', name, e); }
        }
      });
    }

    function readColours() {
      var c = {};
      ['--gv-helix', '--gv-strand', '--gv-loop', '--gv-hero', '--gv-rib-a', '--gv-rib-b', '--gv-density',
       '--gv-half-a', '--gv-half-b', '--gv-gvpc', '--gv-phobic', '--gv-neutral', '--gv-philic',
       '--gv-outline', '--gv-highlight', '--gv-bg'].forEach(function (n) { c[n] = colour(host, n); });
      return c;
    }

    function styleCanvas() {
      var p = plugin.canvas3d.props;
      plugin.canvas3d.setProps({
        transparentBackground: true,
        checkeredTransparentBackground: false,
        camera: Object.assign({}, p.camera, {
          // the viewer places the camera itself: Mol* would otherwise re-aim it whenever
          // nothing visible is in front of it, e.g. as the vesicle starts growing
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
            scale: 1, threshold: 0.33, color: COLORS['--gv-outline'], includeTransparent: true } }   // true: the outline stays while a part fades
        }),
        cameraFog: { name: 'on', params: { intensity: 12 } }
      });
    }

    /* ---- components: built once, on demand, then only shown or hidden ---- */
    var LOOK = { celShaded: true, material: { metalness: 0, roughness: 1, bumpiness: 0 } };
    function uniform(c) { return { color: 'uniform', colorParams: { value: c } }; }
    function cartoon(extra) { return Object.assign({ type: 'cartoon', typeParams: Object.assign({ sizeFactor: 0.28 }, LOOK) }, extra); }
    function blob(extra) {
      return Object.assign({ type: 'gaussian-surface',
        typeParams: Object.assign({ quality: 'custom', resolution: 1.6, radiusOffset: 0.8, smoothness: 1.6 }, LOOK) }, extra);
    }

    // B-factor column = Kyte–Doolittle hydropathy per residue (prep_data.py)
    function hydro() {
      return { type: 'molecular-surface', typeParams: Object.assign({ quality: 'high', probeRadius: 1.4 }, LOOK),
        color: 'uncertainty', colorParams: { domain: [-4.5, 4.5],
          // Mol*'s uncertainty theme runs its scale backwards (low value = last colour), hence phobic first
          list: { kind: 'interpolate', colors: [COLORS['--gv-phobic'], COLORS['--gv-neutral'], COLORS['--gv-philic']] } } };
    }

    var COMPONENTS = {
      heroSS: function () {
        var ss = { alphaHelix: COLORS['--gv-helix'], threeTenHelix: COLORS['--gv-helix'], piHelix: COLORS['--gv-helix'],
                   betaStrand: COLORS['--gv-strand'], betaTurn: COLORS['--gv-loop'], coil: COLORS['--gv-loop'],
                   bend: COLORS['--gv-loop'], turn: COLORS['--gv-loop'], dna: COLORS['--gv-loop'], rna: COLORS['--gv-loop'],
                   carbohydrate: COLORS['--gv-loop'] };
        return rep('hero', { type: 'cartoon', typeParams: Object.assign({ sizeFactor: 0.32 }, LOOK),
          color: 'secondary-structure', colorParams: { saturation: 0, lightness: 0, colors: { name: 'custom', params: ss } } });
      },
      heroHydro: function () { return rep('hero', hydro()); },
      heroFlat: function () {               // a little fatter than its neighbours, so it is found at once
        return rep('hero', Object.assign({ type: 'cartoon', typeParams: Object.assign({ sizeFactor: 0.42 }, LOOK) },
          uniform(COLORS['--gv-hero'])));
      },
      heroBlob: function () { return rep('hero', blob(uniform(COLORS['--gv-hero']))); },
      // step 4: the step-1 GvpA back in its half's colour, so orange means GvpC alone
      heroPlain: function () { return rep('hero', blob(uniform(COLORS[manifest.heroHalf === 2 ? '--gv-half-b' : '--gv-half-a']))); },
      gvpc: function () { return rep('gvpc', blob(uniform(COLORS['--gv-gvpc']))); },
      density: densityRep
    };

    /** a component not listed above: one patch copy (p-NN) or one rib of the shell (hN-k) */
    function made(name) {
      if (/^hy-/.test(name)) return rep(name.slice(3), hydro());   // hy-<assembly id>: the preview's surfaces
      if (/^p-/.test(name)) {
        var rib = manifest.patch.filter(function (p) { return p.id === name; })[0].rib;
        return rep(name, cartoon(uniform(COLORS[rib % 2 ? '--gv-rib-b' : '--gv-rib-a'])));   // ribs alternate
      }
      var colour = uniform(COLORS[name.charAt(1) === '1' ? '--gv-half-a' : '--gv-half-b']);
      return rep(name.charAt(0) === 'q' ? 'mono' : name, blob(colour));   // q…: pool; h…: tip or block
    }

    var structuresP = {};
    function structure(id) {
      if (!structuresP[id]) {
        structuresP[id] = plugin.builders.structure.createStructure(model, { name: 'assembly', params: { id: id } });
      }
      return structuresP[id];
    }

    async function rep(id, props) {
      var sel = await plugin.builders.structure.representation.addRepresentation(
        await structure(id), props, { initialState: { isHidden: true } });
      return sel.ref;
    }

    function volumeReps() {
      var VR = lib.plugin.PluginStateObject.Volume.Representation3D, refs = [];
      plugin.state.data.cells.forEach(function (cell) { if (VR.is(cell.obj)) refs.push(cell.transform.ref); });
      return refs;
    }

    async function densityRep() {
      var before = volumeReps();
      await viewer.loadVolumeFromUrl({ url: base + 'emd-14238.mrc', format: 'ccp4', isBinary: true },
        [{ type: 'absolute', value: manifest.isovalue, color: COLORS['--gv-density'], alpha: 1 }]);
      var ref = volumeReps().filter(function (r) { return before.indexOf(r) < 0; })[0];
      if (!ref) throw new Error('density representation not found');
      // a pale, glassy skin: the fitted model stays readable inside it
      await plugin.build().to(ref).update(function (old) {
        old.type.params.alpha = 0.55;
        old.type.params.material = LOOK.material;
        old.type.params.xrayShaded = true;
        return old;
      }).commit();
      return ref;
    }

    var pending = {};
    function ensure(name) {
      if (built[name]) return Promise.resolve(built[name]);
      if (!pending[name]) {
        pending[name] = (COMPONENTS[name] ? COMPONENTS[name]() : made(name)).then(function (ref) {
          built[name] = ref;
          // whatever it was created as, it shows only if the current step wants it
          if (state.step >= 0) setVisible(name, wanted(state.step).indexOf(name) >= 0 && !state.revealing);
          return ref;
        });
      }
      return pending[name];
    }

    function setVisible(name, on) {
      if (!built[name]) return;
      var cell = plugin.state.data.cells.get(built[name]);
      if (cell && !!cell.state.isHidden === on) plugin.state.data.updateCellState(built[name], { isHidden: !on });
    }

    /** the component names a step shows, with its toggle applied */
    function group(name) {
      if (!manifest) return [];
      if (name === '@hydro4') return manifest.hydro4.ids.map(function (id) { return 'hy-' + id; });
      if (name === '@patch') return manifest.patch.map(function (p) { return p.id; });
      var g = manifest.growth, out = [];
      if (name === '@shell') {
        Object.keys(g.halves).forEach(function (h) { out = out.concat([g.halves[h].tip], g.halves[h].blocks); });
        return out;
      }
      if (name === '@pool') {                       // single GvpA, placed one by one while a block fills
        Object.keys(g.halves).forEach(function (h) {
          for (var q = 0; q < g.block; q++) out.push('q' + h + '-' + (q < 10 ? '0' : '') + q);
        });
        return out;
      }
      return [name];
    }

    function wanted(i) {
      var s = STEPS[i], set = [];
      s.show.forEach(function (n) { set = set.concat(n.charAt(0) === '@' ? group(n) : [n]); });
      if (s.toggle) {
        var t = s.toggle;
        if (state.toggles[t.key]) {
          set = set.filter(function (n) { return t.hide.indexOf(n) < 0; }).concat(t.show);
        }
      }
      return set;
    }

    /** what must be built before step i runs: the vesicle's growth also draws from the pool,
        which no step shows at rest, so `wanted` doesn't list it; unbuilt, its units would just not appear */
    function needs(i) { return STEPS[i].key === 'vesicle' ? wanted(i).concat(group('@pool')) : wanted(i); }

    /* ---- camera ---- */
    function frame(step, ms) {
      if (!plugin || !manifest) return;
      var LA = lib.math.LinearAlgebra, Vec3 = LA.Vec3;
      var f = manifest.frames[step.frame] || manifest[step.frame];
      var cam = plugin.canvas3d.camera;
      var dir = Vec3.create.apply(null, f.dir || manifest.view.dir);
      var up = Vec3.create.apply(null, f.up || manifest.view.up);
      var target = Vec3.create.apply(null, f.target);
      var w = cam.viewport.width, h = cam.viewport.height, aspect = w / Math.max(h, 1);
      var radius = f.radius;
      // from 768 up the step card sits over the scene's top-left: aim a little up and left
      // of the model, so it is drawn right of and below the card
      var card = !preview && getComputedStyle(cardEl).position === 'absolute';
      var right = Vec3.cross(Vec3(), dir, up), shiftX = 0, shiftY = 0;
      if (f.halfWidth) {
        // the vesicle: long axis (z) across the wider side of the canvas
        if (aspect >= 1) radius = Math.max(f.halfHeight, f.halfWidth / aspect);
        else { up = Vec3.create(0, 0, 1); radius = Math.max(f.halfHeight, f.halfWidth * aspect); }
        radius *= card ? 1.12 : 1.04;
        if (card) shiftY = 0.14;
      } else if (card) {
        shiftX = 0.32; shiftY = 0.1;
      }
      Vec3.scaleAndAdd(target, target, right, -shiftX * radius);
      Vec3.scaleAndAdd(target, target, up, shiftY * radius);
      var snap = cam.getInvariantFocus(target, radius, up, dir);
      // with manualReset, Mol* no longer sizes the far clipping plane to the scene: keep
      // the whole vesicle inside it
      var v = manifest.frames.vesicle;
      snap.radiusMax = Math.hypot(v.halfWidth, v.halfHeight) * 1.5;
      // straight to the camera: Mol*'s camera manager drops a move requested while
      // nothing is visible, which is how the vesicle's growth starts
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

      var t = s.toggle, toggled = t && state.toggles[t.key];
      card.innerHTML = (STEPS.length > 1 ? '<span class="gvv-eyebrow">Step ' + (i + 1) + ' of ' + STEPS.length + '</span>' : '') +
        '<strong class="gvv-title">' + s.title + '</strong>' +
        '<span class="gvv-text">' + (toggled && t.text ? t.text : s.text) + '</span>';

      chips.innerHTML = '';
      if (t) {
        var c = el('button', 'gvv-chip', icon(t.icon) + '<span>' + t.label + '</span>');
        c.type = 'button';
        c.setAttribute('aria-pressed', toggled ? 'true' : 'false');
        c.addEventListener('click', function () {
          state.toggles[t.key] = !state.toggles[t.key];
          paintStep(state.step < 0 ? 0 : state.step);
          toggleFade();
        });
        chips.appendChild(c);
      }

      legend.innerHTML = '';
      if (s.legendAlt === 'hydro' && (toggled || !t)) {
        legend.innerHTML = '<li class="gvv-key gvv-key--scale"><span class="gvv-key-label">hydrophilic</span>' +
          '<span class="gvv-swatch gvv-swatch--scale"></span><span class="gvv-key-label">hydrophobic</span></li>';
      } else {
        s.legend.forEach(function (k) {
          if (k[0] === '--gv-density' && t && !state.toggles[t.key]) return;
          var li = el('li', 'gvv-key',
            '<span class="gvv-swatch' + (k[2] ? ' gvv-swatch--' + k[2] : '') + '" style="--gvv-swatch: var(' + k[0] + ')"></span>' +
            '<span class="gvv-key-label">' + k[1] + '</span>');
          legend.appendChild(li);
        });
      }
    }

    /* ---- transitions (forward only; going back just switches) ----
       Each moving part is its own representation and moves by one matrix (Mol*'s
       per-copy "explode" state does not reach instanced visuals in this build):
       - step 2: the 44 GvpA copies drop onto the wall one by one, rib by rib;
       - step 3: the vesicle grows the way vesicles do: two cone tips meet as a small
         bicone, then ribs are inserted at the central seam and push the halves apart;
         the step-1 GvpA rides along with its half, so it always sits on the shell;
       - step 4: the GvpC cage starts wider and turned, then tightens onto the shell. */
    var M4, V3;
    var moved = {};                                   // components not at home right now

    function setMatrix(name, m) {
      var cell = built[name] && plugin.state.data.cells.get(built[name]);
      if (!cell || !cell.obj) return;
      moved[name] = !!m;
      cell.obj.data.repr.setState({ transform: m || M4.identity() });
      plugin.canvas3d.update(cell.obj.data.repr);
    }
    function home(name) { if (moved[name]) setMatrix(name, null); }
    function shift(x, y, z) { return M4.fromTranslation(M4(), V3.create(x, y, z)); }

    function isVisible(name) {
      var cell = built[name] && plugin.state.data.cells.get(built[name]);
      return !!cell && !cell.state.isHidden;
    }

    function ease(k) { return k < .5 ? 4 * k * k * k : 1 - Math.pow(-2 * k + 2, 3) / 2; }   // in-out cubic
    function settleIn(k) { var c = 1.4, j = k - 1; return 1 + (c + 1) * j * j * j + c * j * j; } // out, a little past
    function nextFrame() { return new Promise(function (r) { requestAnimationFrame(function () { r(); }); }); }

    /* ---- the chips (Hydropathy, Density) cross-fade: what leaves fades out as what enters fades in ---- */
    var cur = {};                                     // the alpha each component was last set to
    function baseAlpha(n) { return n === 'density' ? 0.55 : 1; }         // the density skin is glassy by design
    function alphaOf(n) { return cur[n] == null ? baseAlpha(n) : cur[n]; }
    // never exactly 0: Mol* leaves a fully transparent part out of the visible bounds, which can clamp the camera
    function setAlpha(n, a) {
      a = Math.max(a, 0.01);
      if (!built[n] || alphaOf(n) === a) return Promise.resolve();
      cur[n] = a;
      return plugin.build().to(built[n]).update(function (old) { old.type.params.alpha = a; return old; }).commit();
    }
    async function tween(jobs, ms, token) {
      jobs.forEach(function (j) { j.from = alphaOf(j.n); });
      var t0 = performance.now(), p;
      do {
        p = Math.min(1, (performance.now() - t0) / ms);
        var e = p * p * (3 - 2 * p);                  // smoothstep
        await Promise.all(jobs.map(function (j) { return setAlpha(j.n, j.from + (j.to - j.from) * e); }));
        if (token !== state.token) return false;
        if (p < 1) await nextFrame();
      } while (p < 1);
      return true;
    }
    async function toggleFade() {
      if (!plugin) return;
      if (stillMotion()) return apply(state.step, false);
      var want = wanted(state.step);
      await Promise.all(want.map(ensure));
      var token = ++state.token;
      var names = Object.keys(built);
      var leaving = names.filter(function (n) { return want.indexOf(n) < 0 && isVisible(n); });
      var entering = want.filter(function (n) { return !isVisible(n); });
      await Promise.all(entering.map(function (n) { return setAlpha(n, 0); }));
      if (token !== state.token) return;
      entering.forEach(function (n) { setVisible(n, true); });
      var jobs = entering.map(function (n) { return { n: n, to: baseAlpha(n) }; })
        .concat(leaving.map(function (n) { return { n: n, to: 0 }; }));
      if (await tween(jobs, 520, token)) leaving.forEach(function (n) { setVisible(n, false); setAlpha(n, baseAlpha(n)); });
    }

    var PLANS = {
      rib: {
        names: function (want) { return want.filter(function (n) { return /^p-/.test(n) || n === 'density'; }); },
        tracks: function (want) {
          // from in front of the wall and a little above: it falls onto its place
          var d = manifest.view.dir, u = manifest.view.up;
          var o = [-d[0] * 55 + u[0] * 18, -d[1] * 55 + u[1] * 18, -d[2] * 55 + u[2] * 18];
          var list = manifest.patch.map(function (p) { return p.id; });
          var tracks = list.map(function (id, j) {
            return { delay: 150 + j * 45, ms: 320,
              start: function () { setMatrix(id, shift(o[0], o[1], o[2])); setVisible(id, true); },
              tick: function (k) {
                var v = 1 - settleIn(k);
                setMatrix(id, k < 1 ? shift(o[0] * v, o[1] * v, o[2] * v) : null);
              } };
          });
          if (want.indexOf('density') >= 0) {
            tracks.push({ delay: 150 + list.length * 45 + 400, ms: 0, start: function () { setVisible('density', true); } });
          }
          return tracks;
        }
      },
      vesicle: {
        names: function (want) {
          return want.filter(function (n) { return /^h\d-/.test(n) || n === 'heroBlob'; });
        },
        tracks: function () {
          // Each half is one helical filament that grows at the seam, a GvpA at a time.
          // Every unit added there pushes the whole filament one step along its helix:
          // it turns about the long axis and moves out. So with n units grown, the built
          // part (units 0..n-1, tip first) is turned and shifted until its newest unit
          // sits where the last one belongs, at the seam. Finished blocks move as one;
          // the block being filled is drawn unit by unit from the pool.
          var g = manifest.growth, K = g.block, ax = g.axis, halves = Object.keys(g.halves);
          var geo = {};
          halves.forEach(function (h) {
            var ops = shell[h], phi = [], z = [], prev = null, acc = 0;
            ops.forEach(function (o) {
              var a = Math.atan2(o[10] - ax[1], o[9] - ax[0]);
              if (prev !== null) {
                var d = a - prev;
                d -= 2 * Math.PI * Math.round(d / (2 * Math.PI));   // unwrapped: turns add up
                acc += d;
              }
              prev = a;
              phi.push(acc);
              z.push(o[11]);
            });
            geo[h] = { phi: phi, z: z, n: ops.length };
          });
          function at(arr, j) {
            var i = Math.min(Math.floor(j), arr.length - 2), f = j - i;
            return arr[i] + (arr[i + 1] - arr[i]) * f;
          }
          // the screw that takes unit j (fractional) to where the seam end belongs
          function screw(h, j) {
            var G = geo[h], last = G.n - 1;
            if (j >= last) return null;
            var ang = G.phi[last] - at(G.phi, j), dz = G.z[last] - at(G.z, j);
            var m = shift(ax[0], ax[1], dz);
            M4.mul(m, m, M4.fromRotation(M4(), ang, V3.create(0, 0, 1)));
            return M4.mul(m, m, shift(-ax[0], -ax[1], 0));
          }
          function placed(T, o) {                              // T x (this unit's own placement)
            var S = M4.identity();
            S[0] = o[0]; S[4] = o[1]; S[8] = o[2];
            S[1] = o[3]; S[5] = o[4]; S[9] = o[5];
            S[2] = o[6]; S[6] = o[7]; S[10] = o[8];
            S[12] = o[9]; S[13] = o[10]; S[14] = o[11];
            return T ? M4.mul(M4(), T, S) : S;
          }
          return [{ delay: 0, ms: 5600, tick: function (k) {
            var e = k < 1 ? ease(k) : 1;
            halves.forEach(function (h) {
              var H = g.halves[h], N = geo[h].n;
              var nf = 1 + (N - 1) * e, n = Math.floor(nf);
              var T = screw(h, nf - 1);
              setMatrix(H.tip, T); setVisible(H.tip, true);
              H.blocks.forEach(function (id, b) {
                var whole = Math.min((b + 1) * K, N) <= n;
                if (whole) setMatrix(id, T);
                setVisible(id, whole);
              });
              var from = Math.floor(n / K) * K;              // the block being filled
              for (var q = 0; q < K; q++) {
                var pos = from + q, id = 'q' + h + '-' + (q < 10 ? '0' : '') + q;
                var show = pos < n && pos < N && !(h === g.hero.half && pos === g.hero.pos);
                if (show) setMatrix(id, placed(T, shell[h][pos]));
                setVisible(id, show);
              }
              if (h === g.hero.half) {
                var heroIn = n > g.hero.pos;
                if (heroIn) setMatrix('heroBlob', T);
                setVisible('heroBlob', heroIn);
              }
            });
          } }];
        }
      },
      gvpc: {
        names: function (want) { return want.indexOf('gvpc') >= 0 ? ['gvpc'] : []; },
        tracks: function () {
          // about the vesicle's long axis (z) through its centre: wider, turned, then home
          var c = manifest.frames.vesicle.target;
          function wrap(v) {
            if (!v) return null;
            var s = 1 + 0.45 * v, m = shift(c[0], c[1], c[2]);
            M4.mul(m, m, M4.fromRotation(M4(), 0.7 * v, V3.create(0, 0, 1)));
            M4.mul(m, m, M4.fromScaling(M4(), V3.create(s, s, s)));
            return M4.mul(m, m, shift(-c[0], -c[1], -c[2]));
          }
          return [{ delay: 100, ms: 1500,
            start: function () { setMatrix('gvpc', wrap(1)); setVisible('gvpc', true); },
            tick: function (k) { setMatrix('gvpc', wrap(1 - ease(k))); } }];
        }
      }
    };

    /** run tracks {delay, ms, start(), tick(k)} in one frame loop; false if overtaken */
    async function run(tracks, token) {
      var t0 = performance.now(), done = false;
      while (!done) {
        if (token !== state.token) return false;
        var now = performance.now() - t0;
        done = true;
        tracks.forEach(function (t) {
          if (now < t.delay) { done = false; return; }
          if (t.start) { t.start(); t.start = null; }
          if (t.finished) return;
          var k = t.ms ? Math.min((now - t.delay) / t.ms, 1) : 1;
          if (t.tick) t.tick(k);
          if (k < 1) done = false; else t.finished = true;
        });
        plugin.canvas3d.requestDraw();
        if (!done) await nextFrame();
      }
      return true;
    }

    /** show exactly what step i wants; with `animate`, the step's plan moves its parts in.
        `ready` runs once the scene is set up, just before anything moves: the camera
        starts there, so nothing that is still being prepared cuts its flight short. */
    async function apply(i, animate, ready) {
      if (!plugin) return;
      var token = ++state.token;
      var want = wanted(i);
      await Promise.all(needs(i).map(ensure));
      if (token !== state.token) return;
      M4 = lib.math.LinearAlgebra.Mat4; V3 = lib.math.LinearAlgebra.Vec3;
      Object.keys(cur).forEach(function (n) { setAlpha(n, baseAlpha(n)); });         // a fade cut short
      var plan = animate && !stillMotion() ? PLANS[STEPS[i].key] : null;
      var managed = plan ? plan.names(want) : [];
      Object.keys(built).forEach(function (n) {
        if (want.indexOf(n) < 0 || managed.indexOf(n) >= 0) setVisible(n, false);
        else { home(n); setVisible(n, true); }
      });
      if (ready) ready();
      if (plan && !(await run(plan.tracks(want), token))) return;
      want.forEach(home);
    }

    function go(i, first) {
      if (i === state.step) return state.busy;
      var from = state.step;
      state.step = i;
      paintStep(i);
      if (!plugin) return state.busy;
      host.classList.add('is-busy');
      state.busy = (async function () {
        var animate = !first && from >= 0 && i > from;
        state.revealing = animate && !stillMotion();   // parts built now stay hidden until they fly in
        var needsBuild = needs(i).some(function (n) { return !built[n]; });
        if (needsBuild && !first) veilBuild(true);
        try {
          await Promise.all(needs(i).map(ensure));
        } finally { veilBuild(false); }
        if (state.step !== i) return;
        await apply(i, animate, function () {
          // the vesicle grows in its own view: cut there instead of flying
          frame(STEPS[i], first || (STEPS[i].key === 'vesicle' && animate) ? 0 : 1100);
        });
      })().finally(function () {
        if (state.step === i) { host.classList.remove('is-busy'); state.revealing = false; }
      });
      return state.busy;
    }

    function veilBuild(on) {
      if (on) {
        veil.querySelector('.gvv-veil-text').textContent = 'Building the model…';
        veil.classList.add('is-soft');
        veil.hidden = false;
      } else if (host.classList.contains('is-ready')) {
        veil.hidden = true;
        veil.classList.remove('is-soft');
      }
    }
  }

  function init() {
    Array.prototype.forEach.call(document.querySelectorAll('[data-gv-viewer]'), function (h) {
      if (!h.__gvv) { h.__gvv = true; mount(h); }
    });
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
  else init();
})();
