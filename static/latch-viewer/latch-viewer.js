/* Latch (biosensor) viewer: the interactive 3D figure that shows how the GvpA binder construct senses a target.
 *
 * Mounted onto every [data-latch-viewer] element (embed block: proteins/latch_animation/README.md). Built like
 * the gas-vesicle viewer (gv-viewer.js, same card, step buttons and tools; load gv-viewer.css first): Mol* 5.13
 * (MIT, vendored in ../molstar/) is only fetched once the figure scrolls near the viewport.
 *
 * Data (built by proteins/latch_animation/tools/prep_viewer.py, then prep_wall.py):
 *   latch-assembly.cif   every part as its own mmCIF assembly: binder, latch, tbd (target-binding domain, with the
 *                        linker), target; and the polypeptide in unbroken chains: cont2 (binder + latch),
 *                        cont3 (+ linker + domain), arm (latch + linker + domain from the hinge)
 *   latch-wall.ply.gz    the vesicle wall, one coarse surface mesh (prep_wall.py): no atoms, not a structure
 *   latch-viewer.json    hinge, fly-in vector, dock offset, camera frames
 * They load from this script's folder, or from data-base="…/" on the mount element.
 *
 * Every part is a rigid body, so every move is one matrix per representation:
 *   step 2  the latch swings in from its open pose onto the binding site (fading in);
 *   step 3  the target-binding domain fades in at the latch's free end;
 *   step 4  the latch opens first: latch + domain swing off the site as ONE block (block_hinge of the Drylab film:
 *           axis, angle and centre are in latch-viewer.json); then the target flies in and binds the open domain.
 *           The chip "Bind with latch closed" shuts the block again with the target on it, to look at the clash;
 *   step 5  the wall fades in and the construct, which has stood off the wall until now, docks onto GvpA.
 *
 * One cartoon: Mol* ends a cartoon at the end of its chain, so parts that are chains of their own show a seam where they
 * join. Where the parts sit together (steps 2 and 3 at rest) the polypeptide is drawn as one unbroken chain (cont2, cont3);
 * the pieces (binder, latch, tbd, arm) are only used while one part moves against another, and the unbroken chain fades in
 * over them when the move is done. Across the hinge itself (open state) the cartoon is, truly, in two pieces.
 * Colours are --lv-* properties in latch-viewer.css, read back at start-up.
 */
(function () {
  'use strict';

  var HERE = (document.currentScript && document.currentScript.src || '').replace(/[^/]*$/, '');
  var MOLSTAR = HERE + '../molstar/';

  /* ---------------------------------------------------------------- icons -- */
  /* 24 x 24, drawn with currentColor strokes so they take the button's colour. */
  var ICON = {
    // a helical bundle: three capsules
    binder: '<rect x="3.5" y="5" width="4.5" height="14" rx="2.25" /><rect x="9.75" y="5" width="4.5" height="14" rx="2.25" />' +
            '<rect x="16" y="5" width="4.5" height="14" rx="2.25" />',
    // the bundle with a strand folded over it
    latch: '<rect x="3.5" y="10" width="4.5" height="10" rx="2.25" /><rect x="9.75" y="10" width="4.5" height="10" rx="2.25" />' +
           '<path d="M5 6.5h9a4.5 4.5 0 0 1 4.5 4.5v5" stroke-width="2.4" />',
    // ... which ends in a small domain
    tbd: '<rect x="3.5" y="11" width="4.5" height="9" rx="2.25" /><rect x="9.75" y="11" width="4.5" height="9" rx="2.25" />' +
         '<path d="M5 7.5h6" stroke-width="2.4" /><circle cx="16.5" cy="7.5" r="3.5" />',
    // a domain and the thing that binds it
    target: '<circle cx="7.5" cy="12" r="4" /><path d="M12.5 8.5a4 4 0 0 1 0 7" /><circle cx="18" cy="12" r="2.5" />',
    // the whole vesicle: cylinder with conical tips, ribs across it
    vesicle: '<path d="M2 12 6.5 6.5h11L22 12l-4.5 5.5h-11Z" /><path d="M9 6.5v11M12 6.5v11M15 6.5v11" />',
    // two shapes overlapping
    clash: '<circle cx="8.5" cy="12" r="5" /><circle cx="15.5" cy="12" r="5" />',
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

  /* ---------------------------------------------------------------- steps -- */
  /* `show` names components (= assemblies of latch-assembly.cif). The block is what swings in step 4. */
  var BLOCK = ['latch', 'tbd', 'arm', 'target', 'targetX'];           // what swings (targetX: the clash view's PD-L1)
  var CONSTRUCT = ['binder', 'binderX', 'cont2', 'cont3'].concat(BLOCK);
  var CLASH_STEP = 3;                                                // the step that has the clash view

  var STEPS = [
    { key: 'binder', label: 'GvpA binder', icon: 'binder', frame: 'binder',
      title: 'A designed GvpA binder',
      text: 'The construct starts from a protein binder designed to recognise GvpA, the building block of the gas vesicle ' +
            'wall. Its GvpA-binding site (light blue) lies on one face. Bound to the wall, it stabilises the vesicle.',
      show: ['binder'],
      legend: [['--lv-binder', 'GvpA binder'], ['--lv-site', 'GvpA-binding site']] },
    { key: 'latch', label: '+ Latch', icon: 'latch', frame: 'latch',
      title: 'The latch blocks the site',
      text: 'A latch is fused to the end of the binder. It folds back over the binding site and covers it, so the ' +
            'construct cannot bind GvpA yet. The latch is a designed backbone: it has no amino-acid sequence yet.',
      show: ['cont2'], anim: ['binder', 'latch'],                // at rest one chain; the latch swings in as a piece
      legend: [['--lv-binder', 'GvpA binder'], ['--lv-site', 'GvpA-binding site'], ['--lv-latch', 'latch']] },
    { key: 'tbd', label: '+ Binding domain', icon: 'tbd', frame: 'tbd',
      title: 'A target-binding domain at the latch’s free end',
      text: 'A second domain is joined to the free end of the latch by a short flexible linker and recognises the ' +
            'analyte. Here it is a binder against PD-L1, which stands in for any target: swap the domain and the ' +
            'sensor senses something else.',
      show: ['cont3'], anim: ['cont2'],                          // the chain with the domain fades in over the one without
      legend: [['--lv-binder', 'GvpA binder'], ['--lv-latch', 'latch'], ['--lv-tbd', 'target-binding domain']] },
    { key: 'target', label: '+ Target', icon: 'target', frame: 'mechanism',
      title: 'The latch opens, then the target binds',
      text: 'With the latch closed, PD-L1 cannot sit on its domain: it would run into the GvpA binder. So the latch ' +
            'opens first. Latch and domain swing off the binding site as one block, a hinge motion placed so that ' +
            'nothing collides (design intent, not a simulation). Then PD-L1 binds the open state and holds it open.',
      show: ['binder', 'arm', 'target'],                         // open: the binder and the block are two pieces
      toggle: { key: 'clash', label: 'Bind with latch closed', icon: 'clash',
                text: 'Latch shut, PD-L1 on its domain: PD-L1 runs into the GvpA binder, and the residues of both that ' +
                      'overlap are marked in red-orange. This is the state the design rules out, so PD-L1 can only ' +
                      'bind once the latch is open. (The binder model has no side chains, so the real overlap is larger.) ' +
                      'Drag to turn it.' },
      legend: [['--lv-binder', 'GvpA binder'], ['--lv-latch', 'latch'], ['--lv-tbd', 'target-binding domain'],
               ['--lv-target', 'target protein']],
      legendClash: [['--lv-clash', 'residues in clash']] },
    { key: 'vesicle', label: '+ Gas vesicle', icon: 'vesicle', frame: 'dock',
      title: 'With the site free, the binder docks onto GvpA',
      text: 'The open construct binds GvpA in the vesicle wall and stabilises it: the readout is a shift in collapse ' +
            'pressure. Shown is a patch of the wall of a Bacillus megaterium vesicle (Huber et al., Cell 2023; Zenodo ' +
            '10.5281/zenodo.6458345); the binder sits on the five GvpA it was designed against.',
      show: ['binder', 'arm', 'target', 'wall'],
      legend: [['--lv-binder', 'GvpA binder'], ['--lv-latch', 'latch'], ['--lv-tbd', 'target-binding domain'],
               ['--lv-target', 'target protein'], ['--lv-wall', 'vesicle wall']] }
  ];

  /* What each step does when it is entered going forward; `enter` fades in, `move` is moved by the plan. */
  var RESTING = function (i, clash) { return { rot: i >= 3 && !clash ? 1 : 0, away: i >= 4 ? 0 : 1 }; };

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

  /* ================================================================ viewer == */
  function mount(host) {
    var base = host.getAttribute('data-base') || HERE;
    var state = { step: -1, busy: Promise.resolve(), token: 0, clash: false };

    var queue = [];
    var idle = window.requestIdleCallback
      ? function (f, o) { window.requestIdleCallback(f, o); }
      : function (f) { setTimeout(f, 50); };
    /** background builds (see boot), and: the reader is about to look at step i, so build its parts first */
    function prefer(i, clash) {
      var first = [];
      needed(i).concat(clash ? ['binderX', 'targetX', 'latch', 'tbd'] : []).forEach(function (n) { if (!built[n]) first.push(n); });
      queue = first.concat(queue.filter(function (n) { return first.indexOf(n) < 0; }));
    }

    /* ---- markup ---- */
    host.innerHTML = '';
    host.classList.add('gvv', 'lv');
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

    var steps = el('ol', 'gvv-steps');
    steps.setAttribute('aria-label', 'Sensing steps');
    STEPS.forEach(function (s, i) {
      var li = el('li', 'gvv-step');
      var b = el('button', 'gvv-step-btn',
        '<span class="gvv-ring">' + icon(s.icon) + '</span>' +
        '<span class="gvv-step-label"><span class="gvv-step-n">' + (i + 1) + '</span>' + s.label + '</span>');
      b.type = 'button';
      b.title = s.title;
      b.addEventListener('click', function () { go(i); });
      b.addEventListener('pointerenter', function () { prefer(i); });
      b.addEventListener('focus', function () { prefer(i); });
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
        console.error('[latch-viewer]', err);
        fail('The 3D model could not be loaded. Reload the page to try again.');
      });
    }
    if ('IntersectionObserver' in window) {
      var io = new IntersectionObserver(function (es) {
        if (es.some(function (x) { return x.isIntersecting; })) { io.disconnect(); start(); }
      }, { rootMargin: '150% 0px' });
      io.observe(host);
    } else start();

    /* once the page itself has loaded, fetch Mol* and the data at low priority, so the download is done before the
       reader gets here (prefetch reuses the HTTP cache; browsers without it simply load on boot, as before) */
    function prefetch() {
      if (started) return;                       // already loading for real
      [MOLSTAR + 'molstar.js', MOLSTAR + 'molstar.css', base + 'latch-viewer.json', base + 'latch-assembly.cif']
        .forEach(function (href) {
          var l = document.createElement('link');
          l.rel = 'prefetch'; l.href = href;
          document.head.appendChild(l);
        });
    }
    if (document.readyState === 'complete') idle(prefetch, { timeout: 3000 });
    else window.addEventListener('load', function () { idle(prefetch, { timeout: 3000 }); });

    function fail(msg) {
      veil.classList.add('is-error');
      veil.querySelector('.gvv-veil-text').textContent = msg;
      veil.hidden = false;
    }

    var viewer, plugin, lib, manifest, model, built = {}, COLORS, M4, V3;

    async function boot() {
      loadStyle(MOLSTAR + 'molstar.css');
      var manifestP = fetch(base + 'latch-viewer.json').then(function (r) {
        if (!r.ok) throw new Error('latch-viewer.json: HTTP ' + r.status);
        return r.json();
      });
      await loadScript(MOLSTAR + 'molstar.js');
      lib = window.molstar.lib;
      M4 = lib.math.LinearAlgebra.Mat4;
      V3 = lib.math.LinearAlgebra.Vec3;
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
      COLORS = readColours();
      styleCanvas();

      var data = await plugin.builders.data.download(
        { url: base + 'latch-assembly.cif', isBinary: false, label: 'Latch construct' }, { state: { isGhost: true } });
      var traj = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
      model = await plugin.builders.structure.createModel(traj);

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

      // Mol* draws nothing while the scene is still; the one thing that keeps it drawing is the spin, so that stops
      // while the figure is scrolled away and goes on when it is back. (Not canvas3d.pause(): its resume() can leave
      // two frame loops running.)
      if ('IntersectionObserver' in window) {
        new IntersectionObserver(function (es) {
          state.offscreen = !es[es.length - 1].isIntersecting;
          if (state.spin) applySpin();
        }).observe(stage);
      }

      /* build what the later steps need while the reader looks at step 1, so a click is instant, in the order the steps
         need them (the clash view's parts last: they only matter once its chip is pressed), each in an idle moment.
         Every build takes ~30 ms (the wall mesh a good deal more, so it goes last), so none of them is felt. Pointing at (or tabbing to) a
         step button or the clash chip moves its parts to the front. */
      queue = ['latch', 'cont2', 'cont3', 'arm', 'target', 'tbd', 'binderX', 'targetX', 'wall'];
      (async function () {
        await new Promise(function (ok) { idle(ok, { timeout: 250 }); });
        while (queue.length) {
          var name = queue.shift();
          if (built[name]) continue;
          await new Promise(function (ok) { idle(ok, { timeout: 250 }); });
          try { await ensure(name); } catch (e) { console.warn('[latch-viewer] prebuild', name, e); }
        }
      })();
    }

    function readColours() {
      var c = {};
      ['--lv-binder', '--lv-site', '--lv-latch', '--lv-tbd', '--lv-target', '--lv-wall', '--lv-clash',
       '--gv-outline', '--gv-highlight', '--gv-bg'].forEach(function (n) { c[n] = colour(host, n); });
      return c;
    }

    function styleCanvas() {
      var p = plugin.canvas3d.props;
      plugin.canvas3d.setProps({
        transparentBackground: true,
        checkeredTransparentBackground: false,
        camera: Object.assign({}, p.camera, {
          // the viewer places the camera itself: Mol* would otherwise re-aim it whenever nothing visible is in front of it
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
            scale: 1, threshold: 0.33, color: COLORS['--gv-outline'], includeTransparent: false } }  // false: no outline on the wall mesh (alpha < 1); parts lose theirs while they fade
        }),
        cameraFog: { name: 'on', params: { intensity: 12 } }
      });
    }

    /* ---- components: built once, on demand, then only shown, hidden or moved ---- */
    var LOOK = { celShaded: true, material: { metalness: 0, roughness: 1, bumpiness: 0 } };
    function uniform(c) { return { color: 'uniform', colorParams: { value: c } }; }
    function cartoon(c, size) {
      return Object.assign({ type: 'cartoon', typeParams: Object.assign({ sizeFactor: size || 0.3 }, LOOK) }, uniform(c));
    }

    // The prep script stores the part of each residue in the B-factor column: 0 binder, 1 GvpA-binding site, 2 clash with
    // PD-L1 (the clash view), 3 latch, 4 linker and target-binding domain. The theme runs high to low along its list, so the
    // highest value comes first. One thickness for every chain, so that an unbroken chain has no step where one part meets
    // the next (the pieces are drawn the same, so the hand-over between them does not show).
    function chain(clash) {
      return { type: 'cartoon', typeParams: Object.assign({ sizeFactor: 0.34 }, LOOK),
        color: 'uncertainty', colorParams: { domain: [0, 4],
          list: { kind: 'interpolate', colors: [COLORS['--lv-tbd'], COLORS['--lv-latch'],
            COLORS[clash ? '--lv-clash' : '--lv-binder'], COLORS['--lv-site'], COLORS['--lv-binder']] } } };
    }

    var COMPONENTS = {
      binder: function () { return rep('binder', chain(false)); },
      // the clash view keeps the cartoon and paints the residues that overlap PD-L1 in the clash colour
      binderX: function () { return rep('binder', chain(true)); },
      targetX: function () {
        return rep('target', Object.assign({ type: 'cartoon', typeParams: Object.assign({ sizeFactor: 0.3 }, LOOK),
          color: 'uncertainty', colorParams: { domain: [0, 1],
            list: { kind: 'interpolate', colors: [COLORS['--lv-clash'], COLORS['--lv-target']] } } }));
      },
      latch: function () { return rep('latch', chain(false)); },          // binder's last residue + latch
      tbd: function () { return rep('tbd', chain(false)); },              // latch's last residue + linker + domain
      arm: function () { return rep('arm', chain(false)); },              // latch + linker + domain, from the hinge
      cont2: function () { return rep('cont2', chain(false)); },          // binder + latch, unbroken
      cont3: function () { return rep('cont3', chain(false)); },          // + linker + domain, unbroken
      target: function () { return rep('target', cartoon(COLORS['--lv-target'])); },
      wall: wall,
    };
    // the wall is a mesh (latch-wall.ply.gz, made by prep_wall.py), not a structure: no atoms, nothing to pick
    function wall() {
      var T = lib.plugin.StateTransforms, ghost = { state: { isGhost: true } };
      return fetch(base + manifest.wallMesh).then(function (res) {
        if (!res.ok) throw new Error(manifest.wallMesh + ': HTTP ' + res.status);
        return new Response(res.body.pipeThrough(new DecompressionStream('gzip'))).text();
      }).then(function (ply) {
        var b = plugin.state.data.build().toRoot()
          .apply(T.Data.RawData, { data: ply, label: 'Vesicle wall' }, ghost)
          .apply(T.Data.ParsePly, {}, ghost)
          .apply(T.Shape.ShapeFromPly, { label: 'Vesicle wall' }, ghost)
          .apply(T.Representation.ShapeRepresentation3D,
            // Smooth shading and a shade darker than the legend swatch, so that the pale mesh shows against the page. alpha just
            // under 1 puts it among the transparent parts, which the outline leaves out (see the outline in styleCanvas):
            // on a mesh the outline only draws loose dashes along the rib crests at the wall's rim.
            { coloring: { name: 'uniform', params: { color: COLORS['--lv-wall'], saturation: 0, lightness: -0.5 } },
              celShaded: false, material: LOOK.material, alpha: 0.99 },
            { initialState: { isHidden: true } });
        return b.commit().then(function () {
          plugin.state.data.cells.get(b.selector.ref).obj.data.repr.setState({ pickable: false });
          return b.selector.ref;
        });
      });
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

    var pending = {};
    function ensure(name) {
      if (built[name]) return Promise.resolve(built[name]);
      if (!pending[name]) {
        pending[name] = COMPONENTS[name]().then(function (ref) {
          built[name] = ref;
          // whatever it was created as, it shows only if the current step wants it, and sits where it rests
          if (state.step >= 0) {
            setVisible(name, wanted(state.step).indexOf(name) >= 0 && !state.revealing);
            if (!state.playing) setRest(name, state.step);
          }
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

    /** the component names a step shows at rest; the clash view swaps binder and PD-L1 for their painted versions,
        and the block (arm) for its pieces (latch, tbd) */
    function wanted(i) {
      var set = STEPS[i].show.slice();
      if (i === CLASH_STEP && state.clash) {
        set = [].concat.apply([], set.map(function (n) {
          return n === 'binder' ? ['binderX'] : n === 'arm' ? ['latch', 'tbd'] : n === 'target' ? ['targetX'] : [n];
        }));
      }
      return set;
    }

    /** ... and what has to be built to get there: also the pieces its move runs on before the chain takes over */
    function needed(i) {
      var set = wanted(i);
      (STEPS[i].anim || []).forEach(function (n) { if (set.indexOf(n) < 0) set.push(n); });
      return set;
    }

    /* ---- camera ---- */
    function frame(step, ms) {
      if (!plugin || !manifest) return;
      var Vec3 = V3;
      var f = manifest.frames[step.frame];
      var cam = plugin.canvas3d.camera;
      var dir = Vec3.create.apply(null, f.dir);
      var up = Vec3.create.apply(null, f.up);
      var target = Vec3.create.apply(null, f.target);
      // from 768 up the step card sits over the scene's top-left: aim a little up and left
      // of the model, so it is drawn right of and below the card
      var card = getComputedStyle(cardEl).position === 'absolute';
      var right = Vec3.cross(Vec3(), dir, up), shiftX = 0, shiftY = 0;
      if (card) { shiftX = 0.32; shiftY = 0.1; }
      Vec3.scaleAndAdd(target, target, right, -shiftX * f.radius);
      Vec3.scaleAndAdd(target, target, up, shiftY * f.radius);
      var snap = cam.getInvariantFocus(target, f.radius, up, dir);
      // with manualReset, Mol* no longer sizes the far clipping plane to the scene: keep the wall inside it
      snap.radiusMax = 700;
      cam.setState(snap, stillMotion() ? 0 : ms);
      plugin.canvas3d.requestDraw();
    }

    function setSpin(on) {
      state.spin = on;
      bSpin.setAttribute('aria-pressed', on ? 'true' : 'false');
      applySpin();
    }

    function applySpin() {
      if (!plugin) return;
      var on = state.spin && !state.offscreen;
      var tb = plugin.canvas3d.props.trackball;
      plugin.canvas3d.setProps({ trackball: Object.assign({}, tb, {
        animate: on ? { name: 'spin', params: { speed: 0.08, axis: V3.create(0, -1, 0) } } : { name: 'off', params: {} } }) });
    }

    /* ---- the step UI (works before and after the 3D part loads) ---- */
    function paintStep(i) {
      var s = STEPS[i];
      Array.prototype.forEach.call(steps.querySelectorAll('.gvv-step-btn'), function (b, j) {
        if (j === i) b.setAttribute('aria-current', 'step'); else b.removeAttribute('aria-current');
        b.parentNode.classList.toggle('is-done', j < i);
      });
      steps.style.setProperty('--gvv-progress', String(i / (STEPS.length - 1)));
      var t = s.toggle, toggled = t && state.clash;
      card.innerHTML = '<span class="gvv-eyebrow">Step ' + (i + 1) + ' of ' + STEPS.length + '</span>' +
        '<strong class="gvv-title">' + s.title + '</strong>' +
        '<span class="gvv-text">' + (toggled && t.text ? t.text : s.text) + '</span>';
      chips.innerHTML = '';
      if (t) {
        var c = el('button', 'gvv-chip', icon(t.icon) + '<span>' + t.label + '</span>');
        c.type = 'button';
        c.setAttribute('aria-pressed', toggled ? 'true' : 'false');
        c.addEventListener('click', function () { toggleClash(); });
        c.addEventListener('pointerenter', function () { prefer(i, true); });
        c.addEventListener('focus', function () { prefer(i, true); });
        chips.appendChild(c);
      }
      legend.innerHTML = '';
      (s.legend.concat(toggled && s.legendClash ? s.legendClash : [])).forEach(function (k) {
        legend.appendChild(el('li', 'gvv-key',
          '<span class="gvv-swatch" style="--gvv-swatch: var(' + k[0] + ')"></span>' +
          '<span class="gvv-key-label">' + k[1] + '</span>'));
      });
    }

    /* ---- poses: one matrix per part ----
       rot  0 closed .. 1 open (the block's swing about the hinge)
       fly  1 far out .. 0 bound (the target's approach, along the closed frame)
       away 1 standing off the wall .. 0 docked
       Order, right to left: fly (in the closed frame), then the swing, then the offset. */
    function shift(x, y, z) { return M4.fromTranslation(M4(), V3.create(x, y, z)); }

    function pose(name, rot, fly, away) {
      var a = manifest.away;
      var m = shift(a[0] * away, a[1] * away, a[2] * away);
      if (rot && BLOCK.indexOf(name) >= 0) {
        var h = manifest.hinge, c = h.center;
        var r = shift(c[0], c[1], c[2]);
        M4.mul(r, r, M4.fromRotation(M4(), rot * h.deg * Math.PI / 180, V3.create(h.axis[0], h.axis[1], h.axis[2])));
        M4.mul(r, r, shift(-c[0], -c[1], -c[2]));
        M4.mul(m, m, r);
      }
      if (fly && (name === 'target' || name === 'targetX')) {
        var f = manifest.fly;
        M4.mul(m, m, shift(f.dir[0] * f.dist * fly, f.dir[1] * f.dist * fly, f.dir[2] * f.dist * fly));
      }
      return m;
    }

    function setPose(name, rot, fly, away) {
      var cell = built[name] && plugin.state.data.cells.get(built[name]);
      if (!cell || !cell.obj) return;
      var repr = cell.obj.data.repr, m = pose(name, rot, fly, away);
      repr.setState({ transform: m });
      plugin.canvas3d.update(repr);
      // Mol* frustum-culls each mesh by its bounding sphere, which does not follow `transform`: zoomed in, a part whose
      // sphere had slid out of view was dropped while it was still on screen. The sphere moves with the part.
      repr.renderObjects.forEach(function (o) {
        V3.transformMat4(o.values.boundingSphere.ref.value.center, o.values.invariantBoundingSphere.ref.value.center, m);
      });
    }

    function setRest(name, i) {
      if (CONSTRUCT.indexOf(name) < 0) return;
      var r = RESTING(i, state.clash);
      setPose(name, r.rot, 0, r.away);
    }

    function ease(k) { return k < .5 ? 4 * k * k * k : 1 - Math.pow(-2 * k + 2, 3) / 2; }   // in-out cubic
    function nextFrame() { return new Promise(function (r) { requestAnimationFrame(function () { r(); }); }); }

    /* ---- fading: a part enters by its opacity going from 0 to 1 ----
       Set as the representation's alphaFactor: Mol* draws with alpha x alphaFactor in every pass, so this is the same
       picture as changing the `alpha` parameter, without a state transaction (and its rebuild check) on every frame. */
    var cur = {};                                     // the opacity each component was last set to
    // never exactly 0: Mol* leaves a fully transparent part out of the visible bounds, which can clamp the camera
    function setAlpha(n, a) {
      a = Math.max(a, 0.01);
      if (!built[n] || (cur[n] == null ? 1 : cur[n]) === a) return Promise.resolve();
      var cell = plugin.state.data.cells.get(built[n]);
      if (!cell || !cell.obj) return Promise.resolve();
      cur[n] = a;
      cell.obj.data.repr.setState({ alphaFactor: a });
      plugin.canvas3d.update(cell.obj.data.repr);
      return Promise.resolve();
    }
    function fade(names, a) { return Promise.all(names.map(function (n) { return setAlpha(n, a); })); }

    /* ---- transitions (forward only; going back just switches) ---- */
    /* `enter` starts invisible and is shown by a track; `keep` is visible while the move runs and gone after it; `move` is
       posed by the tracks. What the step shows at rest (STEPS[i].show) is what is left when the move is done. */
    var PLANS = {
      latch: {
        enter: ['latch', 'cont2'], keep: ['binder'], move: [],
        tracks: function () {
          var SWING = 1800;
          return [
            // it swings in from the open pose onto the site, fading in as it goes ...
            { ms: SWING,
              start: function () { setPose('latch', 1, 0, 1); setVisible('latch', true); },
              tick: function (k) { setPose('latch', 1 - ease(k), 0, 1); return setAlpha('latch', Math.min(1, k * 3)); } },
            // ... and, once it sits, the two pieces hand over to the unbroken chain
            { delay: SWING, ms: 350,
              start: function () { setPose('cont2', 0, 0, 1); setVisible('cont2', true); },
              tick: function (k) { return setAlpha('cont2', ease(k)); } }
          ];
        }
      },
      tbd: {
        enter: ['cont3'], keep: ['cont2'], move: [],
        tracks: function () {
          // the chain that goes on into linker and domain fades in over the one that ended at the latch
          return [{ ms: 1000,
            start: function () { setPose('cont3', 0, 0, 1); setVisible('cont3', true); },
            tick: function (k) { return setAlpha('cont3', ease(k)); } }];
        }
      },
      target: {
        enter: ['target'], move: ['arm'],
        tracks: function () {
          var SWING = 2000, PAUSE = 350, FLY = 1900;
          return [
            // 1. the latch opens: latch + linker + domain swing off the site as one block, PD-L1 is not there yet
            { ms: SWING,
              start: function () { setPose('arm', 0, 0, 1); setPose('target', 1, 1, 1); },
              tick: function (k) { setPose('arm', ease(k), 0, 1); } },
            // 2. PD-L1 comes in to the open domain (the way in is the closed frame's, swung with the block)
            { delay: SWING + PAUSE, ms: FLY,
              start: function () { setPose('target', 1, 1, 1); setVisible('target', true); },
              tick: function (k) { setPose('target', 1, 1 - ease(k), 1); return setAlpha('target', Math.min(1, k * 4)); } }
          ];
        }
      },
      vesicle: {
        enter: ['wall'], move: CONSTRUCT,
        tracks: function () {
          return [
            { ms: 900,
              start: function () {
                CONSTRUCT.forEach(function (n) { setPose(n, 1, 0, 1); });
                setVisible('wall', true);
              },
              tick: function (k) { return fade(['wall'], ease(k)); } },
            { delay: 700, ms: 2600,
              tick: function (k) { CONSTRUCT.forEach(function (n) { setPose(n, 1, 0, 1 - ease(k)); }); } }
          ];
        }
      }
    };

    /* Mol*'s temporal anti-aliasing adds a jittered sample per frame while the camera holds still. During an animation
       the scene changes every frame, so it restarts every frame and only its first, unjittered sample is ever shown: the
       same picture as without it, at a third more GPU time. Off while things move, back on the moment they stop, so the
       still image refines exactly as before. */
    var calm = 0, aaMode = null;
    function hush(on) {
      var c3 = plugin.canvas3d, ms = c3.props.multiSample;
      if (on && calm++ === 0) {
        aaMode = ms.mode;
        c3.setProps({ multiSample: Object.assign({}, ms, { mode: 'off' }) });
      } else if (!on && --calm === 0) {
        c3.setProps({ multiSample: Object.assign({}, ms, { mode: aaMode }) });
        c3.requestDraw();
      }
    }

    /** run tracks {delay, ms, start(), tick(k)} in one frame loop; false if overtaken */
    async function run(tracks, token) {
      hush(true);
      try { return await frames(tracks, token); } finally { hush(false); }
    }

    async function frames(tracks, token) {
      var t0 = performance.now(), done = false;
      while (!done) {
        if (token !== state.token) return false;
        var now = performance.now() - t0, jobs = [];
        done = true;
        tracks.forEach(function (t) {
          var delay = t.delay || 0;
          if (now < delay) { done = false; return; }
          if (t.start) { t.start(); t.start = null; }
          if (t.finished) return;
          var k = t.ms ? Math.min((now - delay) / t.ms, 1) : 1;
          if (t.tick) jobs.push(t.tick(k));
          if (k < 1) done = false; else t.finished = true;
        });
        await Promise.all(jobs);
        plugin.canvas3d.requestDraw();
        if (!done) await nextFrame();
      }
      return true;
    }

    /** show exactly what step i wants; with `animate`, the step's plan moves its parts in.
        `ready` runs once the scene is set up, just before anything moves: the camera starts there. */
    async function apply(i, animate, ready) {
      if (!plugin) return;
      var token = ++state.token;
      var want = wanted(i);
      await Promise.all(needed(i).map(ensure));
      if (token !== state.token) return;
      await fade(Object.keys(cur), 1);                                 // a fade cut short
      var plan = animate && !stillMotion() ? PLANS[STEPS[i].key] : null;
      var enter = plan ? plan.enter : [], move = plan ? plan.move : [], keep = plan && plan.keep || [];
      if (enter.length) await fade(enter, 0);                          // invisible first, then they fade in
      if (token !== state.token) return;
      Object.keys(built).forEach(function (n) {
        setVisible(n, (want.indexOf(n) >= 0 && enter.indexOf(n) < 0) || keep.indexOf(n) >= 0);
        if (enter.indexOf(n) < 0 && move.indexOf(n) < 0) setRest(n, i);
      });
      if (ready) ready();
      state.playing = !!plan;
      try {
        if (plan && !(await run(plan.tracks(), token))) return;
      } finally { state.playing = false; }
      Object.keys(built).forEach(function (n) {
        setVisible(n, want.indexOf(n) >= 0);                           // the pieces the chain took over from go
        setRest(n, i);
      });
    }

    /** the clash view: the block swings shut with PD-L1 on its domain, binder and PD-L1 turn into spheres; and back */
    async function toggleClash() {
      state.clash = !state.clash;
      paintStep(state.step);
      if (!plugin) return;
      var i = state.step, from = state.clash ? 1 : 0, to = 1 - from;
      var token = ++state.token, want = wanted(i);
      host.classList.add('is-busy');
      await Promise.all(want.map(ensure));
      if (token !== state.token) return;
      await fade(Object.keys(cur), 1);
      Object.keys(built).forEach(function (n) {
        setVisible(n, want.indexOf(n) >= 0);
        if (CONSTRUCT.indexOf(n) >= 0) setPose(n, from, 0, 1);
      });
      frame({ frame: state.clash ? 'clash' : STEPS[i].frame }, 1100);
      state.playing = !stillMotion();
      try {
        if (stillMotion()) { CONSTRUCT.forEach(function (n) { setPose(n, to, 0, 1); }); }
        else await run([{ ms: 1700, tick: function (k) {
          CONSTRUCT.forEach(function (n) { setPose(n, from + (to - from) * ease(k), 0, 1); });
        } }], token);
      } finally { state.playing = false; host.classList.remove('is-busy'); }
      if (token !== state.token) return;
      CONSTRUCT.forEach(function (n) { setRest(n, i); });
    }

    function go(i, first) {
      if (i === state.step) return state.busy;
      var from = state.step;
      state.step = i;
      state.clash = false;
      paintStep(i);
      if (!plugin) return state.busy;
      host.classList.add('is-busy');
      state.busy = (async function () {
        var animate = !first && from >= 0 && i > from;
        state.revealing = animate && !stillMotion();   // parts built now stay hidden until they fly in
        var needsBuild = needed(i).some(function (n) { return !built[n]; });
        if (needsBuild && !first) veilBuild(true);
        try {
          await Promise.all(needed(i).map(ensure));
        } finally { veilBuild(false); }
        if (state.step !== i) return;
        await apply(i, animate, function () { frame(STEPS[i], first ? 0 : 1100); });
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
    Array.prototype.forEach.call(document.querySelectorAll('[data-latch-viewer]'), function (h) {
      if (!h.__lv) { h.__lv = true; mount(h); }
    });
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
  else init();
})();
