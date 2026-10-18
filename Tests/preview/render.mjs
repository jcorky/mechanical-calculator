// Renders preview PNGs of the exported keyboard with three.js in headless Chromium (Playwright).
//
//   build/mcalc_export_obj <dir>
//   node Tests/preview/render.mjs <dir> Docs/images [view ...]
//
// Views: hero (3/4 view of the board on a desk), closeup (a few keys and the OLED, one cap pulled),
// exploded (one switch and its keycap pulled apart), profile (side view of the sculpt).
// three.js is fetched once from jsDelivr into the OS temp directory. Playwright is taken from
// PLAYWRIGHT_PATH or the global install at /opt/node22/lib/node_modules/playwright.
import { createRequire } from 'node:module';
import { execFileSync } from 'node:child_process';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';

const require = createRequire(import.meta.url);
const { chromium } = require(process.env.PLAYWRIGHT_PATH || '/opt/node22/lib/node_modules/playwright');

const [inDir = '.', outDir = 'Docs/images', ...wanted] = process.argv.slice(2);
const views = wanted.length ? wanted : ['hero', 'closeup', 'exploded', 'profile'];

const THREE_VERSION = '0.160.0';
const cacheDir = path.join(os.tmpdir(), `mcalc-three-${THREE_VERSION}`);
const libs = {
	'three.module.js': 'build/three.module.js',
	'RoomEnvironment.js': 'examples/jsm/environments/RoomEnvironment.js',
};

function fetchLibs() {
	fs.mkdirSync(cacheDir, { recursive: true });
	for (const [name, rel] of Object.entries(libs)) {
		const file = path.join(cacheDir, name);
		if (!fs.existsSync(file) || fs.statSync(file).size < 1000) {
			execFileSync('curl', ['-sSfL', '-o', file, `https://cdn.jsdelivr.net/npm/three@${THREE_VERSION}/${rel}`]);
		}
	}
}

// OBJ (as written by mgeo::WriteObj: each group's vertices directly precede its faces) to one
// binary blob: a JSON header followed by per-group Float32 positions/normals/uvs and Uint32 indices.
function packObj(file) {
	const text = fs.readFileSync(file, 'utf8');
	const P = [], N = [], T = [];
	const groups = [];
	let cur = null;
	let pos = 0;
	while (pos < text.length) {
		let end = text.indexOf('\n', pos);
		if (end < 0) end = text.length;
		const c0 = text.charCodeAt(pos), c1 = text.charCodeAt(pos + 1);
		if (c0 === 118) {
			const a = text.slice(pos, end).split(' ');
			if (c1 === 32) P.push(+a[1], +a[2], +a[3]);
			else if (c1 === 110) N.push(+a[1], +a[2], +a[3]);
			else if (c1 === 116) T.push(+a[1], +a[2]);
		} else if (c0 === 102) {
			const a = text.slice(pos, end).split(' ');
			cur.idx.push(parseInt(a[1], 10) - 1, parseInt(a[2], 10) - 1, parseInt(a[3], 10) - 1);
		} else if (c0 === 103) {
			cur = { name: text.slice(pos + 2, end), idx: [] };
			groups.push(cur);
		}
		pos = end + 1;
	}
	const header = [];
	const chunks = [];
	let offset = 0;
	for (const g of groups) {
		if (!g.idx.length) continue;
		let lo = Infinity, hi = -Infinity;
		for (const i of g.idx) { if (i < lo) lo = i; if (i > hi) hi = i; }
		const count = hi - lo + 1;
		const parts = [
			new Float32Array(P.slice(lo * 3, (hi + 1) * 3)),
			new Float32Array(N.slice(lo * 3, (hi + 1) * 3)),
			new Float32Array(T.slice(lo * 2, (hi + 1) * 2)),
			Uint32Array.from(g.idx, (i) => i - lo),
		];
		header.push({ name: g.name, vertices: count, indices: g.idx.length, offset });
		for (const p of parts) { chunks.push(Buffer.from(p.buffer)); offset += p.byteLength; }
	}
	const head = Buffer.from(JSON.stringify(header));
	const pad = Buffer.alloc((4 - ((4 + head.length) % 4)) % 4, 32);
	const len = Buffer.alloc(4);
	len.writeUInt32LE(head.length + pad.length);
	return Buffer.concat([len, head, pad, ...chunks]);
}

// Runs in the page.
async function pageMain(THREE, RoomEnvironment) {
	const scene = await (await fetch('/data/scene.json')).json();
	async function load(name) {
		const buf = await (await fetch(`/data/${name}.bin`)).arrayBuffer();
		const headLen = new DataView(buf).getUint32(0, true);
		const header = JSON.parse(new TextDecoder().decode(new Uint8Array(buf, 4, headLen)));
		const base = 4 + headLen;
		return header.map((g) => {
			const geo = new THREE.BufferGeometry();
			let o = base + g.offset;
			const take = (Ctor, n) => { const a = new Ctor(buf, o, n); o += a.byteLength; return a; };
			geo.setAttribute('position', new THREE.BufferAttribute(take(Float32Array, g.vertices * 3), 3));
			geo.setAttribute('normal', new THREE.BufferAttribute(take(Float32Array, g.vertices * 3), 3));
			geo.setAttribute('uv', new THREE.BufferAttribute(take(Float32Array, g.vertices * 2), 2));
			geo.setIndex(new THREE.BufferAttribute(take(Uint32Array, g.indices), 1));
			return { name: g.name, geo };
		});
	}

	const capColors = {
		alpha: [0xe8e2d2, '#2b2d31'], numpad: [0xe8e2d2, '#2b2d31'], modifier: [0x3f4247, '#e8e2d2'],
		navigation: [0x3f4247, '#e8e2d2'], function: [0x3f4247, '#e8e2d2'], space: [0x3f4247, '#e8e2d2'], accent: [0xd8642e, '#fbf3e6'],
	};
	const phys = (o) => new THREE.MeshPhysicalMaterial(o);
	const materials = {
		aluminium: phys({ color: 0xaab0b8, metalness: 0.85, roughness: 0.42, envMapIntensity: 1.4 }),
		badge: phys({ color: 0xc8a25a, metalness: 1.0, roughness: 0.22 }),
		rubber: phys({ color: 0x141414, roughness: 0.95 }),
		plate: phys({ color: 0xb9a36a, metalness: 1.0, roughness: 0.32 }),
		pcb: phys({ color: 0x1e6a3b, roughness: 0.45, clearcoat: 0.6, clearcoatRoughness: 0.4 }),
		pcb_parts: phys({ color: 0x1a1a1a, roughness: 0.4, metalness: 0.3 }),
		switch_housing_bottom: phys({ color: 0x1c1c1e, roughness: 0.55 }),
		switch_housing_top: phys({ color: 0xdfe6ee, roughness: 0.08, transparent: true, opacity: 0.42, depthWrite: false }),
		switch_stem: phys({ color: 0x7a4a2c, roughness: 0.45 }),
		switch_spring: phys({ color: 0xd6b25e, metalness: 1.0, roughness: 0.25 }),
		switch_contacts: phys({ color: 0xe8b070, metalness: 0.9, roughness: 0.3, envMapIntensity: 1.6 }),
		switch_pins: phys({ color: 0xe8b070, metalness: 0.9, roughness: 0.3, envMapIntensity: 1.6 }),
		stab_housing: phys({ color: 0x202022, roughness: 0.5 }),
		stab_stem: phys({ color: 0x2c2c2e, roughness: 0.5 }),
		stab_wire: phys({ color: 0xc9c9c9, metalness: 1.0, roughness: 0.25 }),
		oled_bezel: phys({ color: 0x101012, roughness: 0.3, metalness: 0.2 }),
	};

	function legendTexture(key, ink) {
		const aspect = key.topW / key.topD;
		const H = 256, W = Math.round(Math.min(4, Math.max(1, aspect)) * H);
		const c = document.createElement('canvas');
		c.width = W; c.height = H;
		const g = c.getContext('2d');
		g.fillStyle = ink;
		g.textBaseline = 'top';
		const pad = 0.16 * H;
		const font = (px, weight = 500) => `${weight} ${px}px Inter, sans-serif`;
		const main = key.legend;
		if (key.role === 'space' || !main) {
			// blank
		} else if (key.shift) {
			g.font = font(0.30 * H);
			g.fillText(key.shift, pad, pad);
			g.fillText(main, pad, 0.52 * H);
		} else if (main.length === 1 || /^[←-↓]$/.test(main)) {
			g.font = font(0.40 * H);
			if (key.role === 'navigation' || key.role === 'numpad') {
				g.textAlign = 'center';
				g.fillText(main, W / 2, 0.28 * H);
			} else {
				g.fillText(main, pad, pad);
			}
		} else {
			g.font = font(0.19 * H, 600);
			const words = main.split(' ');
			words.forEach((w, i) => g.fillText(w, pad, pad + i * 0.24 * H));
		}
		const t = new THREE.CanvasTexture(c);
		t.colorSpace = THREE.SRGBColorSpace;
		t.anisotropy = 8;
		t.wrapS = t.wrapT = THREE.ClampToEdgeWrapping;
		return t;
	}

	// Legend printed where uv.x >= 0 (the cap's top); everything else keeps the cap colour.
	function capMaterial(color, tex) {
		const m = phys({ color, roughness: 0.62, sheen: 0.2 });
		m.onBeforeCompile = (sh) => {
			sh.uniforms.legendMap = { value: tex };
			sh.fragmentShader = sh.fragmentShader
				.replace('#include <common>', '#include <common>\nuniform sampler2D legendMap;')
				.replace('#include <map_fragment>', `#include <map_fragment>
					if (vLegendUv.x >= 0.0) { vec4 lg = texture2D(legendMap, vLegendUv); diffuseColor.rgb = mix(diffuseColor.rgb, lg.rgb, lg.a); }`)
				.replace('#include <common>', '#include <common>\nvarying vec2 vLegendUv;');
			sh.vertexShader = sh.vertexShader
				.replace('#include <common>', '#include <common>\nvarying vec2 vLegendUv;')
				.replace('#include <uv_vertex>', '#include <uv_vertex>\nvLegendUv = uv;');
		};
		m.customProgramCacheKey = () => 'cap';
		return m;
	}

	function frontDecal(key, ink) {
		const w = key.topW + 1.8, h = 3.0;
		const c = document.createElement('canvas');
		c.width = Math.round(64 * w / h); c.height = 64;
		const g = c.getContext('2d');
		g.fillStyle = ink;
		g.font = '600 40px Inter, sans-serif';
		g.textAlign = 'center';
		g.textBaseline = 'middle';
		g.fillText(key.front, c.width / 2, 34);
		const t = new THREE.CanvasTexture(c);
		t.colorSpace = THREE.SRGBColorSpace;
		const m = new THREE.MeshStandardMaterial({ map: t, transparent: true, roughness: 0.6, depthWrite: false, polygonOffset: true, polygonOffsetFactor: -2 });
		const mesh = new THREE.Mesh(new THREE.PlaneGeometry(w, h), m);
		mesh.rotation.x = Math.PI / 2 - key.frontTilt;
		mesh.position.set(key.x, key.y + key.frontY - 0.12, scene.stemTopZ + key.frontZ);
		return mesh;
	}

	function displayTexture() {
		const W = 1024, H = Math.round(W * scene.display.h / scene.display.w);
		const c = document.createElement('canvas');
		c.width = W; c.height = H;
		const g = c.getContext('2d');
		g.fillStyle = '#000'; g.fillRect(0, 0, W, H);
		g.fillStyle = '#62e6ff';
		g.font = `500 ${Math.round(H * 0.17)}px "DejaVu Sans Mono", monospace`;
		g.textBaseline = 'top';
		g.fillText('RAD  FIX4   LIB ORBITS', W * 0.03, H * 0.08);
		g.fillRect(W * 0.03, H * 0.31, W * 0.94, 2);
		g.font = `500 ${Math.round(H * 0.22)}px "DejaVu Sans Mono", monospace`;
		g.fillText('sqrt(2)*sin(pi/4)+2^3', W * 0.03, H * 0.40);
		g.textAlign = 'right';
		g.font = `600 ${Math.round(H * 0.26)}px "DejaVu Sans Mono", monospace`;
		g.fillText('= 9.0000', W * 0.97, H * 0.68);
		const t = new THREE.CanvasTexture(c);
		t.colorSpace = THREE.SRGBColorSpace;
		t.flipY = false; // UV (0,0) is the screen's top-left corner
		return t;
	}

	const renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
	renderer.setPixelRatio(1);
	renderer.toneMapping = THREE.ACESFilmicToneMapping;
	renderer.toneMappingExposure = 1.0;
	renderer.outputColorSpace = THREE.SRGBColorSpace;
	renderer.shadowMap.enabled = true;
	renderer.shadowMap.type = THREE.PCFSoftShadowMap;
	document.body.appendChild(renderer.domElement);
	const pmrem = new THREE.PMREMGenerator(renderer);
	const envMap = pmrem.fromScene(new RoomEnvironment(), 0.04).texture;

	const keysById = Object.fromEntries(scene.keys.map((k) => [k.id, k]));
	function buildModel(parts, { decals = true, pulled = {} } = {}) {
		const model = new THREE.Group();
		for (const { name, geo } of parts) {
			let mat;
			if (name.startsWith('cap:')) {
				const [, role, id] = name.split(':');
				const key = keysById[id] || { id, role, legend: id, shift: '', front: '', topW: 12.6, topD: 13.9 };
				const [color, ink] = capColors[role] || capColors.alpha;
				mat = capMaterial(color, legendTexture(key, ink));
				if (decals && key.front && keysById[id]) {
					const d = frontDecal(key, role === 'alpha' || role === 'numpad' ? '#c4521f' : '#f0a070');
					d.position.z += pulled[id] || 0;
					model.add(d);
				}
			} else if (name === 'oled_glass') {
				mat = phys({ color: 0x050607, roughness: 0.05, emissive: 0xffffff, emissiveMap: displayTexture(), emissiveIntensity: 1.6, clearcoat: 1.0 });
			} else {
				mat = materials[name] || phys({ color: 0xff00ff });
			}
			const mesh = new THREE.Mesh(geo, mat);
			mesh.castShadow = name !== 'switch_housing_top';
			mesh.receiveShadow = true;
			if (name.startsWith('cap:')) mesh.position.z += pulled[name.split(':')[2]] || 0;
			model.add(mesh);
		}
		return model;
	}

	// Board frame (mm, z up, y away from the typist) to three.js (y up, z toward the viewer).
	function place(model, rest) {
		const tilt = new THREE.Group();
		if (rest) {
			tilt.rotation.x = scene.restTiltRad;
			tilt.position.z = scene.restLiftMm;
		}
		tilt.add(model);
		const world = new THREE.Group();
		world.rotation.x = -Math.PI / 2;
		world.add(tilt);
		return world;
	}

	function stage(width, height, background, shadowOnly = false) {
		renderer.setSize(width, height);
		const s = new THREE.Scene();
		s.background = new THREE.Color(background);
		s.fog = new THREE.Fog(background, 900, 2600);
		s.environment = envMap;
		const deskMat = shadowOnly ? new THREE.ShadowMaterial({ opacity: 0.35 }) : new THREE.MeshStandardMaterial({ color: 0x55585e, roughness: 0.8 });
		const desk = new THREE.Mesh(new THREE.PlaneGeometry(8000, 8000), deskMat);
		desk.rotation.x = -Math.PI / 2;
		desk.receiveShadow = true;
		s.add(desk);
		const sun = new THREE.DirectionalLight(0xfff4e6, 2.8);
		sun.position.set(-260, 520, 240);
		sun.castShadow = true;
		sun.shadow.mapSize.set(4096, 4096);
		Object.assign(sun.shadow.camera, { left: -320, right: 320, top: 240, bottom: -240, near: 10, far: 1500 });
		sun.shadow.bias = -0.0004;
		sun.shadow.normalBias = 0.3;
		s.add(sun);
		s.add(new THREE.HemisphereLight(0xdfe8ff, 0x30302c, 0.35));
		return { s, sun };
	}

	const cache = {};
	const parts = async (name) => (cache[name] ||= await load(name));

	window.renderView = async (view) => {
		let width = 1800, height = 1100, cam, s;
		if (view === 'exploded') {
			width = 1100; height = 1400;
			({ s } = stage(width, height, 0x6a6e75, true));
			const model = place(buildModel(await parts('switch_exploded'), { decals: false }), false);
			model.position.y = 9;
			s.add(model);
			cam = new THREE.PerspectiveCamera(31, width / height, 1, 2000);
			cam.position.set(112, 52, 150);
			cam.lookAt(0, 44, 0);
		} else {
			({ s } = stage(width, height, 0x55585e));
			const pulled = view === 'closeup' ? { Num8: 15 } : {};
			s.add(place(buildModel(await parts('board'), { pulled }), true));
			if (view === 'hero') {
				cam = new THREE.PerspectiveCamera(22, width / height, 10, 5000);
				cam.position.set(-360, 440, 640);
				cam.lookAt(12, 4, -6);
			} else if (view === 'closeup') {
				width = 1600; height = 1000;
				renderer.setSize(width, height);
				cam = new THREE.PerspectiveCamera(30, width / height, 5, 3000);
				cam.position.set(92, 100, 128);
				cam.lookAt(165, 18, -12);
			} else if (view.startsWith('cam:')) {
				// cam:x,y,z:tx,ty,tz:fov in three.js coordinates (mm, y up, z toward the typist)
				const [, p, t, fov] = view.split(':');
				const [px, py, pz] = p.split(',').map(Number), [tx, ty, tz] = t.split(',').map(Number);
				cam = new THREE.PerspectiveCamera(Number(fov) || 30, width / height, 1, 5000);
				cam.position.set(px, py, pz);
				cam.lookAt(tx, ty, tz);
			} else if (view === 'back') {
				cam = new THREE.PerspectiveCamera(24, width / height, 10, 5000);
				cam.position.set(140, 160, -520);
				cam.lookAt(40, 10, 0);
			} else {
				width = 1800; height = 800;
				renderer.setSize(width, height);
				cam = new THREE.PerspectiveCamera(8, width / height, 10, 6000);
				cam.position.set(-1100, 45, -6);
				cam.lookAt(-150, 24, -6);
			}
		}
		renderer.render(s, cam);
		return renderer.domElement.toDataURL('image/png');
	};
	window.ready = true;
}

fetchLibs();
const blobs = {};
const need = new Set(views.map((v) => (v === 'exploded' ? 'switch_exploded' : 'board')));
for (const name of need) {
	const t0 = Date.now();
	blobs[name] = packObj(path.join(inDir, `${name}.obj`));
	console.log(`packed ${name}.obj (${(blobs[name].length / 1e6).toFixed(1)} MB) in ${Date.now() - t0} ms`);
}
const page_html = `<!doctype html><html><head><meta charset="utf-8">
<script type="importmap">{"imports": {"three": "/three.module.js"}}</script>
<style>html,body{margin:0;background:#000}</style></head><body>
<script type="module">
import * as THREE from 'three';
import { RoomEnvironment } from '/RoomEnvironment.js';
(${pageMain.toString()})(THREE, RoomEnvironment).catch((e) => { window.failed = String(e && e.stack || e); });
</script></body></html>`;

const browser = await chromium.launch({ args: ['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--ignore-gpu-blocklist'] });
const page = await browser.newPage({ viewport: { width: 1800, height: 1400 } });
page.on('console', (m) => { if (m.type() === 'error') console.log('page:', m.text()); });
await page.route('http://preview.local/**', async (route) => {
	const url = new URL(route.request().url());
	const p = url.pathname;
	if (p === '/') return route.fulfill({ contentType: 'text/html', body: page_html });
	if (p === '/three.module.js' || p === '/RoomEnvironment.js') {
		return route.fulfill({ contentType: 'text/javascript', body: fs.readFileSync(path.join(cacheDir, p.slice(1))) });
	}
	if (p === '/data/scene.json') return route.fulfill({ contentType: 'application/json', body: fs.readFileSync(path.join(inDir, 'scene.json')) });
	const m = p.match(/^\/data\/(.+)\.bin$/);
	if (m && blobs[m[1]]) return route.fulfill({ contentType: 'application/octet-stream', body: blobs[m[1]] });
	return route.fulfill({ status: 404, body: '' });
});
await page.goto('http://preview.local/');
await page.waitForFunction(() => window.ready || window.failed, null, { timeout: 120000 });
const failed = await page.evaluate(() => window.failed);
if (failed) throw new Error(failed);
fs.mkdirSync(outDir, { recursive: true });
for (const view of views) {
	const t0 = Date.now();
	const url = await page.evaluate((v) => window.renderView(v), view);
	const file = path.join(outDir, `${view.startsWith('cam:') ? 'debug-' + views.indexOf(view) : view}.png`);
	fs.writeFileSync(file, Buffer.from(url.split(',')[1], 'base64'));
	console.log(`${file} (${Date.now() - t0} ms)`);
}
await browser.close();
