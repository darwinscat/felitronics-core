<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### tools · wasm — `fc_probe` ships as an ES module too, for a module worker

**`tools/wasm/build.sh` now builds `fcprobe.web.mjs`**: the `fcprobe.web.js` line plus `-sEXPORT_ES6=1`, as
`fcmaster.web.mjs` has been built all along. A module worker cannot load the classic glue; it can `import` this
one, whose default export is `createFcProbe`. `fcprobe.web.js` stays, because `probe.html` loads it with a plain
`<script>`.

**One wasm under both glues.** Both web builds write `fcprobe.web.wasm`, so the node module is built first and
the script compares the web module's hash with `fcprobe.node.wasm` after EACH web build, not once after both,
since a single check at the end would only see the second. The `.mjs` glue goes through `check-no-threads.mjs`
and appears in the size table. Nothing in `fc_probe.cpp` changed, and the `fcprobe.web.js` line is the one it
was.
