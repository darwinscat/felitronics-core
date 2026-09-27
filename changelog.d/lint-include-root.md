### tools — the libm audit's satellite run takes `--include-root`: a third library's headers resolve, and are audited as the satellite's own

`tools/lint/check-det-math.mjs --satellite` resolved `#include <felitronics/...>` only in core's and the satellite's
`modules/*/include`. A satellite that builds against a third library — felitronics-mastering-core's session reads its
config with felitronics-toml, whose `<felitronics/toml/Toml.h>`, `Schema.h` and `Embedded.h` live in neither — failed
with CLOSURE violations for headers that exist: rule 4 refusing, rightly, a header it could not find.

`--include-root <dir>`, repeatable and only with `--satellite`, names such a library's include directory, absolute or
relative to the satellite's root. A CI step reads it from the build, as it already reads core's checkout (for
felitronics-mastering-core, `FELITRONICS_MASTERING_TOML_SOURCE_DIR` from the CMake cache, plus `/include`).

- Every file under the root joins the satellite's pass, under every rule — zone, carrier, manifest, closure and the
  rot checks — and is walked whole, as a module's include directory is, not only where the closure reaches. The same
  headers staged under the satellite's `modules/` get the same verdict, clean and with a libm call planted.
- The satellite's lists name such a file by its include spelling in angle brackets — `zone <felitronics/toml/Toml.h>
  <why>` in `det-math-zone.txt`, `<felitronics/toml/Toml.h>  [...]  <disposition>  <why>` in the manifest — and never
  by a path, which differs between a CI build tree and a working checkout.
- Refused before anything is read, exit 2: a root that does not exist or is not a directory, one given twice, one that
  overlaps another root or a directory either repository already scans as its own, a flag with no directory after it,
  and the flag without `--satellite`. A header that two roots provide, or a root and a repository, is refused like one
  two repositories provide, and two roots holding one spelling are refused even where nothing includes it. A root that
  no include of the parity closure resolves through is the new `ROOT-ROT`.
- A command-line option and not a zone-file line, on purpose: the zone file holds claims that read the same on every
  machine, and where a checkout lives is a fact of the build.

A run without the flag behaves as before — core's own pass prints the same bytes, `--report` included — except that the
unresolved-include message now names the flag. `--self-test` gains 19 cases (67 in all) that build a satellite and a
library in a temporary directory and run the same passes the gate runs. On felitronics-mastering-core with
felitronics-toml's include directory the satellite run is clean; without it, it fails with the same seven CLOSURE
violations as before. The long-double and allocation-counter lints resolve no includes and are unchanged.
