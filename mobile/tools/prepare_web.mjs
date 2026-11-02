// Copies RDKit's WebAssembly build into www/rdkit, so the app works offline.   npm run prepare-web
import { copyFileSync, mkdirSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const here = dirname(dirname(fileURLToPath(import.meta.url)));
const out = join(here, 'www', 'rdkit');
mkdirSync(out, { recursive: true });
for (const f of ['RDKit_minimal.js', 'RDKit_minimal.wasm'])
  copyFileSync(join(here, 'node_modules', '@rdkit', 'rdkit', 'dist', f), join(out, f));
console.log('www/rdkit ready');
