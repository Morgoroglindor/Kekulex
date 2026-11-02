// Serves www/ on http://localhost:8766 to try the app in a PC browser (npm run serve). Optional port argument.
import { createReadStream, statSync } from 'node:fs';
import { createServer } from 'node:http';
import { dirname, extname, join, normalize } from 'node:path';
import { fileURLToPath } from 'node:url';

const www = join(dirname(dirname(fileURLToPath(import.meta.url))), 'www');
const port = Number(process.argv[2] || 8766);
const TYPES = { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.css': 'text/css',
  '.json': 'application/json', '.wasm': 'application/wasm', '.svg': 'image/svg+xml', '.png': 'image/png' };

createServer((req, res) => {
  const p = normalize(decodeURIComponent(new URL(req.url, 'http://x').pathname)).replace(/^([/\\])+/, '');
  let f = join(www, p);
  if (!f.startsWith(www)) { res.writeHead(403).end(); return; }
  try { if (statSync(f).isDirectory()) f = join(f, 'index.html'); statSync(f); }
  catch { res.writeHead(404).end('not found'); return; }
  res.writeHead(200, { 'Content-Type': TYPES[extname(f).toLowerCase()] || 'application/octet-stream', 'Cache-Control': 'no-cache' });
  createReadStream(f).pipe(res);
}).listen(port, '127.0.0.1', () => console.log(`Kekulex on http://localhost:${port}`));
