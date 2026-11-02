// Runs the Android build (android/gradlew <task>) with JAVA_HOME / ANDROID_HOME found on this PC, then copies the APK to
// dist/Kekulex-<version>[-debug].apk.   node tools/gradle.mjs assembleDebug | assembleRelease
import { execFileSync } from 'node:child_process';
import { copyFileSync, existsSync, mkdirSync, readdirSync, readFileSync, writeFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const here = dirname(dirname(fileURLToPath(import.meta.url)));
const task = process.argv[2] || 'assembleDebug';
const env = { ...process.env };

if (!env.JAVA_HOME) {
  const base = 'C:\\Program Files\\Microsoft';
  const jdk = existsSync(base) && readdirSync(base).filter(d => /^jdk-(2[1-9])/.test(d)).sort().pop();
  if (jdk) env.JAVA_HOME = join(base, jdk);
}
if (!env.ANDROID_HOME) env.ANDROID_HOME = join(env.LOCALAPPDATA || '', 'Android', 'Sdk');
if (!env.JAVA_HOME || !existsSync(env.ANDROID_HOME)) {
  console.error('Need JDK 21+ (JAVA_HOME) and the Android SDK (ANDROID_HOME); see mobile/README.txt'); process.exit(1);
}
const android = join(here, 'android');
writeFileSync(join(android, 'local.properties'), `sdk.dir=${env.ANDROID_HOME.replace(/\\/g, '\\\\')}\n`);

const win = process.platform === 'win32';
execFileSync(win ? 'cmd.exe' : './gradlew', win ? ['/c', join(android, 'gradlew.bat'), task, '--console=plain'] : [task, '--console=plain'],
  { cwd: android, env, stdio: 'inherit' });

const kind = task.toLowerCase().includes('release') ? 'release' : 'debug';
const dir = join(android, 'app', 'build', 'outputs', 'apk', kind);
const apk = existsSync(dir) && readdirSync(dir).find(f => f.endsWith('.apk'));
if (apk) {
  const ver = JSON.parse(readFileSync(join(here, 'package.json'), 'utf8')).version;
  mkdirSync(join(here, 'dist'), { recursive: true });
  const out = join(here, 'dist', `Kekulex-${ver}${kind === 'debug' ? '-debug' : ''}.apk`);
  copyFileSync(join(dir, apk), out);
  console.log(`APK: ${out}`);
}
