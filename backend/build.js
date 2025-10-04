import * as ESBuild from 'esbuild';

console.log(await ESBuild.build({
    entryPoints: ['ts_src/index.ts'],
    bundle: true,
    platform: 'node',
    treeShaking: true,
    target: [ 'node20' ],
    minify: true,
    format: 'cjs',
    outfile: 'build/index.js',
    external: ['node:*']
}))
