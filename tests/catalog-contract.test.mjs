import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { readFile, readdir, stat } from 'node:fs/promises';
import test from 'node:test';
import { gunzipSync } from 'node:zlib';

const ROOT = new URL('../', import.meta.url);
const EXPORTS = new URL('../exports/', import.meta.url);
const EXPECTED = Object.freeze({
  'synthetic-cell-types-2d': Object.freeze({
    cells: 72,
    genes: 8,
    dimension: 2,
    edges: 144,
    vector: false,
  }),
  'synthetic-development-3d': Object.freeze({
    cells: 96,
    genes: 10,
    dimension: 3,
    edges: 187,
    vector: true,
  }),
  'synthetic-trajectory-1d': Object.freeze({
    cells: 48,
    genes: 6,
    dimension: 1,
    edges: 93,
    vector: true,
  }),
});

function exactKeys(value, keys, label) {
  assert.ok(
    value !== null && typeof value === 'object' && !Array.isArray(value),
    `${label} must be an object`
  );
  assert.deepEqual(Object.keys(value).sort(), [...keys].sort(), `${label} keys`);
}

async function readJson(url) {
  return JSON.parse(await readFile(url, 'utf8'));
}

function payload(bytes, relativePath) {
  return relativePath.endsWith('.gz') ? gunzipSync(bytes) : bytes;
}

async function walk(url, prefix = '') {
  const output = [];
  for (const entry of await readdir(url, { withFileTypes: true })) {
    const relative = `${prefix}${entry.name}`;
    const child = new URL(
      `${encodeURIComponent(entry.name)}${entry.isDirectory() ? '/' : ''}`,
      url
    );
    if (entry.isDirectory()) {
      output.push(...await walk(child, `${relative}/`));
    } else if (entry.isFile()) {
      output.push(relative);
    } else {
      assert.fail(`${relative} must be a regular file or directory`);
    }
  }
  return output.sort();
}

test('catalog is the exact three-dataset current contract', async () => {
  const catalog = await readJson(new URL('datasets.json', EXPORTS));
  exactKeys(catalog, ['version', 'default', 'datasets'], 'catalog');
  assert.equal(catalog.version, 1);
  assert.equal(catalog.default, 'synthetic-cell-types-2d');
  assert.deepEqual(
    catalog.datasets.map(dataset => dataset.id),
    Object.keys(EXPECTED).sort()
  );
  for (const entry of catalog.datasets) {
    exactKeys(
      entry,
      ['id', 'path', 'name', 'description', 'n_cells', 'n_genes'],
      `catalog dataset ${entry.id}`
    );
    assert.equal(entry.path, `${entry.id}/`);
    assert.equal(entry.n_cells, EXPECTED[entry.id].cells);
    assert.equal(entry.n_genes, EXPECTED[entry.id].genes);
    assert.ok(entry.name.length > 0);
    assert.ok(entry.description.length > 0);
  }
});

test('landing page references only committed local runtime assets', async () => {
  const index = await readFile(new URL('index.html', ROOT), 'utf8');
  const favicon = await readFile(new URL('favicon.svg', ROOT), 'utf8');
  assert.match(index, /<link rel="icon" href="\.\/favicon\.svg"/);
  assert.match(index, /fetch\('\.\/exports\/datasets\.json'/);
  assert.doesNotMatch(index, /<(?:link|script)\b[^>]+(?:href|src)="https?:/i);
  assert.match(favicon, /^<svg xmlns="http:\/\/www\.w3\.org\/2000\/svg"/);
  assert.ok(favicon.endsWith('\n'));
});

for (const [datasetId, expected] of Object.entries(EXPECTED)) {
  test(`${datasetId} identity and binary axes agree exactly`, async () => {
    const root = new URL(`${datasetId}/`, EXPORTS);
    const identity = await readJson(new URL('dataset_identity.json', root));
    const identityKeys = [
      'version',
      'id',
      'name',
      'description',
      'created_at',
      'cellucid_data_version',
      'stats',
      'embeddings',
      'obs_fields',
      'export_settings',
      'source',
      ...(expected.vector ? ['vector_fields'] : []),
    ];
    exactKeys(identity, identityKeys, `${datasetId} identity`);
    assert.equal(identity.version, 2);
    assert.equal(identity.id, datasetId);
    assert.equal(identity.created_at, '2026-07-27T00:00:00Z');
    assert.equal(identity.cellucid_data_version, '0.9.1');
    assert.equal(identity.stats.n_cells, expected.cells);
    assert.equal(identity.stats.n_genes, expected.genes);
    assert.equal(identity.stats.n_edges, expected.edges);
    assert.deepEqual(identity.embeddings.available_dimensions, [expected.dimension]);
    assert.equal(identity.embeddings.default_dimension, expected.dimension);

    const pointPath = identity.embeddings.files[`${expected.dimension}d`];
    const pointBytes = payload(await readFile(new URL(pointPath, root)), pointPath);
    assert.equal(pointBytes.length, expected.cells * expected.dimension * 4);

    const obs = await readJson(new URL('obs_manifest.json', root));
    const variable = await readJson(new URL('var_manifest.json', root));
    assert.equal(obs._format, 'compact_v1');
    assert.equal(variable._format, 'compact_v1');
    assert.equal(obs.n_points, expected.cells);
    assert.equal(variable.n_points, expected.cells);
    assert.equal(variable.fields.length, expected.genes);

    const connectivity = await readJson(
      new URL('connectivity_manifest.json', root)
    );
    assert.equal(connectivity.format, 'edge_pairs');
    assert.equal(connectivity.n_cells, expected.cells);
    assert.equal(connectivity.n_edges, expected.edges);
    assert.equal(connectivity.index_dtype, 'uint16');
    assert.equal(connectivity.weight_dtype, 'float64');
    const sources = payload(
      await readFile(new URL(connectivity.sourcesPath, root)),
      connectivity.sourcesPath
    );
    const destinations = payload(
      await readFile(new URL(connectivity.destinationsPath, root)),
      connectivity.destinationsPath
    );
    const weights = payload(
      await readFile(new URL(connectivity.weightsPath, root)),
      connectivity.weightsPath
    );
    assert.equal(sources.length, expected.edges * 2);
    assert.equal(destinations.length, expected.edges * 2);
    assert.equal(weights.length, expected.edges * 8);
    let previousSource = -1;
    let previousDestination = -1;
    for (let index = 0; index < expected.edges; index++) {
      const source = sources.readUInt16LE(index * 2);
      const destination = destinations.readUInt16LE(index * 2);
      assert.ok(source < destination && destination < expected.cells);
      assert.ok(
        source > previousSource ||
          (source === previousSource && destination > previousDestination)
      );
      const weight = weights.readDoubleLE(index * 8);
      assert.ok(Number.isFinite(weight) && weight > 0);
      previousSource = source;
      previousDestination = destination;
    }

    if (expected.vector) {
      const vectorField = identity.vector_fields.fields.velocity_umap;
      assert.deepEqual(vectorField.available_dimensions, [expected.dimension]);
      const vectorPath = vectorField.files[`${expected.dimension}d`];
      const vectorBytes = payload(
        await readFile(new URL(vectorPath, root)),
        vectorPath
      );
      assert.equal(vectorBytes.length, expected.cells * expected.dimension * 4);
    }
  });
}

test('SHA256SUMS owns every export file and every digest', async () => {
  const checksumText = await readFile(new URL('SHA256SUMS', EXPORTS), 'utf8');
  assert.ok(checksumText.endsWith('\n'));
  const records = checksumText.trimEnd().split('\n').map(line => {
    assert.match(line, /^[0-9a-f]{64}  [A-Za-z0-9_.\/-]+$/);
    return { digest: line.slice(0, 64), path: line.slice(66) };
  });
  assert.deepEqual(
    records.map(record => record.path),
    records.map(record => record.path).sort()
  );
  const actualPaths = (await walk(EXPORTS)).filter(path => path !== 'SHA256SUMS');
  assert.deepEqual(records.map(record => record.path), actualPaths);
  for (const record of records) {
    const bytes = await readFile(new URL(record.path, EXPORTS));
    assert.equal(createHash('sha256').update(bytes).digest('hex'), record.digest);
  }
});

test('repository text and generated bytes have portable ownership', async () => {
  assert.equal(await readFile(new URL('.gitattributes', ROOT), 'utf8'),
    '* text=auto eol=lf\n*.bin binary\n*.gz binary\n'
  );
  assert.equal((await stat(new URL('exports/', ROOT))).isDirectory(), true);
  const ignore = await readFile(new URL('.gitignore', ROOT), 'utf8');
  assert.match(ignore, /^\.ruff_cache\/$/m);
  assert.match(ignore, /^__pycache__\/$/m);
  assert.match(ignore, /^\*\.py\[cod\]$/m);
});
