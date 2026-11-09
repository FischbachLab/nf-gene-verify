"""Validated inputs shared by the EC2 pipeline and independent plotter."""
import csv
import gzip
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

from Bio import SeqIO, Phylo

FASTA_SUFFIXES = ('.fna', '.fasta', '.fa', '.fas')
MISSING = {'', 'NA', 'N/A', 'NAN', 'NONE', '?', '.', '-'}


def sample_id(value):
    name = Path(str(value).strip()).name
    if name.lower().endswith('.gz'):
        name = name[:-3]
    for suffix in FASTA_SUFFIXES:
        if name.lower().endswith(suffix):
            name = name[:-len(suffix)]
            break
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,63}', name):
        raise ValueError(f'Invalid sample identifier {value!r}. Use 1-64 letters, digits, underscores, dots or hyphens.')
    return name


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def save_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')
    temporary.replace(path)


def read_tsv(path, required):
    with Path(path).open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise ValueError(f'Duplicated column in {path}')
        if not set(required).issubset(fields):
            raise ValueError(f'{path} must have tab-separated columns: {", ".join(required)}')
        rows = []
        for line, row in enumerate(reader, 2):
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f'Wrong number of fields in {path}, line {line}')
            rows.append({k: v.strip() for k, v in row.items()})
    return fields, rows


def read_genes(path, samples, allow_missing=False):
    _, rows = read_tsv(path, ['sample', 'gene', 'present'])
    calls, genes = {}, []
    for row in rows:
        sample, gene, value = sample_id(row['sample']), row['gene'], row['present']
        if sample not in samples:
            raise ValueError(f'Gene table has an unknown sample: {row["sample"]}')
        if not gene or any(ord(c) < 32 for c in gene):
            raise ValueError('Gene labels must be nonempty and contain no control characters')
        key = sample, gene
        if key in calls:
            raise ValueError(f'Duplicated sample/gene call: {sample}, {gene}')
        if value in ('0', '1'):
            calls[key] = int(value)
        elif value.upper() in MISSING:
            calls[key] = None
        else:
            raise ValueError(f'Invalid present value {value!r} for {sample}/{gene}; use 0, 1, or NA')
        if gene not in genes:
            genes.append(gene)
    if not genes:
        raise ValueError('Gene table has no marker rows')
    absent_rows = [(s, g) for s in samples for g in genes if (s, g) not in calls]
    if absent_rows and not allow_missing:
        examples = ', '.join(f'{s}/{g}' for s, g in absent_rows[:5])
        raise ValueError(f'{len(absent_rows)} sample/marker rows are missing ({examples}). Add explicit calls or use --allow-missing to display them as unknown.')
    return calls, genes, len(absent_rows)


def normalize_st(value):
    if value.upper() in MISSING:
        return ''
    value = re.sub(r'^ST\s*', '', value, flags=re.IGNORECASE)
    if re.fullmatch(r'\d+\.0', value):
        value = value[:-2]
    if not re.fullmatch(r'\d+', value):
        raise ValueError(f'Invalid sequence type: {value!r}; use an ST number or NA')
    return str(int(value))


def read_metadata(path, samples):
    metadata = {sample: {'ST': ''} for sample in samples}
    if path is None:
        return metadata, ['ST']
    fields, rows = read_tsv(path, [])
    id_columns = ('sample','genome','strainID')
    key = next((column for column in id_columns if column in fields),None)
    if key is None:
        raise ValueError('Metadata needs a sample, genome or strainID column')
    columns = [f for f in fields if f not in id_columns]
    if 'ST' not in columns:
        columns.insert(0, 'ST')
    seen = set()
    for row in rows:
        sample = sample_id(row[key])
        if sample not in samples:
            raise ValueError(f'Metadata has an unknown sample: {row[key]}')
        if sample in seen:
            raise ValueError(f'Duplicated metadata sample: {sample}')
        seen.add(sample)
        metadata[sample] = {column: row.get(column, '') for column in columns}
        metadata[sample]['ST'] = normalize_st(metadata[sample]['ST'])
    return metadata, columns


def write_metadata(path, metadata, columns):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['sample'] + columns, delimiter='\t')
        writer.writeheader()
        for sample, row in metadata.items():
            writer.writerow(dict(sample=sample, **{c: row.get(c, '') for c in columns}))


def discover_genomes(directory):
    directory = Path(directory).resolve()
    if not directory.is_dir():
        raise ValueError(f'Genome folder does not exist: {directory}')
    genomes = {}
    for path in sorted(directory.iterdir()):
        lower = path.name.lower()
        if path.is_file() and any(lower.endswith(ext) or lower.endswith(ext+'.gz') for ext in FASTA_SUFFIXES):
            sample = sample_id(path.name)
            if sample in genomes:
                raise ValueError(f'Two genome files produce the same identifier: {sample}')
            genomes[sample] = path
    if not genomes:
        raise ValueError('No .fna/.fasta/.fa/.fas genome files found directly in the supplied folder')
    return genomes


def read_fasta(path):
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    with opener(path, 'rt') as stream:
        records = list(SeqIO.parse(stream, 'fasta'))
    if not records or len({r.id for r in records}) != len(records):
        raise ValueError(f'Empty FASTA or duplicated contig identifier: {path}')
    for record in records:
        sequence = str(record.seq).upper()
        if not sequence or set(sequence) - set('ACGTRYSWKMBDHVN'):
            raise ValueError(f'Empty contig or invalid DNA alphabet in {path}: {record.id}')
        record.seq = record.seq.upper()
    return records


def assembly_stats(sample, path):
    records = read_fasta(path)
    lengths = sorted([len(r.seq) for r in records], reverse=True)
    counts = Counter(''.join(str(r.seq) for r in records))
    canonical = sum(counts[b] for b in 'ACGT')
    if not canonical:
        raise ValueError(f'Genome contains no A/C/G/T bases: {path}')
    total, cumulative, n50 = sum(lengths), 0, None
    for length in lengths:
        cumulative += length
        if cumulative >= total/2:
            n50 = length
            break
    return dict(sample=sample, source=str(path), sha256=sha256(path), contigs=len(lengths),
                bases=total, N50=n50, GC_percent=round(100*(counts['G']+counts['C'])/canonical, 5),
                ambiguous_bases=total-canonical)


def read_tree(path, samples=None):
    tree = Phylo.read(path, 'newick')
    tips = tree.get_terminals()
    if len(tips) < 2:
        raise ValueError('Tree must have at least two tips')
    for tip in tips:
        if not tip.name:
            raise ValueError('Tree contains an unnamed tip')
        tip.name = sample_id(tip.name)
    labels = [tip.name for tip in tips]
    if len(set(labels)) != len(labels):
        raise ValueError('Tree has duplicated sample identifiers after filename normalization')
    for clade in tree.find_clades():
        if clade.branch_length is not None:
            import math
            if not math.isfinite(clade.branch_length) or clade.branch_length < 0:
                raise ValueError('Tree branch lengths must be finite and nonnegative')
    if samples is not None and set(labels) != set(samples):
        missing, extra = sorted(set(samples)-set(labels)), sorted(set(labels)-set(samples))
        raise ValueError(f'Tree/sample mismatch. Missing tips: {missing}; extra tips: {extra}')
    return tree
