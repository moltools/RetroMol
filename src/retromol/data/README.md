# RetroMol rules

This directory contains three YAML configuration files:

| File | Purpose |
| --- | --- |
| [rxn.yml](rxn.yml) | Reaction SMARTS that reverse bonds or modifications and generate candidate fragments. |
| [mxn.yml](mxn.yml) | SMILES that identify whole fragments as monomers, modified units, or tailoring groups. |
| [pmp.yml](pmp.yml) | Maps predictions from BGC annotations/models onto matching-rule names or pseudonyms. It does not control compound disconnections. |

Settings that affect parsing or sequence readout belong at the top level of each rule. `props` is for descriptive metadata such as references, notes, and example structures.

## How compound parsing uses the rules

1. RetroMol standardizes the input and gives its atoms positive isotope tags to track their origin. These internal tags are separate from reaction SMARTS atom-map numbers such as `:1`.
2. For each fragment, it tries the matching rules. Matching requires a structural match and equal atom and bond counts: a matching rule identifies the entire fragment, not just a substructure inside it.
3. A terminal match stops expansion of that fragment. Otherwise, reaction rules generate alternative sets of products, which are processed in turn. Eligible uncontested bulk reactions run before the ordinary alternative-reaction search.
4. The resulting reaction graph records the alternatives. An assembly selects compatible fragments from those alternatives, subject to terminal-boundary protection and atom ownership.

Coverage measures the original input heavy atoms assigned to identified fragments in one compatible assembly. Reaction-added atoms do not increase coverage. A complete assembly is a structural explanation under the rules; it is not, by itself, proof of the biological pathway or the loaded substrates.

## Reaction rules: `rxn.yml`

The file is a YAML list. Each `- name:` begins a rule.

| Field | Type / default | Meaning |
| --- | --- | --- |
| `name` | String; required | Human-readable reaction name, also recorded in reaction provenance. Use a descriptive, distinct name. |
| `smarts` | String; required | RDKit reaction SMARTS in the direction RetroMol should apply it, usually mature structure → precursor fragments. |
| `allowed_in_bulk` | Boolean; `false` | Allows the rule to participate in eager preprocessing when its reactive atoms do not overlap competing bulk reactions. |
| `product_max_matches` | Mapping of SMARTS strings to nonnegative integers; `{}` | Maximum number of matches allowed in each individual product of this rule. If any product exceeds any limit, the entire candidate product set is rejected. |
| `props` | Mapping; `{}` | Optional metadata. Values such as `note`, `reference`, `positive_example`, and `negative_example` document the rule; examples here are not automatically executed as tests. |

For example, this is the saturated PK cleavage rule with its product constraints:

```yaml
- name: reversed polyketide synthesis (saturated)
  smarts: '[C,c:1][C;!R:2]-[C;!R:3]-[C:4](=[O:5])[OH:6]>>[C,c:1]C(=O)[OH].[OH][S][C:2]-[C:3]-[C:4](=[O:5])[OH:6]'
  allowed_in_bulk: false
  product_max_matches:
    '[0S]-[0OH]': 1
    '[0C](=[0O])-[0OH]': 1
  props:
    note: Keep products from accumulating competing artificial chain ends.
```

### Reaction SMARTS and atom tracking

`>>` separates reactant and product templates. RetroMol supplies one reactant molecule to a rule. Map atoms that must survive using the same map number on both sides, for example `[C:1]`. Unmapped atoms introduced on the product side are artificial additions, such as the `SO` cap or a continuation carboxyl group in the example above.

Each output molecule must be connected and chemically sanitizable. For intermolecular cleavage, separate product templates with `.`. For an intramolecular ring opening where both ends remain connected elsewhere, the existing rules use a grouped product such as `>>([fragment1].[fragment2])` so the disconnected template portions describe one output molecule. Follow an existing rule of the same kind and check the actual products.

The engine sanitizes products, restores tracked stereochemistry where possible, rejects products that share original atom tags, and deduplicates equivalent product sets. These checks do not guarantee conservation of every input atom. If a removed atom should contribute to the assembly, retain its map in an appropriate byproduct and provide a matching rule for that byproduct. Check explicit hydrogens and valence as well as connectivity when changing reaction products.

### What `product_max_matches` does

The keys are ordinary SMARTS queries; this is a general product filter, not a list of hardcoded chemical groups. All specified limits must pass. An empty mapping imposes no limits, and a limit of `0` forbids the pattern entirely. Counts use RDKit's unique substructure matches.

In the PK example:

- `[0S]-[0OH]` counts an S–OH group whose two atoms both have isotope zero: an artificial cap in the normal tagged parsing workflow.
- `[0C](=[0O])-[0OH]` counts a carboxyl group whose carbon and both oxygens were introduced by reactions.
- Each product may contain at most one of each. Two separate products can each have one cap and pass. Two caps in the same product fail.
- Original sulfur groups and original carboxyl groups carry positive atom tags and do not match these queries. A carboxyl group containing a mixture of original and added atoms also does not match the entirely artificial carboxyl query.

The acid limit prevents a second artificial continuation end from accumulating in the same fragment. The cap limit prevents accumulating multiple artificial unit ends. Together they reject certain routes that cut a PK chain from both ends and collide in the middle. They do not select a preferred direction: alternative directions can still be explored if their products pass the limits.

Limits apply only to products of the rule carrying them; they are not global constraints automatically applied to every other reaction. They are also conservative: some legitimate complex intermediates may need different constraints, so inspect rejected routes before relaxing them.

### When to use `allowed_in_bulk`

Use `true` for preprocessing transformations that can be applied eagerly when uncontested. A typical example is an independently cleavable linkage. Keep it `false` when competing reaction choices need to remain explicit alternatives. This setting can affect which routes are explored, not just their speed.

## Matching rules: `mxn.yml`

This file is also a YAML list.

| Field | Type / default | Meaning |
| --- | --- | --- |
| `name` | String; required | Identity assigned to the matched fragment and used in assembly graphs. Multiple structural entries can share a name, for example keto/enol representations of one unit. |
| `smiles` | String; required | The whole fragment structure to recognize, including any artificial caps expected after disconnection. This field takes SMILES, not a SMARTS query. |
| `requires_smarts` | SMARTS string or `null`; `null` | Additional substructure condition on the tagged fragment. Both this query and the whole-fragment `smiles` match must succeed. Compiled and validated when loading the rule. |
| `cases` | List; `[]` | Conditional `name` and/or `sequence_role` overrides after the structural match. Each case requires `requires_smarts`; the first matching case wins, otherwise the rule's defaults apply. |
| `pseudonyms` | List of strings; `[]` | Alternative names or broader groups, such as `PK_A` and `PK`, used by downstream representations and BGC matching. They do not change the structural match. |
| `terminal` | Boolean; `true` | Stop expanding a fragment identified by this rule. With `false`, keep its identity but also allow further reactions and finer assemblies. |
| `stereochemistry` | Boolean; `false` | Marks a stereo-specific rule. These rules are prioritized when global stereochemical matching is enabled and skipped when it is disabled. |
| `preserve_across_routes` | Boolean; `true` | For terminal matches, also protect the original atom region against subdivision through alternative routes. See below. |
| `display_smiles` | String or `null`; `null` | Optional alternative structure for depiction. It does not replace `smiles` during matching or change coverage. |
| `sequence_role` | `backbone`, `tailoring`, or `null`; `null` | Role of this monomer in linear assembly readouts. It does not change matching, expansion, terminal protection, or coverage. See below. |
| `props` | Mapping; `{}` | Descriptive metadata, such as a reference or an explanation of a modified residue. |

For example:

```yaml
- name: A1
  smiles: 'O=C(O)CC(=O)SO'
  display_smiles: 'O=C(CC([1*])=O)S[2*]'
  pseudonyms: [PK_A, PK]
  terminal: true
  stereochemistry: false
  preserve_across_routes: false
  props:
    note: Formal capped beta-keto PK unit.
```

### Conditional names and sequence roles

Define the structure once, then use `cases` to assign a name and role based on
the matched fragment's existing atom tags. The default PK cleavage rules add a
complete carboxyl cap to the remaining chain. Reversing thiazole formation instead
retains the original acyl carbon and restores its carbonyl oxygen; subsequent
amide cleavage adds the acid hydroxyl. Both starter representations have two
reaction-added acid oxygens (tag zero). Ordinary O- or N-acetyl substituents
retain their original carbonyl oxygen (a positive tag):

```yaml
- name: acetic acid
  smiles: 'CC(=O)O'
  sequence_role: tailoring
  cases:
    - requires_smarts: '[C](=[0O])[0OH]'
      name: acetic acid starter
      sequence_role: backbone
```

### Backbone and tailoring roles

Set `sequence_role` on a monomer matching rule when its role in the sequence view
is known:

```yaml
- name: methylation
  smiles: C
  sequence_role: tailoring
```

`tailoring` means a separately represented modification that should be read
outside the backbone path, even if attached through C–C or C–N. Connected tailoring
units form their own paths, including multi-unit moieties. A modification already
contained in a matched monomer, such as N-methylalanine, stays part of that monomer;
its name alone does not make the whole unit tailoring. This is a representation
policy, not evidence that a reaction happened after backbone synthesis.

Connections between two explicit `backbone` units are initially retained regardless
of bond elements. The readout then distinguishes backbone continuation from side
attachments using connection sites in the fragments: an incoming carbon marked
by an artificial S–OH cap, or an alpha-amino acid's backbone nitrogen. If that
entry is already connected to a preceding unit's carboxyl side, an extra
ester/thioester connection to the unit's O/S is omitted
from traversal. This covers ring closures and ester attachments between separate
backbones; both chains remain backbone sequences. Extra ether links to the already
connected incoming carbon are also excluded, so they cannot redirect the path.
Esters without that competing entry evidence remain eligible, including ordinary
ester-linked backbone steps. All original bonds remain in the assembly graph.
For example, enterobactin's three N-acylated serines yield three identical
`['2,3-dihydroxybenzoic acid', 'serine']` paths; their side-chain ester links do
not merge those paths. Repeated parts are retained as separate occurrences.

Where either non-tailoring endpoint has no role, the readout falls
back to C–C/C–N connections. A `tailoring` endpoint always separates from an untagged
or `backbone` endpoint. The defaults explicitly mark PK units and conditional
acid starters as backbone and modification identities as tailoring. Other units
can remain untagged; no role is inferred from size or degree.
Role choices can depend on the intended compound class, so use custom matching
rules where a motif plays a different role.

Roles are saved with matching rules in result JSON. Older results without the
field remain readable and use the fallback; reparse with annotated rules to
apply the new defaults to those compounds.

`assembly_to_primary_sequence(assembly)` returns the resulting paths as
`list[list[str]]`. It retains branches in additional paths and does not label an
unannotated leftover branch as tailoring. Each identified monomer occurs once;
unresolved fragments and the synthetic unassigned node are omitted. Unresolved
regions break paths. Assembly coverage measures how much of the compound these
named units represent: 100% means all submitted heavy atoms have a monomer
identity. See the [Python example](../../../README.md#linear-sequences-from-sampled-assemblies).

### `terminal` versus `preserve_across_routes`

`terminal` controls expansion of the matched fragment itself. `preserve_across_routes` controls whether that terminal match also restricts how **other derivations** divide the same original atoms.

| Settings | Effect |
| --- | --- |
| `terminal: true`, `preserve_across_routes: true` | Stop at this match and protect its atom region from finer subdivision through other routes. This is the default. |
| `terminal: true`, `preserve_across_routes: false` | Stop at this match in its own route, while allowing independently reached routes to divide the region differently. |
| `terminal: false` | The recognized fragment can be retained or expanded further. It creates no terminal-region protection, regardless of `preserve_across_routes`. |

For example, if a terminal match covers atoms `A+B`, preservation normally excludes an alternative assembly that splits them into separate `A` and `B` units. Setting preservation to `false` allows that subdivision if another route reaches it; it does not allow expanding the terminal `AB` node itself.

Preservation is not a “first route wins” switch. Crossing terminal assignments can compete: `AB + C` and `A + BC` may both be valid, even when both `AB` and `BC` are protected. The fully subdivided `A + B + C` is excluded. Keeping the protected region within a larger selected fragment is also allowed.

The default PK units use `preserve_across_routes: false` because a unit found in one decomposition should not veto a different valid decomposition reached by reversing other cyclizations or tailoring steps.

### Matching order and stereochemistry

The parser records the first eligible structural match for a given fragment. It does not collect every matching rule for that same fragment. This is separate from retaining multiple reaction routes and assemblies.

When loading files, stereo-specific matching rules are placed first, preserving file order within the stereo-specific and generic groups. With global `match_stereochemistry=False`, stereo-specific rules are skipped. With it enabled, matching uses chirality and tries the stereo-specific entries before the generic fallbacks. Put more specific or preferred representations before competing entries within each group.

Suffixes such as `^R`, `^S`, `^L`, and `^D` are naming conventions; they do not encode chemistry on their own. The `smiles` and stereochemistry settings determine the structural match. Pseudonyms do not create alternative structural assignments.

## Loading, serialization, and compatibility

```python
from retromol.model.rules import RuleSet

rules = RuleSet.load_from_files(
    reaction_rules_path="my_rxn.yml",
    matching_rules_path="my_mxn.yml",
    match_stereochemistry=False,
)
```

Omit either path, or pass `None`, to use that packaged default file. Supplying a file replaces that category of rules; it does not append to the defaults. The CLI exposes the paths as `--rxn-rules` and `--mxn-rules`.

Use YAML booleans `true` and `false`, not quoted strings. Quote SMILES/SMARTS; single quotes are convenient for backslashes. The product limits must be nonnegative integers, and `preserve_across_routes` must be an actual boolean. Unknown top-level keys are currently ignored, so check spelling carefully.

Older YAML files and saved results with `props.preserve_across_routes` or `props.product_max_matches` are still accepted through `from_dict`. Loading moves these values out of metadata. If both locations occur, the explicit top-level value wins, including `true` or an empty `{}`. Serialization always writes the new top-level fields. When constructing Python rule objects directly, pass these settings as named constructor arguments.

Reaction rule IDs hash their SMARTS. Matching rule IDs hash their SMILES and,
when present, the top-level `requires_smarts`. Cases choose output identities
without changing this structural ID, so variants of one rule share an ID.
Matching rules without a top-level condition retain their SMILES-only IDs.
Ruleset hashes include the ordered serialized definitions, including behavior
settings, names, cases, and metadata. Adding serialized fields therefore changes
ruleset hashes even when effective behavior is unchanged; old saved provenance
remains historical. Treat loaded rules and their nested mappings as immutable.

## Prediction mapping: `pmp.yml`

This file is a mapping with two top-level sections, rather than a list of chemical rules:

- `predictors.<group>.<axis>` selects a method name.
- `methods.<group>.<axis>.<method_name>` defines that method.

The configured axes are `pks.extender`, `pks.starter`, `pks.kr_stereochemistry`, and `nrps.a_domain`. For example, the default extender method maps an antiSMASH `mmal` call to the digit `"2"`, which combines with the module's reduction state to produce `A2`, `B2`, `C2`, or `D2`.

| Method field | Type / default | Meaning |
| --- | --- | --- |
| `source` | `qualifier` or `model`; `qualifier` | Read a GenBank qualifier or a registered inference model's top-scoring domain result. |
| `feature` | String; required for `qualifier` | Declares the expected feature type, normally `aSDomain`. |
| `domain` | String; required for `qualifier` | Declares the expected domain type, such as `PKS_AT`, `PKS_KR`, or `AMP-binding`. |
| `qualifier` | String; required for `qualifier` | Qualifier to read, normally `specificity`. |
| `key` | String; required for `qualifier` | Text before `: ` identifying a particular qualifier value, such as `substrate consensus`. |
| `model` | String; required for `model` | Registered inference model name, such as `paras_cli`. Adding a name here does not implement a new model. |
| `combine_with` | String or `null`; `null` | Token convention: `reduction_letter` for a PK extender digit, `stereo_suffix` for a KR beta-hydroxyl stereo suffix, or `null` for a complete monomer name. |
| `mapping` | String-to-string mapping; `{}` | Reword predicted values into names, pseudonyms, or token fragments. Quote numeric-looking values such as `"2"`. |

An empty or incomplete ordinary mapping leaves unmapped labels unchanged. NRPS model predictions can fall back to structural matching when they provide a SMILES and their label cannot be resolved; qualifier predictions carry no such structure. KR `stereo_suffix` mapping is different: an absent mapping entry means no suffix, not passthrough of raw labels such as `A1`. Current KR mapping only assigns the beta-hydroxyl stereo of B-series units, not the alpha stereocenter.

Current implementation limits matter when editing this file: module collection chooses the relevant domains itself, so `feature` and `domain` are required declarations, not general configurable feature filters. PK extender/KR collection currently consumes qualifier methods; NRPS A-domain collection supports both sources. `combine_with` is not a general transformation pipeline. The `pks.starter` axis is present in the configuration but is not currently consumed by module collection. NRPS epimerization and N-methylation are derived from module anatomy, outside these predictor settings.

Load a custom mapping with `PredictionConfig.load_from_file(path)` from `retromol_antismash.predictions` and pass that configuration to the BGC processing code. See the inline comments in [pmp.yml](pmp.yml) for the supplied predictor mappings.

## Adding or changing rules

Use a reaction rule when a general chemical transformation is missing, and a matching rule when an appropriate generated fragment lacks an identity. Often both are needed. Adding a matching rule for an entire unresolved natural product can raise coverage without explaining its assembly, so prefer reusable units supported by the intended chemistry.

For chemical regression coverage, use the existing files:

1. Add the compound's name and SMILES to [tests/data/integration_retromol.py](../../../tests/data/integration_retromol.py).
2. Add one or more independently checked expected assembly graphs under the same name in [tests/data/integration_retromol_assemblies.yml](../../../tests/data/integration_retromol_assemblies.yml). Nodes specify identity and original-atom count; edges specify node indices and original connecting bond types. Include expected coverage.
3. Every graph listed for that compound must occur among the parsed assemblies. If you list two of five valid assemblies, both must match; extra parsed assemblies are allowed. Do not blindly regenerate expectations from the implementation being tested.
4. Check a nearby negative example for broad reaction patterns, atom conservation, and unintended changes in existing compounds. Behavioral settings are covered in the existing parsing and assembly-sampling tests.

From the repository root, with the RetroMol environment active:

```sh
python -m pytest tests/test_integration_retromol.py tests/test_parsing.py tests/test_assembly_sampling.py
python -m pytest tests
```

The implementation lives in [rules.py](../model/rules.py), [parsing.py](../pipelines/parsing.py), and [assembly_sampling.py](../model/assembly_sampling.py).
