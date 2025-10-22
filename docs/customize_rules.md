# Custom rule file formats

Below are the specifications for custom reaction and matching rule files, as well as wave configuration file, used by RetroMol. These files are written in YAML format. Reaction and matching rule files are dumps of lists of rules, while the configuration file structures them into a workflow.

## Wave configuration

### Creating a wave_config.yml

RetroMol processes molecules in waves: ordered stages that (1) apply selected reaction rules to split or transform the current “frontier” structures, and then (2) match the resulting fragments to motif classes. The behavior of each stage is controlled by a YAML file, wave_config.yml.

If you don’t pass a config, RetroMol uses `retromol/data/default_wave_config.yml` (the same structure as the example below). You can point the CLI at a custom file with --wave-config (see `src/retromol/cli.py`) or pass it via the API (api.run_retromol(..., wave_configs=...)).

### File structure

A wave config is a YAML list, where each item is one wave. Waves run top-to-bottom. Each wave supports the keys below.

Supported keys (per wave):
- `wave_name` (`string`, required): Human-readable label stored on nodes produced in that wave (see `resolve_mol` in `src/retromol/apply.py`).

- `reaction_groups` (`list[string]`, required): Names of reaction-rule groups to apply in this wave. These must correspond to groups defined in your reaction rules YAML (loaded via `rules.load_rules_from_files(...)`). Examples: preprocessing, linearization, NRP disassembly, PK disassembly, etc.

- `matching_groups` (`list[string]`, optional): If provided, only matching rules whose groups intersect this list are considered when assigning node identities in this wave (see filtering in `apply.resolve_mol`). If omitted, all matching rules are eligible. Names correspond to groups defined in your matching rules YAML (loaded via `rules.load_rules_from_files(...)`). Examples: amino acid, polyketide building block, glycosylation, etc.

- `only_leaf_nodes` (`bool`, optional; default: true): After applying the reaction rules in this wave, RetroMol selects which reaction-graph nodes to turn into motif graph nodes:
  * `true`: use leaf nodes created in the previous wave only (preferred; avoids duplicating intermediate steps).
  * `false`: use all nodes created in the previous wave.

  Note: If a “leafs only” wave produces no leaves (e.g., an uncontested rule returns itself), RetroMol falls back to all nodes for that wave.

- `parse_identified_nodes` (`bool`, optional; default: `false`):
Controls eligibility of nodes for re-parsing at the current deepest nesting level (see `find_eligible_nodes` in `src/retromol/api.py`):
  * `false`: skip nodes that already have an identity.
  * `true`: allow re-parsing nodes even if they already carry an identity. This is useful when one wave (e.g., linearization) changes structures in a way that reveals new tailoring motifs for the next wave.


### Execution model (what happens under the hood)

- Wave 1 operates on the input molecule to produce the initial motif graph (see `api.run_retromol` > first call to `apply.resolve_mol`).

- Subsequent waves operate on the current frontier: nodes at the deepest nesting level whose graph is `None` (and, unless `parse_identified_nodes: true`, also have no identity).

- In each wave:
  * Reaction rules are chosen from `reaction_groups`.
  * The reaction graph is computed; RetroMol selects leaf or all nodes according to `only_leaf_nodes`.
  * Selected nodes are matched to motifs using the (optionally filtered) `matching_groups`.
  * Each selected node becomes a node in the motif graph and is annotated with `wave_name`, `identity`, `props`, `smiles`, `smiles_no_tags`, and `tags`.
  * For nodes that should be further expanded in later waves, a nested graph is attached.

Note: Stereochemistry matching is controlled separately (CLI `--matchstereochem` or API `match_stereochemistry`), not in the wave config.

### Minimal example

A minimal file needs just a name and at least one reaction group:

```yaml
- wave_name: preprocessing
  reaction_groups:
    - preprocessing
```

This will apply the preprocessing reaction rules and match against all matching rules.

## Reaction rules

The default reaction rules can be found at `retromol/data/default_reaction_rules.yml`.

A minimal reaction rules file looks like:

```yaml
# reactions.yml
- rid: reverse O-methylation
  smarts: "[O;D2:1][CH3:2]>>[O:1].[C:2]"
  groups: [preprocessing]
  props: {} # optional

- rid: break ester bond (intermolecular)
  smarts: "[C:1][C;!R:2](=[O:3])[O;!R:4][C:5]>>[C:1][C:2](=[O:3])[OH].[OH:4][C:5]"
  groups: [linearization]
```

```yaml
- rid: <string>                  # unique human-readable identifier
  smarts: <reaction SMARTS>      # RDKit reaction SMARTS (LHS>>RHS)
  groups: [<string>, ...]        # arbitrary grouping labels
  props:                         # optional metadata and global conditions
    conditions:
      # pre-conditions (whole molecule)
      reactant:
        requires_any:  ["<SMARTS>", ...]
        requires_all:  ["<SMARTS>", ...]
        forbids_any:   ["<SMARTS>", ...]
        min_counts:    {"<SMARTS>": <int>, ...}
        max_counts:    {"<SMARTS>": <int>, ...}
        ring_count:    {min: <int>, max: <int>}
        atom_count:    {min: <int>, max: <int>}
        total_charge:  {min: <int>, max: <int>}
        custom_props:  {has_metal: <bool>, is_macrocycle: <bool>}
      # post-conditions (applied to each product)
      product:
        requires_any:  ["<SMARTS>", ...]
        requires_all:  ["<SMARTS>", ...]
        forbids_any:   ["<SMARTS>", ...]
        min_counts:    {"<SMARTS>": <int>, ...}
        max_counts:    {"<SMARTS>": <int>, ...}
        ring_count:    {min: <int>, max: <int>}
        atom_count:    {min: <int>, max: <int>}
        total_charge:  {min: <int>, max: <int>}
        custom_props:  {has_metal: <bool>, is_macrocycle: <bool>}
```

What each field means:
* `rid`: A stable, unique identifier (used for logging/debugging).
* `smarts`: Reaction SMARTS using RDKit’s format, including mapped atoms on both sides. Element identities and multiplicities for each map number must match (the loader enforces this).
* `groups`: Tags used to select rule subsets (preprocessing, linearization, NRP disassembly, etc.).
* `props.conditions.reactant`: “Global pre-filter” — the reaction only runs if the whole reactant satisfies these.
* `props.conditions.product`: “Global post-filter” — each product must satisfy these (result dropped otherwise).

Tip: Use local constraints inside the SMARTS for atom-level logic (e.g., “:1 must not be acyl”), and global conditions for whole-molecule predicates (e.g., “forbid nitro anywhere”).

## Matching rules

The default matching rules can be found at `retromol/data/default_matching_rules.yml`.

A minimal matching rules file looks like:

```yaml
- rid: valine
  mol: "CC(C)C(N)C(=O)O"
  groups: [amino_acids]
  props: {}
```

Matching uses exact topology equality (same atom/bond counts) plus substructure match.

If you need stereochemistry matching, it’s supported by adding ste_mols in code (optional advanced usage).

```yaml
- rid: <string>
  mol: "<SMILES>"  # exact structure to match (full match, not subgraph)
  groups: [<string>, ...]
  props: {}
```
