# Convex Dice Design (convex_dice)

Compact reference for the convex dice optimization tool (convex_dice executable).
Examples of section 5.1 of the [paper](https://hbaktash.github.io/projects/putting-rigid-bodies-to-rest/index.html) are provided as an example in ./examples .

convex_dice runs an optimization that adjusts a convex shape (vertex positions) to match a target distribution of resting orientations (a "normal → probability" assignment).

## Build

Requirements (same as main README)
- C++20, CMake ≥ 3.16, Eigen3, Qhull
- (Optional) Polyscope for the GUI

Typical build (with GUI):
```
cmake -S . -B build -DRESTING_RIGID_BODIES_USE_POLYSCOPE=ON -DRESTING_RIGID_BODIES_ENABLE_INVERSE_DESIGN=ON
cmake --build build -j
```

Note: Polyscope and inverse-design features must be enabled in CMake to build convex_dice.

## Usage

Binary: build/convex_dice

Basic CLI options:
- --mesh, -m <path>  
  Path to input mesh.

- --params <file>  
  Load optimization parameters and regularized coefficients from JSON (see Parameters section). Missing entries use defaults.

- --normal_probs <file>  
  Load a normal→probability assignment (JSON). The assignment is used as the initial goal for resting probabilities of the provided orientations. 

- --out, -o <path>  
  Output path for the optimized OBJ. Parent directory will be created if missing. Several auxiliary files are written to the same directory (params.json, G_*.txt, normal-probs JSONs).

Upon execution a polyscope window opens, and the "dice energy opt" button will start the optimization. When done, the final optimized shape will be shown and can be saved (along with the mentiond log files) to the output directory provided.

Example:

```
# Load params and normal probs, write optimized mesh to a custom folder
./convex_dice -m ../opt_release/2D6_hendecahedron_V2/hendecahedron.obj --out ../opt_release/2D6_hendecahedron_V2/optimized_dice.obj --normal_probs ../opt_release/2D6_hendecahedron_V2/initial_normal_probs.json --params ../opt_release/2D6_hendecahedron_V2/params.json
```

## Parameters (params.json)

Parameters are grouped in JSON under `optimization` and `sobolev`; the sobolev section is only used for large convex hulls (kitten, armadillo.. models). The code expects the same structure produced by the UI save function. Example minimal structure:
```json
{
  "optimization": {
    "bary_reg": 0.1,
    "coplanar_reg": 0.0,
    "cluster_distance_reg": 0.0,
    "unstable_attraction_thresh": 0.1,
    "dice_energy_step": 0.01,
    "dice_search_decay": 0.98,
    "DE_step_count": 40,
    "fair_sides_count": 6,
    "frozen_G": false
  },
  "sobolev": {
    "do_sobolev_dice_grads": false,
    "sobolev_lambda": 2.0,
    "sobolev_lambda_decay": 0.8,
    "sobolev_p": 2
  }
}
```
Use `--params` to load this JSON. The GUI "save optimized hull" also writes `params.json` to the output directory.

## Normal→Probability JSON

The normal-probability JSON format used by load/save utilities contains:
- "policy": optional string describing the policy used to produce the assignment; 
so far the only supported policy is `manualCluster` which is the clustering scheme mentioned in section 3 of the [paper](https://hbaktash.github.io/projects/putting-rigid-bodies-to-rest/index.html).
- "pairs": array of objects { "normal": [x,y,z], "probability": p }

Example:
```json
{
  "type": "normal_probability_assignment",
  "policy": "manualCluster cube",
  "pairs": [
    { "normal": [ 0, 0, -1 ], "probability": 0.1667 },
    { "normal": [ 0, 0,  1 ], "probability": 0.1667 }
  ]
}
```
When loaded via `--normal_probs`, the assignment is used directly as the optimization target.

## GUI Controls (Polyscope)

When built with Polyscope, a small UI is available:
- Sliders for iterations, step size, fair sides, and regularizers
- Checkboxes for Sobolev smoothing, frozen G (Center of Mass)
- Buttons:
  - "dice energy opt" — run optimization with current settings
  - "save optimized hull" — write optimized OBJ + params.json + G file + normal-probs JSONs to the output directory

## Outputs

When you press "save optimized hull" or specify `--out`, the program writes:
- optimized OBJ (as specified by `--out`)
- params.json (current optimization parameters) in same directory
- G_<basename>.txt (center of mass) if `frozen_G` is true
- final_normal_probs.json and initial_normal_probs.json (assignments used)

Directories are created automatically when needed.

## Notes & Tips

- If you provide a normal-probability file, its internal "policy" string may override the CLI `--policy` interpretation.
- Parameters not present in a loaded JSON will use defaults.
- For reproducible runs save the params.json and the normal-probs JSON alongside the optimized output.
- If you encounter missing includes (e.g. stan/math.hpp), ensure inverse-design dependencies are available and the inverse-design CMake option is enabled.

License and citation follow the project README.